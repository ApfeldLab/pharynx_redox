from typing import List, Dict, Union

import numpy as np
import tqdm
import xarray as xr
from scipy import ndimage as ndi
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
from scipy.spatial.distance import cdist
from skimage import measure, transform
from skimage.measure import label


def center_and_rotate_pharynxes(fl_images, seg_images, reference_wavelength='410') -> (np.ndarray, np.ndarray):
    """
    Given a fluorescence stack and a pharyngeal mask stack, center and rotate each frame of both the FL and mask such
    that the pharynx is in the center of the image, with its anterior on the left

    Parameters
    ----------
    fl_images
    seg_images
    reference_wavelength

    Returns
    -------

    """
    img_center_y, img_center_x = (fl_images.y.size // 2, fl_images.x.size // 2)  # (y, x)

    fl_rotated_stack = fl_images.copy()
    seg_rotated_stack = seg_images.copy()

    blurred_seg = fl_images.copy()
    blurred_seg_data = ndi.gaussian_filter(fl_images, sigma=(0, 0, 0, 6, 6))
    blurred_seg_data = blurred_seg_data > 1000
    blurred_seg.data = blurred_seg_data

    for img_idx in tqdm.trange(fl_images.strain.size):
        for wvl in fl_images.wavelength.data:
            for pair in fl_images.pair.data:
                # Optimization potential here...
                # this recalculates all region properties for the reference each time
                reference_seg = blurred_seg.isel(strain=img_idx).sel(wavelength=reference_wavelength, pair=pair)
                img = fl_images.isel(strain=img_idx).sel(wavelength=wvl, pair=pair)
                seg = seg_rotated_stack.isel(strain=img_idx).sel(wavelength=wvl, pair=pair)

                props = measure.regionprops(measure.label(reference_seg), coordinates='rc')[0]

                # pharynx_center_y, pharynx_center_x = props.centroid
                pharynx_center_y, pharynx_center_x = np.mean(np.nonzero(reference_seg), axis=1)
                pharynx_orientation = props.orientation

                translation_matrix = transform.EuclideanTransform(
                    translation=(-(img_center_x - pharynx_center_x), -(img_center_y - pharynx_center_y))
                )

                rotated_img = rotate(img.data, translation_matrix, pharynx_orientation)
                rotated_seg = rotate(seg.data, translation_matrix, pharynx_orientation)

                fl_rotated_stack.loc[dict(wavelength=wvl, pair=pair)][img_idx] = rotated_img
                seg_rotated_stack.loc[dict(wavelength=wvl, pair=pair)][img_idx] = rotated_seg

    return fl_rotated_stack, seg_rotated_stack


def extract_largest_binary_object(bin_img):
    labels = label(bin_img)
    if labels.max() == 0:
        # No connected components (TL images)
        return bin_img
    return labels == np.argmax(np.bincount(labels.flat)[1:]) + 1


# noinspection PyUnresolvedReferences
def segment_pharynxes(fl_stack: xr.DataArray, threshold=2000):
    """

    Parameters
    ----------
    threshold
    fl_stack

    Returns
    -------

    """
    seg = fl_stack > threshold
    for img_idx in range(fl_stack.strain.size):
        for wvl_idx in range(fl_stack.wavelength.size):
            for pair in range(fl_stack.pair.size):
                seg_img = seg[dict(strain=img_idx, wavelength=wvl_idx, pair=pair)]
                seg[dict(strain=img_idx, wavelength=wvl_idx, pair=pair)] = extract_largest_binary_object(seg_img)
    return seg


def get_centroids(fl_stack: xr.DataArray, threshold=1000, gaussian_sigma=6):
    image_data = fl_stack.copy()
    image_data.data = ndi.gaussian_filter(image_data.data, sigma=(0, 0, 0, gaussian_sigma, gaussian_sigma))
    image_data.data[image_data.data < threshold] = 0
    image_data.data[image_data.data > threshold] = 1


def rotate(data, tform, orientation):
    """

    Parameters
    ----------
    data
    tform
    orientation

    Returns
    -------

    """
    # noinspection PyTypeChecker
    return transform.rotate(
        transform.warp(
            data, tform, preserve_range=True, mode='wrap'),
        np.degrees(np.pi / 2 - orientation), mode='edge')


def calculate_midlines(rot_seg_stack: xr.DataArray, rot_fl_stack: xr.DataArray,
                       s=1e8, ext=0) -> List[Dict[str, List[UnivariateSpline]]]:
    """Calculate the midlines for the given stack. Only calculates midlines for NON-TL images
    Parameters
    ----------
    rot_fl_stack
    rot_seg_stack
        The rotated mask with which midlines should be calculated. It should have the following dimensions::

            (strain, wavelength, pair, height, width)
    s
        smoothing constraint
    ext
        extrapolation

    Returns
    -------
    list
        A list of dictionaries with the following structure::

            [
                {
                    wavelength0: [midline_pair_0, midline_pair_1, ...],
                    wavelength1: [midline_pair_0, midline_pair_1, ...]
                }
            ]
    """
    return [
        {
            wvl:
                [calculate_midline(rot_seg_stack.isel(strain=img_idx).sel(wavelength=wvl, pair=pair),
                                   rot_fl_stack.isel(strain=img_idx).sel(wavelength=wvl, pair=pair), s, ext) for pair
                 in rot_seg_stack.pair.data]
            for wvl in rot_seg_stack.wavelength.data if 'tl' not in wvl.lower()
        }
        for img_idx in tqdm.trange(rot_seg_stack.strain.size)
    ]


def calculate_midline(rot_seg_img: Union[np.ndarray, xr.DataArray], rot_fl: Union[np.ndarray, xr.DataArray], s=1e8,
                      ext=0):
    """Calculate a the midline for a single image. Right now this only works for images that have been centered and aligned
    with their anterior-posterior along the horizontal.

    Parameters
    ----------
    rot_fl
    rot_seg_img: Union[np.ndarray, xr.DataArray]
        The rotated masked pharynx image
    s
        The smoothing constraint for the spline
    ext
        The extrapolation method

    Returns
    -------
    UnivariateSpline
        An approximation of the midline of the pharynx for this frame

    """
    # TODO center segmentation before fitting splines
    seg_coords = measure.regionprops(measure.label(rot_seg_img))[0].coords

    # Adding tiny amounts of random noise to the coordinates because the spline fitting function expects that all
    # x-coordinates are unique.
    # TODO: change this; I think this could fail if two random numbers happen to be the same... could use a while loop checking for uniqueness
    seg_coords = seg_coords + (np.random.random(seg_coords.shape) / 100)
    seg_coords = seg_coords[np.argsort(seg_coords, axis=0)[:, 1]]

    xs = seg_coords[:, 1]
    ys = seg_coords[:, 0]

    try:
        return UnivariateSpline(xs, ys, s=s, ext=ext)
    except:
        return center_of_mass_midline(rot_fl.data, s=s, ext=ext)


def measure_under_midline(fl: xr.DataArray, mid: UnivariateSpline, xs: np.ndarray) -> np.ndarray:
    """
    Measure the intensity profile of the given image under the given midline at the given x-coordinates.

    Parameters
    ----------
    fl: np.ndarray
        The fluorescence image to measure
    mid: UnivariateSpline
        The midline under which to measure
    xs: np.ndarray
        The x-coordinates to evaluate the midline at

    Returns
    -------
    zs: np.ndarray
        The intensity profile of the image measured under the midline at the given x-coordinates.

    """
    # old way
    # print('measuring')
    ys = xr.DataArray(mid(xs), dims='z')
    xs = xr.DataArray(xs, dims='z')
    return fl.interp(x=xs, y=ys).data.T

    # # Use line thickness!
    # print('measuring')
    # ys = mid(xs)
    # der = mid.derivative()
    # normal_slopes = -1 / der(xs)
    # normal_thetas = np.arctan(normal_slopes)
    #
    # mag = 1.5
    # x0 = np.cos(normal_thetas) * mag
    # y0 = np.sin(normal_thetas) * mag
    #
    # x1 = np.cos(normal_thetas) * -mag
    # y1 = np.sin(normal_thetas) * -mag
    #
    # xs0 = xs + x0
    # xs1 = xs + x1
    # ys0 = ys + y0
    # ys1 = ys + y1
    #
    # prof = []
    # for i in range(len(xs)):
    #     line = measure.profile._line_profile_coordinates((xs0[i], ys0[i]), (xs1[i], ys1[i]))[:,:,0]
    #     line_xs = xr.DataArray(line[0], dims='z')
    #     line_ys = xr.DataArray(line[1], dims='z')
    #     prof.append(np.mean(fl.interp(x=line_xs, y=line_ys).data, 0))
    # # TODO: optionally use gaussian weight for the profile
    # return prof


def measure_under_midlines(fl_stack: xr.DataArray, midlines: List[Dict[str, List[UnivariateSpline]]],
                           x_range: tuple, n_points: int) -> xr.DataArray:
    """
    Parameters
    ----------
    fl_stack
        The fluorescence stack under
    midlines: dict
        A list of dictionaries of midlines with the following structure:
        ```
            [
                {'410': [UnivariateSpline, UnivariateSpline, ...]},
                {'470': [UnivariateSpline, UnivariateSpline, ...]},
                ...
            ]
        ```
    x_range: tuple
        the range at which to evaluate the midlines
    n_points: int
        the number of points to sample under the midline

    Returns
    -------

    """
    non_tl_wvls = list(filter(lambda x: x != 'TL', fl_stack.wavelength.data))
    raw_intensity_data = xr.DataArray(
        np.zeros(
            (fl_stack.strain.size, len(non_tl_wvls), fl_stack.pair.size, n_points)
        ),
        dims=['strain', 'wavelength', 'pair', 'position'],
        coords={'strain': fl_stack.strain, 'wavelength': non_tl_wvls, 'pair': fl_stack.pair}
    )

    xs = np.linspace(x_range[0], x_range[1], n_points)

    for img_idx in tqdm.trange(fl_stack.strain.size):
        for wvl_idx, wvl in enumerate(raw_intensity_data.wavelength.data):
            for pair in range(fl_stack.pair.size):
                img = fl_stack.sel(wavelength=wvl, pair=pair).isel(strain=img_idx)
                raw_intensity_data[img_idx, wvl_idx, pair, :] = \
                    measure_under_midline(img, midlines[img_idx][wvl][pair], xs)

    return raw_intensity_data


def align_pa(intensity_data, reference_wavelength='410', reference_pair=0):
    """

    Parameters
    ----------
    intensity_data
    reference_pair
    reference_wavelength

    Returns
    -------

    """
    data = trim_profiles(intensity_data, threshold=2000, new_length=100)

    ref_data = data.sel(wavelength=reference_wavelength, pair=reference_pair)
    ref_profile = ref_data.isel(strain=0).data

    ref_vecs = np.tile(ref_profile, (data.strain.size, 1))
    unflipped = data.sel(wavelength=reference_wavelength, pair=reference_pair).data
    flipped = np.fliplr(unflipped)

    should_flip = cdist(ref_vecs, unflipped, 'cosine')[0, :] > cdist(ref_vecs, flipped, 'cosine')[0, :]

    intensity_data[should_flip] = np.flip(intensity_data[should_flip], axis=3)

    mean_intensity = trim_profile(
        np.mean(intensity_data.sel(wavelength=reference_wavelength, pair=reference_pair), axis=0).data, threshold=2000,
        new_length=100)

    peaks, _ = find_peaks(mean_intensity, distance=.2 * len(mean_intensity), prominence=200, wlen=10)

    if len(peaks) < 2:
        return intensity_data

    if peaks[0] < len(mean_intensity) - peaks[1]:
        intensity_data = np.flip(intensity_data, axis=3)

    return intensity_data


def trim_profile(profile, threshold, new_length):
    """
    Parameters
    ----------
    profile
    threshold
    new_length

    Returns
    -------

    """
    first = np.argmax(profile > threshold)
    last = len(profile) - np.argmax(np.flip(profile > threshold))

    trimmed = profile[first:last + 1]
    new_xs = np.linspace(0, len(trimmed), new_length)
    old_xs = np.arange(0, len(trimmed))

    return np.interp(new_xs, old_xs, trimmed)


def get_trim_boundaries(data, ref_wvl='410', thresh=2000):
    prof_len = data.position.size
    l_bound = np.argmax(data.sel(wavelength=ref_wvl) >= thresh, axis=2).data - 1
    r_bound = prof_len - np.argmax(np.flip(data.sel(wavelength=ref_wvl), axis=2) >= thresh, axis=2).data
    return l_bound, r_bound


def trim_profiles(intensity_data, threshold, new_length, ref_wvl='410'):
    """
    Parameters
    ----------
    ref_wvl
    intensity_data
    threshold
    new_length

    Returns
    -------

    """
    trimmed_intensity_data = xr.DataArray(
        np.zeros(
            (intensity_data.strain.size, intensity_data.wavelength.size, intensity_data.pair.size, new_length)
        ),
        dims=['strain', 'wavelength', 'pair', 'position'],
        coords={'strain': intensity_data.strain, 'wavelength': intensity_data.wavelength, 'pair': intensity_data.pair}
    )

    l, r = get_trim_boundaries(intensity_data, ref_wvl=ref_wvl, thresh=threshold)

    for img_idx in range(intensity_data.strain.size):
        for wvl_idx in range(intensity_data.wavelength.size):
            wvl = intensity_data.wavelength.data[wvl_idx]
            if 'tl' not in wvl.lower():
                for pair in range(intensity_data.pair.size):
                    data = intensity_data.sel(wavelength=wvl, pair=pair).isel(strain=img_idx).data

                    trimmed = data[l[img_idx, pair]:r[img_idx, pair]]
                    new_xs = np.linspace(0, len(trimmed), new_length)
                    old_xs = np.arange(0, len(trimmed))
                    resized = np.interp(new_xs, old_xs, trimmed)

                    trimmed_intensity_data[img_idx, wvl_idx, pair, :] = resized

    return trimmed_intensity_data


def r_to_oxd(r, r_min=0.852, r_max=6.65, instrument_factor=0.171):
    """

    Parameters
    ----------
    r
    r_min
    r_max
    instrument_factor

    Returns
    -------

    """
    return (r - r_min) / ((r - r_min) + instrument_factor * (r_max - r))


def oxd_to_redox_potential(oxd, midpoint_potential=-265, z=2, temperature=22):
    """Convert OxD to redox potential

    NOTE: may return NaN

    Parameters
    ----------
    oxd
    midpoint_potential
    z
    temperature

    Returns
    -------

    """
    return midpoint_potential - (8314.462 * (273.15 + temperature) / (z * 96485.3415)) * np.log((1 - oxd) / oxd)


def center_of_mass_midline(rot_fl, s, ext):
    """ Calculate the midline using the Center of Mass method

    Parameters
    ----------
    ext
    s
    rot_fl

    Returns
    -------

    """
    ys = np.arange(rot_fl.shape[0])
    midline_ys = []
    xs = np.arange(rot_fl.shape[1])
    for i in xs:
        midline_ys.append(np.average(ys, weights=rot_fl[:, i].data))
    return UnivariateSpline(xs, np.array(midline_ys), s=s, ext=ext)
