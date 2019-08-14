from typing import List, Dict, Union

import numpy as np
import tqdm
import xarray as xr
from numpy.polynomial.polynomial import Polynomial
from scipy import ndimage as ndi
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
from scipy.spatial.distance import cdist
from skimage import measure, transform
from skimage.measure import label
from skimage.transform import warp, AffineTransform

from pharynx_analysis import profile_processing


def get_lr_bounds(rot_seg_stack: xr.DataArray, pad: int = 0, ref_wvl: str = '410', ref_pair: int = 0) -> np.ndarray:
    """
    Get the Left and Right boundaries of the rotated pharynxes
    Parameters
    ----------
    rot_seg_stack
        the rotated segmented pharynxes
    pad
        the amount of padding on the left/right of the  bounds
    ref_wvl
        the wavelength to use as reference
    ref_pair
        the pair to use as reference

    Returns
    -------
    bounds
        An (m, 2) array where m = number of animals,  the  first column is the left bound and the second column is the right
        bound
    """
    imgs = rot_seg_stack.sel(wavelength=ref_wvl, pair=ref_pair)
    bounds = np.zeros((imgs.strain.size, 2))  # (0, (x, y))
    for i, img in enumerate(imgs):
        l, _, _, r = measure.regionprops(measure.label(img))[0].bbox
        bounds[i, :] = [l - pad, r + pad]
    return bounds


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


def calculate_midlines(rot_seg_stack: xr.DataArray, degree: int = 4) -> List[Dict[str, List[Polynomial]]]:
    """Calculate the midlines for the given stack. Only calculates midlines for NON-TL images
    Parameters
    ----------
    rot_seg_stack
        The rotated mask with which midlines should be calculated. It should have the following dimensions::

            (strain, wavelength, pair, height, width)
    degree
        The degree of the polynomial fit

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
            wvl: [
                calculate_midline(rot_seg_stack.isel(strain=img_idx).sel(wavelength=wvl, pair=pair), degree=degree)
                for pair in rot_seg_stack.pair.data
            ]
            for wvl in rot_seg_stack.wavelength.data if 'tl' not in wvl.lower()
        }
        for img_idx in tqdm.trange(rot_seg_stack.strain.size)
    ]


def calculate_midline(rot_seg_img: Union[np.ndarray, xr.DataArray], degree: int = 4, pad: int = 5) -> Polynomial:
    """Calculate a the midline for a single image. Right now this only works for images that have been centered and aligned
    with their anterior-posterior along the horizontal.

    Parameters
    ----------
    pad
    degree
    rot_seg_img: Union[np.ndarray, xr.DataArray]
        The rotated masked pharynx image

    Returns
    -------
    UnivariateSpline
        An approximation of the midline of the pharynx for this frame

    """
    rp = measure.regionprops(measure.label(rot_seg_img))[0]
    xs, ys = rp.coords[:, 1], rp.coords[:, 0]
    left_bound, _, _, right_bound = rp.bbox
    # noinspection PyTypeChecker
    return Polynomial.fit(xs, ys, degree, domain=[left_bound - pad, right_bound + pad])


def measure_under_midline(fl: xr.DataArray, mid: Polynomial, xs: np.ndarray, thickness: float = 0.0) -> np.ndarray:
    """
    Measure the intensity profile of the given image under the given midline at the given x-coordinates.

    Parameters
    ----------
    fl
        The fluorescence image to measure
    mid
        The midline under which to measure
    xs
        The x-coordinates to evaluate the midline at
    thickness
        The thickness of the line to measure under. WARNING: this is a lot slower right now

    Returns
    -------
    zs: np.ndarray
        The intensity profile of the image measured under the midline at the given x-coordinates.

    """
    if thickness == 0:
        ys = xr.DataArray(mid(xs), dims='z')
        xs = xr.DataArray(xs, dims='z')
        return fl.interp(x=xs, y=ys).data.T
    else:
        ys = mid(xs)
        der = mid.deriv()
        normal_slopes = -1 / der(xs)
        normal_thetas = np.arctan(normal_slopes)

        mag = 1.5
        x0 = np.cos(normal_thetas) * mag
        y0 = np.sin(normal_thetas) * mag

        x1 = np.cos(normal_thetas) * -mag
        y1 = np.sin(normal_thetas) * -mag

        xs0 = xs + x0
        xs1 = xs + x1
        ys0 = ys + y0
        ys1 = ys + y1

        prof = []
        for i in range(len(xs)):
            line = measure.profile._line_profile_coordinates((xs0[i], ys0[i]), (xs1[i], ys1[i]))[:, :, 0]
            line_xs = xr.DataArray(line[0], dims='z')
            line_ys = xr.DataArray(line[1], dims='z')
            prof.append(np.mean(fl.interp(x=line_xs, y=line_ys).data, 0))
        # TODO: optionally use gaussian weight for the profile
        return np.array(prof)


def measure_under_midlines(fl_stack: xr.DataArray, midlines: List[Dict[str, List[Polynomial]]], n_points: int = 300,
                           frame_specific: bool = False) -> xr.DataArray:
    """
    Parameters
    ----------
    fl_stack
        The fluorescence stack under which to measure
    midlines: dict
        A list of dictionaries of midlines with the following structure:
        ```
            [
                {'410': [UnivariateSpline, UnivariateSpline, ...]},
                {'470': [UnivariateSpline, UnivariateSpline, ...]},
                ...
            ]
        ```
    n_points: int
        the number of points to sample under the midline
    frame_specific
        whether to use a different midline for each frame. if False, a single midline will be used within all
        wavelengths in a pair

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

    ref_wvl = '410'

    for img_idx in tqdm.trange(fl_stack.strain.size):
        for pair in range(fl_stack.pair.size):
            for wvl_idx, wvl in enumerate(raw_intensity_data.wavelength.data):
                img = fl_stack.sel(wavelength=wvl, pair=pair).isel(strain=img_idx)

                if frame_specific:
                    mid = midlines[img_idx][wvl][pair]
                else:
                    mid = midlines[img_idx][ref_wvl][pair]

                # Why do we use the bounds for the 410 midlines?
                ref_mid = midlines[img_idx]['410'][pair]
                xs = np.linspace(ref_mid.domain[0], ref_mid.domain[1], n_points)

                raw_intensity_data[img_idx, wvl_idx, pair, :] = measure_under_midline(img, mid, xs)

    raw_intensity_data.values = np.nan_to_num(raw_intensity_data.values)
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
    data = intensity_data

    ref_data = data.sel(wavelength=reference_wavelength, pair=reference_pair)
    ref_profile = ref_data.isel(strain=0).data

    ref_vecs = np.tile(ref_profile, (data.strain.size, 1))
    unflipped = data.sel(wavelength=reference_wavelength, pair=reference_pair).data
    flipped = np.fliplr(unflipped)

    should_flip = cdist(ref_vecs, unflipped, 'cosine')[0, :] > cdist(ref_vecs, flipped, 'cosine')[0, :]

    intensity_data[should_flip] = np.flip(intensity_data[should_flip], axis=3)

    mean_intensity = profile_processing.trim_profile(
        np.mean(intensity_data.sel(wavelength=reference_wavelength, pair=reference_pair), axis=0).data, threshold=2000,
        new_length=100)

    peaks, _ = find_peaks(mean_intensity, distance=.2 * len(mean_intensity), prominence=200, wlen=10)

    if len(peaks) < 2:
        return intensity_data

    if peaks[0] < len(mean_intensity) - peaks[1]:
        intensity_data = np.flip(intensity_data, axis=3)

    return intensity_data


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


def shift(image, vector):
    transform = AffineTransform(translation=vector)
    shifted = warp(image, transform, mode='wrap', preserve_range=True)

    shifted = shifted.astype(image.dtype)
    return shifted
