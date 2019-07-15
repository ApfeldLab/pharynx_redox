from typing import Union, List, Iterable

import numpy as np
import xarray as xr
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
from skimage import measure, transform


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

    for img_idx in range(fl_images.strain.size):
        for wvl in fl_images.wavelength.data:
            for pair in fl_images.pair.data:
                # Optimization potential here...
                # this recalculates all region properties for the reference each time
                reference_seg = seg_images.isel(strain=img_idx).sel(wavelength=reference_wavelength, pair=pair)
                img = fl_images.isel(strain=img_idx).sel(wavelength=wvl, pair=pair)
                seg = seg_images.isel(strain=img_idx).sel(wavelength=wvl, pair=pair)

                props = measure.regionprops(measure.label(reference_seg), coordinates='rc')[0]

                pharynx_center_y, pharynx_center_x = props.centroid
                pharynx_orientation = props.orientation

                translation_matrix = transform.EuclideanTransform(
                    translation=(-(img_center_x - pharynx_center_x), -(img_center_y - pharynx_center_y))
                )

                rotated_img = rotate(img.data, translation_matrix, pharynx_orientation)
                rotated_seg = rotate(seg.data, translation_matrix, pharynx_orientation)

                fl_rotated_stack.loc[dict(wavelength=wvl, pair=pair)][img_idx] = rotated_img
                seg_rotated_stack.loc[dict(wavelength=wvl, pair=pair)][img_idx] = rotated_seg

    return fl_rotated_stack, seg_rotated_stack


def segment_pharynxes(fl_stack: Union[np.ndarray, xr.DataArray]):
    """

    Parameters
    ----------
    fl_stack

    Returns
    -------

    """
    return fl_stack > 2000


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
        np.degrees(np.pi / 2 - orientation), mode='wrap')


def calculate_midlines(rot_seg_stack, s=1e8, ext=0) -> List[dict]:
    """Calculate the midlines for the given stack. Only calculates midlines for NON-TL images
    Parameters
    ----------
    rot_seg_stack: xr.DataArray
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
                [calculate_midline(rot_seg_stack.isel(strain=img_idx).sel(wavelength=wvl, pair=pair), s, ext) for pair
                 in rot_seg_stack.pair.data]
            for wvl in rot_seg_stack.wavelength.data if 'tl' not in wvl.lower()
        }
        for img_idx in range(rot_seg_stack.strain.size)
    ]


def calculate_midline(rot_seg_img, s=1e8, ext=0):
    """Calculate a the midline for a single image. Right now this only works for images that have been centered and aligned
    with their anterior-posterior along the horizontal.

    Parameters
    ----------
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
    return UnivariateSpline(xs, ys, s=s, ext=ext)


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

    ys = xr.DataArray(mid(xs), dims='z')
    xs = xr.DataArray(xs, dims='z')
    return fl.interp(x=xs, y=ys).data.T


def measure_under_midlines(fl_stack: xr.DataArray, midlines: Iterable, x_range: tuple, n_points: int) -> xr.DataArray:
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

    raw_intensity_data = xr.DataArray(
        np.zeros(
            (fl_stack.strain.size, fl_stack.wavelength.size, fl_stack.pair.size, n_points)
        ),
        dims=['strain', 'wavelength', 'pair', 'position'],
        coords={'strain': fl_stack.strain, 'wavelength': fl_stack.wavelength, 'pair': fl_stack.pair}
    )

    xs = np.linspace(x_range[0], x_range[1], n_points)

    # TODO: implement non-frame-specific midlines
    for img_idx in range(fl_stack.strain.size):
        for wvl_idx in range(fl_stack.wavelength.size):
            wvl = fl_stack.wavelength.data[wvl_idx]
            if 'tl' not in wvl.lower():
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

    Returns
    -------

    """
    # peaks_all = xr.DataArray(
    #     np.full((intensity_data.strain.size, intensity_data.wavelength.size, intensity_data.pair.size, 2), np.nan),
    #     dims=['strain', 'wavelength', 'pair', 'peaks'],
    #     coords={'strain': intensity_data.strain, 'wavelength': intensity_data.wavelength, 'pair': intensity_data.pair}
    # )
    # for img_idx in range(intensity_data.strain.size):
    #     for wvl_idx in range(intensity_data.wavelength.size):
    #         if intensity_data.wavelength[wvl_idx] != 'TL':
    #             for pair_idx in range(intensity_data.pair.size):
    #                 ys = intensity_data.isel(strain=img_idx, wavelength=wvl_idx, pair=pair_idx).data
    #                 p, peak_props = find_peaks(ys, distance=.3 * len(ys), prominence=100, wlen=10)
    #                 # Get 2 largest peaks
    #                 p = np.argpartition(p, len(p) - 2)[-2:]
    #                 peaks_all[dict(strain=img_idx, wavelength=wvl_idx, pair=pair_idx)] = p
    #
    #                 should_flip = False
    #                 if len(p) == 1:
    #                     # Usually if there is only 1 peak, it's the anterior peak
    #                     if p[0] < len(ys) // 2:
    #                         should_flip = True
    #                 if len(p) == 2:
    #                     proms = peak_props['prominences']
    #                     if proms[0] > proms[1]:
    #                         should_flip = True
    #                     if p[0] < len(ys) - p[1]:
    #                         should_flip = True
    #
    #                 if should_flip:
    #                     intensity_data[dict(strain=img_idx, wavelength=wvl_idx, pair=pair_idx)] = np.flip(
    #                         intensity_data[dict(strain=img_idx, wavelength=wvl_idx, pair=pair_idx)])

    data = trim_profiles(intensity_data, threshold=2000, new_length=100)
    ref_data = data.sel(wavelength=reference_wavelength, pair=reference_pair)
    ref_profile = ref_data.isel(strain=0)

    unflipped_mse = np.sum(np.power(ref_data - ref_profile, 2), axis=1).data
    flipped_mse = np.sum(np.power(np.flip(ref_data) - ref_profile, 2), axis=1).data

    intensity_data[unflipped_mse > flipped_mse] = np.flip(intensity_data[unflipped_mse > flipped_mse], axis=3)

    mean_intensity = trim_profile(np.mean(intensity_data.sel(wavelength=reference_wavelength, pair=reference_pair), axis=0).data, threshold=2000, new_length=100)

    peaks, _ = find_peaks(mean_intensity, distance=.2 * len(mean_intensity), prominence=200, wlen=10)

    if peaks[0] < len(mean_intensity) - peaks[1]:
        intensity_data = np.flip(intensity_data, axis=3)

    return intensity_data


def trim_profile(profile, threshold, new_length):
    """
    TODO: Documentation
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

    # TODO: also return first and last idx... ummm why?
    return np.interp(new_xs, old_xs, trimmed)


def trim_profiles(intensity_data, threshold, new_length):
    """
    TODO: Documentation
    Parameters
    ----------
    intensity_data
    threshold
    new_length

    Returns
    -------

    """
    # TODO: use trim boundaries from 410 to trim 470
    trimmed_intensity_data = xr.DataArray(
        np.zeros(
            (intensity_data.strain.size, intensity_data.wavelength.size, intensity_data.pair.size, new_length)
        ),
        dims=['strain', 'wavelength', 'pair', 'position'],
        coords={'strain': intensity_data.strain, 'wavelength': intensity_data.wavelength, 'pair': intensity_data.pair}
    )
    for img_idx in range(intensity_data.strain.size):
        for wvl_idx in range(intensity_data.wavelength.size):
            wvl = intensity_data.wavelength.data[wvl_idx]
            if 'tl' not in wvl.lower():
                for pair in range(intensity_data.pair.size):
                    untrimmed_profile = intensity_data.sel(wavelength=wvl, pair=pair).isel(strain=img_idx).data
                    trimmed_intensity_data[img_idx, wvl_idx, pair, :] = \
                        trim_profile(untrimmed_profile, threshold, new_length)

    return trimmed_intensity_data


def r_to_oxd(r, r_min=0.852, r_max=6.65, instrument_factor=0.171):
    return (r - r_min) / ((r - r_min) + instrument_factor * (r_max - r))


def oxd_to_redox_potential(oxd, midpoint_potential=-265, z=2, temperature=22):
    # TODO: returns NaN sometimes?
    return midpoint_potential - (8314.462 * (273.15 + temperature) / (z * 96485.3415)) * np.log((1 - oxd) / oxd)


def center_of_mass_midline(rot_fl):
    """ Calculate the midline using the Center of Mass method

    Parameters
    ----------
    rot_fl

    Returns
    -------

    """
    ys = np.arange(rot_fl.shape[0])
    midline_ys = []
    for i in range(rot_fl.shape[1]):
        midline_ys.append(np.average(ys, weights=rot_fl[:, i].data))
    return np.array(midline_ys)
