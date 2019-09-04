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
from skimage.transform import warp, AffineTransform

from pharynx_analysis import profile_processing


def get_lr_bounds(
    rot_seg_stack: xr.DataArray, pad: int = 0, ref_wvl: str = "410", ref_pair: int = 0
) -> np.ndarray:
    """
    Get the left and right boundaries of the rotated pharynxes
    Parameters
    ----------
    rot_seg_stack
        the rotated segmented pharynxes
    pad
        the amount of padding on the left/right of the  bounds
    ref_wvl
        the wavelength to use for calculating bounds
    ref_pair
        the pair to use for calculating bounds

    Returns
    -------
    bounds
        An (m, 2) array where m = number of animals, the first column is the left bound
        and the second column is the right bound
    """
    imgs = rot_seg_stack.sel(wavelength=ref_wvl, pair=ref_pair)
    bounds = np.zeros((imgs.strain.size, 2))  # (0, (x, y))
    for i, img in enumerate(imgs):
        _, l, _, r = measure.regionprops(measure.label(img))[0].bbox
        bounds[i, :] = [l - pad, r + pad - 1]
    return bounds.astype(np.int)


def center_and_rotate_pharynxes(
    fl_images, seg_images, reference_wavelength="410"
) -> (np.ndarray, np.ndarray):
    """
    Given a fluorescence stack and a pharyngeal mask stack, center and rotate each frame of both the FL and mask such
    that the pharynx is in the center of the image, with its anterior on the left.

    Parameters
    ----------
    fl_images
        The fluorescence images to rotate and align
    seg_images
        The segmented images to rotate and align
    reference_wavelength
        The wavelength to use for calculating center of mass and angle of orientation

    Returns
    -------
    (rotated_fl_stack, rotated_seg_stack)
        A 2-tuple where the first item is the rotated fluorescence stack and the second is the rotated mask stack

    Notes
    -----
    This function uses the a reference wavelength to calculate the center of mass and angle of orientation, then applies
    the according translation to all wavelengths for that animal/pair.

    The current implementation uses a gaussian blur on the fluorescence images and segments that to calculate the
    centroid and the orientation angle.
    """
    img_center_y, img_center_x = (
        fl_images.y.size // 2,
        fl_images.x.size // 2,
    )  # (y, x)

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
                reference_seg = blurred_seg.isel(strain=img_idx).sel(
                    wavelength=reference_wavelength, pair=pair
                )
                img = fl_images.isel(strain=img_idx).sel(wavelength=wvl, pair=pair)
                seg = seg_rotated_stack.isel(strain=img_idx).sel(
                    wavelength=wvl, pair=pair
                )

                props = measure.regionprops(
                    measure.label(reference_seg), coordinates="rc"
                )[0]

                # pharynx_center_y, pharynx_center_x = props.centroid
                pharynx_center_y, pharynx_center_x = np.mean(
                    np.nonzero(reference_seg), axis=1
                )
                pharynx_orientation = props.orientation

                translation_matrix = transform.EuclideanTransform(
                    translation=(
                        -(img_center_x - pharynx_center_x),
                        -(img_center_y - pharynx_center_y),
                    )
                )

                rotated_img = rotate(img.data, translation_matrix, pharynx_orientation)
                rotated_seg = rotate(seg.data, translation_matrix, pharynx_orientation)

                fl_rotated_stack.loc[dict(wavelength=wvl, pair=pair)][
                    img_idx
                ] = rotated_img
                seg_rotated_stack.loc[dict(wavelength=wvl, pair=pair)][
                    img_idx
                ] = rotated_seg

    return fl_rotated_stack, seg_rotated_stack


def extract_largest_binary_object(
    bin_img: Union[xr.DataArray, np.ndarray]
) -> Union[xr.DataArray, np.ndarray]:
    """
    Extracts the largest binary object from the given binary image

    Parameters
    ----------
    bin_img
        The binary image to process

    Returns
    -------
    bin_img
        The binary image containing only the largest binary object from the input

    """
    labels = measure.label(bin_img)
    if labels.max() == 0:
        # If there are no objects in the image... simply return the image
        return bin_img
    return labels == np.argmax(np.bincount(labels.flat)[1:]) + 1


# noinspection PyUnresolvedReferences
def segment_pharynxes(fl_stack: xr.DataArray, threshold: int = 2000) -> xr.DataArray:
    """
    Segment the pharynxes in the given fluorescence image stack

    Parameters
    ----------
    fl_stack
        the images to segment
    threshold
        pixels with brightness above this intensity are considered to be in the pharynx

    Returns
    -------
    seg
        the image stack containing the segmented masks of the pharynxes in the input fl_stack

    Notes
    -----
    This function currently uses a static threshold to segment, then extracts the largest binary object. More
    sophisticated segmentation strategies should be tried in the future.
    """
    seg = fl_stack > threshold
    for img_idx in range(fl_stack.strain.size):
        for wvl_idx in range(fl_stack.wavelength.size):
            for pair in range(fl_stack.pair.size):
                seg_img = seg[dict(strain=img_idx, wavelength=wvl_idx, pair=pair)]
                seg[
                    dict(strain=img_idx, wavelength=wvl_idx, pair=pair)
                ] = extract_largest_binary_object(seg_img)
    return seg


def get_centroids(fl_stack: xr.DataArray, threshold=1000, gaussian_sigma=6):
    """
    Obtain the centers-of-mass for each pharynx in the given fluorescence image stack

    Parameters
    ----------
    fl_stack
        the fluorescence image stack to measure
    threshold
        the segmentation threshold for the blurred pharynx images
    gaussian_sigma
        the degree to blur the pharynxes before segmentation

    Returns
    -------
    centroids
        the centers-of-mass of each pharynx in the given stack

    """
    image_data = fl_stack.copy()
    image_data.data = ndi.gaussian_filter(
        image_data.data, sigma=(0, 0, 0, gaussian_sigma, gaussian_sigma)
    )
    image_data.data[image_data.data < threshold] = 0
    image_data.data[image_data.data > threshold] = 1

    # TODO finish implementation


def rotate(img: Union[np.ndarray, xr.DataArray], tform, orientation):
    """
    Rotate the

    Parameters
    ----------
    img
        the image to rotate
    tform
        the translation matrix to apply
    orientation
        the angle of orientation (in radians)

    Returns
    -------
    rotated
        the translated and rotated image

    """
    # noinspection PyTypeChecker
    return transform.rotate(
        transform.warp(img, tform, preserve_range=True, mode="wrap"),
        np.degrees(np.pi / 2 - orientation),
        mode="edge",
    )


def calculate_midlines(
    rot_seg_stack: xr.DataArray, degree: int = 4
) -> List[Dict[str, List[Polynomial]]]:
    """
    Calculate a midline for each animal in the given stack

    Parameters
    ----------
    rot_seg_stack
        The rotated mask with which midlines should be calculated.

    degree
        The degree of the polynomial fit

    Returns
    -------
    list of dict
        A list of dictionaries with the following structure::

            [
                {
                    wavelength0: [midline_pair_0, midline_pair_1, ...],
                    wavelength1: [midline_pair_0, midline_pair_1, ...]
                },
                ...
            ]

    See Also
    --------
    calculate_midline
    """
    return [
        {
            wvl: [
                calculate_midline(
                    rot_seg_stack.isel(strain=img_idx).sel(wavelength=wvl, pair=pair),
                    degree=degree,
                )
                for pair in rot_seg_stack.pair.data
            ]
            for wvl in rot_seg_stack.wavelength.data
            if "tl" not in wvl.lower()
        }
        for img_idx in tqdm.trange(rot_seg_stack.strain.size)
    ]


def calculate_midline(
    rot_seg_img: Union[np.ndarray, xr.DataArray], degree: int = 4, pad: int = 5
) -> Polynomial:
    """
    Calculate a the midline for a single image by fitting a polynomial to the segmented
    pharynx

    Parameters
    ----------
    rot_seg_img: Union[np.ndarray, xr.DataArray]
        The rotated masked pharynx image
    degree
        the degree of the polynomial
    pad
        the number of pixels to "pad" the domain of the midline with respect to the
        boundaries of the segmentation mask

    Returns
    -------
    Polynomial
        the estimated midline

    Notes
    -----
    Right now this only works for images that have been centered and aligned with their
    anterior-posterior along the horizontal.
    """
    rp = measure.regionprops(measure.label(rot_seg_img))[0]
    xs, ys = rp.coords[:, 1], rp.coords[:, 0]
    left_bound, _, _, right_bound = rp.bbox
    # noinspection PyTypeChecker
    return Polynomial.fit(xs, ys, degree, domain=[left_bound - pad, right_bound + pad])


def measure_under_midline(
    fl: xr.DataArray, mid: Polynomial, n_points: int = 100, thickness: float = 0.0
) -> np.ndarray:
    """
    Measure the intensity profile of the given image under the given midline at the given x-coordinates.

    Parameters
    ----------
    fl
        The fluorescence image to measure
    mid
        The midline under which to measure
    n_points
        The number of points to measure under
    thickness
        The thickness of the line to measure under. WARNING: this is a lot slower right now

    Returns
    -------
    zs: np.ndarray
        The intensity profile of the image measured under the midline at the given x-coordinates.

    """
    if thickness == 0:
        xs, ys = mid.linspace(n=n_points)
        fl = np.asarray(fl)
        return ndi.map_coordinates(fl, np.stack([ys, xs]), order=1)
    else:
        xs, ys = mid.linspace()
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
            line = measure.profile._line_profile_coordinates(
                (xs0[i], ys0[i]), (xs1[i], ys1[i])
            )[:, :, 0]
            line_xs = xr.DataArray(line[0], dims="z")
            line_ys = xr.DataArray(line[1], dims="z")
            prof.append(np.mean(fl.interp(x=line_xs, y=line_ys).data, 0))
        # TODO: optionally use gaussian weight for the profile
        return np.array(prof)


def measure_under_midlines(
    fl_stack: xr.DataArray,
    midlines: List[Dict[str, List[Polynomial]]],
    n_points: int = 300,
    frame_specific: bool = False,
) -> xr.DataArray:
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
    non_tl_wvls = list(filter(lambda x: x != "TL", fl_stack.wavelength.data))
    raw_intensity_data = xr.DataArray(
        np.zeros(
            (fl_stack.strain.size, len(non_tl_wvls), fl_stack.pair.size, n_points)
        ),
        dims=["strain", "wavelength", "pair", "position"],
        coords={
            "strain": fl_stack.strain,
            "wavelength": non_tl_wvls,
            "pair": fl_stack.pair,
        },
    )

    ref_wvl = "410"

    for img_idx in tqdm.trange(fl_stack.strain.size):
        for pair in range(fl_stack.pair.size):
            for wvl_idx, wvl in enumerate(raw_intensity_data.wavelength.data):
                img = fl_stack.sel(wavelength=wvl, pair=pair).isel(strain=img_idx)

                if frame_specific:
                    mid = midlines[img_idx][wvl][pair]
                else:
                    mid = midlines[img_idx][ref_wvl][pair]

                raw_intensity_data[img_idx, wvl_idx, pair, :] = measure_under_midline(
                    img, mid, n_points
                )

    raw_intensity_data.values = np.nan_to_num(raw_intensity_data.values)
    return raw_intensity_data


def align_pa(
    intensity_data: xr.DataArray,
    reference_wavelength: str = "410",
    reference_pair: int = 0,
) -> xr.DataArray:
    """
    Given intensity profile data, flip each animal along their anterior-posterior axis
    if necessary, so that all face the same direction

    Parameters
    ----------
    intensity_data
        the data to align
    reference_pair: optional
        the pair to calculate the alignment for
    reference_wavelength: optional
        the wavelength to calculate the alignment for

    Returns
    -------
    aligned_intensity_data
        the PA-aligned intensity data

    Notes
    -----
    The alignments are calculated for a single wavelength and pair for each animal, then
    applied to all wavelengths and pairs for that animal.

    The algorithm works as follows:

        - take the derivative of the (trimmed) intensity profiles (this accounts for
          differences in absolute intensity between animals)
        - use the first animal in the stack as the reference profile
        - for all animals:

           - compare a forward and reverse profile to the reference profile (using the
             cosine-similarity metric)
           - keep either the forward or reverse profile accordingly

        - finally, determine the location of the peaks in the *average* profile

            - reverse all profiles if necessary (this will be necessary if the first
              animal happens to be reversed)

    """
    data = intensity_data

    ref_data = data.sel(wavelength=reference_wavelength, pair=reference_pair)
    ref_profile = ref_data.isel(strain=0).data

    ref_vecs = np.tile(ref_profile, (data.strain.size, 1))
    unflipped = data.sel(wavelength=reference_wavelength, pair=reference_pair).data
    flipped = np.fliplr(unflipped)

    should_flip = (
        cdist(ref_vecs, unflipped, "cosine")[0, :]
        > cdist(ref_vecs, flipped, "cosine")[0, :]
    )

    intensity_data[should_flip] = np.flip(intensity_data[should_flip], axis=3)

    mean_intensity = profile_processing.trim_profile(
        np.mean(
            intensity_data.sel(wavelength=reference_wavelength, pair=reference_pair),
            axis=0,
        ).data,
        threshold=2000,
        new_length=100,
    )

    peaks, _ = find_peaks(
        mean_intensity, distance=0.2 * len(mean_intensity), prominence=200, wlen=10
    )

    if len(peaks) < 2:
        return intensity_data

    if peaks[0] < len(mean_intensity) - peaks[1]:
        intensity_data = np.flip(intensity_data, axis=3)

    return intensity_data


def center_of_mass_midline(rot_fl: xr.DataArray, s: float, ext: str):
    """
    Calculate the midline using the Center of Mass method

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


def shift(image: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """
    Translate the image according to the given movement vector

    Parameters
    ----------
    image
        the image to translate
    vector :
        translation parameters ``(dx, dy)``

    Returns
    -------
    img: np.ndarray
        the translated image

    """
    tform = AffineTransform(translation=vector)
    shifted = warp(image, tform, mode="wrap", preserve_range=True)
    shifted = shifted.astype(image.dtype)
    return shifted
