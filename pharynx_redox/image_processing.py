from typing import List, Dict, Union

import numpy as np
import numpy.ma as ma
<<<<<<< HEAD
=======
from tqdm.auto import tqdm
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
import xarray as xr
from numpy.polynomial.polynomial import Polynomial
from scipy import ndimage as ndi
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
from scipy.spatial.distance import cdist
from skimage import measure, transform
from skimage.transform import warp, AffineTransform
from skimage.external import tifffile
import matplotlib as mpl
import matplotlib.pyplot as plt
<<<<<<< HEAD
from scipy.stats import norm
=======
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80

from pharynx_redox import profile_processing


def subtract_medians(
    imgs: Union[np.ndarray, xr.DataArray], img_dims=(-2, -1)
) -> Union[np.ndarray, xr.DataArray]:
    """
    Subtract the median from each image.

    Parameters
    ----------
    imgs
        the images to subtract the median from. May be a high-dimensional array.
    img_dims
        the dimensions that the images are stored in the `imgs` array. 
    """

    return imgs - np.median(imgs, axis=(-2, -1), keepdims=True).astype(imgs.dtype)


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
    bounds = np.zeros((imgs.spec.size, 2))  # (0, (x, y))
    for i, img in enumerate(imgs):
        _, l, _, r = measure.regionprops(measure.label(img))[0].bbox
        bounds[i, :] = [l - pad, r + pad - 1]
    return bounds.astype(np.int)


def center_and_rotate_pharynxes(
    fl_images: xr.DataArray, seg_images: xr.DataArray, reference_wavelength: str = "410"
) -> (xr.DataArray, xr.DataArray):
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

<<<<<<< HEAD
    for img_idx in range(fl_images.spec.size):
=======
    for img_idx in tqdm(
        range(fl_images.spec.size), leave=False, desc="aliging pharynxes"
    ):
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
        for wvl in fl_images.wavelength.data:
            for pair in fl_images.pair.data:
                # Optimization potential here...
                # this recalculates all region properties for the reference each time
                reference_seg = blurred_seg.isel(spec=img_idx).sel(
                    wavelength=reference_wavelength, pair=pair
                )
                img = fl_images.isel(spec=img_idx).sel(wavelength=wvl, pair=pair)
                seg = seg_rotated_stack.isel(spec=img_idx).sel(
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

    return fl_rotated_stack.astype(fl_images.dtype), seg_rotated_stack


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
    for img_idx in range(fl_stack.spec.size):
        for wvl_idx in range(fl_stack.wavelength.size):
            for pair in range(fl_stack.pair.size):
                seg_img = seg[dict(spec=img_idx, wavelength=wvl_idx, pair=pair)]
                seg[
                    dict(spec=img_idx, wavelength=wvl_idx, pair=pair)
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
        
        accessed like so::

            midlines[img_idx][wvl][pair]

    See Also
    --------
    calculate_midline
    """
    return [
        {
            wvl: [
                calculate_midline(
                    rot_seg_stack.isel(spec=img_idx).sel(wavelength=wvl, pair=pair),
                    degree=degree,
                )
                for pair in rot_seg_stack.pair.data
            ]
            for wvl in rot_seg_stack.wavelength.data
            if "tl" not in wvl.lower()
        }
<<<<<<< HEAD
        for img_idx in range(rot_seg_stack.spec.size)
=======
        for img_idx in tqdm(
            range(rot_seg_stack.spec.size), leave=False, desc="calculating midlines"
        )
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
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
<<<<<<< HEAD
    Right now this only works for images that this have been centered this and aligned this with their
=======
    Right now this only works for images that have been centered and aligned with their
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
    anterior-posterior along the horizontal.
    """
    rp = measure.regionprops(measure.label(rot_seg_img))[0]
    xs, ys = rp.coords[:, 1], rp.coords[:, 0]
    left_bound, _, _, right_bound = rp.bbox
    # noinspection PyTypeChecker
    return Polynomial.fit(xs, ys, degree, domain=[left_bound - pad, right_bound + pad])


def measure_under_midline(
<<<<<<< HEAD
    fl: xr.DataArray,
    mid: Polynomial,
    n_points: int = 100,
    thickness: float = 0.0,
    order=1,
    norm_scale=1,
    flatten=True
=======
    fl: xr.DataArray, mid: Polynomial, n_points: int = 100, thickness: float = 0.0, order=1
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
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
<<<<<<< HEAD
        The thickness of the line to measure under. 
    
    Notes
    -----
    Using thickness is 5-8 times slower, depending on the amount of thickness 

    On my machine (2GHz Intel Core i5), as of 12/4/19:
        0-thickness:
            492 µs ± 16.6 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
        2-thickness:
            1.99 ms ± 65.8 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
        10-thickness:
            3.89 ms ± 92.1 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
=======
        The thickness of the line to measure under. WARNING: this is a lot slower right now
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80

    Returns
    -------
    zs: np.ndarray
        The intensity profile of the image measured under the midline at the given x-coordinates.

    """
<<<<<<< HEAD

    if thickness == 0:
        xs, ys = mid.linspace(n=n_points)
        fl = np.asarray(fl)
        return ndi.map_coordinates(fl, np.stack([ys, xs]), order=1)
    else:
        # Gets a bit wonky, but makes sense
        # TODO: maybe refactor this stuff out into individual functions?
        
        # We need to get the normal lines from each point in the midline
        # then measure under those lines.

        # First, get the coordinates of the midline
        xs, ys = mid.linspace(n=n_points)
        
        # Now, we get the angles of each normal vector
=======
    if thickness == 0:
        xs, ys = mid.linspace(n=n_points)
        fl = np.asarray(fl)
        return ndi.map_coordinates(fl, np.stack([ys, xs]), order=order)
    else:
        xs, ys = mid.linspace()
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
        der = mid.deriv()
        normal_slopes = -1 / der(xs)
        normal_thetas = np.arctan(normal_slopes)

<<<<<<< HEAD
        # We get the x and y components of the start/end of the normal vectors
        mag = thickness / 2
=======
        mag = 1.5
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
        x0 = np.cos(normal_thetas) * mag
        y0 = np.sin(normal_thetas) * mag

        x1 = np.cos(normal_thetas) * -mag
        y1 = np.sin(normal_thetas) * -mag

<<<<<<< HEAD
        # These are the actual coordinates of the starts/ends of the normal vectors as they move
        # from (x,y) coordinates in the midline
=======
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
        xs0 = xs + x0
        xs1 = xs + x1
        ys0 = ys + y0
        ys1 = ys + y1

<<<<<<< HEAD
        # This is kinda weird... But we need to do it.
        # TODO: Could be made faster w/ vectorization? Too tired to figure that out when I wrote it though

        # We need to measure in a consistent direction along the normal line
        # y0 < y1, we're going to be measuring in an opposite direction along the line... so we need flip the coordinates
        for y0, y1, x0, x1, i in zip(ys0, ys1, xs0, xs1, range(len(xs0))):
            if y0 < y1:
                tx = xs0[i]
                xs0[i] = xs1[i]
                xs1[i] = tx

                ty = ys0[i]
                ys0[i] = ys1[i]
                ys1[i] = ty


        n_line_pts = thickness

        all_xs = np.linspace(xs0, xs1, n_line_pts)
        all_ys = np.linspace(ys0, ys1, n_line_pts)

        straightened = ndi.map_coordinates(fl, [all_ys, all_xs])
        
        if flatten:
            # Create a normal distribution centered around 0 with the given scale (see scipy.norm.pdf)
            # the distribution is then tiled to be the same shape as the straightened pharynx
            # then, this resultant matrix is the weights for averaging
            w = np.tile(norm.pdf(np.linspace(-1, 1, n_line_pts), scale=norm_scale), (n_points,1)).T
            profile = np.average(straightened, axis=0, weights=w)

            return profile
        else:
            return straightened
=======
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
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80


def measure_under_midlines(
    fl_stack: xr.DataArray,
    midlines: List[Dict[str, List[Polynomial]]],
    n_points: int = 300,
    frame_specific: bool = False,
    order=1,
<<<<<<< HEAD
    thickness=0,
    flatten=True
=======
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
) -> xr.DataArray:
    """
    Measure under all midlines in stack

    Parameters
    ----------
    fl_stack
        The fluorescence stack under which to measure
    midlines: dict
        A list of dictionaries of midlines with the following structure::

            [
                {'410': [UnivariateSpline, UnivariateSpline, ...]},
                {'470': [UnivariateSpline, UnivariateSpline, ...]},
                ...
            ]
    n_points: int
        the number of points to sample under the midline
    frame_specific: bool
        whether to use a different midline for each frame. if False, a single midline
        will be used within all wavelengths in a pair
<<<<<<< HEAD
    thickness: float
        the thickness of the midline to measure under
=======
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80

    Returns
    -------
    profile_data: xr.DataArray
        the intensity profiles for each image in the stack
    """
    non_tl_wvls = list(filter(lambda x: x != "TL", fl_stack.wavelength.data))
    raw_intensity_data = xr.DataArray(
        0,
        dims=["spec", "wavelength", "pair", "position"],
        coords={
            "spec": fl_stack.spec,
            "wavelength": non_tl_wvls,
            "pair": fl_stack.pair,
            "position": np.arange(n_points),
            "strain": ("spec", fl_stack.strain),
        },
    )

    ref_wvl = "410"

<<<<<<< HEAD
    for img_idx in range(fl_stack.spec.size):
=======
    for img_idx in tqdm(
        range(fl_stack.spec.size), leave=False, desc="measuring under midlines"
    ):
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
        for pair in range(fl_stack.pair.size):
            for wvl_idx, wvl in enumerate(raw_intensity_data.wavelength.data):
                img = fl_stack.sel(wavelength=wvl, pair=pair).isel(spec=img_idx)

                if frame_specific:
                    mid = midlines[img_idx][wvl][pair]
                else:
                    mid = midlines[img_idx][ref_wvl][pair]

                raw_intensity_data[img_idx, wvl_idx, pair, :] = measure_under_midline(
<<<<<<< HEAD
                    img, mid, n_points, order=order, thickness=thickness, flatten=True
=======
                    img, mid, n_points, order=order
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
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
    ref_profile = ref_data.isel(spec=0).data

    ref_vecs = np.tile(ref_profile, (data.spec.size, 1))
    unflipped = data.sel(wavelength=reference_wavelength, pair=reference_pair).data
    flipped = np.fliplr(unflipped)

    # do the actual cosine-similarity measurements
    should_flip = (
        cdist(ref_vecs, unflipped, "cosine")[0, :]
        > cdist(ref_vecs, flipped, "cosine")[0, :]
    )

    # position needs to be reindexed, otherwise xarray freaks out
    intensity_data[should_flip] = np.flip(intensity_data[should_flip], axis=3).reindex(
        position=np.arange(intensity_data.position.size)
    )

    mean_intensity = profile_processing.trim_profile(
        np.mean(
            intensity_data.sel(wavelength=reference_wavelength, pair=reference_pair),
            axis=0,
        ).data,
        threshold=2000,
        new_length=100,
    )

    # parameters found experimentally
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


def normalize_images_by_wvl_pair(
    fl_imgs: xr.DataArray, profiles: xr.DataArray, percent_to_clip: float = 2.0
):
    """
    Normalize images by subtracting mean profile then min-max rescaling to [0, 1]
    Parameters
    ----------
    fl_imgs
        the images to normalize
    profiles
        the intensity profiles corresponding to the images
    percent_to_clip
        how much to clip the profile when calculating mean/min/max, expressed as a percentage of the length of the profile

    Returns
    -------
    xr.DataArray
        the normalized images
    """
    idx_to_clip = int(profiles.shape[-1] * percent_to_clip / 100)
    profiles = profiles[:, idx_to_clip:-idx_to_clip]

    norm_fl = fl_imgs.copy().astype(np.float)
    for pair in fl_imgs.pair:
        for wvl in fl_imgs.wavelength.values:
            if wvl not in profiles.wavelength.values:
                continue
            for animal in range(fl_imgs.spec.size):
                prof = profiles.sel(wavelength=wvl, pair=pair).isel(spec=animal)
                img = fl_imgs.sel(wavelength=wvl, pair=pair)[animal].astype(np.float)

                # First, center according to mean
                img = img - prof.mean()

                # Then rescale to [0, 1]
                img = (img - prof.min()) / (prof.max() - prof.min())

                norm_fl.loc[dict(wavelength=wvl, pair=pair)][animal] = img

    return norm_fl


def normalize_images_single_wvl(
    fl_imgs: Union[np.ndarray, xr.DataArray],
    profiles: Union[np.ndarray, xr.DataArray],
    percent_to_clip: float = 2.0,
) -> Union[np.ndarray, xr.DataArray]:
    """
    Normalize single wavelength image stack by subtracting the mean of the corresponding
    intensity profile, then min-max rescaling to [0, 1]

    Parameters
    ----------
    fl_imgs
        an array-like structure of shape (frame, row, col)
    profiles
        an array-like structure of shape (frame, position_along_midline)
    percent_to_clip
        how much to clip the profile when calculating mean/min/max, expressed as a percentage of the length of the profile

    Returns
    -------
    Union[np.ndarray, xr.DataArray]
        normalized images
    """

    if fl_imgs.ndim != 3:
        raise ValueError("images must have shape (frame, row, col)")
    if profiles.ndim != 2:
        raise ValueError("profiles must have shape (frame, position_along_midline)")

    normed_imgs = fl_imgs.copy().astype(np.float32)

    idx_to_clip = int(profiles.shape[-1] * percent_to_clip / 100)
    profiles = profiles[:, idx_to_clip:-idx_to_clip]

    prof_means = np.mean(profiles, axis=1)

    profiles = profiles - prof_means
    normed_imgs = normed_imgs - prof_means

    prof_mins = np.min(profiles, axis=1)
    prof_maxs = np.max(profiles, axis=1)

    normed_imgs = (normed_imgs - prof_mins) / (prof_maxs - prof_mins)

    return normed_imgs

<<<<<<< HEAD

=======
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
def z_normalize_with_masks(imgs, masks):
    """
    Perform z-normalization [0] on the entire image (relative to the content within the masks).

<<<<<<< HEAD
    That is to say, we center the pixels (within the mask) such that their mean is 0, and ensure their standard deviation is ~1.
=======
    That is to say, we center the pixels within the mask such that their mean is 0, and ensure their standard deviation is ~1.
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80

    This allows us to see spatial patterns within the masked region (even if pixels outside of the masked region
    fall very far above or below those inside) by setting the colormap center around 0.

    [0] - https://jmotif.github.io/sax-vsm_site/morea/algorithm/znorm.html
    """
    masked = ma.masked_array(imgs, np.logical_not(masks))
    mu = np.mean(masked, axis=(-2, -1), keepdims=True)
    sigma = np.std(masked, axis=(-2, -1), keepdims=True)
<<<<<<< HEAD

    return (imgs - mu) / sigma


def create_normed_rgb_ratio_stack(
    r_imgs, seg_imgs, vmin=-7, vmax=7, cmap="coolwarm", output_filename=None
):
=======
    
    return (imgs - mu) / sigma 

def create_normed_rgb_ratio_stack(r_imgs, seg_imgs, vmin=-7, vmax=7, cmap='coolwarm', output_filename=None):
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
    """
    Z-normalize the images (relative to the masks), then transform them into RGB with the given colormap
    """
    r_znormed = z_normalize_with_masks(r_imgs, seg_imgs)
<<<<<<< HEAD

=======
    
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
    normalizer = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    # TODO generalize dtype? for now, 32-bit only
<<<<<<< HEAD
    rgb_img = cmap(normalizer(r_znormed))[:, :, :, :3].astype(np.float16)

    if output_filename is not None:
        with open(output_filename, "wb") as f:
            tifffile.imsave(f, rgb_img)

    return rgb_img
=======
    rgb_img = cmap(normalizer(r_znormed))[:,:,:,:3].astype(np.float16)
    
    if output_filename is not None:
        with open(output_filename, 'wb') as f:
            tifffile.imsave(f, rgb_img)
    
    return rgb_img
>>>>>>> 9e81c2636da482517a0c2c4d99584fa8716aed80
