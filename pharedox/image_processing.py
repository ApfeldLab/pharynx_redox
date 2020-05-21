import logging
from typing import Union

import SimpleITK as sITK
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from numpy.polynomial.polynomial import Polynomial
from scipy import ndimage as ndi
from scipy.stats import norm, zscore
from skimage import io
from skimage import measure, transform
from skimage.measure import label, regionprops
from skimage.transform import AffineTransform, warp


def measure_under_labels(
    imgs: xr.DataArray,
    masks: xr.DataArray,
    ref_wvl: str = "410",
    ratio_numerator="410",
    ratio_denominator="470",
):
    """Measure the intensities of each channel under the label image"""
    df = []
    imgs = imgs.where(imgs.wavelength != "TL", drop=True)

    for a in imgs.animal.values:
        for tp in imgs.timepoint.values:
            for p in imgs.pair.values:
                for wvl in imgs.wavelength.values:
                    img_selector = dict(animal=a, timepoint=tp, pair=p, wavelength=wvl)

                    if "wavelength" in masks.dims:
                        seg_frame = masks.sel(
                            animal=a, timepoint=tp, pair=p, wavelength=ref_wvl
                        )
                    else:
                        # single wavelength was passed
                        seg_frame = masks.sel(animal=a, timepoint=tp, pair=p)

                    labels = measure.label(seg_frame)

                    sub_df = pd.DataFrame(
                        measure.regionprops_table(
                            labels,
                            intensity_image=imgs.sel(**img_selector).values,
                            properties=["label", "mean_intensity", "area"],
                        )
                    )

                    sub_df["animal"] = a
                    sub_df["timepoint"] = tp
                    sub_df["pair"] = p
                    sub_df["wavelength"] = wvl
                    sub_df["strain"] = imgs.sel(**img_selector).strain.values

                    df.append(sub_df)

    df = pd.concat(df)
    df = df.set_index(["animal", "timepoint", "pair", "wavelength", "label"]).unstack(
        "wavelength"
    )
    df[("mean_intensity", "r")] = (
        df["mean_intensity"][ratio_numerator] / df["mean_intensity"][ratio_denominator]
    )
    df[("area", "r")] = df[("area", ratio_numerator)]
    df[("strain", "r")] = df[("strain", ratio_numerator)]

    df = df.stack("wavelength")

    return df


def subtract_medians(imgs: xr.DataArray) -> xr.DataArray:
    """
    Subtract the median from each image, keeping the datatype the same.

    Parameters
    ----------
    imgs
        the images to subtract the median from. May be a high-dimensional array.
    """

    submed = imgs.copy()
    submed.values = np.maximum(imgs - imgs.median(dim=["x", "y"]), 0)
    submed.loc[dict(wavelength="TL")] = imgs.sel(wavelength="TL")
    submed = submed.astype(imgs.dtype)
    return submed


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
    bounds = np.zeros((imgs.animal.size, 2))  # (animal, (l, r))
    for i, img in enumerate(imgs):
        _, l, _, r = measure.regionprops(measure.label(img))[0].bbox
        bounds[i, :] = [l - pad, r + pad - 1]
    return bounds.astype(np.int)


def center_and_rotate_pharynxes(
    fl_images: xr.DataArray, seg_images: xr.DataArray
) -> (xr.DataArray, xr.DataArray):
    """
    Given a fluorescence stack and a pharyngeal mask stack, center and rotate each frame
    of both the FL and mask such that the pharynx is in the center of the image, with 
    its anterior on the left.

    Parameters
    ----------
    fl_images
        The fluorescence images to rotate and align
    seg_images
        The segmented images to rotate and align

    Returns
    -------
    (rotated_fl_stack, rotated_seg_stack)
        A 2-tuple where the first item is the rotated fluorescence stack and the second 
        is the rotated mask stack
    """
    img_center_y, img_center_x = (
        fl_images.y.size // 2,
        fl_images.x.size // 2,
    )

    fl_rotated_stack = fl_images.copy()
    seg_rotated_stack = seg_images.copy()

    # STACK_ITERATION
    for img_idx in range(fl_images.animal.size):
        for wvl in fl_images.wavelength.data:
            for pair in fl_images.pair.data:
                for tp in fl_images.timepoint.values:
                    # Optimization potential here...
                    # this recalculates all region properties for the reference each time
                    img = fl_images.isel(animal=img_idx).sel(
                        wavelength=wvl, pair=pair, timepoint=tp
                    )
                    ref_seg = seg_images.isel(animal=img_idx).sel(
                        pair=pair, timepoint=tp
                    )

                    try:
                        props = measure.regionprops(measure.label(ref_seg))[0]
                    except IndexError:
                        raise ValueError(
                            f"No binary objects found in image @ [idx={img_idx} ; wvl={wvl} ; pair={pair}]"
                        )

                    # pharynx_center_y, pharynx_center_x = props.centroid
                    pharynx_center_y, pharynx_center_x = np.mean(
                        np.nonzero(ref_seg), axis=1
                    )
                    pharynx_orientation = props.orientation

                    translation_matrix = transform.EuclideanTransform(
                        translation=(
                            -(img_center_x - pharynx_center_x),
                            -(img_center_y - pharynx_center_y),
                        )
                    )

                    rotated_img = rotate(
                        img.data, translation_matrix, pharynx_orientation
                    )
                    rotated_seg = rotate(
                        ref_seg.data, translation_matrix, pharynx_orientation, order=1
                    )

                    fl_rotated_stack.loc[dict(wavelength=wvl, pair=pair, timepoint=tp)][
                        img_idx
                    ] = rotated_img

                    seg_rotated_stack.loc[dict(pair=pair, timepoint=tp)][
                        img_idx
                    ] = rotated_seg

    fl_rotated_stack.values = fl_rotated_stack.values.astype(fl_images.dtype)
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


def get_area_of_largest_object(mask: np.ndarray) -> int:
    """Returns the area (px) of the largest object in a binary image

    Parameters
    ----------
    mask : np.ndarray
        the binary image

    Returns
    -------
    int
        the area of the largest object 
    """
    try:
        return measure.regionprops(measure.label(mask))[0].area
    except IndexError:
        return 0


def segment_pharynx(
    fl_img: xr.DataArray, target_area: int = 450, area_range: int = 100
) -> xr.DataArray:
    """Generate a mask for the given image containing a pharynx.

    Parameters
    ----------
    fl_img : xr.DataArray
        a fluorescent image containing a single pharynx 
    target_area : int, optional
        the presumptive area (in px) of a pharynx, by default 450
    area_range : int, optional
        the acceptable range (in px) above/below the target_area, by default 100

    Returns
    -------
    xr.DataArray
        an image containing the segmented pharynx (dtype: np.uint8). Pixels of value=1
        indicate the pharynx, pixels of value=0 indicate the background.
    """
    # target_area = 450  # experimentally derived
    # area_range = 100
    min_area = target_area - area_range
    max_area = target_area + area_range

    max_iter = 300

    p = 0.15
    t = fl_img.max() * p
    mask = fl_img > t

    area = get_area_of_largest_object(mask)

    i = 0
    while (min_area > area) or (area > max_area):
        if i >= max_iter:
            return mask
        area = get_area_of_largest_object(mask)

        logging.debug(f"Setting p={p}")
        if area > max_area:
            p = p + 0.01
        if area < min_area:
            p = p - 0.01
        i = i + 1

        t = fl_img.max() * p
        mask = fl_img > t

        if p < 0:
            # break out if loop gets stuck w/ sensible default
            logging.warning("Caught infinite loop")
            return fl_img > (fl_img.max() * 0.15)
        if p > 0.9:
            logging.warning("Caught infinite loop")
            return fl_img > (fl_img.max() * 0.15)

    mask = extract_largest_binary_object(mask).astype(np.uint8)

    return mask


def segment_pharynxes(
    fl_stack: xr.DataArray,
    wvl: str = "410",
    target_area: int = 450,
    area_range: int = 100,
) -> xr.DataArray:
    """Segment a hyperstack of pharynxes

    Parameters
    ----------
    fl_stack : xr.DataArray
        the fluorescent images to segment
    wvl : str, optional
        the wavelength to segment, by default "410"
    target_area : int, optional
        the presumptive area of a pharynx, in pixels, by default 450
    area_range : int, optional
        the acceptable range of pharyngeal areas, by default 100

    Returns
    -------
    xr.DataArray
        the masks for the specified wavelength
    """

    to_segment = fl_stack.sel(wavelength=wvl)
    seg = xr.apply_ufunc(
        segment_pharynx,
        to_segment,
        input_core_dims=[["y", "x"]],
        output_core_dims=[["y", "x"]],
        vectorize=True,
        kwargs={"target_area": target_area, "area_range": area_range},
    )
    return seg


def rotate(img: Union[np.ndarray, xr.DataArray], tform, orientation, order=1):
    """
    Rotate the given image with the given translation matrix and orientation angle

    Parameters
    ----------
    img
        the image to rotate
    tform
        the translation matrix to apply
    orientation
        the angle of orientation (radians)
    order
        the order of the interpolation

    Returns
    -------
    rotated
        the translated and rotated image

    """

    return transform.rotate(
        transform.warp(img, tform, preserve_range=True, mode="wrap", order=order),
        np.degrees(np.pi / 2 - orientation),
        mode="edge",
        order=order,
    )


def calculate_midlines(rot_seg_stack: xr.DataArray, degree: int = 4) -> xr.DataArray:
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
    midlines
        a DataArray containing the midline objects

    See Also
    --------
    calculate_midline
    """
    return xr.apply_ufunc(
        calculate_midline,
        rot_seg_stack,
        input_core_dims=[["y", "x"]],
        vectorize=True,
        keep_attrs=True,
        kwargs={"degree": degree},
    )


def calculate_midline(
    rot_seg_img: Union[np.ndarray, xr.DataArray], degree: int = 4, pad: int = 10
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
    Right now this only works for images that this have been centered and aligned with their
    anterior-posterior along the horizontal.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            rp = measure.regionprops(measure.label(rot_seg_img))[0]
            xs, ys = rp.coords[:, 1], rp.coords[:, 0]
            left_bound, _, _, right_bound = rp.bbox

            return Polynomial.fit(
                xs, ys, degree, domain=[left_bound - pad, right_bound + pad]
            )
        except IndexError:
            # Indicates trying to measure on TL for example
            return None


def measure_under_midline(
    fl: xr.DataArray,
    mid: Polynomial,
    n_points: int = 100,
    thickness: float = 0.0,
    order=1,
    norm_scale=1,
    flatten=True,
) -> np.ndarray:
    """
    Measure the intensity profile of the given image under the given midline at the given x-coordinates.

    Parameters
    ----------
    flatten
    norm_scale
    order
        the interpolation order
    fl
        The fluorescence image to measure
    mid
        The midline under which to measure
    n_points
        The number of points to measure under
    thickness
        The thickness of the line to measure under. 
    
    Notes
    -----
    Using thickness is slower, depending on the amount of thickness 

    On my machine (2GHz Intel Core i5), as of 12/4/19:
        0-thickness:
            492 µs ± 16.6 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
        2-thickness:
            1.99 ms ± 65.8 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)
        10-thickness:
            3.89 ms ± 92.1 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)

    Returns
    -------
    zs: np.ndarray
        The intensity profile of the image measured under the midline at the given 
        x-coordinates.

    """
    # Make sure the image orientation matches with the expected order of map_coordinates
    try:
        if thickness == 0:
            xs, ys = mid.linspace(n=n_points)
            fl = np.asarray(fl)
            return ndi.map_coordinates(fl, np.stack([xs, ys]), order=order)
        else:
            # Gets a bit wonky, but makes sense

            # We need to get the normal lines from each point in the midline
            # then measure under those lines.

            # First, get the coordinates of the midline
            xs, ys = mid.linspace(n=n_points)

            # Now, we get the angles of each normal vector
            der = mid.deriv()
            normal_slopes = -1 / der(xs)
            normal_thetas = np.arctan(normal_slopes)

            # We get the x and y components of the start/end of the normal vectors
            mag = thickness / 2
            x0 = np.cos(normal_thetas) * mag
            y0 = np.sin(normal_thetas) * mag

            x1 = np.cos(normal_thetas) * -mag
            y1 = np.sin(normal_thetas) * -mag

            # These are the actual coordinates of the starts/ends of the normal vectors as they move
            # from (x,y) coordinates in the midline
            xs0 = xs + x0
            xs1 = xs + x1
            ys0 = ys + y0
            ys1 = ys + y1

            # We need to measure in a consistent direction along the normal line
            # if y0 < y1, we're going to be measuring in an opposite direction along the line... so we need flip the coordinates
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

            straightened = ndi.map_coordinates(fl, [all_xs, all_ys], order=order)

            if flatten:
                # Create a normal distribution centered around 0 with the given scale (see scipy.norm.pdf)
                # the distribution is then tiled to be the same shape as the straightened pharynx
                # then, this resultant matrix is the weights for averaging
                w = np.tile(
                    norm.pdf(np.linspace(-1, 1, n_line_pts), scale=norm_scale),
                    (n_points, 1),
                ).T
                profile = np.average(straightened, axis=0, weights=w)

                return profile
            else:
                return straightened
    except AttributeError:
        # This happens if the image is TL. Then it will have `None` instead of
        # a midline object
        pass
    except Exception as e:
        # Here, something actually went wrong
        logging.warning(f"measuring under midline failed with error {e}")

    return np.zeros((1, n_points))


def measure_under_midlines(
    fl_stack: xr.DataArray,
    midlines: xr.DataArray,
    n_points: int = 300,
    order=1,
    thickness=0,
) -> xr.DataArray:
    """
    Measure under all midlines in stack

    Parameters
    ----------
    order
    fl_stack
        The fluorescence stack under which to measure
    midlines: dict
        A DataArray containing the midlines 
    n_points: int
        the number of points to sample under the midline
    thickness: float
        the thickness of the midline to measure under

    Returns
    -------
    profile_data: xr.DataArray
        the intensity profiles for each image in the stack
    """
    measurements = xr.apply_ufunc(
        measure_under_midline,
        fl_stack,
        midlines,
        input_core_dims=[["x", "y"], []],
        output_core_dims=[["position"]],
        vectorize=True,
        keep_attrs=True,
        kwargs={
            "n_points": n_points,
            "thickness": thickness,
            "order": order,
            "flatten": True,
        },
    )

    measurements = measurements.assign_coords(
        {"position": np.linspace(0, 1, measurements.position.size)},
    )
    try:
        measurements = measurements.assign_coords(time=fl_stack.time)
    except AttributeError:
        pass

    return measurements


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
            for animal in range(fl_imgs.animal.size):
                prof = profiles.sel(wavelength=wvl, pair=pair).isel(animal=animal)
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


def z_normalize_with_masks(imgs, masks):
    """
    Perform z-normalization [0] on the entire image (relative to the content within the
    masks).

    That is to say, we center the pixels (within the mask) such that their mean is 0,
    and ensure their standard deviation is ~1.

    This allows us to see spatial patterns within the masked region (even if pixels
    outside of the masked region fall very far above or below those inside) by setting
    the colormap center around 0.

    [0] - https://jmotif.github.io/sax-vsm_site/morea/algorithm/znorm.html
    """
    masked = ma.masked_array(imgs, np.logical_not(masks))
    mu = np.mean(masked, axis=(-2, -1), keepdims=True)
    sigma = np.std(masked, axis=(-2, -1), keepdims=True)

    return (imgs - mu) / sigma


def create_normed_rgb_ratio_stack(
    r_imgs, seg_imgs, vmin=-7, vmax=7, cmap="coolwarm", output_filename=None
):
    """
    Z-normalize the images (relative to the masks), then transform them into RGB with
    the given colormap
    """
    r_znormed = z_normalize_with_masks(r_imgs, seg_imgs)
    # noinspection PyUnresolvedReferences
    normalizer = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    # TODO generalize dtype? for now, 32-bit only
    rgb_img = cmap(normalizer(r_znormed))[:, :, :, :3].astype(np.float16)

    if output_filename is not None:
        io.imsave(output_filename, rgb_img)

    return rgb_img


def get_bbox(m, pad=5):
    try:
        y_min, x_min, y_max, x_max = np.array(regionprops(label(m))[0].bbox)

        y_min = max(int(y_min - (pad / 2)), 0)
        x_min = max(int(x_min - (pad / 2)), 0)
        y_max = min(int(y_max + (pad / 2)), m.shape[0])
        x_max = min(int(x_max + (pad / 2)), m.shape[1])

        return np.array([y_min, x_min, y_max, x_max]).astype(np.float)
    except IndexError:
        return [np.nan, np.nan, np.nan, np.nan]


def bspline_intra_modal_registration(
    fixed_image, moving_image, fixed_image_mask=None, point_width=5.0,
):

    registration_method = sITK.ImageRegistrationMethod()

    # Determine the number of BSpline control points using the physical spacing we want
    # for the control grid.
    grid_physical_spacing = [
        point_width,
        point_width,
        point_width,
    ]
    image_physical_size = [
        size * spacing
        for size, spacing in zip(fixed_image.GetSize(), fixed_image.GetSpacing())
    ]
    mesh_size = [
        int(image_size / grid_spacing + 0.5)
        for image_size, grid_spacing in zip(image_physical_size, grid_physical_spacing)
    ]

    initial_transform = sITK.BSplineTransformInitializer(
        image1=fixed_image, transformDomainMeshSize=mesh_size, order=2
    )
    registration_method.SetInitialTransform(initial_transform)

    registration_method.SetMetricAsMeanSquares()

    if fixed_image_mask:
        registration_method.SetMetricFixedMask(fixed_image_mask)

    registration_method.SetInterpolator(sITK.sitkLinear)
    registration_method.SetOptimizerAsLBFGSB(
        gradientConvergenceTolerance=1e-5, numberOfIterations=10
    )

    return registration_method.Execute(fixed_image, moving_image)


def register_image(fixed, moving, mask=None, point_width=5.0):
    z_fixed = zscore(fixed.values)
    z_moving = zscore(moving.values)

    if mask is not None:
        mask = sITK.GetImageFromArray(mask * 255)

    tx = bspline_intra_modal_registration(
        sITK.GetImageFromArray(z_fixed),
        sITK.GetImageFromArray(z_moving),
        fixed_image_mask=mask,
        point_width=point_width,
    )

    reg_moving = sITK.GetArrayFromImage(
        sITK.Resample(
            sITK.GetImageFromArray(moving),
            sITK.GetImageFromArray(fixed),
            tx,
            sITK.sitkLinear,
        )
    )

    return reg_moving


def crop(img, bbox):
    y_min, x_min, y_max, x_max = bbox.values.astype(np.int)

    return img[y_min:y_max, x_min:x_max]


def register_all_images(
    imgs,
    masks,
    bbox_pad=10,
    point_width=6.0,
    fixed_wvl="410",
    moving_wvl="470",
    mask_wvl="410",
):
    bboxes = xr.apply_ufunc(
        get_bbox,
        masks,
        input_core_dims=[["y", "x"]],
        output_core_dims=[["pos"]],
        vectorize=True,
        kwargs={"pad": bbox_pad},
    ).assign_coords({"pos": ["min_row", "max_row", "min_col", "max_col"]})

    reg_imgs = imgs.copy()

    for animal in imgs.animal:
        for pair in imgs.pair:
            for timepoint in imgs.timepoint:
                fixed = imgs.sel(
                    animal=animal, pair=pair, timepoint=timepoint, wavelength=fixed_wvl
                )
                moving = imgs.sel(
                    animal=animal, pair=pair, timepoint=timepoint, wavelength=moving_wvl
                )
                mask = masks.sel(
                    animal=animal, pair=pair, timepoint=timepoint, wavelength=mask_wvl
                )
                bbox = bboxes.sel(
                    animal=animal, pair=pair, timepoint=timepoint, wavelength=mask_wvl
                )

                # crop image
                crop_fixed = crop(fixed, bbox)
                crop_moving = crop(moving, bbox)
                crop_mask = crop(mask, bbox)

                # register image
                reg_moving = register_image(
                    crop_fixed, crop_moving, mask=crop_mask, point_width=point_width
                )

                # paste cropped images back into correct location (from bbox)
                y_min, x_min, y_max, x_max = bbox.values.astype(np.int)
                reg_imgs.loc[
                    dict(
                        animal=animal,
                        pair=pair,
                        timepoint=timepoint,
                        wavelength=moving_wvl,
                    )
                ][y_min:y_max, x_min:x_max] = reg_moving

    return reg_imgs
