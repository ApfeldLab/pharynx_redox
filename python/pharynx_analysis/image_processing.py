from typing import Union

import numpy as np
import xarray as xr
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
from skimage import measure, transform


def center_and_rotate_pharynxes(fl_stack, seg_stack, crop_width=None, crop_height=None) -> (np.ndarray, np.ndarray):
    """ TODO: Documentation
    :param fl_stack:
    :param seg_stack:
    :param crop_width:
    :param crop_height:
    :return:
    """
    all_props = [measure.regionprops(x, coordinates='rc')[0] for x in measure.label(seg_stack)]

    rotated_fl = fl_stack.copy()
    rotated_seg = seg_stack.copy()
    rotated_seg['wavelength'] = rotated_fl['wavelength']

    img_center = ((fl_stack.shape[2] / 2) - .5, (fl_stack.shape[1] / 2) - .5)
    for i, p in enumerate(all_props):
        translation_matrix = transform.EuclideanTransform(
            translation=(-(img_center[0] - p.centroid[1]), -(img_center[1] - p.centroid[0])))
        rotated_fl[i, :, :] = rotate(fl_stack[i].data, translation_matrix, p.orientation)
        rotated_seg[i, :, :] = rotate(seg_stack[i].data, translation_matrix, p.orientation)

    # go through the rotated images to check that the orientation is correct
    mean_intensity_ap = np.mean(rotated_fl, axis=1)

    peak_dists = np.abs(np.vstack(
        [find_peaks(ys, distance=20, height=500)[0] for ys in mean_intensity_ap]
    ) - rotated_fl.shape[2] / 2)

    should_flip_arr = peak_dists[:, 0] > peak_dists[:, 1]

    # noinspection PyTypeChecker
    for i, should_flip in enumerate(should_flip_arr):
        if should_flip:
            rotated_fl[i] = np.fliplr(rotated_fl[i])
            rotated_seg[i] = np.fliplr(rotated_seg[i])

    center_x = rotated_fl.shape[2] // 2
    center_y = rotated_fl.shape[1] // 2
    if crop_width:
        rotated_fl = rotated_fl[:, :, center_x - crop_width // 2:center_x + crop_width // 2]
        rotated_seg = rotated_seg[:, :, center_x - crop_width // 2:center_x + crop_width // 2]
    if crop_height:
        rotated_fl = rotated_fl[:, center_y - crop_height // 2:center_y + crop_height // 2, :]
        rotated_seg = rotated_seg[:, :, center_y - crop_height // 2:center_y + crop_height // 2]

    return rotated_fl, rotated_seg


def segment_pharynxes(fl_stack: Union[np.ndarray, xr.DataArray]):
    """
    TODO: Documentation
    :param fl_stack:
    :return:
    """
    return fl_stack > 2000


def rotate(data, tform, orientation):
    """
    TODO: Documentation
    :param data:
    :param tform:
    :param orientation:
    :return:
    """
    # noinspection PyTypeChecker
    return transform.rotate(transform.warp(
        data, tform, preserve_range=True, mode='wrap'),
        np.degrees(np.pi / 2 - orientation), mode='wrap')


def calculate_midlines(rot_seg_stack, s=1e8, ext=0):
    """
    TODO: Documentation
    :param rot_seg_stack:
    :param s:
    :param ext:
    :return:
    """
    return [calculate_midline(x, s, ext) for x in rot_seg_stack]


def calculate_midline(rot_seg_img, s=1e8, ext=0):
    # TODO center segmentation before fitting splines
    """
    TODO: Documentation
    :param rot_seg_img:
    :param s:
    :param ext:
    :return:
    """
    seg_coords = measure.regionprops(measure.label(rot_seg_img))[0].coords

    # Adding tiny amounts of random noise to the coordinates because the spline fitting function expects that all
    # x-coordinates are unique.
    # TODO: change this; I think this could fail if two random numbers happen to be the same
    seg_coords = seg_coords + (np.random.random(seg_coords.shape) / 100)
    seg_coords = seg_coords[np.argsort(seg_coords, axis=0)[:, 1]]

    xs = seg_coords[:, 1]
    ys = seg_coords[:, 0]
    return UnivariateSpline(xs, ys, s=s, ext=ext)


def measure_under_midline(fl: np.ndarray, mid: UnivariateSpline, xs: np.ndarray) -> np.ndarray:
    """
    TODO: Documentation
    :param fl:
    :param mid:
    :param xs:
    :return:
    """
    ys = mid(xs)
    zs = np.zeros((1, len(xs)), dtype=np.uint16)
    for i, (x, y) in enumerate(zip(np.int_(xs), np.int_(ys))):
        zs[i] = fl[y, x].data
    return zs


def measure_under_midlines(fl_stack: np.ndarray, midlines: [UnivariateSpline],
                           x_range: tuple, n_points: int) -> np.ndarray:
    """
    TODO: Documentation
    :param fl_stack:
    :param midlines:
    :param x_range:
    :param n_points:
    :return:
    """
    xs = np.linspace(x_range[0], x_range[1], n_points)
    return np.asarray([
        measure_under_midline(fl_stack[i], midlines[i], xs) for i in range(fl_stack.shape[0])
    ])


def trim_profile(profile, threshold, new_length):
    """
    TODO: Documentation
    :param profile:
    :param threshold:
    :param new_length:
    :return:
    """
    first = np.argmax(profile > threshold)
    last = len(profile) - np.argmax(np.flip(profile > threshold))

    trimmed = profile[first:last]
    new_xs = np.linspace(0, len(trimmed), new_length)
    old_xs = np.arange(0, len(trimmed))

    return np.interp(new_xs, old_xs, trimmed)


def trim_profiles(intensity_data_stack, threshold, new_length):
    """
    TODO: Documentation
    :param intensity_data_stack:
    :param threshold:
    :param new_length:
    :return:
    """
    trimmed_data = np.ndarray(
        shape=(len(intensity_data_stack.wavelength), intensity_data_stack.strain.shape[0], new_length))
    for i, wvl in enumerate(intensity_data_stack.wavelength):
        for j in range(intensity_data_stack.strain.shape[0]):
            untrimmed_profile = intensity_data_stack.sel(wavelength=wvl).isel(strain=j).data
            trimmed_data[i, j, :] = trim_profile(untrimmed_profile, threshold, new_length)
    return xr.DataArray(trimmed_data, dims=['wavelength', 'strain', 'position'],
                        coords={'wavelength': intensity_data_stack.wavelength, 'strain': intensity_data_stack.strain})
