from scipy.signal import find_peaks
from skimage import measure, transform
from scipy.interpolate import UnivariateSpline
import numpy as np


def center_and_rotate_pharynxes(fl_stack: np.ndarray, seg_stack: np.ndarray) -> (np.ndarray, np.ndarray):
    all_props = [measure.regionprops(x, coordinates='rc')[0] for x in measure.label(seg_stack)]

    rotated_fl = fl_stack.copy()
    rotated_seg = seg_stack.copy()

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

    return rotated_fl, rotated_seg


def segment_pharynxes(fl_stack: np.ndarray) -> np.ndarray:
    return fl_stack > 2000


def rotate(data, tform, orientation):
    # noinspection PyTypeChecker
    return transform.rotate(transform.warp(
        data, tform, preserve_range=True, mode='wrap'),
        np.degrees(np.pi / 2 - orientation), mode='wrap')


def calculate_midlines(rot_seg_stack, s=1e4, ext=0):
    return [calculate_midline(x, s, ext) for x in rot_seg_stack]


def calculate_midline(rot_seg_img, s=1e4, ext=0):
    seg_coords = measure.regionprops(measure.label(rot_seg_img))[0].coords

    # Adding tiny amounts of random noise to the coordinates because the spline fitting function expects that all
    # x-coordinates are unique. TODO: figure out if it's possible to not do this
    seg_coords = seg_coords + (np.random.random(seg_coords.shape) / 100)
    seg_coords = seg_coords[np.argsort(seg_coords, axis=0)[:, 1]]

    xs = seg_coords[:, 1]
    ys = seg_coords[:, 0]
    return UnivariateSpline(xs, ys, s=s, ext=ext)
