from collections import Counter

import numpy as np
import pandas as pd
import typing
from skimage.measure import regionprops, label


def jaccard(im1, im2):
    """
    Computes the Jaccard metric, a measure of set similarity.

    Parameters
    ----------
    im1 : array-like, bool
        Any array of arbitrary size. If not boolean, will be converted.
    im2 : array-like, bool
        Any other array of identical size. If not boolean, will be converted.

    Returns
    -------
    jaccard : float
        Jaccard metric returned is a float on range [0,1].
        Maximum similarity = 1
        No similarity = 0
    """
    im1 = np.asarray(im1).astype(np.bool)
    im2 = np.asarray(im2).astype(np.bool)

    intersection = np.logical_and(im1, im2)
    union = np.logical_or(im1, im2)

    return intersection.sum() / float(union.sum())


def create_occurrence_count_tuples(l: iter) -> [(object, int)]:
    # TODO: test
    """ Given a list of things, return a list of tuples (item, nth_occurrence)

    The n_th occurrence is 0-indexed

    :param l: an iterable
    :return: [(item, nth_occurrence), ...]
    """
    count_tuples = []
    c = Counter()
    for item in l:
        c.update([item])
        count_tuples.append((item, c[item] - 1))
    return count_tuples


def figure_to_np_array(fig):
    """Given a figure, return a 3D numpy array, where the dimensions are (height, width, RGB)"""
    fig.tight_layout(pad=0)
    fig.canvas.draw()
    data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return data


def calc_max_bbox(rot_seg_stack, ref_pair=0, ref_wvl="410"):
    b_boxes = []

    for i in range(rot_seg_stack.strain.size):
        props = regionprops(
            label(rot_seg_stack.sel(pair=ref_pair, wavelength=ref_wvl).isel(strain=i))
        )[0]
        b_boxes.append(props.bbox)

    b_boxes = np.vstack(b_boxes)
    min_row = np.min(b_boxes[:, 0])
    min_col = np.min(b_boxes[:, 1])
    max_row = np.max(b_boxes[:, 2])
    max_col = np.max(b_boxes[:, 3])

    return min_row, min_col, max_row, max_col


def get_mvmt_pair_i(mvmt, pair):
    return mvmt.loc[pd.IndexSlice[:, pair], :]
