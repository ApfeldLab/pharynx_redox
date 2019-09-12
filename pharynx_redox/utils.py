import re
from collections import Counter

import numpy as np
import pandas as pd
import xarray as xr
import typing
from skimage.measure import regionprops, label


def scale_region_boundaries(regions: dict, profile_length: int):
    """
    Scale the region boundaries

    the length of the profiles may be changed as a parameter of the experiment. As such,
    region boundaries are expressed through percentages of profile length, then scaled
    to match the length of the profile vectors.

    Examples
    --------

        >>> regions = {'a': [0, 0.1], 'b': [.2, .4]}
        >>> scale_region_boundaries(regions, 10)
        {'a': [0, 1], 'b': [2, 4]}


    Parameters
    ----------
    regions
        a dictionary mapping region names to [left_bound, right_bound] expressed as
        percentages of the profile
    profile_length
        the length of the profile that the resultant scaled regions will be used on

    Returns
    -------
    dict
        the scaled region boundary map

    """
    return {
        region: np.int_(profile_length * np.asarray(regions[region]))
        for region in regions.keys()
    }


def rmse(u, v, axis=0):
    """
    Calculate the root mean squared error (RMSE) between the two given vectors

    Parameters
    ----------
    u
        the first vector
    v
        the second vector
    axis
        the axis over which the mean should be taken

    Returns
    -------
    float
        RMSE between the two vectors

    """
    return np.sqrt(np.mean((u - v) ** 2, axis=axis))


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


def create_occurrence_count_tuples(l: typing.Iterable) -> [(typing.Any, int)]:
    # TODO: test
    """
    Given a list of things, return a list of tuples ``(item, nth_occurrence)``


    Parameters
    ----------
    l
        the list of things

    Returns
    -------
    occurrence_count_list
        a list of tuples ``(item, nth_occurrence)``

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


def calc_max_bbox(
    rot_seg_stack: xr.DataArray, ref_pair: int = 0, ref_wvl: str = "410"
) -> (int, int, int, int):
    """
    Calculate the smallest bounding-box that encapsulates all of the the pharynxes in
    the given set of rotated, segmented pharynxes.


    Parameters
    ----------
    rot_seg_stack
        the rotated and segmented stack of pharynx images
    ref_pair
        the pair to consider when calculating bounding boxes
    ref_wvl
        the wavelength to consider when calculating bounding boxes

    Returns
    -------
    a tuple of ``(min_row, min_col, max_row, max_col)``

    """
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


def get_valid_filename(s):
    """
    Return the given string converted to a string that can be used for a clean
    filename. Remove leading and trailing spaces; convert other spaces to
    underscores; and remove anything that is not an alphanumeric, dash,
    underscore, or dot.
    >>> get_valid_filename("john's portrait in 2004.jpg")
    'johns_portrait_in_2004.jpg'

    From https://github.com/django/django/blob/master/django/utils/text.py
    """
    s = str(s).strip().replace(" ", "_")
    return re.sub(r"(?u)[^-\w.]", "", s)


def z_transform(fd):
    """
    Z-standardize a functional data grid object

    computed by taking (x-mean(x))/std(x)

    Parameters
    ----------
    fd

    Returns
    -------

    """
    data = np.squeeze(fd.data_matrix)
    means = np.mean(data, axis=1)
    stds = np.std(data, axis=1)
    fd.data_matrix = ((data.T - means) / stds).T
    return fd
