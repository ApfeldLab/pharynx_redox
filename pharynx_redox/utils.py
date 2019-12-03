import re
from collections import Counter
from pharynx_redox import experiment

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage as ndi
import typing
from skimage.measure import regionprops, label
from pathlib import Path
from tqdm import tqdm
import argparse


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
    """
    Given a list of things, return a list of tuples ``(item, nth_occurrence)``

    .. todo::
        test

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

    for i in range(rot_seg_stack.spec.size):
        props = regionprops(
            label(rot_seg_stack.sel(pair=ref_pair, wavelength=ref_wvl).isel(spec=i))
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

    From: https://github.com/django/django/blob/master/django/utils/text.py

    Examples
    --------

        >>> get_valid_filename("john's portrait in 2004.jpg")
        'johns_portrait_in_2004.jpg'

    Parameters
    ----------
    s
        a string to convert to valid filename

    Returns
    -------
    str
        a valid filename

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
    try:
        data = np.squeeze(fd.data_matrix)
    except AttributeError:
        data = fd
    means = np.mean(data, axis=1)
    stds = np.std(data, axis=1)
    try:
        fd.data_matrix = ((data.T - means) / stds).T
    except AttributeError:
        return ((data.T - means) / stds).T
    return fd


def measure_shifted_midlines(
    experiment,
    shift_range: typing.Iterable[float],
    shift_steps: int,
    n_points: int = 200,
) -> (xr.DataArray, typing.List[float], typing.List[int]):
    """
    Measure under shifted midlines for use in synthetic movement analysis

    First, **only non-moving animals** are selected for analysis (this guarantees that
    for 0-shifts, there will be no errors).

    shifted midlines are then generated according to the following pseudocode::

        for each pair:
            for each animal:
                for each shift (dx):
                    measure under 410 and 470 using the 410 midline, storing the resultant
                    measurements as pair 0

                    shift the 470 midline by dx

                    measure under 410 and 470 using 410-dx midline, storing the resultant
                    measurements as pair 1

    Parameters
    ----------
    experiment
        The experiment to analyze
    shift_range
        a tuple of ``(min, max)`` indicating the bounds of the shifts
    shift_steps
        the number of shifts to use within the bounds of ``shift_range``
    n_points
        the number of points to measure under each midline

    Returns
    -------
    (xr.DataArray, typing.List[float], typing.List[int])
        a tuple containing:

            - the measurements
            - a list of shifts (where the index corresponds to the first dimension of
              the measurements (the animal))
            - a list of original indices (i.e. the experimental indices) for each
              synthetic index

    """
    shifts = np.linspace(*shift_range, shift_steps)

    ex_meas = experiment.trimmed_profiles

    df = experiment.movement
    df = pd.DataFrame(df.to_records())

    # indexed by pair
    stationary_animals = np.array(
        [
            df[
                (df.pair == pair) & (df.anterior == 0) & (df.posterior == 0)
            ].animal.values
            for pair in df.pair.unique()
        ]
    )

    measurements = xr.DataArray(
        np.zeros(
            (
                len(np.concatenate(stationary_animals)) * len(shifts),
                ex_meas.wavelength.size,
                ex_meas.pair.size,
                n_points,
            )
        ),
        dims=["spec", "wavelength", "pair", "position"],
        coords={"wavelength": ex_meas.wavelength, "pair": ex_meas.pair},
    )

    all_shifts = []
    orig_idx = []
    new_animal_idx = 0
    for pair in range(stationary_animals.shape[0]):
        for orig_animal_idx in stationary_animals[pair]:
            midline = experiment.midlines[orig_animal_idx]["410"][pair]

            # TODO: generalize over wavelengths
            i410 = experiment.rot_fl.sel(wavelength="410", pair=pair)[
                orig_animal_idx
            ].values
            i470 = experiment.rot_fl.sel(wavelength="470", pair=pair)[
                orig_animal_idx
            ].values

            for shift in shifts:
                all_shifts.append(shift)
                orig_idx.append(orig_animal_idx)

                xs, ys = midline.linspace(n=n_points)

                # First, measure under 410 and 470 like normal, (will be pair 0)
                unshifted_410 = ndi.map_coordinates(i410, np.stack([ys, xs]), order=1)
                unshifted_470 = ndi.map_coordinates(i470, np.stack([ys, xs]), order=1)

                # now, shift the coordinates of the midline and measure under 470
                # (will be pair 1)
                shifted_470 = ndi.map_coordinates(
                    i470, np.stack([ys, xs + shift]), order=1
                )

                measurements[new_animal_idx].loc[
                    {"pair": 0, "wavelength": "410"}
                ] = unshifted_410
                measurements[new_animal_idx].loc[
                    {"pair": 0, "wavelength": "470"}
                ] = unshifted_470

                measurements[new_animal_idx].loc[
                    {"pair": 1, "wavelength": "410"}
                ] = unshifted_410
                measurements[new_animal_idx].loc[
                    {"pair": 1, "wavelength": "470"}
                ] = shifted_470

                new_animal_idx += 1

    return measurements, all_shifts, orig_idx


def run_all_analyses(meta_dir: str, imaging_scheme: str, **kwargs):

    meta_dir = Path(meta_dir)

    exps = list(filter(lambda x: x.is_dir(), meta_dir.iterdir()))
    for exp_dir in tqdm(exps):
        experiment.PairExperiment(
            experiment_dir=exp_dir,
            imaging_scheme=imaging_scheme,
            strategy="_".join(f"{k}={v}" for k, v in kwargs.items()),
        ).full_pipeline()


def cli_run_all_analyses():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "meta_dir",
        metavar="D",
        type=str,
        help="the directory containing all experiment directories",
    )
    parser.add_argument(
        "imaging_scheme",
        metavar="S",
        type=str,
        help="the imaging scheme (e.g. TL/470/410/470/410)",
    )

    args = parser.parse_args()
    run_all_analyses(meta_dir=args.meta_dir, imaging_scheme=args.imaging_scheme)


if __name__ == "__main__":
    cli_run_all_analyses()
