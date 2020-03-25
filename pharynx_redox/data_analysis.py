from pathlib import Path
from typing import Union, List

import numpy as np
import pandas as pd
import typing
import xarray as xr
from scipy import ndimage as ndi
from functools import reduce

from pandas.core.indexes.base import InvalidIndexError

from pharynx_redox import profile_processing
from pharynx_redox import io as pio


def fold_v_point_table(
    data: xr.DataArray, regions: dict, moving_regions: typing.List[str] = None, **kwargs
) -> pd.DataFrame:
    """
    summarize the given data with errors using both the point-wise and region-wise
    formulations.
    
    Parameters
    ----------
    data : xr.DataArray
        the data to summarize. expected to have animal, wavelength, pair, position
        dimensions.
    regions: dict
        a region dictionary. should map `[region -> [left_bound, right_bound]]` where
        `region` is a string and `left_bound` and `right_bound` are floats from 
        `[0, 1]`. See the `pharynx_redox.constants` module
    
    Returns
    -------
    pd.DataFrame
        a summary DataFrame with 
    """
    df = []
    for wvl, val_name in [
        ("410", "I410"),
        ("470", "I470"),
        ("r", "R_point"),
    ]:
        sub_df = []
        for pair in [0, 1]:
            pair_df = profile_processing.summarize_over_regions(
                data.sel(wavelength=wvl, pair=pair), regions, value_name=val_name,
            )
            pair_df["pair"] = pair
            sub_df.append(pair_df)
        sub_df = pd.concat(sub_df, sort=False)
        df.append(sub_df)

    df[0]["I470"] = df[1]["I470"]
    df[0]["R_point"] = df[2]["R_point"]
    df = df[0]
    df["R_region"] = df["I410"] / df["I470"]

    ## Now add the fold_errors
    # we have to do this separately and massage the data a little bit since fold_error
    # gets ride of pair information... so we are essentially duplicating the errors
    # across pairs, getting rid of the movement annotations, then merging into the
    # big table

    fold_df_pairs = []
    for pair in [0, 1]:
        fold_error_df = profile_processing.summarize_over_regions(
            fold_error(data),
            regions=regions,
            rescale=True,
            add_attrs=False,
            value_name="fold_error_point",
        )
        fold_error_df["pair"] = pair
        fold_df_pairs.append(fold_error_df)
    fold_error_df = pd.concat(fold_df_pairs)
    fold_error_df = fold_error_df[
        fold_error_df.columns.drop(list(fold_error_df.filter(regex="^mvmt-")))
    ]

    df["fold_error_point"] = fold_error_df["fold_error_point"]

    ## Now we stack the pair data

    df = df.set_index(["region", "pair"], append=True).unstack()

    fold_error_region = fold_error_pairs(
        df["R_region"][0].values, df["R_region"][1].values
    )
    df[("fold_error_region", 0)] = fold_error_region
    df[("fold_error_region", 1)] = fold_error_region

    for k, v in kwargs.items():
        df[k] = v

    if moving_regions is not None:
        mvmt_cols = [f"mvmt-{region}" for region in moving_regions]
        pair_0_mvmt = (df[mvmt_cols].xs(0, level=1, axis=1).sum(axis=1) > 0).values
        pair_1_mvmt = (df[mvmt_cols].xs(1, level=1, axis=1).sum(axis=1) > 0).values

        st_bools = (pair_0_mvmt == False) & (pair_1_mvmt == False)
        mv_bools = ((pair_0_mvmt == False) & (pair_1_mvmt == True)) | (
            (pair_0_mvmt == True) & (pair_1_mvmt == False)
        )
        neither_bools = pair_0_mvmt & pair_1_mvmt

        df["moving-AP"] = 0
        df.loc[mv_bools, "moving-AP"] = 1
        df.loc[neither_bools, "moving-AP"] = -1

    return df


def get_mvmt_codes(data: xr.DataArray, regions: typing.Union[str, typing.List[str]]):
    mvmts = np.repeat(-1, len(data))
    mv, st = get_moving_idx(data, regions)
    mvmts[mv] = 1
    mvmts[st] = 0

    return mvmts


def get_moving_idx(
    data: xr.DataArray,
    regions: typing.Union[str, typing.List[str]],
    mvmt_thresh: int = 2,
) -> typing.Tuple[xr.DataArray, xr.DataArray]:
    """
    Returns a boolean indexing array to select the moving animals from the data array.

    The resultant array will be the same length as the `animal` index of the input
    array, allowing one to use it to select the moving (or stationary) animals from the
    original array by using traditional numpy boolean indexing.
    
    Parameters
    ----------
    data : xr.DataArray
        The data
    regions : typing.Union[str, typing.List[str]]
        which region(s) should be considered for movement stratification
    
    Returns
    -------
    moving: xr.DataArray
        A boolean array where `True` indicates that the animal was moving at that index
        and `False` indicates that the animal was either stationary or unacceptable (as
        in the case where the animal moves in both the 0th and 1st pair)
    stationary: xr.DataArray
        A boolean array where `True` indicates that the animal was stationary at that index
        and `False` indicates that the animal was either moving or unacceptable (as
        in the case where the animal moves in both the 0th and 1st pair)
    """
    if type(regions) == str:
        regions = [regions]

    re_str = f"mvmt-({'|'.join(regions)})$"
    mvmt_df = (
        data.coords.to_dataset()
        .to_dataframe()
        .filter(regex=re_str, axis=1)
        .groupby(["animal", "timepoint", "pair"])
        .max()
        .applymap(lambda x: x >= 2)
    )
    mvmt_df = mvmt_df.any("columns").unstack("pair")

    mv_df = np.logical_xor(mvmt_df[0], mvmt_df[1])
    mv_idx = {
        k: np.unique(v)
        for k, v in mv_df[mv_df].index.to_frame().to_dict(orient="list").items()
    }

    st_df = (
        data.coords.to_dataset()
        .to_dataframe()
        .filter(regex=re_str, axis=1)
        .groupby(["animal", "timepoint", "pair"])
        .max()
        .applymap(lambda x: x == 0)
        .all("columns")
        .unstack("pair")
        .all("columns")
    )
    st_idx = {
        k: np.unique(v)
        for k, v in st_df[st_df].index.to_frame().to_dict(orient="list").items()
    }

    return mv_idx, st_idx  # , mv_both_idx, mv_first_idx


def resample_moving(
    data: xr.DataArray,
    p_moving: float,
    n: float,
    regions: typing.Union[str, typing.List[str]],
) -> xr.DataArray:
    """
    Sample from moving and stationary populations in proportion to a hypothetical
    scenario wherein ``p`` percent animals are moving.
    
    Parameters
    ----------
    data : xr.DataArray
        the Profile Data to resampel
    p_moving : float
        the percent moving for the hypothetical population
    n : float
        the hypothetical populationsize
    regions : typing.Union[str, typing.List[str]]
        which region(s) should be considered for movement stratification

    Returns
    -------
    xr.DataArray
        the resampled hypothetical population
    """
    mv_idx, st_idx = get_moving_idx(data, regions)

    # Create the probability of selection map - each idx indicates the probability of
    # selection during re-sampling... initialize with 0s for all idx
    p_selected = np.zeros((data.shape[0],))

    # we need to divide by the number of animals in each movement group
    p_selected[mv_idx] = p_moving / np.sum(mv_idx)
    p_selected[st_idx] = (1 - p_moving) / np.sum(st_idx)

    # Finally, normalize the entire probability array so that it sums to 1, which is
    # required by np.random.choice
    p_selected = p_selected / p_selected.sum(0)

    # The idx that are False for mv_idx AND st_idx were never assigned probabilities,
    # so they stay at 0.0
    # This is what we want, where if the animal moved in both frames, we are not
    # interested in it

    idx = np.random.choice(np.arange(data.shape[0]), size=n, replace=True, p=p_selected)

    return data[idx]


def select_by_mvmt(
    prof_data: xr.DataArray, regions: typing.Union[str, typing.List[str]]
) -> typing.Union[xr.DataArray, xr.DataArray]:
    """
    Return the data in the given DataArray, filtered according to movement in the 
    specified region.

    Stationary animals are those that are labelled as stationary in both pair 0 and pair
    1.

    Moving animals are those that are labelled as stationary in pair 0 and moving in 
    pair 1.
    
    Parameters
    ----------
    prof_data : xr.DataArray
        the profile data to filter
    region: str
        the region to filter movement on
    moving : bool
        if True, returns moving animals. If False, returns stationary animals.
    
    Returns
    -------
    moving: xr.DataArray
        [description]
    stationary: xr.DataArray
    """
    mv_idx, st_idx = get_moving_idx(prof_data, regions)

    return prof_data[mv_idx], prof_data[st_idx]


def load_all_cached_profile_data(meta_dir, glob_pattern):
    return xr.concat(
        (pio.load_profile_data(p) for p in sorted(meta_dir.glob(glob_pattern))),
        dim="animal",
    )


def get_resid_rr_pairs(
    pair0, pair1, summarize=False, **summarize_kwargs
) -> Union[xr.DataArray, typing.Tuple[xr.DataArray, pd.DataFrame]]:
    with np.errstate(divide="ignore"):
        # prof_data = np.power(np.e, np.abs(np.log((pair0 / pair1)))) - 1
        prof_data = 1 - (pair0 / pair1)
    if summarize:
        summary_table = profile_processing.summarize_over_regions(
            prof_data, **summarize_kwargs
        )
        return prof_data, summary_table
    else:
        return prof_data


def get_resid_rr(
    data: xr.DataArray, summarize=False, **summarize_kwargs
) -> Union[xr.DataArray, typing.Tuple[xr.DataArray, pd.DataFrame]]:
    try:
        pair0 = data.sel(wavelength="r", pair=0)
        pair1 = data.sel(wavelength="r", pair=1)
    except KeyError:
        pair0 = data.sel(wavelength="410", pair=0) / data.sel(wavelength="470", pair=0)
        pair1 = data.sel(wavelength="410", pair=1) / data.sel(wavelength="470", pair=1)

    return get_resid_rr_pairs(pair0, pair1, summarize=summarize, **summarize_kwargs)


def relative_error(data):
    try:
        pair0 = data.sel(wavelength="r", pair=0)
        pair1 = data.sel(wavelength="r", pair=1)
    except KeyError:
        pair0 = data.sel(wavelength="410", pair=0) / data.sel(wavelength="470", pair=0)
        pair1 = data.sel(wavelength="410", pair=1) / data.sel(wavelength="470", pair=1)
    return relative_error_pairs(pair0, pair1)


def relative_error_pairs(pair0, pair1):
    err = 1 - (pair1 / pair0)

    err = err.assign_attrs(pair0.attrs).assign_attrs(pair1.attrs)
    for region in ["posterior", "anterior", "sides_of_tip", "tip"]:
        mvmt_coords = {
            f"mvmt-{region}": (
                ("animal",),
                profile_processing.get_mvmt(
                    pair0[f"mvmt-{region}"].values, pair1[f"mvmt-{region}"].values
                ),
            )
        }
        err = err.assign_coords(mvmt_coords)
    return err


def fold_error(data):
    try:
        pair0 = data.sel(wavelength="r", pair=0)
        pair1 = data.sel(wavelength="r", pair=1)
    except KeyError:
        pair0 = data.sel(wavelength="410", pair=0) / data.sel(wavelength="470", pair=0)
        pair1 = data.sel(wavelength="410", pair=1) / data.sel(wavelength="470", pair=1)

    return fold_error_pairs(pair0, pair1)


def fold_error_pairs(pair0: xr.DataArray, pair1: xr.DataArray):
    with np.errstate(divide="ignore"):
        prof_data = np.power(np.e, np.abs(np.log((pair0 / pair1)))) - 1

    try:
        prof_data = prof_data.assign_attrs(pair0.attrs).assign_attrs(pair1.attrs)

        try:
            for region in ["posterior", "anterior", "sides_of_tip", "tip"]:
                mvmt_coords = {
                    f"mvmt-{region}": (
                        ("animal",),
                        profile_processing.get_mvmt(
                            pair0[f"mvmt-{region}"].values,
                            pair1[f"mvmt-{region}"].values,
                        ),
                    )
                }
                prof_data = prof_data.assign_coords(mvmt_coords)
        except:
            pass
    except AttributeError:
        # this is if the pairs are numpy arrays not dataarrays
        return prof_data

    return prof_data


def filter_only_moving_roi(df, pair, roi):
    """
    [summary]
    
    Parameters
    ----------
    df : [type]
        [description]
    pair : [type]
        [description]
    roi : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """
    other_pair = (pair - 1) % 2
    return df[
        (df["movement"].loc[:, pair][roi] > 0)
        & (df["movement"].loc[:, pair].drop(roi, axis=1).sum(axis=1) == 0)
        & (df["movement"].loc[:, other_pair].sum(axis=1) == 0)
    ]


def mvmt_long_to_wide(mvmt_df):
    """
    Pivot the given dataframe from the following structure::

                       experiment  animal  movement  pair    region
        0  2017_02_22-HD233_SAY47       0         0     0  anterior
        1  2017_02_22-HD233_SAY47       1         0     0  anterior
        2  2017_02_22-HD233_SAY47       2         0     0  anterior
        3  2017_02_22-HD233_SAY47       3         0     0  anterior
        4  2017_02_22-HD233_SAY47       4         0     0  anterior

    to the following structure::

                                    movement                                       
        region                      anterior      posterior       sides_of_tip    tip   
        pair                                 0  1         0  1            0  1   0  1
        experiment             animal                                                
        2017_02_22-HD233_SAY47 0             0  0         0  0            0  0   1  0
                               1             0  1         0  1            1  1   0  1
                               2             0  0         0  0            0  0   0  0
                               3             0  0         0  0            0  0   0  0
                               4             0  0         0  0            0  0   0  0
    """

    return mvmt_df.pivot_table(
        index=["experiment", "animal"], columns=["region", "pair"], values=["movement"]
    )


def synthetic_shift(
    rot_fl: xr.DataArray, midlines: xr.DataArray, shifts: List[float], n_points=200
) -> xr.DataArray:
    """
    Move the midlines by `shifts` and measure under them, simulating movement in the
    X- and Y- direction.

    Parameters
    ----------
    rot_fl : xr.DataArray
        the images to measure under
    midlines : xr.DataArray
        the midlines to measure with
    shifts : List[float]
        the shift magnitudes to use
    n_points : int, optional
        the number of points to measure under the midline with, by default 200

    Returns
    -------
    xr.DataArray
        the shifted data, dimensions: 
            `(animal, pair, wavelength, shift, direction, position)`
    """
    shift_data = xr.DataArray(
        0.0,
        coords={
            "animal": rot_fl.animal.values,
            "pair": rot_fl.pair.values,
            "wavelength": ["410", "470", "r"],
            "shift": shifts,
            "direction": ["x", "y"],
            "position": np.arange(n_points),
        },
        dims=["animal", "pair", "wavelength", "shift", "direction", "position"],
    )

    for pair in range(rot_fl.pair.size):
        for shift in shifts:
            for direction in ["x", "y"]:
                i410_ = np.zeros((rot_fl.animal.size, n_points), dtype=rot_fl.dtype)
                i470_ = i410_.copy()

                for img_idx in range(rot_fl.animal.size):
                    if direction == "x":
                        dx, dy = shift, 0
                    if direction == "y":
                        dx, dy = 0, shift

                    # select image
                    mid_xs, mid_ys = (
                        midlines.sel(wavelength="410", pair=pair)[img_idx]
                        .values[()]
                        .linspace(n=n_points)
                    )
                    im410 = rot_fl.sel(wavelength="410", pair=pair)[img_idx]
                    im470 = rot_fl.sel(wavelength="410", pair=pair)[img_idx]

                    # measure under midline
                    i410_[img_idx, :] = ndi.map_coordinates(
                        im410, np.stack([mid_ys, mid_xs]), order=1
                    )
                    i470_[img_idx, :] = ndi.map_coordinates(
                        im470, np.stack([mid_ys + dy, mid_xs + dx]), order=1
                    )

                # save
                shift_data.loc[
                    dict(pair=pair, wavelength="410", shift=shift, direction=direction)
                ] = i410_
                shift_data.loc[
                    dict(pair=pair, wavelength="470", shift=shift, direction=direction)
                ] = i470_
                shift_data.loc[
                    dict(pair=pair, wavelength="r", shift=shift, direction=direction)
                ] = (i410_ / i470_)
    return shift_data
