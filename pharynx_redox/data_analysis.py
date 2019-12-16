from pathlib import Path
from typing import Union, List

import numpy as np
import pandas as pd
import typing
import xarray as xr
from scipy import ndimage as ndi

from pandas.core.indexes.base import InvalidIndexError

from . import profile_processing
from . import pharynx_io as pio


def load_all_cached_profile_data(meta_dir, glob_pattern):
    try:
        return xr.concat(
            (pio.load_profile_data(p) for p in sorted(meta_dir.glob(glob_pattern))),
            dim="animal",
        )
    except InvalidIndexError:
        return xr.concat(
            (pio.load_profile_data(p) for p in sorted(meta_dir.glob(glob_pattern))),
            dim="strain",
        )


def load_all_summaries(meta_dir: Union[Path, str]) -> pd.DataFrame:
    if isinstance(meta_dir, str):
        meta_dir = Path(meta_dir)

    return pd.concat(
        (pd.read_csv(x) for x in sorted(meta_dir.glob("**/*summary*csv"))), sort=False
    )


def load_all_movement(meta_dir: Union[Path, str]) -> pd.DataFrame:
    if isinstance(meta_dir, str):
        meta_dir = Path(meta_dir)

    return pd.concat(pd.read_csv(x) for x in sorted(meta_dir.glob("**/*-mvmt.csv")))


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


def relative_error(data, summarize=False, **summarize_kwargs):
    try:
        pair0 = data.sel(wavelength="r", pair=0)
        pair1 = data.sel(wavelength="r", pair=1)
    except KeyError:
        pair0 = data.sel(wavelength="410", pair=0) / data.sel(wavelength="470", pair=0)
        pair1 = data.sel(wavelength="410", pair=1) / data.sel(wavelength="470", pair=1)

    return fold_error_pairs(pair0, pair1, summarize=summarize, **summarize_kwargs)


def fold_error(data, summarize=False, **summarize_kwargs):
    try:
        pair0 = data.sel(wavelength="r", pair=0)
        pair1 = data.sel(wavelength="r", pair=1)
    except KeyError:
        pair0 = data.sel(wavelength="410", pair=0) / data.sel(wavelength="470", pair=0)
        pair1 = data.sel(wavelength="410", pair=1) / data.sel(wavelength="470", pair=1)

    return fold_error_pairs(pair0, pair1, summarize=summarize, **summarize_kwargs)


def fold_error_pairs(pair0, pair1, summarize=False, **summarize_kwargs):
    with np.errstate(divide="ignore"):
        prof_data = np.power(np.e, np.abs(np.log((pair0 / pair1)))) - 1
        # prof_data = np.abs(1 - (pair0 / pair1))
    if summarize:
        summary_table = profile_processing.summarize_over_regions(
            prof_data, **summarize_kwargs
        )
        return prof_data, summary_table
    else:
        return prof_data


def filter_only_moving_roi(df, pair, roi):
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


def split_by_movement_types(df, roi, t=0):
    # TODO: write function that generates the required format for this function
    """ Return a set of filtered DataFrames, according to the following scheme:

    m_0_0:
        no movement in the 0th or 1st pair
    m_0_1:
        no movement in the 0th pair, movement in the 1st pair
    m_1_0:
        movement in the 0th pair, no movement in the 1st pair
    m_1_1:
        movement in both pairs

    In each case, movement is classified as such if the movement call *within the given
    ROI* is greater than the specified threshold (default=0).

    Requires a pandas DataFrame in the following format::

                                    movement                                       
        region                      anterior      posterior       sides_of_tip    tip   
        pair                                 0  1         0  1            0  1   0  1
        experiment             animal                                                
        2017_02_22-HD233_SAY47 0             0  0         0  0            0  0   1  0
                               1             0  1         0  1            1  1   0  1
                               2             0  0         0  0            0  0   0  0
                               3             0  0         0  0            0  0   0  0
                               4             0  0         0  0            0  0   0  0

    See Also
    --------
    mvmt_long_to_wide
        for generating the DataFrame with the required format

    returned in the following order::

        (m_0_0, m_0_1, m_1_0, m_1_1)
    """

    m_0_0 = df[(df["movement"][roi][0] <= t) & (df["movement"][roi][1] <= t)]
    m_0_1 = df[(df["movement"][roi][0] <= t) & (df["movement"][roi][1] > t)]
    m_1_0 = df[(df["movement"][roi][0] > t) & (df["movement"][roi][1] <= t)]
    m_1_1 = df[(df["movement"][roi][0] > t) & (df["movement"][roi][1] > t)]
    return m_0_0, m_0_1, m_1_0, m_1_1


def label_by_movement_types(df, t=0):
    pass


def synthetic_shift(
    rot_fl: xr.DataArray, shifts: List[float], n_points=300
) -> xr.DataArray:
    shift_data = xr.DataArray(
        0.0,
        coords={
            "animal": rot_fl.animal.values,
            "pair": rot_fl.pair.values,
            "wavelengths": ["410", "470", "r"],
            "shifts": shifts,
            "direction": ["x", "y"],
            "position": np.arange(n_points),
        },
        dims=["strain", "pair", "wavelength", "shift", "direction", "position"],
    )

    for pair in range(rot_fl.pair.size):
        for shift in shifts:
            for direction in ["x", "y"]:
                i410_ = np.zeros((rot_fl.strain.size, n_points), dtype=rot_fl.dtype)
                i470_ = i410_.copy()

                for img_idx in range(rot_fl.strain.size):
                    if direction == "x":
                        dx, dy = shift, 0
                    if direction == "y":
                        dx, dy = 0, shift

                    # select image
                    mid_xs, mid_ys = midlines[img_idx]["410"][pair].linspace(n=200)
                    im410 = rot_fl.sel(wavelength="410", pair=pair)[img_idx]
                    im470 = rot_fl.sel(wavelength="470", pair=pair)[img_idx]

                    # measure under midline
                    i410_[img_idx, :] = ndi.map_coordinates(
                        im410, np.stack([m_ys, m_xs]), order=1
                    )
                    i470_[img_idx, :] = ndi.map_coordinates(
                        im470, np.stack([m_ys + dy, m_xs + dx]), order=1
                    )

                # smooth / resample
                new_xs = np.linspace(0, 1, n_points)
                i410_sm = np.squeeze(
                    profile_processing.smooth_profile_data(i410_)[0](new_xs)
                )
                i470_sm = np.squeeze(
                    profile_processing.smooth_profile_data(i470_)[0](new_xs)
                )

                # save
                shift_data.loc[
                    dict(pair=pair, wavelength="410", shift=shift, direction=direction)
                ] = i410_sm
                shift_data.loc[
                    dict(pair=pair, wavelength="470", shift=shift, direction=direction)
                ] = i470_sm
                shift_data.loc[
                    dict(pair=pair, wavelength="r", shift=shift, direction=direction)
                ] = (i410_sm / i470_sm)
