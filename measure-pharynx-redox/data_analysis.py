from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd


def load_all_summaries(meta_dir: Union[Path, str]) -> pd.DataFrame:
    if isinstance(meta_dir, str):
        meta_dir = Path(meta_dir)

    return pd.concat(pd.read_csv(x) for x in sorted(meta_dir.glob("**/*summary*csv")))


def load_all_movement(meta_dir: Union[Path, str]) -> pd.DataFrame:
    if isinstance(meta_dir, str):
        meta_dir = Path(meta_dir)

    return pd.concat(pd.read_csv(x) for x in sorted(meta_dir.glob("**/*mvmt.csv")))


def get_resid_rr_pairs(pair1, pair2):
    return np.power(np.e, np.abs(np.log(pair1 / pair2) - 1))


def get_resid_rr(data):
    return (
        np.power(
            np.e,
            np.abs(
                np.log(
                    data.sel(wavelength="r", pair=0) / data.sel(wavelength="r", pair=1)
                )
            ),
        )
        - 1
    )


def pandas_df_to_markdown_table(df):
    fmt = ["---" for i in range(len(df.columns))]
    df_fmt = pd.DataFrame([fmt], columns=df.columns)
    df_formatted = pd.concat([df_fmt, df])
    return df_formatted.to_csv(sep="|", index=False)


def filter_only_moving_roi(df, pair, roi):
    return df[
        (df["movement"].loc[:, pair][roi] > 0)
        & (df["movement"].loc[:, pair].drop(roi, axis=1).sum(axis=1) == 0)
    ]


def split_by_movement_types(df, roi, t=0):
    """ Return a set of filtered DataFrames, according to the following scheme:

    m_0_0:
        no movement in the 0th or 1st pair
    m_0_1:
        no movement in the 0th pair, movement in the 1st pair
    m_1_0:
        movment in the 0th pair, no movement in the 1st pair
    m_1_1:
        movment in both pairs

    In each case, movement is classified as such if the movement call *within the given ROI* is greater than
    the specified threshold (default=0).

    returned in the following order:
        (m_0_0, m_0_1, m_1_0, m_1_1)


    """
    m_0_0 = df[
        ((df["movement"].loc[:, 0][roi] <= t) & (df["movement"].loc[:, 1][roi] <= t))
    ]
    m_0_1 = df[
        ((df["movement"].loc[:, 0][roi] <= t) & (df["movement"].loc[:, 1][roi] > t))
    ]
    m_1_0 = df[
        ((df["movement"].loc[:, 0][roi] > t) & (df["movement"].loc[:, 1][roi] <= t))
    ]
    m_1_1 = df[
        ((df["movement"].loc[:, 0][roi] > t) & (df["movement"].loc[:, 1][roi] > t))
    ]

    return m_0_0, m_0_1, m_1_0, m_1_1
