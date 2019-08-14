from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd


def load_all_summaries(meta_dir: Union[Path, str]) -> pd.DataFrame:
    if isinstance(meta_dir, str):
        meta_dir = Path(meta_dir)

    return pd.concat(pd.read_csv(x) for x in sorted(meta_dir.glob('**/*summary*csv')))


def load_all_movement(meta_dir: Union[Path, str]) -> pd.DataFrame:
    if isinstance(meta_dir, str):
        meta_dir = Path(meta_dir)

    return pd.concat(pd.read_csv(x) for x in sorted(meta_dir.glob('**/*mvmt.csv')))


def get_resid_rr(data):
    return np.power(
        np.e,
        np.abs(np.log(data.sel(wavelength='r', pair=0) / data.sel(wavelength='r', pair=1)))) - 1
