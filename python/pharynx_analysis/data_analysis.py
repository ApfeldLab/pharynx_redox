from typing import Union

import pandas as pd
from pathlib import Path


def load_all_summaries(meta_dir: Union[Path, str]) -> pd.DataFrame:
    if isinstance(meta_dir, str):
        meta_dir = Path(meta_dir)

    return pd.concat(pd.read_csv(x) for x in meta_dir.glob('**/*summary*csv'))
