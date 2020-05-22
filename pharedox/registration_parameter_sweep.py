from sklearn.model_selection import ParameterGrid

from pharedox import io as pio
from pharedox import profile_processing as pp
import numpy as np
import logging
from pathlib import Path
import xarray as xr
from matlab import engine


LOG_FILE = "/Users/sean/Desktop/registration_log.log"
META_DIR = "/Users/sean/code/pharedox/data/paired_ratio/"
OUTPUT_FOLDER = "~/Desktop/reg_sweep/"

logging.basicConfig(
    filename=LOG_FILE,
    filemode="w",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s:%(message)s",
)

meta_dir = Path(META_DIR)
prof_raw = xr.concat(
    [
        xr.load_dataarray(x)
        for x in meta_dir.glob("**/2020-01-09_unregistered/*untrimmed*.nc")
    ],
    dim="animal",
)

# n_derivs = [1.0, 2.0]
# smooth_lambdas = [float(x) for x in np.power(10, np.linspace(0, 2, 15))]
warp_lambdas = [1.0, 10.0, 100.0, 1000.0, 10000.0]
warp_n_basis = np.linspace(10, 50, 15)
smooth_lambda = np.linspace(0.01, 1, 10)


all_params = list(
    ParameterGrid(
        {
            "n_deriv": [1.0,],
            "rough_lambda": [1e-2,],
            "rough_n_breaks": [200.0,],
            "rough_order": [4.0,],
            "smooth_lambda": smooth_lambda,
            "smooth_n_breaks": [50.0,],
            "smooth_order": [4.0,],
            "warp_lambda": warp_lambdas,
            "warp_n_basis": warp_n_basis,
            "warp_order": [4.0,],
        }
    )
)

eng = engine.connect_matlab()
i = 0
for reg_params in all_params:
    param_string = "-".join([f"{k}={v}" for k, v in reg_params.items()])
    logging.info(f"Starting Registration ({i}/{len(all_params)})")
    try:
        reg_data, _ = pp.channel_register(prof_raw, eng=eng, **reg_params)
        pio.save_profile_data(reg_data, OUTPUT_FOLDER + param_string + ".nc")
    except Exception as e:
        logging.error(f"Failed to register {param_string}")
        logging.error(str(e))
    finally:
        i += 1
