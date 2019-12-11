import pharynx_io as pio
import image_processing as ip
import experiment
import plots
import profile_processing as pp
import utils
import data_analysis as da
import numpy as np
import matlab
import logging
from pathlib import Path

logging.basicConfig(
    filename="/Users/sean/Desktop/registration_log.log",
    filemode="w",
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s:%(message)s",
)

prof_raw = pio.load_profile_data(
    "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/analyses/2019-10-11_unreg/2017_02_22-HD233_SAY47-untrimmed_profile_data.nc"
)

meta_dir = Path(
    "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47"
)
prof_raw = da.load_all_cached_profile_data(meta_dir, "**/*11_unreg/*-untrimmed*.nc")

eng = matlab.engine.connect_matlab()

n_derivs = [1.0, 2.0]
rough_lambdas = [1.0e-3, 1.0e-2, 1.0e-1, 1]
smooth_lambdas = [1.0e-1, 1.0, 10.0, 100.0]
warp_lambdas = [1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e3, 1.0e4]

all_params = []
for n_deriv in n_derivs:
    for sl in smooth_lambdas:
        for wl in warp_lambdas:
            for rl in rough_lambdas:
                all_params.append(
                    {
                        "n_deriv": n_deriv,
                        "rough_lambda": rl,
                        "rough_n_breaks": 300.0,
                        "rough_order": 6.0,
                        "smooth_lambda": sl,
                        "smooth_n_breaks": 100.0,
                        "smooth_order": 6.0,
                        "warp_lambda": wl,
                        "warp_n_basis": 10.0,
                        "warp_order": 4.0,
                    }
                )

i = 0
for reg_params in all_params:
    param_string = "_".join([f"{k}={v}" for k, v in reg_params.items()])
    logging.info(f"Starting Registration ({i}/{len(all_params)})")
    try:
        reg_data = pp.register_profiles_pairs(prof_raw, eng=eng, **reg_params)
        pio.save_profile_data(reg_data, f"/tmp/{param_string}.nc")
        logging.info(f"Registered {param_string}")
    except Exception as e:
        logging.error(f"Failed to register {param_string}")
        logging.error(str(e))
    finally:
        i += 1
