from pharynx_redox import (
    io as pio,
    image_processing as ip,
    experiment,
    plots,
    profile_processing as pp,
    utils,
    data_analysis as da,
    constants,
)

import numpy as np
import pandas as pd
from pathlib import Path
import logging
import xarray as xr

logging.basicConfig(
    format="%(asctime)s %(levelname)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%I:%M:%S",
)

meta_dir = Path("/Users/sean/code/pharynx_redox/data/paired_ratio")

prof_raw = xr.concat(
    [
        pio.load_profile_data(x)
        for x in sorted(meta_dir.glob("**/2020-01-09_unregistered/*untrimmed*.nc"))
    ],
    dim="animal",
)
prof_raw = prof_raw.assign_coords(
    {"position": np.linspace(0, 1, prof_raw.position.size)}
)

rot_fl = xr.load_dataarray(
    "/Users/sean/code/pharynx_redox/data/paired_ratio/all_rot_fl.nc"
).rename({"spec": "animal"})
rot_seg = xr.load_dataarray(
    "/Users/sean/code/pharynx_redox/data/paired_ratio/all_rot_seg.nc"
).rename({"spec": "animal"})

midlines = ip.calculate_midlines(rot_seg)
shifts = np.arange(-1, 1.25, step=0.25)

shift_fp = "~/Desktop/shift_data.nc"
try:
    logging.info("trying to load shift data")
    shifted_data = xr.load_dataarray(shift_fp)
except IOError:
    logging.info("no shift data found")
    logging.info("shifting data")
    shifted_data = da.synthetic_shift(rot_fl, midlines, shifts=shifts, n_points=200)
    logging.info(f"saving shift data to {shift_fp}")
    shifted_data.to_netcdf(shift_fp)
    logging.info("saved shift data")

data = shifted_data.sel(shift=0.25, direction="x")

logging.info("starting tabularizing")
da.fold_v_point_table(data, constants.untrimmed_regions_with_medial)
