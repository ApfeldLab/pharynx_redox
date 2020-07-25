import logging
import os
from pathlib import Path
import numpy as np
import pytest
import requests
from tqdm import tqdm

from pharedox import profile_processing as pp, utils, pio

test_data_path = Path(os.path.join(os.path.dirname(__file__), "test_data"))


class TestUtils:
    @pytest.fixture(scope="function")
    def paired_imgs(self, shared_datadir):

        return pio.load_tiff_as_hyperstack(
            (
                shared_datadir
                / "experiments"
                / "2017_02_22-HD233_SAY47"
                / "2017_02_22-HD233_SAY47.tif"
            ),
            manual_metadata=(
                shared_datadir
                / "experiments"
                / "2017_02_22-HD233_SAY47"
                / "2017_02_22-HD233_SAY47-frame_map.csv"
            ),
            mvmt_metadata=(
                shared_datadir
                / "experiments"
                / "2017_02_22-HD233_SAY47"
                / "2017_02_22-HD233_SAY47-mvmt.csv"
            ),
        )

    def test_expand_dimension(self, paired_imgs):
        data = paired_imgs
        r = data.sel(wavelength="410") / data.sel(wavelength="470")
        oxd = pp.r_to_oxd(r)

        data = utils.expand_dimension(
            data, dim="wavelength", new_coords={"r": r, "oxd": oxd}
        )

        assert np.array_equal(data.sel(wavelength="r").values, r.values)
        assert np.array_equal(data.sel(wavelength="oxd").values, oxd.values)

    def test_add_derived_wavelengths(self, paired_imgs):
        data = paired_imgs

        r = data.sel(wavelength="410") / data.sel(wavelength="470")
        oxd = pp.r_to_oxd(r)
        e = pp.oxd_to_redox_potential(oxd)

        data = utils.add_derived_wavelengths(data)

        assert np.allclose(data.sel(wavelength="r").values, r.values, equal_nan=True)
        assert np.allclose(
            data.sel(wavelength="oxd").values, oxd.values, equal_nan=True
        )
        assert np.allclose(data.sel(wavelength="e").values, e.values, equal_nan=True)
