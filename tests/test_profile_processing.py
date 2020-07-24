import logging

import pytest

from pharedox import pio
from pharedox import profile_processing as pp


class TestProfileProcessing:

    redox_params = {
        "ratio_numerator": "410",
        "ratio_denominator": "470",
        "r_min": 0.852,
        "r_max": 6.65,
        "instrument_factor": 0.171,
        "midpoint_potential": -265.0,
        "z": 2,
        "temperature": 22.0,
    }

    reg_d2_parameters = {
        "n_deriv": 2.0,
        "warp_n_basis": 30.0,
        "warp_order": 4.0,
        "warp_lambda": 5000.0,
        "smooth_lambda": 100.0,
        "smooth_n_breaks": 100.0,
        "smooth_order": 4.0,
        "rough_lambda": 0.01,
        "rough_n_breaks": 300.0,
        "rough_order": 4.0,
    }

    @pytest.fixture(scope="module")
    def matlab_engine(self):
        try:
            import matlab.engine

            return matlab.engine.start_matlab()
        except ModuleNotFoundError:
            logging.warning("no matlab installation found, returning None as engine")
            return None

    @pytest.mark.matlab
    def test_channel_registration(self, matlab_engine, shared_datadir):
        raw_data = pio.load_profile_data(
            shared_datadir
            / "experiments"
            / "2017_02_23-HD233_HD236"
            / "analyses"
            / "2020-05-22"
            / "2017_02_23-HD233_HD236-untrimmed_profile_data.nc"
        ).isel(animal=[0, 1, 2])

        reg_data, warp_data = pp.channel_register(
            raw_data,
            redox_params=self.redox_params,
            eng=matlab_engine,
            reg_params=self.reg_d2_parameters,
        )

        assert raw_data.animal.size == reg_data.animal.size

    @pytest.mark.matlab
    def test_standardization(self, matlab_engine, shared_datadir):
        raw_data = pio.load_profile_data(
            shared_datadir
            / "experiments"
            / "2017_02_23-HD233_HD236"
            / "analyses"
            / "2020-05-22"
            / "2017_02_23-HD233_HD236-untrimmed_profile_data.nc"
        ).isel(animal=[0, 1, 2])

        std_data, std_warp = pp.standardize_profiles(
            raw_data,
            redox_params=self.redox_params,
            eng=matlab_engine,
            **self.reg_d2_parameters,
        )

        assert raw_data.animal.size == std_data.animal.size
