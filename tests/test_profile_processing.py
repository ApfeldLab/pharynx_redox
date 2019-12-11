import pytest
import os
from pathlib import Path
import numpy as np

from pharynx_redox import profile_processing as pp, pharynx_io as pio

test_data_path = Path(__file__).resolve().parent.joinpath("test_data")


class TestProfileProcessing:
    unreg_profile_data = pio.load_profile_data(
        test_data_path.joinpath("2017_02_22-HD233_SAY47-untrimmed_profile_data.nc")
    )

    def test_registration_single_profile_single_pair(self):
        data_to_register = self.unreg_profile_data.sel(pair=0)[0]
        registered_data = pp.register_profiles(data_to_register)

        assert data_to_register.shape == registered_data.shape
        assert data_to_register.coords == registered_data.coords

    def test_registeration_multiple_profiles(self):
        pytest.fail()

    def test_registration_multiple_pairs(self):
        pytest.fail()
