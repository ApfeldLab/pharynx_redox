import logging

import pytest

from pharedox import pio
from pharedox import profile_processing as pp


class TestImageProcessing:
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

    def test_measure_under_labels(self):
        pass

    def test_subtract_medians(self):
        pass
