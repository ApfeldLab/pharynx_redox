import os
from unittest.mock import patch

import numpy as np
import pytest
import xarray as xr

from pharedox import experiment
from pharedox import image_processing as ip
from pharedox import pio


class TestExperiment:
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

    @pytest.mark.slow
    def test_full_experiment_run_paired_single_timepoint(self, shared_datadir):
        exp = experiment.Experiment(
            shared_datadir / "experiments" / "2017_02_22-HD233_SAY47"
        )
        exp.full_pipeline()

    @patch.object(ip, "segment_pharynxes")
    def test_experiment_loads_masks_if_they_exist(
        self, mock_segment_pharynxes, shared_datadir
    ):
        exp = experiment.Experiment(
            shared_datadir / "experiments" / "2017_02_23-HD233_HD236"
        )
        exp.segment_pharynxes()
        mock_segment_pharynxes.assert_not_called()

    @patch.object(ip, "segment_pharynxes")
    def test_experiment_doesnt_loads_masks_if_they_dont_exist(
        self, mock_segment_pharynxes, shared_datadir
    ):
        exp = experiment.Experiment(
            shared_datadir / "experiments" / "2017_02_22-HD233_SAY47"
        )
        exp.segment_pharynxes()
        mock_segment_pharynxes.assert_called()
