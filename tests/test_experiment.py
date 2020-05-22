import os
import pytest
from unittest.mock import patch
import xarray as xr
import numpy as np
from pharedox import io, experiment, image_processing as ip


class TestExperiment:
    @pytest.fixture(scope="function")
    def paired_imgs(self, shared_datadir):
        return io.load_tiff_as_hyperstack(
            (
                shared_datadir
                / "experiments"
                / "2017_02_22-HD233_SAY47"
                / "2017_02_22-HD233_SAY47.tif"
            ),
            indexer_path=(
                shared_datadir
                / "experiments"
                / "2017_02_22-HD233_SAY47"
                / "2017_02_22-HD233_SAY47-indexer.csv"
            ),
            channel_order=["TL", "470", "410", "470", "410"],
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
