import os
from pathlib import Path
import numpy as np
import pytest

from pharynx_redox import image_processing as ip, io as pio, data_analysis as da

test_data_path = Path(os.path.join(os.path.dirname(__file__), "test_data"))

untrimmed_profile_data_path = test_data_path.joinpath(
    "2017_02_22-HD233_SAY47-untrimmed_profile_data-new.nc"
)
untrimmed_profile_data = pio.load_profile_data(untrimmed_profile_data_path)


def test_moving_idx_one_region_moving_correct():
    mv, st = da.get_moving_idx(untrimmed_profile_data, "posterior")

    assert mv[1] == 1
    assert st[1] == 0

    assert mv[10] == 0
    assert st[10] == 0


def test_moving_idx_one_region_stationary_correct():
    mv, st = da.get_moving_idx(untrimmed_profile_data, "posterior")

    assert mv[0] == 0
    assert st[0] == 1


def test_moving_idx_one_region_neither_correct():
    mv, st = da.get_moving_idx(untrimmed_profile_data, "posterior")

    assert mv[10] == 0
    assert st[10] == 0


def test_moving_idx_one_region_correct_length():
    mv, st = da.get_moving_idx(untrimmed_profile_data, "posterior")

    assert len(mv) == untrimmed_profile_data.animal.size
    assert len(st) == untrimmed_profile_data.animal.size


def test_moving_idx_multiple_regions_stationary_correct():
    regions = ["posterior", "anterior"]
    mv, st = da.get_moving_idx(untrimmed_profile_data, regions)

    assert mv[0] == 0
    assert st[0] == 1


def test_moving_idx_multiple_regions_neither_correct():
    # "Unacceptable" in the anterior
    regions = "anterior"
    mv, st = da.get_moving_idx(untrimmed_profile_data, regions)

    assert mv[5] == 0
    assert st[5] == 0

    # Stationary in the posterior
    regions = "posterior"
    mv, st = da.get_moving_idx(untrimmed_profile_data, regions)

    assert mv[5] == 0
    assert st[5] == 1

    # "Unacceptable" when considering both regions
    regions = ["anterior", "posterior"]
    mv, st = da.get_moving_idx(untrimmed_profile_data, regions)

    assert mv[5] == 0
    assert st[5] == 0


def test_resample_moving_correct_percent_moving():
    percent_moving = 0.5
    resampled = da.resample_moving(
        untrimmed_profile_data, percent_moving, 3000, "posterior"
    )

    mv, st = da.get_moving_idx(resampled, "posterior")

    np.testing.assert_almost_equal(
        (np.sum(mv) / len(resampled)), percent_moving, decimal=1
    )


def test_resample_moving_correct_percent_moving():
    percent_moving = 1.0
    resampled = da.resample_moving(
        untrimmed_profile_data, percent_moving, 3000, "posterior"
    )

    mv, st = da.get_moving_idx(resampled, "posterior")

    np.testing.assert_almost_equal(
        (np.sum(mv) / len(resampled)), percent_moving, decimal=1
    )
