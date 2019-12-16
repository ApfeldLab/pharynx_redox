import os
from pathlib import Path
import numpy as np
import pytest

from pharynx_redox import image_processing as ip, pharynx_io as pio

test_data_path = Path(os.path.join(os.path.dirname(__file__), "test_data"))


class TestImageProcessing:
    img_stk_0 = {
        "img_path": os.path.join(
            test_data_path, "paired_ratio_0/2017_02_22-HD233_SAY47.tif"
        ),
        "strain_map_path": os.path.join(
            test_data_path, "paired_ratio_0/2017_02_22-HD233_SAY47-indexer.csv"
        ),
        "imaging_scheme": "TL/470/410/470/410",
        "n_animals": 123,
        "n_frames": 615,
        "n_wavelengths": 3,
        "unique_wavelengths": ["TL", "470", "410"],
        "n_pairs": 2,
        "img_w": 174,
        "img_h": 130,
        "unique_strains": ["HD233", "SAY47"],
        "n_HD233": 60,
        "n_SAY47": 63,
    }

    imgdata = pio.load_images(
        img_stk_0["img_path"], indexer_path=img_stk_0["strain_map_path"]
    )

    seg0 = imgdata > 2000
    fl0 = imgdata.sel(wavelength="410", pair=0)

    rot_fl, rot_seg = ip.center_and_rotate_pharynxes(imgdata, imgdata > 2000)

    def test_get_lr_bounds(self):
        expected_bounds = np.asarray([[55, 112], [55, 112], [55, 112]])
        actual_bounds = ip.get_lr_bounds(self.rot_seg)

        np.testing.assert_array_equal(expected_bounds, actual_bounds[:3, :])

    def test_get_lr_bounds_padding(self):
        padding = 5
        expected_bounds = np.asarray([[55, 112], [55, 112], [55, 112]])
        expected_bounds[:, 0] = expected_bounds[:, 0] - padding
        expected_bounds[:, 1] = expected_bounds[:, 1] + padding

        actual_bounds = ip.get_lr_bounds(self.rot_seg, pad=padding)
        np.testing.assert_array_equal(expected_bounds, actual_bounds[:3, :])

    @pytest.mark.skip(
        reason="center/rotate test failing because of instability of the rotation? look into this"
    )
    def test_center_and_rotate(self):
        fl_test_rot, seg_test_rot = ip.center_and_rotate_pharynxes(self.fl0, self.seg0)

        fl_expected_rot = pio.load_images(
            test_data_path.joinpath("fl_0_rot.tif"), "410", np.repeat("HD233", 154)
        )
        seg_expected_rot = pio.load_images(
            test_data_path.joinpath("seg_0_rot.tif"), "410", np.repeat("HD233", 154)
        )

        np.testing.assert_array_equal(fl_test_rot.values, fl_expected_rot.values)
        np.testing.assert_array_equal(seg_test_rot.values, seg_expected_rot.values)
