import os
import numpy as np
from pathlib import Path

from pharynx_redox import pharynx_io as pio

test_data_path = Path(os.path.join(os.path.dirname(__file__), "test_data"))


class TestPharynxIO:
    img_stk_0 = {
        "img_path": test_data_path.joinpath(
            "paired_ratio_0/2017_02_22-HD233_SAY47.tif"
        ),
        "strain_map_path": test_data_path.joinpath(
            "paired_ratio_0/2017_02_22-HD233_SAY47-indexer.csv"
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

    def test_save_load_profile_data(self):
        # TODO
        pass

    def test_parse_illum_setting_numerical(self):
        assert pio._parse_illum_setting("pH Sensitive GFP 470") == "470"
        assert pio._parse_illum_setting("pH Sensitive GFP 410") == "410"

    def test_parse_illum_setting_transmitted_light(self):
        assert pio._parse_illum_setting("Transmitted Light") == "TL"
        assert pio._parse_illum_setting("a Transmitted Light") == "TL"

    def test_parse_illum_setting_unkown(self):
        assert pio._parse_illum_setting("something else") == "something else"

    def test_get_metadata_from_tiff_correct_length(self):
        img_stack = pio.load_tiff_from_disk(self.img_stk_0["img_path"])
        metadata = pio.get_metadata_from_tiff(self.img_stk_0["img_path"])

        assert len(metadata) == img_stack.shape[0]

    def test_load_images_loads_metadata_if_required(self):
        img_stack, img_metadata = pio.load_tiff_from_disk(
            self.img_stk_0["img_path"], return_metadata=True
        )

        # now call it without asking for metadata
        _ = pio.load_tiff_from_disk(self.img_stk_0["img_path"])

        assert len(img_metadata) == img_stack.shape[0]

    def test_load_tiff_from_disk_shape(self):
        img_stack = pio.load_tiff_from_disk(self.img_stk_0["img_path"])

        assert img_stack.shape[0] == self.img_stk_0["n_frames"]
        assert img_stack.shape[1] == self.img_stk_0["img_h"]
        assert img_stack.shape[2] == self.img_stk_0["img_w"]

    def test_load_tiff_from_disk_dtype(self):
        img_stack = pio.load_tiff_from_disk(self.img_stk_0["img_path"])
        assert img_stack.dtype == "uint16"

    def test_load_images_shape(self):
        img_stack = pio.load_images(
            self.img_stk_0["img_path"],
            pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"]),
        )

        assert img_stack.shape[0] == self.img_stk_0["n_animals"]
        assert img_stack.shape[1] == self.img_stk_0["n_pairs"]
        assert img_stack.shape[2] == self.img_stk_0["n_wavelengths"]
        assert img_stack.shape[3] == self.img_stk_0["img_h"]
        assert img_stack.shape[4] == self.img_stk_0["img_w"]

    def test_load_images_dimension_ordering(self):
        img_stack = pio.load_images(
            self.img_stk_0["img_path"],
            pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"]),
        )

        assert img_stack.dims == ("animal", "pair", "wavelength", "y", "x")

    def test_load_strain_map_shape(self):
        strains = pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"])

        assert (self.img_stk_0["n_animals"],) == strains.shape

    def test_load_strain_map_strains(self):
        strains = pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"])

        assert sorted(np.unique(strains)) == sorted(self.img_stk_0["unique_strains"])

    def test_load_strain_map_length(self):
        strains = pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"])

        assert len(strains) == self.img_stk_0["n_animals"]

    def test_load_strain_map_n_strains(self):
        strains = pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"])

        assert len(strains[strains == "HD233"]) == self.img_stk_0["n_HD233"]
        assert len(strains[strains == "SAY47"]) == self.img_stk_0["n_SAY47"]
