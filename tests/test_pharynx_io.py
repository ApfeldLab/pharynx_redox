import os
import numpy as np

from pharynx_redox import pharynx_io as pio

test_data_path = os.path.join(os.path.dirname(__file__), "test_data")
print(test_data_path)


class TestPharynxAnalysis:
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
            self.img_stk_0["imaging_scheme"],
            pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"]),
        )

        assert img_stack.shape[0] == self.img_stk_0["n_animals"]
        assert img_stack.shape[1] == self.img_stk_0["n_wavelengths"]
        assert img_stack.shape[2] == self.img_stk_0["n_pairs"]
        assert img_stack.shape[3] == self.img_stk_0["img_h"]
        assert img_stack.shape[4] == self.img_stk_0["img_w"]

    def test_load_images_dimension_ordering(self):
        img_stack = pio.load_images(
            self.img_stk_0["img_path"],
            self.img_stk_0["imaging_scheme"],
            pio.load_strain_map_from_disk(self.img_stk_0["strain_map_path"]),
        )

        assert img_stack.dims == ("strain", "wavelength", "pair", "y", "x")

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