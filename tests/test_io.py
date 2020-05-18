import os
import pytest
from unittest.mock import MagicMock
from unittest.mock import patch
import xarray as xr
import numpy as np
from pharedox import io


class TestIO:
    @pytest.fixture(scope="function")
    def paired_imgs(self, shared_datadir):
        return io.load_images(
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

    def test_load_images_returns_dataarray(self, paired_imgs):
        assert type(paired_imgs) == xr.DataArray

    def test_load_images_dimensions(self, paired_imgs):
        assert paired_imgs.animal.size == 123
        assert paired_imgs.timepoint.size == 1
        assert paired_imgs.pair.size == 2
        assert paired_imgs.y.size == 130
        assert paired_imgs.x.size == 174

        # datatype
        assert paired_imgs.dtype == np.uint16

    def test_load_images_no_metadata(self, shared_datadir):
        paired_imgs = io.load_images(
            (shared_datadir / "misc_imgs" / "processed_2017_02_22-HD233_SAY47.tif"),
            indexer_path=(
                shared_datadir
                / "experiments"
                / "2017_02_22-HD233_SAY47"
                / "2017_02_22-HD233_SAY47-indexer.csv"
            ),
            channel_order=["TL", "470", "410", "470", "410"],
        )

        assert type(paired_imgs) == xr.DataArray

        # Dimensions
        assert paired_imgs.animal.size == 123
        assert paired_imgs.timepoint.size == 1
        assert paired_imgs.pair.size == 2
        assert paired_imgs.y.size == 130
        assert paired_imgs.x.size == 174

        # datatype
        assert paired_imgs.dtype == np.uint16

    def test_save_load_from_disk_identical(self, paired_imgs, shared_datadir):
        path = shared_datadir / "paired_imgs.nc"

        io.save_profile_data(paired_imgs, path)
        imgs_from_disk = io.load_profile_data(path)

        # For now, we drop the time coordinates because xarray loses millisecond precision
        # when round-tripping a DataArray to/from disk
        # https://github.com/pydata/xarray/issues/4045
        assert paired_imgs.reset_coords("time", drop=True).equals(
            imgs_from_disk.reset_coords("time", drop=True)
        )

    def test_parse_illum_setting(self):
        assert io._parse_illum_setting("pH Sensitive GFP 470") == "470"
        assert io._parse_illum_setting("pH Sensitive GFP 410") == "410"
        assert io._parse_illum_setting("Transmitted Light") == "TL"

    def test_get_metadata_from_tiff_on_metamorph_file(self, shared_datadir):
        """Tests that a file straight from metamorph successfully loads metadata"""
        md = io.get_metadata_from_tiff(
            shared_datadir
            / "experiments"
            / "2017_02_22-HD233_SAY47"
            / "2017_02_22-HD233_SAY47.tif"
        )

        assert len(md["data"]) == 615

    def test_get_metadata_from_tiff_on_file_no_metadata_raises(self, shared_datadir):
        """Tests that files without metadata (e.g. after having been processed by FIJI) raises an AttributeError"""
        with pytest.raises(AttributeError):
            io.get_metadata_from_tiff(
                shared_datadir / "misc_imgs" / "processed_2017_02_22-HD233_SAY47.tif"
            )

    def test_load_tiff_from_disk_returns_array_and_metadata_for_metamorph_file(
        self, shared_datadir
    ):
        img_data, metadata = io.load_tiff_from_disk(
            shared_datadir
            / "experiments"
            / "2017_02_22-HD233_SAY47"
            / "2017_02_22-HD233_SAY47.tif"
        )

        assert type(img_data) == np.ndarray
        assert type(metadata) == dict

    def test_load_tiff_from_disk_returns_array_and_none_for_non_metamorph(
        self, shared_datadir
    ):
        """Files without metadata (e.g. after having been processed by FIJI) should contain `None` as the metadata"""
        img_data, metadata = io.load_tiff_from_disk(
            shared_datadir / "misc_imgs" / "processed_2017_02_22-HD233_SAY47.tif"
        )

        assert type(img_data) == np.ndarray
        assert metadata == None

    def test_load_strain_map_from_disk(self, shared_datadir):

        strain_map = io.load_strain_map_from_disk(
            shared_datadir
            / "experiments"
            / "2017_02_22-HD233_SAY47"
            / "2017_02_22-HD233_SAY47-indexer.csv"
        )

        assert len(strain_map) == 123

    def test_save_images_xarray_to_disk(self, paired_imgs, shared_datadir):
        io.save_images_xarray_to_disk(
            paired_imgs, shared_datadir / "test_save_imgs_dir", z_axis="animal"
        )

        files = list((shared_datadir / "test_save_imgs_dir").glob("*.tiff"))
        assert (
            len(files)
            == paired_imgs.timepoint.size
            * paired_imgs.pair.size
            * paired_imgs.wavelength.size
        )

    def test_save_images_xarray_to_disk_pair_zaxis(self, paired_imgs, shared_datadir):
        io.save_images_xarray_to_disk(
            paired_imgs, shared_datadir / "test_save_imgs_dir", z_axis="pair"
        )

        files = list((shared_datadir / "test_save_imgs_dir").glob("*.tiff"))
        assert (
            len(files)
            == paired_imgs.animal.size
            * paired_imgs.timepoint.size
            * paired_imgs.wavelength.size
        )
