import numpy as np
import pytest
import xarray as xr

from pharedox import pio


class TestIO:
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

    def test_load_images_returns_dataarray(self, paired_imgs):
        assert type(paired_imgs) == xr.DataArray

    def test_load_images_dtype_float(self, paired_imgs):
        assert paired_imgs.dtype == np.float64

    def test_load_images_dimensions(self, paired_imgs):
        assert paired_imgs.animal.size == 123
        assert paired_imgs.timepoint.size == 1
        assert paired_imgs.pair.size == 2
        assert paired_imgs.y.size == 130
        assert paired_imgs.x.size == 174

    def test_load_images_processed_by_imagej(self, shared_datadir):
        paired_imgs = pio.load_tiff_as_hyperstack(
            (
                shared_datadir
                / "misc_data_files"
                / "processed_2017_02_22-HD233_SAY47.tif"
            ),
            manual_metadata=(
                shared_datadir
                / "misc_data_files"
                / "2017_02_22-HD233_SAY47-frame_map.csv"
            ),
        )

        assert type(paired_imgs) == xr.DataArray

    def test_save_load_from_disk_identical(self, paired_imgs, shared_datadir):
        path = shared_datadir / "paired_imgs.nc"

        pio.save_profile_data(paired_imgs, path)
        imgs_from_disk = pio.load_profile_data(path)

        # For now, we drop the time coordinates because xarray loses millisecond precision
        # when round-tripping a DataArray to/from disk
        # https://github.com/pydata/xarray/issues/4045

        # 2020-07-24: I removed adding the time coordinate to the data while messing
        # around with metadata, so I commended this assertion out
        # assert paired_imgs.reset_coords("time", drop=True).equals(
        #     imgs_from_disk.reset_coords("time", drop=True)
        # )

        assert paired_imgs.equals(imgs_from_disk)

    def test_extract_metadata_metamorph_acquire(self, shared_datadir):
        """Tests that a file straight from metamorph successfully loads metadata"""
        md = pio.get_image_metadata_metamorph_acquire(
            shared_datadir
            / "experiments"
            / "2017_02_22-HD233_SAY47"
            / "2017_02_22-HD233_SAY47.tif"
        )

        assert len(md) == 615

    def test_get_metadata_from_tiff_on_file_no_metadata_raises(self, shared_datadir):
        """Tests that files without metadata (e.g. after having been processed by FIJI) raises an AttributeError"""
        with pytest.raises(AttributeError):
            pio.get_image_metadata(
                shared_datadir
                / "misc_data_files"
                / "processed_2017_02_22-HD233_SAY47.tif",
                "acquire",
            )

    def test_save_images_xarray_to_disk(self, paired_imgs, shared_datadir):
        pio.save_images_xarray_to_tiffs(
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
        pio.save_images_xarray_to_tiffs(
            paired_imgs, shared_datadir / "test_save_imgs_dir", z_axis="pair"
        )

        files = list((shared_datadir / "test_save_imgs_dir").glob("*.tiff"))
        assert (
            len(files)
            == paired_imgs.animal.size
            * paired_imgs.timepoint.size
            * paired_imgs.wavelength.size
        )

    def test_load_tiff_as_hyperstack_attaches_movement_correctly(self, paired_imgs):
        # TODO
        pass
