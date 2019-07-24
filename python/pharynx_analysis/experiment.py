import logging
import typing
from pathlib import Path

import attr
import numpy as np
import xarray as xr

from pharynx_analysis import image_processing as ip
from pharynx_analysis import pharynx_io as pio



def get_non_tl(data_array):
    return data_array.where(data_array.wavelength != 'TL', drop=True)


@attr.s(auto_attribs=True)
class PairExperiment:
    """
    This is the paired ratio experiment

    Attributes
    ----------
    raw_image_path
    imaging_scheme
    strain_map
    midline_smoothing
    """
    raw_image_path: str
    imaging_scheme: str
    strain_map: typing.List[str]

    # Pipeline Parameters
    midline_smoothing: int = 1e8
    trimmed_profile_length: int = 100
    n_midline_pts: int = 200
    seg_threshold: int = 2000

    # R -> OxD parameters
    r_min: float = 0.852
    r_max: float = 6.65
    instrument_factor: float = 0.171

    # OxD -> E parameters
    midpoint_potential: float = -265.0
    z: int = 2
    temperature: float = 22.0

    image_display_order: typing.List[str] = [
        '410_1', '470_1', 'r1',
        '410_2', '470_2', 'r2'
    ]

    regions: dict = {
        'pm3': [.07, .28],
        'pm4': [.33, .45],
        'pm5': [.53, .69],
        'pm6': [.78, .83],
        'pm7': [.85, .925],
    }

    scaled_regions: dict = None

    raw_image_data: xr.DataArray = None
    seg_images: xr.DataArray = None

    rot_fl: xr.DataArray = None
    rot_seg: xr.DataArray = None

    midlines: typing.List[dict] = None

    raw_intensity_data: xr.DataArray = None
    trimmed_intensity_data: xr.DataArray = None

    r: xr.DataArray = None
    oxd: xr.DataArray = None
    e: xr.DataArray = None

    def __attrs_post_init__(self):
        self.scale_region_boundaries()
        self.full_pipeline()

    def full_pipeline(self):
        logging.info(f'Starting full pipeline run for {self.raw_image_path}')

        self.load_images()
        self.scale_region_boundaries()
        if self.seg_images is None:
            self.seg_images = self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        self.trim_data()
        self.calculate_redox()

        logging.info(f'Finished full pipeline run for {self.raw_image_path}')

    def trim_data(self):
        logging.info('Trimming intensity data')
        self.trimmed_intensity_data = ip.trim_profiles(
            self.raw_intensity_data, self.seg_threshold, self.trimmed_profile_length
        )

    def calculate_redox(self):
        logging.info('Calculating redox measurements')
        self.r = self.trimmed_intensity_data.sel(wavelength='410') / self.trimmed_intensity_data.sel(wavelength='470')
        self.oxd = ip.r_to_oxd(self.r)
        self.e = ip.oxd_to_redox_potential(self.oxd)

    def align_and_center(self):
        logging.info('Centering and rotating pharynxes')
        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(self.raw_image_data, self.seg_images)

    def calculate_midlines(self):
        logging.info('Calculating midlines')
        self.midlines = ip.calculate_midlines(self.rot_seg, self.rot_fl)

    def measure_under_midlines(self):
        logging.info('Measuring under midlines')
        step = self.raw_image_data.x.size // 4
        self.midline_xs = np.linspace(step, self.raw_image_data.x.size - step, self.n_midline_pts)
        self.raw_intensity_data = ip.measure_under_midlines(
            self.rot_fl, self.midlines, (step, self.raw_image_data.x.size - step), n_points=self.n_midline_pts
        )
        self.raw_intensity_data.data = np.nan_to_num(self.raw_intensity_data.data)
        self.raw_intensity_data = ip.align_pa(self.raw_intensity_data)

    def load_images(self):
        logging.info('Loading Images')
        self.raw_image_data = pio.load_images(self.raw_image_path, self.imaging_scheme, self.strain_map)

    def segment_pharynxes(self):
        logging.info('Segmenting pharynxes')
        return ip.segment_pharynxes(self.raw_image_data, self.seg_threshold)

    def flip_at(self, idx):
        np.fliplr(self.rot_fl[:, idx])
        np.fliplr(self.rot_seg[:, idx])
        np.fliplr(self.raw_intensity_data[:, idx])

    def scale_region_boundaries(self):
        self.scaled_regions = {
            region: np.int_(self.trimmed_profile_length * np.asarray(self.regions[region]))
            for region in self.regions.keys()
        }

    def persist_to_disk(self, parent_dir=None):
        # TODO
        raw_img_path = Path(self.raw_image_path)
        pio.save_images_xarray_to_disk(
            self.seg_images, raw_img_path.parent.joinpath('processed_images'), raw_img_path.stem, 'seg')


if __name__ == '__main__':
    img_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif"
    strain_map_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/indexer.csv"
    strains = pio.load_strain_map_from_disk(strain_map_path)
    ex = PairExperiment(img_path, "TL/470/410/470/410", strains)
