import logging
import typing
from pathlib import Path

import attr
import numpy as np
import pandas as pd
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

    # Required initialization parameters
    experiment_dir: typing.Union[str, Path]
    imaging_scheme: str

    strains: np.ndarray = None
    experiment_id: str = None

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

    region_means: pd.DataFrame = None
    movement: pd.DataFrame = None

    def __attrs_post_init__(self):
        logging.info(f'Starting full pipeline run for {self.experiment_dir}')
        self.experiment_dir = Path(self.experiment_dir)
        self.experiment_id = self.experiment_dir.stem
        self.load_strains()
        self.scale_region_boundaries()
        self.load_movement_annotation()
        self.full_pipeline()

    def full_pipeline(self):

        self.load_images()
        self.scale_region_boundaries()
        if self.seg_images is None:
            self.seg_images = self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        self.trim_data()
        self.calculate_redox()
        self.calc_region_means()

        logging.info(f'Finished full pipeline run for {self.experiment_dir}')

    def trim_data(self):
        logging.info('Trimming intensity data')
        self.trimmed_intensity_data = ip.trim_profiles(
            self.raw_intensity_data, self.seg_threshold, self.trimmed_profile_length
        )

    def calculate_redox(self):
        logging.info('Calculating redox measurements')

        # Expand the trimmed_intensity_data to include new wavelengths
        new_wvls = np.append(self.trimmed_intensity_data.wavelength.data, ['r', 'oxd', 'e'])
        self.trimmed_intensity_data = self.trimmed_intensity_data.reindex(wavelength=new_wvls)

        self.trimmed_intensity_data.loc[dict(wavelength='r')] = \
            self.trimmed_intensity_data.sel(wavelength='410') / self.trimmed_intensity_data.sel(wavelength='470')
        self.trimmed_intensity_data.loc[dict(wavelength='oxd')] = ip.r_to_oxd(
            self.trimmed_intensity_data.loc[dict(wavelength='r')])
        self.trimmed_intensity_data.loc[dict(wavelength='e')] = ip.oxd_to_redox_potential(
            self.trimmed_intensity_data.loc[dict(wavelength='oxd')])

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
        raw_image_path = Path(self.experiment_dir).joinpath(self.experiment_id + '.tif')

        self.raw_image_data = pio.load_images(raw_image_path, self.imaging_scheme, self.strains)

    def load_movement_annotation(self):
        df = pd.read_csv(self.experiment_dir.joinpath(self.experiment_id + '-mvmt.csv'))
        df = df.pivot_table(index='animal', columns=['region', 'pair'], values='movement')
        df = df.stack('pair')
        self.movement = df

    def load_strains(self):
        self.strains = pio.load_strain_map_from_disk(self.experiment_dir.joinpath(self.experiment_id + '-indexer.csv'))

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

    def calc_region_means(self):
        data = self.trimmed_intensity_data
        dfs = []
        for region, bounds in self.scaled_regions.items():
            region_dfs = []
            for wvl in data.wavelength.data:
                sub_df = data[dict(position=range(bounds[0], bounds[1]))] \
                    .mean(dim='position').sel(wavelength=wvl).to_pandas()
                sub_df = sub_df.reset_index()
                sub_df['animal'] = range(len(sub_df))
                sub_df['region'] = region
                sub_df['experiment'] = self.experiment_id
                sub_df = sub_df.melt(
                    value_vars=[0, 1], var_name='pair',
                    id_vars=['animal', 'strain', 'region', 'experiment'], value_name=wvl)
                region_dfs.append(sub_df)
            df_tmp = pd.concat(region_dfs, axis=1)
            df_tmp = df_tmp.loc[:, ~df_tmp.columns.duplicated()]
            dfs.append(df_tmp)
        df = pd.concat(dfs)
        df.reset_index(drop=True, inplace=True)

        self.region_means = df

        if self.movement is not None:
            self.region_means = self.region_means.join(self.movement, on=['animal', 'pair'])

    def persist_to_disk(self, output_dir: str = None):
        if output_dir is None:
            output_dir = self.experiment_dir
        else:
            output_dir = Path(output_dir)

        logging.info(f'Saving {self.experiment_id} inside {output_dir}')

        pio.save_images_xarray_to_disk(self.seg_images, output_dir.joinpath('processed_images'), self.experiment_id,
                                       'seg')

        # Persist the region means
        self.calc_region_means()
        self.region_means.to_csv(output_dir.joinpath(self.experiment_id + '-region_means.csv'), index=False)


if __name__ == '__main__':
    experiment_path = "/Users/sean/code/wormAnalysis/data/paired_ratio/2017_02_22-HD233_SAY47/"
    ex = PairExperiment(experiment_path, "TL/470/410/470/410")
