import datetime
from pathlib import Path
from typing import List, Dict, Union

import attr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from scipy.interpolate import UnivariateSpline
import logging

from pharynx_analysis import (image_processing as ip, profile_processing, pharynx_io as pio, plots)


@attr.s(auto_attribs=True)
class Experiment:
    experiment_dir: Path
    imaging_scheme: str

    strains: np.ndarray = None
    experiment_id: str = None

    raw_image_data: xr.DataArray = None
    seg_images: xr.DataArray = None

    include_idx = None

    strategy: str = None

    # R -> OxD parameters
    r_min: float = 0.852
    r_max: float = 6.65
    instrument_factor: float = 0.171

    # OxD -> E parameters
    midpoint_potential: float = -265.0
    z: int = 2
    temperature: float = 22.0

    regions: dict = {
        'pm3': [.07, .28],
        'pm4': [.33, .45],
        'pm5': [.53, .70],
        'pm6': [.80, .86],
        'pm7': [.88, .96],
    }
    scaled_regions: dict = None

    # Pipeline Parameters
    trimmed_profile_length: int = 300
    n_midline_pts: int = 200
    seg_threshold: int = 2000
    trim_threshold: int = 3000
    reg_lambda: float = 0.01
    frame_specific_midlines: bool = False
    should_register: bool = False

    rot_fl: xr.DataArray = None
    rot_seg: xr.DataArray = None

    midlines: List[Dict[str, List[UnivariateSpline]]] = None

    untrimmed_profiles: xr.DataArray = None
    trimmed_profiles: xr.DataArray = None

    summary_table: pd.DataFrame = None
    movement: pd.DataFrame = None

    load_from_disk: bool = False

    save_summary_plots: bool = False

    def __attrs_post_init__(self):

        self.experiment_dir = Path(self.experiment_dir)
        self.experiment_id = self.experiment_dir.stem

        self.load_strains()
        self.scale_region_boundaries()
        self.load_movement_annotation()
        if self.include_idx is None:
            self.include_idx = np.ones(len(self.strains))
        if self.load_from_disk:
            self.load_exp_from_disk()

    ####################################################################################################################
    # PIPELINE
    ####################################################################################################################
    def segment_pharynxes(self):
        logging.info('Segmenting pharynxes')
        return ip.segment_pharynxes(self.raw_image_data, self.seg_threshold)

    def align_and_center(self):
        logging.info('Centering and rotating pharynxes')
        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(self.raw_image_data, self.seg_images)

    def calculate_midlines(self):
        logging.info('Calculating midlines')
        # noinspection PyTypeChecker
        self.midlines = ip.calculate_midlines(self.rot_seg, degree=4)

    def measure_under_midlines(self):
        logging.info('Measuring under midlines')
        self.untrimmed_profiles = ip.measure_under_midlines(self.rot_fl, self.midlines, n_points=self.n_midline_pts,
                                                            frame_specific=self.frame_specific_midlines)
        self.untrimmed_profiles = ip.align_pa(self.untrimmed_profiles)
        self.untrimmed_profiles = self.add_experiment_metadata_to_data_array(self.untrimmed_profiles)

    def register(self):
        logging.info('Registering profiles')
        self.untrimmed_profiles = profile_processing.register_profiles(self.untrimmed_profiles)

    def trim_data(self):
        logging.info('Trimming intensity data')

        self.trimmed_profiles = self.add_experiment_metadata_to_data_array(profile_processing.trim_profiles(
            self.untrimmed_profiles, self.trim_threshold, self.trimmed_profile_length
        ))

    def calculate_redox(self):
        logging.info('Calculating redox measurements')

        # Expand the trimmed_intensity_data to include new wavelengths
        new_wvls = np.append(self.trimmed_profiles.wavelength.data, ['r', 'oxd', 'e'])

        self.trimmed_profiles = self.trimmed_profiles.reindex(wavelength=new_wvls)

        self.trimmed_profiles.loc[dict(wavelength='r')] = self.trimmed_profiles.sel(
            wavelength='410') / self.trimmed_profiles.sel(wavelength='470')
        self.trimmed_profiles.loc[dict(wavelength='oxd')] = profile_processing.r_to_oxd(
            self.trimmed_profiles.loc[dict(wavelength='r')])
        self.trimmed_profiles.loc[dict(wavelength='e')] = profile_processing.oxd_to_redox_potential(
            self.trimmed_profiles.loc[dict(wavelength='oxd')])

    def generate_summary_table(self):
        dfs = []
        for region, bounds in self.scaled_regions.items():
            region_dfs = []
            for wvl in self.trimmed_profiles.wavelength.data:
                sub_df = self.trimmed_profiles[dict(position=range(bounds[0], bounds[1]))] \
                    .mean(dim='position').sel(wavelength=wvl).to_pandas()
                sub_df = sub_df.reset_index()
                sub_df['animal'] = range(len(sub_df))
                sub_df['region'] = region
                sub_df['experiment'] = self.experiment_id
                sub_df['strategy'] = self.strategy
                sub_df = sub_df.melt(
                    value_vars=self.trimmed_profiles.pair.data, var_name='pair',
                    id_vars=['animal', 'strain', 'region', 'experiment', 'strategy'], value_name=wvl)
                region_dfs.append(sub_df)
            df_tmp = pd.concat(region_dfs, axis=1)
            df_tmp = df_tmp.loc[:, ~df_tmp.columns.duplicated()]
            dfs.append(df_tmp)
        df = pd.concat(dfs, sort=False)
        df.reset_index(drop=True, inplace=True)

        self.summary_table = pd.concat(dfs)

        if self.movement is not None:
            self.summary_table = self.summary_table.join(self.movement, on=['animal', 'pair'])

        return self.summary_table

    ####################################################################################################################
    # LIVE EDITING
    ####################################################################################################################
    def flip_at(self, idx):
        np.fliplr(self.rot_fl[:, idx])
        np.fliplr(self.rot_seg[:, idx])
        np.fliplr(self.untrimmed_profiles[:, idx])
        np.fliplr(self.trimmed_profiles[:, idx])

    def exclude(self, idx: Union[int, np.ndarray]):
        self.include_idx[idx] = 0

    def unexclude(self, idx: Union[int, np.ndarray]):
        self.include_idx[idx] = 1

    def scale_region_boundaries(self):
        self.scaled_regions = {
            region: np.int_(self.trimmed_profile_length * np.asarray(self.regions[region]))
            for region in self.regions.keys()
        }

    ####################################################################################################################
    # PERSISTENCE/IO
    ####################################################################################################################
    def load_images(self):
        logging.info('Loading Images')
        raw_image_path = Path(self.experiment_dir).joinpath(self.experiment_id + '.tif')

        self.raw_image_data = pio.load_images(raw_image_path, self.imaging_scheme, self.strains)

    def load_strains(self):
        self.strains = pio.load_strain_map_from_disk(self.experiment_dir.joinpath(self.experiment_id + '-indexer.csv'))

    def load_movement_annotation(self):
        try:
            df = pd.read_csv(self.experiment_dir.joinpath(self.experiment_id + '-mvmt.csv'))
            df = df.pivot_table(index='animal', columns=['region', 'pair'], values='movement')
            df = df.stack('pair')
            self.movement = df
        except FileNotFoundError:
            pass

    def get_analysis_dir(self):
        date_str = datetime.datetime.now().strftime('%Y-%m-%d')
        strategy_str = self.strategy.replace(' ', '_')
        analysis_dir = self.experiment_dir.joinpath('analyses', f'{date_str}_{strategy_str}')
        analysis_dir.mkdir(parents=True, exist_ok=True)
        return analysis_dir

    def make_fig_dir(self):
        fig_dir = self.get_analysis_dir().joinpath('figs')
        fig_dir.mkdir(parents=True, exist_ok=True)
        return fig_dir

    def persist_profile_data(self):
        analysis_dir = self.get_analysis_dir()

        profile_data_filename = analysis_dir.joinpath(self.experiment_id + '-profile_data.nc')
        logging.info(f'Saving profile data to {profile_data_filename}')
        self.trimmed_profiles.to_netcdf(profile_data_filename)

    def load_profile_data(self):
        analysis_dir = self.get_analysis_dir()
        profile_data_raw_filename = analysis_dir.joinpath(self.experiment_id + '-profile_data.nc')
        logging.info(f'Loading profile data from {profile_data_raw_filename}')
        self.trimmed_profiles = xr.open_dataarray(profile_data_raw_filename)

    def load_summary_table(self):
        analysis_dir = self.get_analysis_dir()
        self.summary_table = pd.read_csv(analysis_dir.joinpath(self.experiment_id + '-summary_table.csv'))

    def save_profile_summary_plots(self, prof_data: xr.DataArray):
        fig_dir = self.make_fig_dir()
        for wvl in prof_data.wavelength.data:
            for pair in prof_data.pair.data:
                fig, ax = plots.plot_profile_avg_by_strain(prof_data.sel(wavelength=wvl, pair=pair),
                                                           ax_title=f'{self.experiment_id}-{wvl}')
                plots.add_regions_to_axis(ax, self.scaled_regions)
                fig.savefig(fig_dir.joinpath(f'{self.experiment_id}-{wvl}-{pair}.pdf'))

    def save_cat_plots(self):
        fig_dir = self.make_fig_dir()
        cat_plot_dir = fig_dir.joinpath('statistical')
        cat_plot_dir.mkdir(parents=True, exist_ok=True)

        for wvl in self.trimmed_profiles.wavelength.data:
            for pair in self.trimmed_profiles.pair.data:
                plt.clf()
                sns.violinplot(
                    y=wvl,
                    x='strain',
                    data=self.summary_table[self.summary_table.pair == pair],
                    inner='quartiles'
                )
                plt.savefig(cat_plot_dir.joinpath(f'violin-{wvl}-{pair}.pdf'))

                plt.clf()
                sns.catplot(
                    y=wvl,
                    x='strain',
                    col='region',
                    col_wrap=2,
                    kind='violin',
                    data=self.summary_table[self.summary_table.pair == pair]
                )
                plt.savefig(cat_plot_dir.joinpath(f'violin_regions-{wvl}-{pair}.pdf'))

    def persist_to_disk(self, summary_plots=False):
        logging.info(f'Saving {self.experiment_id} inside {self.experiment_dir}')

        analysis_dir = self.get_analysis_dir()

        # Persist the region means
        summary_table_filename = analysis_dir.joinpath(self.experiment_id + '-summary_table.csv')
        logging.info(f'Saving region means to {summary_table_filename}')
        self.generate_summary_table()
        self.summary_table.to_csv(summary_table_filename, index=False)

        # Persist the profile data
        self.persist_profile_data()

        # Plots
        if summary_plots:
            self.save_profile_summary_plots(self.trimmed_profiles)
            self.save_cat_plots()

    def load_exp_from_disk(self):
        pass

    ####################################################################################################################
    # MISC / HELPER
    ####################################################################################################################
    def add_experiment_metadata_to_data_array(self, data_array):
        return data_array.assign_attrs(
            r_min=self.r_min,
            r_max=self.r_max,
            instrument_factor=self.instrument_factor,
            midpoint_potential=self.midpoint_potential,
            z=self.z,
            temperature=self.temperature,
            strategy=self.strategy
        )

    def filter_by_exclude_status(self, data):
        return data.loc[dict(strain=np.logical_not(self.include_idx))]


@attr.s(auto_attribs=True)
class PairExperiment(Experiment):
    """
    This is the paired ratio experiment
    """

    strategy: str = "frame specific midlines with registration"

    # Required initialization parameters
    image_display_order: List[str] = [
        '410_1', '470_1', 'r1',
        '410_2', '470_2', 'r2'
    ]

    def full_pipeline(self):
        logging.info(f'Starting full pipeline run for {self.experiment_dir}')
        self.load_images()
        self.scale_region_boundaries()
        if self.seg_images is None:
            self.seg_images = self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        if self.should_register:
            self.register()
        self.trim_data()
        self.calculate_redox()
        self.generate_summary_table()
        self.persist_to_disk(summary_plots=self.save_summary_plots)

        logging.info(f'Finished full pipeline run for {self.experiment_dir}')


@attr.s(auto_attribs=True)
class CataExperiment(Experiment):
    strategy: str = "cata"
    frame_specific_midlines: bool = False

    def full_pipeline(self):
        logging.info(f'Starting full pipeline run for {self.experiment_dir}')
        self.load_images()
        self.scale_region_boundaries()
        if self.seg_images is None:
            self.seg_images = self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        self.trim_data()
        self.calculate_redox()
        self.generate_summary_table()
        self.persist_to_disk(summary_plots=self.save_summary_plots)

        logging.info(f'Finished full Cata pipeline run for {self.experiment_dir}')

    def segment_pharynxes(self):
        ref_wvl = '410'
        seg_images = super().segment_pharynxes()

        # In Cata's pipeline, the fluorescent images are also masked
        for animal_idx in np.arange(seg_images.strain.size):
            for wvl_idx in np.arange(seg_images.wavelength.size):
                for pair in seg_images.pair:
                    ref_seg = seg_images.sel(wavelength=ref_wvl, pair=pair).isel(strain=animal_idx)
                    img = self.raw_image_data.isel(strain=animal_idx, wavelength=wvl_idx, pair=pair)

                    masked = ref_seg * img

                    self.raw_image_data[animal_idx, wvl_idx, pair] = masked

        return seg_images

    def measure_under_midlines(self, ref_wvl='410'):
        logging.info('Measuring under midlines')
        self.untrimmed_profiles = ip.measure_under_midlines(
            self.rot_fl, self.midlines, n_points=self.n_midline_pts, frame_specific=self.frame_specific_midlines)
        self.untrimmed_profiles = ip.align_pa(self.untrimmed_profiles)
        self.untrimmed_profiles = self.add_experiment_metadata_to_data_array(self.untrimmed_profiles)


if __name__ == '__main__':
    import logging

    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', level=logging.DEBUG, datefmt='%I:%M:%S')
    experiment_path = Path("/Users/sean/code/wormAnalysis/data/paired_ratio/2017_02_22-HD233_SAY47/")
    ex = PairExperiment(experiment_path, "TL/470/410/470/410")
    ex.full_pipeline()
