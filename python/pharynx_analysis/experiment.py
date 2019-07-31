import logging
from pathlib import Path
from typing import List, Dict, Union

import attr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tqdm
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from scipy.interpolate import UnivariateSpline

from pharynx_analysis import (image_processing as ip, profile_processing, pharynx_io as pio, plots)


def get_non_tl(data_array):
    return data_array.where(data_array.wavelength != 'TL', drop=True)


@attr.s(auto_attribs=True)
class Experiment:
    experiment_dir: Path
    imaging_scheme: str

    strains: np.ndarray = None
    experiment_id: str = None

    raw_image_data: xr.DataArray = None
    seg_images: xr.DataArray = None

    include_idx = None

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
    midline_smoothing: int = 1e8
    trimmed_profile_length: int = 300
    n_midline_pts: int = 200
    seg_threshold: int = 2000
    trim_threshold: int = 4000
    reg_lambda = 0.05

    rot_fl: xr.DataArray = None
    rot_seg: xr.DataArray = None

    midlines: List[Dict[str, List[UnivariateSpline]]] = None

    raw_profiles: xr.DataArray = None
    reg_profiles: xr.DataArray = None

    trimmed_raw_profiles: xr.DataArray = None
    trimmed_reg_profiles: xr.DataArray = None

    summary_table: pd.DataFrame = None
    movement: pd.DataFrame = None

    def load_images(self):
        logging.info('Loading Images')
        raw_image_path = Path(self.experiment_dir).joinpath(self.experiment_id + '.tif')

        self.raw_image_data = pio.load_images(raw_image_path, self.imaging_scheme, self.strains)

    def add_experiment_metadata_to_data_array(self, data_array):
        return data_array.assign_attrs(
            r_min=self.r_min,
            r_max=self.r_max,
            instrument_factor=self.instrument_factor,
            midpoint_potential=self.midpoint_potential,
            z=self.z,
            temperature=self.temperature
        )

    def flip_at(self, idx):
        np.fliplr(self.rot_fl[:, idx])
        np.fliplr(self.rot_seg[:, idx])
        np.fliplr(self.raw_profiles[:, idx])
        np.fliplr(self.reg_profiles[:, idx])
        np.fliplr(self.trimmed_raw_profiles[:, idx])
        np.fliplr(self.trimmed_reg_profiles[:, idx])

    def scale_region_boundaries(self):
        self.scaled_regions = {
            region: np.int_(self.trimmed_profile_length * np.asarray(self.regions[region]))
            for region in self.regions.keys()
        }

    def make_fig_dir(self):
        fig_dir = self.experiment_dir.joinpath('figs')
        fig_dir.mkdir(parents=True, exist_ok=True)
        return fig_dir

    def exclude(self, idx: Union[int, np.ndarray]):
        self.include_idx[idx] = 0

    def unexclude(self, idx: Union[int, np.ndarray]):
        self.include_idx[idx] = 1

    def register(self):
        logging.info('Registering profiles')
        self.reg_profiles = profile_processing.register_profiles(self.raw_profiles, lam=self.reg_lambda)

    def trim_data(self):
        logging.info('Trimming intensity data')

        self.trimmed_raw_profiles = self.add_experiment_metadata_to_data_array(ip.trim_profiles(
            self.raw_profiles, self.trim_threshold, self.trimmed_profile_length
        ))
        self.trimmed_reg_profiles = self.add_experiment_metadata_to_data_array(ip.trim_profiles(
            self.reg_profiles, self.trim_threshold, self.trimmed_profile_length
        ))

    def calculate_redox(self):
        logging.info('Calculating redox measurements')

        # Expand the trimmed_intensity_data to include new wavelengths
        new_wvls = np.append(self.trimmed_raw_profiles.wavelength.data, ['r', 'oxd', 'e'])
        self.trimmed_raw_profiles = self.trimmed_raw_profiles.reindex(wavelength=new_wvls)
        self.trimmed_reg_profiles = self.trimmed_reg_profiles.reindex(wavelength=new_wvls)

        self.trimmed_raw_profiles.loc[dict(wavelength='r')] = \
            self.trimmed_raw_profiles.sel(wavelength='410') / self.trimmed_raw_profiles.sel(wavelength='470')
        self.trimmed_raw_profiles.loc[dict(wavelength='oxd')] = ip.r_to_oxd(
            self.trimmed_raw_profiles.loc[dict(wavelength='r')])
        self.trimmed_raw_profiles.loc[dict(wavelength='e')] = ip.oxd_to_redox_potential(
            self.trimmed_raw_profiles.loc[dict(wavelength='oxd')])

        self.trimmed_reg_profiles.loc[dict(wavelength='r')] = \
            self.trimmed_reg_profiles.sel(wavelength='410') / self.trimmed_reg_profiles.sel(wavelength='470')
        self.trimmed_reg_profiles.loc[dict(wavelength='oxd')] = ip.r_to_oxd(
            self.trimmed_reg_profiles.loc[dict(wavelength='r')])
        self.trimmed_reg_profiles.loc[dict(wavelength='e')] = ip.oxd_to_redox_potential(
            self.trimmed_reg_profiles.loc[dict(wavelength='oxd')])

    def align_and_center(self):
        logging.info('Centering and rotating pharynxes')
        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(self.raw_image_data, self.seg_images)

    def calculate_midlines(self):
        logging.info('Calculating midlines')
        # noinspection PyTypeChecker
        self.midlines = ip.calculate_midlines(self.rot_seg, self.rot_fl)

    def measure_under_midlines(self):
        logging.info('Measuring under midlines')
        step = self.raw_image_data.x.size // 4
        self.raw_profiles = ip.measure_under_midlines(
            self.rot_fl, self.midlines, (step, self.raw_image_data.x.size - step), n_points=self.n_midline_pts
        )
        self.raw_profiles.data = np.nan_to_num(self.raw_profiles.data)
        self.raw_profiles = ip.align_pa(self.raw_profiles)

        # Add attributes
        self.raw_profiles = self.add_experiment_metadata_to_data_array(self.raw_profiles)

    def load_movement_annotation(self):
        try:
            df = pd.read_csv(self.experiment_dir.joinpath(self.experiment_id + '-mvmt.csv'))
            df = df.pivot_table(index='animal', columns=['region', 'pair'], values='movement')
            df = df.stack('pair')
            self.movement = df
        except FileNotFoundError:
            pass

    def load_strains(self):
        self.strains = pio.load_strain_map_from_disk(self.experiment_dir.joinpath(self.experiment_id + '-indexer.csv'))

    def segment_pharynxes(self):
        logging.info('Segmenting pharynxes')
        return ip.segment_pharynxes(self.raw_image_data, self.seg_threshold)

    def generate_summary_table(self):
        data = self.trimmed_raw_profiles
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
                    value_vars=data.pair.data, var_name='pair',
                    id_vars=['animal', 'strain', 'region', 'experiment'], value_name=wvl)
                region_dfs.append(sub_df)
            df_tmp = pd.concat(region_dfs, axis=1)
            df_tmp = df_tmp.loc[:, ~df_tmp.columns.duplicated()]
            dfs.append(df_tmp)
        df = pd.concat(dfs)
        df.reset_index(drop=True, inplace=True)

        self.summary_table = df

        if self.movement is not None:
            self.summary_table = self.summary_table.join(self.movement, on=['animal', 'pair'])

        return self.summary_table

    def filter_by_exclude_status(self, data):
        return data.loc[dict(strain=np.logical_not(self.include_idx))]

    def persist_to_disk(self):
        output_dir = self.experiment_dir
        logging.info(f'Saving {self.experiment_id} inside {output_dir}')

        # Persist the region means
        summary_table_filename = output_dir.joinpath(self.experiment_id + '-summary_table.csv')
        logging.info(f'Saving region means to {summary_table_filename}')
        self.generate_summary_table()
        self.summary_table.to_csv(summary_table_filename, index=False)

        # Persist the profile data
        self.persist_profile_data(output_dir, output_format='netcdf')

        # Plots
        self.save_profile_summary_plots(self.trimmed_reg_profiles, self.experiment_id)
        self.save_cat_plots()

    def persist_profile_data(self, output_dir: Path, output_format='netcdf'):
        if output_format == 'netcdf':
            profile_data_filename = output_dir.joinpath(self.experiment_id + '-profile_data.nc')
            logging.info(f'Saving profile data to {profile_data_filename}')
            self.trimmed_raw_profiles.to_netcdf(profile_data_filename)

        else:
            raise ValueError('invalid profile data persistence type')

    def save_profile_summary_plots(self, prof_data: xr.DataArray, fname_stem: str):
        fig_dir = self.make_fig_dir()
        for wvl in prof_data.wavelength.data:
            for pair in prof_data.pair.data:
                fig, ax = plots.plot_profile_avg_by_strain(prof_data.sel(wavelength=wvl, pair=pair),
                                                           ax_title=f'{fname_stem}-{wvl}')
                plots.add_regions_to_axis(ax, self.scaled_regions)
                fig.savefig(fig_dir.joinpath(f'{fname_stem}-{wvl}-{pair}.pdf'))

    def save_cat_plots(self):
        fig_dir = self.make_fig_dir()
        cat_plot_dir = fig_dir.joinpath('statistical')
        cat_plot_dir.mkdir(parents=True, exist_ok=True)

        for wvl in self.trimmed_reg_profiles.wavelength.data:
            for pair in self.trimmed_reg_profiles.pair.data:
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

    def generate_per_animal_reports(self):
        fig_dir = self.make_fig_dir()

        f_name = 'per_animal_reports.pdf'
        output_filename = fig_dir.joinpath(f_name)

        logging.info(f'Saving per-animal reports to {output_filename}')
        with PdfPages(str(output_filename)) as pdf:
            for i in tqdm.trange(len(self.strains)):
                fig = self.single_animal_diagnostic_plot(i)
                pdf.savefig(fig)
                plt.close()

    def single_animal_diagnostic_plot(self, i):
        fig = plt.figure(constrained_layout=True, figsize=(15, 15))
        gs = GridSpec(5, 3, figure=fig)
        midline_xs = np.arange(40, 120)

        for pair in range(self.raw_image_data.pair.size):
            i410 = self.rot_fl.sel(wavelength='410', pair=pair).isel(strain=i)
            i470 = self.rot_fl.sel(wavelength='470', pair=pair).isel(strain=i)
            ax = fig.add_subplot(gs[pair, 0])
            ax.imshow(i410)
            ax.plot(midline_xs, self.midlines[i]['410'][pair](midline_xs), color='orange')
            ax.set_title(f'410-{pair}')

            ax = fig.add_subplot(gs[pair, 1])
            ax.imshow(i470)
            ax.plot(midline_xs, self.midlines[i]['470'][pair](midline_xs), color='r')
            ax.set_title(f'470-{pair}')

            ax = fig.add_subplot(gs[pair, 2])
            ax.imshow(i410 / i470)
            ax.plot(midline_xs, self.midlines[i]['410'][pair](midline_xs), color='orange', label='410',
                    alpha=0.5)
            ax.plot(midline_xs, self.midlines[i]['470'][pair](midline_xs), color='r', label='470', alpha=0.5)
            ax.set_title(f'(410/470)-{pair}')
            ax.legend()

        ax = fig.add_subplot(gs[2, :])
        for pair in range(self.raw_image_data.pair.size):
            ax.plot(self.trimmed_raw_profiles.sel(wavelength='410', pair=pair).isel(strain=i),
                    label=f'410-{pair}')
            ax.plot(self.trimmed_raw_profiles.sel(wavelength='470', pair=pair).isel(strain=i),
                    label=f'470-{pair}')
        ax.legend()

        ax = fig.add_subplot(gs[3:, :])
        for pair in range(self.raw_image_data.pair.size):
            ax.plot(self.trimmed_raw_profiles.sel(wavelength='e', pair=pair).isel(strain=i),
                    label=f'E-{pair}')
            ax.set_ylim([np.nanquantile(self.trimmed_raw_profiles.sel(wavelength='e').data, .01),
                         np.nanquantile(self.trimmed_raw_profiles.sel(wavelength='e').data, .99)])
        ax.legend()

        plt.suptitle(f'Animal {i} ({self.strains[i]})')

        return fig


@attr.s(auto_attribs=True)
class PairExperiment(Experiment):
    """
    This is the paired ratio experiment
    """

    # Required initialization parameters
    image_display_order: List[str] = [
        '410_1', '470_1', 'r1',
        '410_2', '470_2', 'r2'
    ]

    def __attrs_post_init__(self):
        logging.info(f'Starting full pipeline run for {self.experiment_dir}')
        self.experiment_dir = Path(self.experiment_dir)
        self.experiment_id = self.experiment_dir.stem
        self.load_strains()
        self.scale_region_boundaries()
        self.load_movement_annotation()
        if self.include_idx is None:
            self.include_idx = np.ones(len(self.strains))
        self.full_pipeline()

    def full_pipeline(self):
        self.load_images()
        self.scale_region_boundaries()
        if self.seg_images is None:
            self.seg_images = self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        self.register()
        self.trim_data()
        self.calculate_redox()
        self.generate_summary_table()

        logging.info(f'Finished full pipeline run for {self.experiment_dir}')


if __name__ == '__main__':
    experiment_path = Path("/Users/sean/code/wormAnalysis/data/paired_ratio/2017_02_22-HD233_SAY47/")
    ex = PairExperiment(experiment_path, "TL/470/410/470/410")
