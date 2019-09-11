"""
This module contains the experiment data structures and classes - the central objects
that manage data and perform the measurement extraction, quantification, and analysis.

``Experiment`` is the base class, on which different types of experiments may be built.
For example, the ``PairExperiment`` class encompasses the acquisition strategy wherein
multiple pairs of images are taken sequentially for each animal.

"""

import datetime
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from cached_property import cached_property
from numpy.polynomial import Polynomial
import logging

from pharynx_redox import (
    image_processing as ip,
    profile_processing,
    pharynx_io as pio,
    plots,
    utils,
)


@dataclass
class Experiment:
    """
    TODO: Documentation
    """

    ###################################################################################
    # REQUIRED PARAMETERS
    ###################################################################################

    experiment_dir: Path
    imaging_scheme: str

    ###################################################################################
    # PIPELINE PARAMETERS
    ###################################################################################

    # R -> OxD parameters
    r_min: float = 0.852
    r_max: float = 6.65
    instrument_factor: float = 0.171

    # OxD -> E parameters
    midpoint_potential: float = -265.0
    z: int = 2
    temperature: float = 22.0

    # Pipeline Parameters
    trimmed_profile_length: int = 300
    n_midline_pts: int = 200
    seg_threshold: int = 2000
    trim_threshold: int = 3000
    reg_lambda: float = 0.01
    frame_specific_midlines: bool = False

    # Registration Parameters
    register: bool = False
    smooth_lambda: int = 1e-5
    rough_lambda: float = 1e-7
    warp_lambda: int = 1e-1
    smooth_nbasis: int = 50
    rough_nbasis: int = 200
    warp_to_mean: bool = False

    regions: dict = field(
        default_factory=lambda: {
            "pm3": [0.07, 0.28],
            "pm4": [0.33, 0.45],
            "pm5": [0.53, 0.70],
            "pm6": [0.80, 0.86],
            "pm7": [0.88, 0.96],
        }
    )
    strategy: str = None

    ###################################################################################
    # COMPUTED PROPERTIES
    ###################################################################################

    _scaled_regions: dict = None
    _movement: pd.DataFrame = None
    _strains: np.ndarray = None
    _raw_image_data: xr.DataArray = None
    _summary_table: pd.DataFrame = None

    seg_images: xr.DataArray = None

    rot_fl: xr.DataArray = None
    rot_seg: xr.DataArray = None

    midlines: List[Dict[str, List[Polynomial]]] = None

    untrimmed_profiles: xr.DataArray = None
    trimmed_profiles: xr.DataArray = None

    save_summary_plots: bool = False

    ####################################################################################
    # PROPERTIES
    ####################################################################################

    @property
    def scaled_regions(self):
        self._scaled_regions = {
            region: [int(self.trimmed_profile_length * x) for x in bounds]
            for region, bounds in self.regions.items()
        }
        return self._scaled_regions

    @cached_property
    def strains(self):
        self._strains = pio.load_strain_map_from_disk(
            self.experiment_dir.joinpath(self.experiment_id + "-indexer.csv")
        )
        return self._strains

    @property
    def experiment_id(self):
        return self.experiment_dir.stem

    @cached_property
    def images(self):
        # TODO: allow '.tiff' as well
        raw_image_path = Path(self.experiment_dir).joinpath(self.experiment_id + ".tif")
        self._raw_image_data = pio.load_images(
            raw_image_path, self.imaging_scheme, self.strains
        )
        return self._raw_image_data

    @cached_property
    def movement(self):
        try:
            df = pd.read_csv(
                self.experiment_dir.joinpath(self.experiment_id + "-mvmt.csv")
            )
            df = df.pivot_table(
                index="animal", columns=["region", "pair"], values="movement"
            )
            df = df.stack("pair")
            self._movement = df
            return self._movement
        except FileNotFoundError:
            return None

    @cached_property
    def analysis_dir(self):
        date_str = datetime.datetime.now().strftime("%Y-%m-%d")
        analysis_dir_ = self.experiment_dir.joinpath(
            "analyses", utils.get_valid_filename(f"{date_str}_{self.strategy}")
        )
        analysis_dir_.mkdir(parents=True, exist_ok=True)
        return analysis_dir_

    @cached_property
    def summary_table(self):
        dfs = []
        for region, bounds in self.scaled_regions.items():
            region_dfs = []
            for wvl in self.trimmed_profiles.wavelength.data:
                sub_df = (
                    self.trimmed_profiles[dict(position=range(bounds[0], bounds[1]))]
                    .mean(dim="position")
                    .sel(wavelength=wvl)
                    .to_pandas()
                )
                sub_df = sub_df.reset_index()
                sub_df["animal"] = range(len(sub_df))
                sub_df["region"] = region
                sub_df["experiment"] = self.experiment_id
                sub_df["strategy"] = self.strategy
                sub_df = sub_df.melt(
                    value_vars=self.trimmed_profiles.pair.data,
                    var_name="pair",
                    id_vars=["animal", "strain", "region", "experiment", "strategy"],
                    value_name=wvl,
                )
                region_dfs.append(sub_df)
            df_tmp = pd.concat(region_dfs, axis=1)
            df_tmp = df_tmp.loc[:, ~df_tmp.columns.duplicated()]
            dfs.append(df_tmp)
        df = pd.concat(dfs, sort=False)
        df.reset_index(drop=True, inplace=True)

        self._summary_table = pd.concat(dfs)

        if self.movement is not None:
            self._summary_table = self.summary_table.join(
                self.movement, on=["animal", "pair"]
            )

        return self._summary_table

    ####################################################################################
    # PIPELINE
    ####################################################################################
    def segment_pharynxes(self):
        logging.info("Segmenting pharynxes")

        return ip.segment_pharynxes(self.images, self.seg_threshold)

    def align_and_center(self):
        logging.info("Centering and rotating pharynxes")
        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(
            self.images, self.seg_images
        )

    def calculate_midlines(self):
        logging.info("Calculating midlines")
        # noinspection PyTypeChecker
        self.midlines = ip.calculate_midlines(self.rot_seg, degree=4)

    def measure_under_midlines(self):
        logging.info("Measuring under midlines")
        self.untrimmed_profiles = ip.measure_under_midlines(
            self.rot_fl,
            self.midlines,
            n_points=self.n_midline_pts,
            frame_specific=self.frame_specific_midlines,
        )
        self.untrimmed_profiles = ip.align_pa(self.untrimmed_profiles)
        self.untrimmed_profiles = self.add_experiment_metadata_to_data_array(
            self.untrimmed_profiles
        )

    def register_profiles(self):
        logging.info("Registering profiles")
        reg_data = profile_processing.register_profiles(
            self.untrimmed_profiles,
            smooth_lambda=self.smooth_lambda,
            rough_lambda=self.rough_lambda,
            warp_lam=self.warp_lambda,
            smooth_nbasis=self.smooth_nbasis,
            rough_nbasis=self.rough_nbasis,
            warp_to_mean=self.warp_to_mean,
        )
        self.untrimmed_profiles = reg_data.reg_data

    def trim_data(self):
        logging.info("Trimming intensity data")

        self.trimmed_profiles = self.add_experiment_metadata_to_data_array(
            profile_processing.trim_profiles(
                self.untrimmed_profiles,
                self.trim_threshold,
                self.trimmed_profile_length,
            )
        )

    def calculate_redox(self):
        logging.info("Calculating redox measurements")

        # Expand the trimmed_intensity_data to include new wavelengths
        new_wvls = np.append(self.trimmed_profiles.wavelength.data, ["r", "oxd", "e"])

        self.trimmed_profiles = self.trimmed_profiles.reindex(wavelength=new_wvls)

        self.trimmed_profiles.loc[dict(wavelength="r")] = self.trimmed_profiles.sel(
            wavelength="410"
        ) / self.trimmed_profiles.sel(wavelength="470")
        self.trimmed_profiles.loc[dict(wavelength="oxd")] = profile_processing.r_to_oxd(
            self.trimmed_profiles.loc[dict(wavelength="r")]
        )
        self.trimmed_profiles.loc[
            dict(wavelength="e")
        ] = profile_processing.oxd_to_redox_potential(
            self.trimmed_profiles.loc[dict(wavelength="oxd")]
        )

    def flip_at(self, idx):
        np.fliplr(self.rot_fl[:, idx])
        np.fliplr(self.rot_seg[:, idx])
        np.fliplr(self.untrimmed_profiles[:, idx])
        np.fliplr(self.trimmed_profiles[:, idx])

    ####################################################################################################################
    # PERSISTENCE / IO
    ####################################################################################################################
    def make_fig_dir(self):
        fig_dir = self.analysis_dir.joinpath("figs")
        fig_dir.mkdir(parents=True, exist_ok=True)
        return fig_dir

    def persist_profile_data(self):
        analysis_dir = self.analysis_dir

        profile_data_filename = analysis_dir.joinpath(
            self.experiment_id + "-profile_data.nc"
        )
        logging.info(f"Saving profile data to {profile_data_filename}")
        self.trimmed_profiles.to_netcdf(profile_data_filename)

    def load_profile_data(self):
        analysis_dir = self.analysis_dir
        profile_data_raw_filename = analysis_dir.joinpath(
            self.experiment_id + "-profile_data.nc"
        )
        logging.info(f"Loading profile data from {profile_data_raw_filename}")
        self.trimmed_profiles = xr.open_dataarray(profile_data_raw_filename)

    def save_profile_summary_plots(self, prof_data: xr.DataArray):
        fig_dir = self.make_fig_dir()
        for wvl in prof_data.wavelength.data:
            for pair in prof_data.pair.data:
                fig, ax = plots.plot_profile_avg_by_strain(
                    prof_data.sel(wavelength=wvl, pair=pair),
                    ax_title=f"{self.experiment_id}-{wvl}",
                )
                plots.add_regions_to_axis(ax, self.scaled_regions)
                fig.savefig(fig_dir.joinpath(f"{self.experiment_id}-{wvl}-{pair}.pdf"))

    def save_cat_plots(self):
        fig_dir = self.make_fig_dir()
        cat_plot_dir = fig_dir.joinpath("statistical")
        cat_plot_dir.mkdir(parents=True, exist_ok=True)

        for wvl in self.trimmed_profiles.wavelength.data:
            for pair in self.trimmed_profiles.pair.data:
                plt.clf()
                sns.violinplot(
                    y=wvl,
                    x="strain",
                    data=self.summary_table[self.summary_table.pair == pair],
                    inner="quartiles",
                )
                plt.savefig(cat_plot_dir.joinpath(f"violin-{wvl}-{pair}.pdf"))

                plt.clf()
                sns.catplot(
                    y=wvl,
                    x="strain",
                    col="region",
                    col_wrap=2,
                    kind="violin",
                    data=self.summary_table[self.summary_table.pair == pair],
                )
                plt.savefig(cat_plot_dir.joinpath(f"violin_regions-{wvl}-{pair}.pdf"))
                plt.close()

    def persist_to_disk(self, summary_plots=False):
        logging.info(f"Saving {self.experiment_id} inside {self.experiment_dir}")

        # Persist the region means
        summary_table_filename = self.analysis_dir.joinpath(
            self.experiment_id + "-summary_table.csv"
        )
        logging.info(f"Saving region means to {summary_table_filename}")
        # self.summary_table.to_csv(summary_table_filename, index=False)

        # Persist the profile data
        self.persist_profile_data()

        # Plots
        if summary_plots:
            self.save_profile_summary_plots(self.trimmed_profiles)
            self.save_cat_plots()

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
            strategy=self.strategy,
        )


@dataclass
class PairExperiment(Experiment):
    """
    TODO: Documentation
    """

    strategy: str = "frame specific midlines with registration"
    register: bool = True

    # Required initialization parameters
    image_display_order: List[str] = field(
        default_factory=lambda: ["410_1", "470_1", "r1", "410_2", "470_2", "r2"]
    )

    def full_pipeline(self):
        logging.info(f"Starting full pipeline run for {self.experiment_dir}")

        if self.seg_images is None:
            self.seg_images = self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        if self.register:
            self.register_profiles()
        self.trim_data()
        self.calculate_redox()
        self.persist_to_disk(summary_plots=self.save_summary_plots)

        logging.info(f"Finished full pipeline run for {self.experiment_dir}")

        return self


@dataclass
class CataExperiment(Experiment):
    """
    TODO: Documentation
    """

    strategy: str = "cata"
    frame_specific_midlines: bool = False

    def full_pipeline(self):
        logging.info(f"Starting full pipeline run for {self.experiment_dir}")
        if self.seg_images is None:
            self.seg_images = self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        self.trim_data()
        self.calculate_redox()
        self.persist_to_disk(summary_plots=self.save_summary_plots)

        logging.info(f"Finished full Cata pipeline run for {self.experiment_dir}")

        return self

    def segment_pharynxes(self):
        ref_wvl = "410"
        seg_images = super().segment_pharynxes()

        # In Cata's pipeline, the fluorescent images are also masked
        for animal_idx in np.arange(seg_images.strain.size):
            for wvl_idx in np.arange(seg_images.wavelength.size):
                for pair in seg_images.pair:
                    ref_seg = seg_images.sel(wavelength=ref_wvl, pair=pair).isel(
                        strain=animal_idx
                    )
                    img = self.images.isel(
                        strain=animal_idx, wavelength=wvl_idx, pair=pair
                    )

                    masked = ref_seg * img

                    self.images[animal_idx, wvl_idx, pair] = masked

        return seg_images

    def measure_under_midlines(self, ref_wvl="410"):
        logging.info("Measuring under midlines")
        self.untrimmed_profiles = ip.measure_under_midlines(
            self.rot_fl,
            self.midlines,
            n_points=self.n_midline_pts,
            frame_specific=self.frame_specific_midlines,
        )
        self.untrimmed_profiles = ip.align_pa(self.untrimmed_profiles)
        self.untrimmed_profiles = self.add_experiment_metadata_to_data_array(
            self.untrimmed_profiles
        )


if __name__ == "__main__":
    import logging

    logging.basicConfig(
        format="%(asctime)s %(levelname)s:%(message)s",
        level=logging.DEBUG,
        datefmt="%I:%M:%S",
    )
    experiment_path = Path(
        "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/"
    )
    ex = PairExperiment(experiment_path, "TL/470/410/470/410")
    ex.full_pipeline()
