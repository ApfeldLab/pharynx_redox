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
from typing import List, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import skfda
import xarray as xr
from cached_property import cached_property
from numpy.polynomial import Polynomial
import logging

from marshmallow import Schema, fields, INCLUDE, pprint, ValidationError

from pharynx_redox import (
    image_processing as ip,
    profile_processing,
    pharynx_io as pio,
    plots,
    utils,
    constants,
)


class ExperimentSchema(Schema):
    experiment_id = fields.Str(required=True)

    strategy = fields.Str()

    r_min = fields.Float(required=True)
    r_max = fields.Float(required=True)
    instrument_factor = fields.Float(required=True)

    midpoint_potential = fields.Float(required=True)
    z = fields.Float(required=True)
    temperature = fields.Float(required=True)


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
    strategy: str = None

    register: bool = False

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
    seg_threshold: int = 500
    trim_threshold: int = 3000
    frame_specific_midlines: bool = False
    smooth_unregistered_data: bool = False
    measurement_order: int = 1
    measure_thickness: float = 0.0

    # Registration Parameters
    register: bool = False

    warp_n_basis: float = 20.0
    warp_order: float = 4.0
    warp_lambda: float = 10.0

    smooth_lambda: float = 0.001
    smooth_n_breaks: float = 50.0
    smooth_order: float = 4.0

    rough_lambda: float = 0.001
    rough_n_breaks: float = 200.0
    rough_order = 4.0

    n_deriv: float = 0.0

    # Summarization parameters
    trimmed_regions: dict = field(default_factory=lambda: constants.trimmed_regions)
    untrimmed_regions: dict = field(default_factory=lambda: constants.untrimmed_regions)
    pointwise_summaries: bool = False
    save_summary_plots: bool = False
    should_save_profile_data: bool = True
    should_save_summary_data: bool = True

    ###################################################################################
    # COMPUTED PROPERTY PLACEHOLDERS
    ###################################################################################

    _scaled_regions: dict = None
    _movement: pd.DataFrame = None
    _strains: np.ndarray = None
    _raw_image_data: xr.DataArray = None
    _image_data: xr.DataArray = None
    _summary_table: pd.DataFrame = None

    seg_images: xr.DataArray = None

    rot_fl: xr.DataArray = None
    rot_seg: xr.DataArray = None

    midlines: List[Dict[str, List[Polynomial]]] = None

    untrimmed_profiles: xr.DataArray = None
    trimmed_profiles: xr.DataArray = None

    warps: List[skfda.FDataGrid] = field(default_factory=list)

    _parameter_dict: dict = None

    ####################################################################################
    # COMPUTED PROPERTIES
    ####################################################################################

    @property
    def parameter_dict(self):
        return {}

    @property
    def scaled_regions(self):
        self._scaled_regions = {
            region: [int(self.trimmed_profile_length * x) for x in bounds]
            for region, bounds in self.trimmed_regions.items()
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
        """
        This returns the median-subtracted images
        """
        # TODO: allow '.tiff' as well
        # self._image_data = ip.subtract_medians(self.raw_images).astype(np.uint16)
        self._image_data = self.raw_images
        return self._image_data

    @cached_property
    def raw_images(self):
        """
        This returns the raw (non-median-subtracted) images
        """
        raw_image_path = Path(self.experiment_dir).joinpath(self.experiment_id + ".tif")
        self._raw_image_data = pio.load_images(
            raw_image_path, self.imaging_scheme, self.strains
        )
        return self._raw_image_data

    @cached_property
    def movement(self):
        movement_filepath = self.experiment_dir.joinpath(
            self.experiment_id + "-mvmt.csv"
        )
        try:
            df = pd.read_csv(movement_filepath)
            df = df.pivot_table(
                index="animal", columns=["region", "pair"], values="movement"
            )
            df = df.stack("pair")
            self._movement = df
            return self._movement
        except FileNotFoundError:
            logging.warning(f"Tried to access {movement_filepath}; file was not found")
            return None

    @cached_property
    def analysis_dir(self):
        date_str = datetime.datetime.now().strftime("%Y-%m-%d")
        analysis_dir_ = self.experiment_dir.joinpath(
            "analyses", utils.get_valid_filename(f"{date_str}_{self.strategy}")
        )
        analysis_dir_.mkdir(parents=True, exist_ok=True)
        return analysis_dir_

    @property
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
                ).reset_index()
                sub_df["animal"] = range(len(sub_df))
                sub_df["region"] = region
                sub_df["experiment"] = self.experiment_id
                sub_df["strategy"] = self.strategy
                sub_df["strain"] = self.strains

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
            self._summary_table = self._summary_table.join(
                self.movement, on=["animal", "pair"]
            )

        if self.pointwise_summaries:
            self._summary_table["pointwise"] = True
        else:
            self._summary_table["pointwise"] = False
            # Calculate R, OxD, and E using region summaries
            self._summary_table["r"] = (
                self._summary_table["410"] / self._summary_table["470"]
            )
            self._summary_table["oxd"] = profile_processing.r_to_oxd(
                self._summary_table["r"],
                r_min=self.r_min,
                r_max=self.r_max,
                instrument_factor=self.instrument_factor,
            )
            self._summary_table["e"] = profile_processing.oxd_to_redox_potential(
                self._summary_table["oxd"],
                midpoint_potential=self.midpoint_potential,
                z=self.z,
                temperature=self.temperature,
            )

        return self._summary_table

    @property
    def untrimmed_profile_data_filepath(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-untrimmed_profile_data.nc"
        )

    @property
    def trimmed_profile_data_filepath(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-trimmed_profile_data.nc"
        )

    @property
    def warp_data_filepath(self):
        return self.analysis_dir.joinpath(self.experiment_id + "-warp_data.npy")

    ####################################################################################
    # PIPELINE
    ####################################################################################

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

    def segment_pharynxes(self):
        logging.info("Segmenting pharynxes")

        return ip.segment_pharynxes(self.images, self.seg_threshold)

    def align_and_center(self):
        logging.info("Centering and rotating pharynxes")
        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(
            self.images, self.seg_images, blur_seg_thresh=self.seg_threshold
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
            order=self.measurement_order,
            thickness=self.measure_thickness
        )
        self.untrimmed_profiles = ip.align_pa(self.untrimmed_profiles)
        self.untrimmed_profiles = self.add_experiment_metadata_to_data_array(
            self.untrimmed_profiles
        )
        if (self.register == False) and (self.smooth_unregistered_data):
            self.untrimmed_profiles = profile_processing.smooth_profile_data(
                self.untrimmed_profiles
            )

    def register_profiles(self):
        logging.info("Registering profiles")

        self.untrimmed_profiles = profile_processing.register_profiles_matlab(
            self.untrimmed_profiles,
            n_deriv=self.n_deriv,
            rough_lambda=self.rough_lambda,
            smooth_lambda=self.smooth_lambda,
        )[0]

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
        self.trimmed_profiles = utils.add_derived_wavelengths(self.trimmed_profiles)
        self.untrimmed_profiles = utils.add_derived_wavelengths(self.untrimmed_profiles)

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
        logging.info(
            f"Saving untrimmed profile data to {self.untrimmed_profile_data_filepath}"
        )
        pio.save_profile_data(
            self.untrimmed_profiles, self.untrimmed_profile_data_filepath
        )

        logging.info(
            f"Saving trimmed profile data to {self.trimmed_profile_data_filepath}"
        )
        pio.save_profile_data(self.trimmed_profiles, self.trimmed_profile_data_filepath)

        logging.info(f"Saving warp data to {self.warp_data_filepath}")
        with open(self.warp_data_filepath, "wb") as f:
            np.save(f, self.warps)

    def load_profile_data(self):
        logging.info(
            f"Loading trimmed profile data from {self.trimmed_profile_data_filepath}"
        )
        self.trimmed_profiles = pio.load_profile_data(
            self.trimmed_profile_data_filepath
        )

        logging.info(
            f"Loading untrimmed profile data from {self.untrimmed_profile_data_filepath}"
        )
        self.untrimmed_profiles = pio.load_profile_data(
            self.untrimmed_profile_data_filepath
        )

    def save_summary_data(self):
        # Persist the region means
        summary_table_filename = self.analysis_dir.joinpath(
            self.experiment_id + "-summary_table.csv"
        )
        logging.info(f"Saving region means to {summary_table_filename}")
        self.summary_table.to_csv(summary_table_filename, index=False)

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

        if self.should_save_summary_data:
            self.save_summary_data()

        if self.should_save_profile_data:
            self.persist_profile_data()

        # Plots
        if summary_plots:
            self.save_profile_summary_plots(self.trimmed_profiles)
            self.save_cat_plots()

    def save_normed_ratio_images(self, vmin=-5, vmax=5):
        for pair in self.images.pair.values:
            r = self.images.sel(wavelength="410", pair=pair) / self.images.sel(
                wavelength="470", pair=pair
            )
            seg = self.seg_images.sel(wavelength="410", pair=pair)

            processed_img_dir = self.experiment_dir.joinpath("processed_images")
            processed_img_dir.mkdir(parents=True, exist_ok=True)
            output_fname = processed_img_dir.joinpath(
                f"{self.experiment_id}_normed-ratio_pair={pair}.tif"
            )
            ip.create_normed_rgb_ratio_stack(r, seg, output_filename=output_fname)

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
            experiment_id=self.experiment_id,
        )


if __name__ == "__main__":
    import logging

    logging.basicConfig(
        format="%(asctime)s %(levelname)s:%(message)s",
        level=logging.DEBUG,
        datefmt="%I:%M:%S",
    )
    experiment_path = Path(
        "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_08_23-HD233_4mm_lev"
    )
    ex = Experiment(experiment_path, "TL/470/410/470/410", register=False)
    ex.full_pipeline()

    utils.measure_shifted_midlines(ex, (-2, 2), 5)
