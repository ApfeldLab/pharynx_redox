"""
This module contains the Experiment class. This class is the object that orchestrates
the analysis pipeline for redox imaging experiments.
"""

import sys

sys.path.append("/Users/sean/code/pharynx_redox/")

import datetime
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List
import os
import sys
import traceback
import yaml
import warnings

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import xarray as xr
from numpy.polynomial import Polynomial
import click
from tqdm import tqdm

from pharynx_redox import constants
from pharynx_redox import image_processing as ip
from pharynx_redox import io as pio
from pharynx_redox import profile_processing, utils, plots


@dataclass
class Experiment:
    """
    This class orchestrates the analysis pipeline for our redox imaging experiments.
    """

    ###################################################################################
    # REQUIRED PARAMETERS
    ###################################################################################
    experiment_dir: Path

    ###################################################################################
    # PIPELINE PARAMETERS
    ###################################################################################

    channel_order: List[str] = None
    strategy: str = ""
    reference_wavelength: str = "410"

    # R Parameters
    fl_wvls = ["410", "470"]
    ratio_numerator: str = "410"
    ratio_denominator: str = "470"

    # R -> OxD parameters
    r_min: float = 0.852
    r_max: float = 6.65
    instrument_factor: float = 0.171

    # OxD -> E parameters
    midpoint_potential: float = -265.0
    z: int = 2
    temperature: float = 22.0

    # Pipeline Parameters
    trimmed_profile_length: int = 200
    untrimmed_profile_length: int = 200
    seg_threshold: int = 2000
    measurement_order: int = 1
    measure_thickness: float = 0.0

    # Registration Parameters
    channel_register: int = 1
    population_register: int = 1

    n_deriv: float = 0.0

    warp_n_basis: float = 20.0
    warp_order: float = 4.0
    warp_lambda: float = 10.0

    smooth_lambda: float = 0.001
    smooth_n_breaks: float = 50.0
    smooth_order: float = 4.0

    rough_lambda: float = 0.001
    rough_n_breaks: float = 200.0
    rough_order = 4.0

    ####################################################################################
    # Summarization parameters
    ####################################################################################
    trimmed_regions: dict = field(
        default_factory=lambda: constants.trimmed_regions_with_medial
    )
    untrimmed_regions: dict = field(
        default_factory=lambda: constants.untrimmed_regions_with_medial
    )

    ####################################################################################
    # Persistence / IO flags
    ####################################################################################
    should_save_plots: bool = False
    should_save_profile_data: bool = True
    should_save_summary_data: bool = True

    ###################################################################################
    # COMPUTED PROPERTY PLACEHOLDERS
    ###################################################################################
    curation_dict: dict = field(default_factory=dict)
    seg_images: xr.DataArray = None
    rot_fl: xr.DataArray = None
    rot_seg: xr.DataArray = None

    midlines: List[Dict[str, List[Polynomial]]] = None

    untrimmed_profiles: xr.DataArray = None
    trimmed_profiles: xr.DataArray = None

    warps: List = field(default_factory=list)

    ####################################################################################
    # COMPUTED PROPERTIES
    ####################################################################################

    def __post_init__(self):
        self.experiment_id = self.experiment_dir.stem
        self.settings_filepath = self.experiment_dir.joinpath("settings.yaml")
        self.try_to_load_from_config_file()

        if self.channel_order is None:
            raise (AttributeError("channel_order not specified"))

        # compute the filenames/paths for this experiment
        self.raw_img_stack_filepath = self.experiment_dir.joinpath(
            self.experiment_id + ".tif"
        )
        self.processed_images_dir = self.experiment_dir.joinpath("processed_images")
        self.rot_seg_dir = self.processed_images_dir.joinpath("rot_seg")
        self.rot_fl_dir = self.processed_images_dir.joinpath("rot_fl")
        self.fl_imgs_dir = self.processed_images_dir.joinpath("fluorescent_images")
        self.analysis_dir = self.get_analysis_dir()
        self.fig_dir = self.analysis_dir.joinpath("figs")

        self.orig_images_filepath = self.processed_images_dir.joinpath("images.nc")
        self.aligned_images_filepath = self.processed_images_dir.joinpath(
            "aligned_images.nc"
        )
        self.seg_images_filepath = self.processed_images_dir.joinpath("seg_images.nc")
        self.aligned_seg_images_filepath = self.processed_images_dir.joinpath(
            "aligned_seg_images.nc"
        )

        self.movement_filepath = self.experiment_dir.joinpath(
            self.experiment_id + "-mvmt.csv"
        )
        self.indexer_filepath = self.experiment_dir.joinpath(
            self.experiment_id + "-indexer.csv"
        )
        self.untrimmed_profile_data_filepath = self.analysis_dir.joinpath(
            self.experiment_id + "-untrimmed_profile_data.nc"
        )
        self.trimmed_profile_data_filepath = self.analysis_dir.joinpath(
            self.experiment_id + "-trimmed_profile_data.nc"
        )
        self.warp_data_filepath = self.analysis_dir.joinpath(
            self.experiment_id + "-warp_data.npy"
        )

        self.untrimmed_profile_data_csv_filepath = self.analysis_dir.joinpath(
            self.experiment_id + "-untrimmed_profile_data.csv"
        )

        self.trimmed_profile_data_csv_filepath = self.analysis_dir.joinpath(
            self.experiment_id + "-trimmed_profile_data.csv"
        )

        self.untrimmed_region_data_filepath = self.analysis_dir.joinpath(
            self.experiment_id + "-untrimmed_region_data.csv"
        )

        self.trimmed_region_data_filepath = self.analysis_dir.joinpath(
            self.experiment_id + "-trimmed_region_data.csv"
        )

        # Other computed properties
        self.strains = pio.load_strain_map_from_disk(self.indexer_filepath)

        self.scaled_regions_trimmed = {
            region: [int(self.trimmed_profile_length * x) for x in bounds]
            for region, bounds in self.trimmed_regions.items()
        }

        self.scaled_regions_untrimmed = {
            region: [int(self.untrimmed_profile_length * x) for x in bounds]
            for region, bounds in self.untrimmed_regions.items()
        }

        # load images
        self.raw_images = self._load_raw_images()
        self.images = ip.subtract_medians(self.raw_images)

    def try_to_load_from_config_file(self):
        try:
            with open(self.settings_filepath, "r") as f:
                config_dict = yaml.safe_load(f)
                try:
                    pipeline_params = config_dict["pipeline"]
                except KeyError as e:
                    logging.error(
                        f"Incorrect settings file format. There must be `pipeline` at the root. See example config file on Github repo."
                    )
                    raise e

                logging.info(
                    f"Loading configuration file from {self.settings_filepath}"
                )

                try:
                    self.curation_dict = config_dict["curation"]
                    logging.info("Setting curation information")
                except KeyError as e:
                    logging.info("No curation information specified in settings file")

                for key, val in pipeline_params.items():
                    try:
                        getattr(self, key)
                        setattr(self, key, val)
                    except AttributeError:
                        logging.warning(f"Parameter not recognized({key}={val})")
        except FileNotFoundError:
            logging.info("No configuration file found. Using defaults.")

    def _load_raw_images(self):
        """
        This returns the raw (non-median-subtracted) images
        """
        logging.info(f"Loading image data from {self.raw_img_stack_filepath}")
        raw_image_data = pio.load_images(
            img_stack_path=self.raw_img_stack_filepath,
            strain_map=self.strains,
            movement_path=self.movement_filepath,
            channel_order=self.channel_order,
        )

        raw_image_data = raw_image_data.assign_coords(
            {
                "experiment_id": (
                    ("animal",),
                    np.repeat(self.experiment_id, raw_image_data.animal.size),
                )
            }
        )

        raw_image_data = self.add_experiment_metadata_to_data_array(raw_image_data)

        return raw_image_data

    def _load_movement(self) -> pd.DataFrame:
        movement_filepath = self.experiment_dir.joinpath(
            self.experiment_id + "-mvmt.csv"
        )
        try:
            df = pd.read_csv(movement_filepath)
            df = df.pivot_table(
                index="animal", columns=["region", "pair"], values="movement"
            )
            df = df.stack("pair")
            return df
        except FileNotFoundError:
            logging.warning(f"Tried to access {movement_filepath}; file was not found")
            return None

    def get_analysis_dir(self) -> Path:
        date_str = datetime.datetime.now().strftime("%Y-%m-%d")
        analysis_dir_ = self.experiment_dir.joinpath(
            "analyses", utils.get_valid_filename(f"{date_str}_{self.strategy}")
        )
        # analysis_dir_.mkdir(parents=True, exist_ok=True)
        return analysis_dir_

    def make_analysis_dir(self) -> None:
        self.get_analysis_dir().mkdir(parents=True, exist_ok=True)

    @property
    def trimmed_summary_table(self):
        df = profile_processing.summarize_over_regions(
            self.trimmed_profiles,
            regions=self.trimmed_regions,
            ratio_numerator=self.ratio_numerator,
            ratio_denominator=self.ratio_denominator,
        )
        return df

    @property
    def untrimmed_summary_table(self):
        df = profile_processing.summarize_over_regions(
            self.untrimmed_profiles,
            regions=self.untrimmed_regions,
            ratio_numerator=self.ratio_numerator,
            ratio_denominator=self.ratio_denominator,
        )
        return df

    ####################################################################################
    # PIPELINE
    ####################################################################################

    def full_pipeline(self):
        logging.info(f"Starting full pipeline run for {self.experiment_dir}")
        logging.info(f"Making analysis directory at {self.analysis_dir}")
        self.make_analysis_dir()

        logging.info(f"Saving fluorescent images to {self.fl_imgs_dir}")
        pio.save_images_xarray_to_disk(
            self.raw_images, self.fl_imgs_dir, prefix=self.experiment_id
        )

        self.segment_pharynxes()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        self.register_profiles()
        self.trim_data()
        self.calculate_redox()
        # self.do_manual_AP_flips()
        self.persist_to_disk()

        logging.info(f"Finished full pipeline run for {self.experiment_dir}")

        return self

    def load_seg_imgs_from_disk(self, filepath):
        logging.info(f"Skipping segmentation; loading from {filepath}")

    def segment_pharynxes(self):
        if self.seg_images is not None:
            logging.info("masks have been specified. skipping mask generation")
            return

        try:
            self.load_masks()
        except IOError:
            logging.warn(f"Failed to load masks from {self.seg_images_filepath}")

            logging.info("Generating masks")
            self.seg_images = ip.segment_pharynxes(self.images, self.seg_threshold)

            self.save_masks()

    def align_and_center(self):
        logging.info("Centering and rotating pharynxes")
        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(
            self.images,
            self.seg_images,
            blur_seg_thresh=self.seg_threshold,
            reference_wavelength=self.reference_wavelength,
        )

        logging.info(f"Saving rotated FL images to {self.aligned_images_filepath}")
        pio.save_profile_data(self.rot_fl, self.aligned_images_filepath)

        logging.info(f"Saving rotated masks to {self.aligned_seg_images_filepath}")
        pio.save_profile_data(self.rot_seg, self.aligned_seg_images_filepath)

    def calculate_midlines(self):
        logging.info("Calculating midlines")
        # noinspection PyTypeChecker
        self.midlines = ip.calculate_midlines(self.rot_seg, degree=4)

    def measure_under_midlines(self):
        logging.info("Measuring under midlines")
        self.untrimmed_profiles = ip.measure_under_midlines(
            self.rot_fl,
            self.midlines,
            n_points=self.untrimmed_profile_length,
            frame_specific=False,
            order=self.measurement_order,
            thickness=self.measure_thickness,
        )
        self.untrimmed_profiles = profile_processing.align_pa(self.untrimmed_profiles)
        self.untrimmed_profiles = self.add_experiment_metadata_to_data_array(
            self.untrimmed_profiles
        )

    def register_profiles(self):
        reg_params = [
            "n_deriv",
            "warp_n_basis",
            "warp_order",
            "warp_lambda",
            "smooth_lambda",
            "smooth_n_breaks",
            "smooth_order",
            "rough_lambda",
            "rough_n_breaks",
            "rough_order",
        ]

        reg_param_dict = {k: getattr(self, k) for k in reg_params}

        if self.population_register:
            logging.info("Standardizing profiles")
            (
                self.untrimmed_profiles,
                self.warps,
            ) = profile_processing.register_profiles_pop(
                self.untrimmed_profiles,
                ratio_numerator=self.ratio_numerator,
                ratio_denominator=self.ratio_denominator,
                **reg_param_dict,
            )

        if self.channel_register:
            logging.info("Channel-Registering profiles")
            (
                self.untrimmed_profiles,
                self.warps,
            ) = profile_processing.channel_register(
                self.untrimmed_profiles,
                ratio_numerator=self.ratio_numerator,
                ratio_denominator=self.ratio_denominator,
                **reg_param_dict,
            )

    def trim_data(self):
        logging.info("Trimming intensity data")
        self.trimmed_profiles = self.add_experiment_metadata_to_data_array(
            profile_processing.trim_profiles(
                self.untrimmed_profiles,
                self.seg_threshold,
                ref_wvl=self.ratio_numerator,
            )
        )

    def calculate_redox(self):
        logging.info("Calculating redox measurements")

        # Images
        self.images = utils.add_derived_wavelengths(
            self.images,
            numerator=self.ratio_denominator,
            denominator=self.ratio_denominator,
        )
        self.rot_fl = utils.add_derived_wavelengths(
            self.rot_fl,
            numerator=self.ratio_denominator,
            denominator=self.ratio_denominator,
        )

        # Expand the trimmed_intensity_data to include new wavelengths
        self.trimmed_profiles = utils.add_derived_wavelengths(
            self.trimmed_profiles,
            numerator=self.ratio_numerator,
            denominator=self.ratio_denominator,
        )
        self.untrimmed_profiles = utils.add_derived_wavelengths(
            self.untrimmed_profiles,
            numerator=self.ratio_numerator,
            denominator=self.ratio_denominator,
        )

    def do_manual_AP_flips(self):
        try:
            to_flip = self.curation_dict["flip"]
            logging.info(f"Flipping animals {to_flip}")
            for idx in to_flip:
                self.flip_at(idx)

            # need to re-save all images after flipping
            logging.info("Re-Saving processed images after flipping")
            logging.info(f"Saving rotated FL images to {self.rot_fl_dir}")
            pio.save_images_xarray_to_disk(
                self.rot_fl, self.rot_fl_dir, prefix=self.experiment_id
            )

            logging.info(f"Saving rotated masks to {self.rot_seg_dir}")
            pio.save_images_xarray_to_disk(
                self.rot_seg, self.rot_seg_dir, prefix=self.experiment_id
            )

        except KeyError:
            logging.info(
                "No manual flips specified in settings file. Skipping AP flips."
            )

    def flip_at(self, idx):
        logging.debug(f"manually flipping animal {idx}")
        np.fliplr(self.rot_fl[:, idx])
        np.fliplr(self.rot_seg[:, idx])
        np.fliplr(self.untrimmed_profiles[:, idx])
        np.fliplr(self.trimmed_profiles[:, idx])

    ####################################################################################
    # PERSISTENCE / IO
    ####################################################################################

    def save_images(self):
        """Save this experiment's images to disk as netCDF4 files"""
        imgs_paths = [
            (self.images, self.orig_images_filepath),
            (self.rot_fl, self.aligned_images_filepath),
            (self.seg_images, self.seg_images_filepath),
            (self.rot_seg, self.aligned_seg_images_filepath),
        ]
        for img, path in imgs_paths:
            if img is not None:
                logging.info(f"Saving images to {path}")
                img.to_netcdf(path)

    # def load_images(self):
    # pass

    def make_fig_dir(self):
        fig_dir = self.analysis_dir.joinpath("figs")
        fig_dir.mkdir(parents=True, exist_ok=True)
        return fig_dir

    def save_plots(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig_dir = self.make_fig_dir()

            # First, save profile data plots
            profile_fig_dir = fig_dir.joinpath("profile_data")

            # need both trimmed and untrimmed data
            for prefix, data in zip(
                ("un", ""), (self.untrimmed_profiles, self.trimmed_profiles)
            ):
                profile_data_fig_dir = profile_fig_dir.joinpath(
                    prefix + "trimmed_profiles"
                )

                # individual data
                individual_data_fig_dir = profile_data_fig_dir.joinpath("individual")
                individual_data_fig_dir.mkdir(exist_ok=True, parents=True)
                for title, fig in plots.generate_wvl_pair_profile_plots(data):
                    title = title.replace(" ", "")
                    fig.savefig(
                        individual_data_fig_dir.joinpath(
                            f"{self.experiment_id}-{title}-individuals.pdf"
                        )
                    )
                    plt.close(fig)

                # avg. data
                avgs_data_fig_dir = profile_data_fig_dir.joinpath("avgs")
                avgs_data_fig_dir.mkdir(exist_ok=True, parents=True)
                for title, fig in plots.generate_avg_wvl_pair_profile_plots(data):
                    title = title.replace(" ", "")
                    fig.savefig(
                        avgs_data_fig_dir.joinpath(
                            f"{self.experiment_id}-{title}-avgs.pdf"
                        )
                    )
                    plt.close(fig)

            # frame-normed Ratio Images
            mvmt_annotation_img_path = self.fig_dir.joinpath(
                f"{self.experiment_id}-movement_annotation_imgs.pdf"
            )
            imgs = utils.add_derived_wavelengths(self.images)
            with PdfPages(mvmt_annotation_img_path) as pdf:
                for i in tqdm(range(self.raw_images.animal.size)):
                    fig = plots.plot_pharynx_R_imgs(imgs[i], mask=self.seg_images[i])
                    fig.suptitle(f"animal = {i}")
                    pdf.savefig(fig)
                    if (i % 20) == 0:
                        plt.close("all")

            # Pop-normed ratio images
            u = self.trimmed_profiles.sel(wavelength="r").mean()
            std = self.trimmed_profiles.sel(wavelength="r").std()

            for pair in self.rot_fl.pair.values:
                for tp in self.rot_fl.timepoint.values:
                    ratio_img_path = self.fig_dir.joinpath(
                        f"{self.experiment_id}-ratio_images-pair={pair};timepoint={tp}.pdf"
                    )
                    with PdfPages(ratio_img_path) as pdf:
                        logging.info(f"Saving ratio images to {ratio_img_path}")
                        for i in tqdm(range(self.rot_fl.animal.size)):
                            fig, ax = plt.subplots(dpi=300)
                            R = (
                                self.rot_fl.sel(
                                    wavelength=self.ratio_numerator,
                                    pair=pair,
                                    timepoint=tp,
                                )
                                / self.rot_fl.sel(
                                    wavelength=self.ratio_denominator,
                                    pair=pair,
                                    timepoint=tp,
                                )
                            )[i]
                            I = self.rot_fl.sel(
                                wavelength=self.ratio_numerator, pair=pair, timepoint=tp
                            )[i]
                            im, cbar = plots.imshow_ratio_normed(
                                R,
                                I,
                                r_min=u - (std * 1.96),
                                r_max=u + (std * 1.96),
                                colorbar=True,
                                i_max=5000,
                                i_min=1000,
                                ax=ax,
                            )
                            ax.plot(
                                *self.midlines.sel(
                                    wavelength=self.ratio_numerator,
                                    pair=pair,
                                    timepoint=tp,
                                )[i]
                                .values[()]
                                .linspace(),
                                color="green",
                                alpha=0.3,
                            )
                            strain = self.rot_fl.strain.values[i]
                            ax.set_title(f"Animal={i} ; Pair={pair} ; Strain={strain}")
                            cax = cbar.ax
                            for j in range(len(self.trimmed_profiles)):
                                cax.axhline(
                                    self.trimmed_profiles.sel(
                                        wavelength="r", pair=pair, timepoint=tp
                                    )[j].mean(),
                                    color="k",
                                    alpha=0.1,
                                )
                            cax.axhline(
                                self.trimmed_profiles.sel(
                                    wavelength="r", pair=pair, timepoint=tp
                                )[i].mean(),
                                color="k",
                            )
                            pdf.savefig()
                            if (i % 20) == 0:
                                plt.close("all")

    def persist_profile_data(self):
        # First, the netCDF4 format
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

        # Warps, if necessary
        if self.channel_register:
            logging.info(f"Saving warp data to {self.warp_data_filepath}")
            with open(self.warp_data_filepath, "wb") as f:
                np.save(f, self.warps)

        # Now, save the profile data in "tidy" format as CSV
        profile_processing.to_dataframe(self.trimmed_profiles, "value").to_csv(
            self.trimmed_profile_data_csv_filepath
        )
        profile_processing.to_dataframe(self.untrimmed_profiles, "value").to_csv(
            self.untrimmed_profile_data_csv_filepath
        )

    def save_summary_data(self):
        # Persist the region means
        logging.info(
            f"Saving untrimmed region means to {self.untrimmed_region_data_filepath}"
        )
        self.untrimmed_summary_table.to_csv(self.untrimmed_region_data_filepath)
        logging.info(
            f"Saving trimmed region means to {self.trimmed_region_data_filepath}"
        )
        self.trimmed_summary_table.to_csv(self.trimmed_region_data_filepath)

    def save_masks(self):
        logging.info(f"writing masks to {self.seg_images_filepath}")
        pio.save_profile_data(self.seg_images, self.seg_images_filepath)
        logging.info(f"Saved masks to {self.seg_images_filepath}")

    def load_masks(self):
        self.seg_images = pio.load_profile_data(self.seg_images_filepath)
        logging.info(f"Loaded masks from {self.seg_images_filepath}")

    def persist_to_disk(self):
        logging.info(f"Saving {self.experiment_id} inside {self.experiment_dir}")

        if self.should_save_summary_data:
            self.save_summary_data()

        if self.should_save_profile_data:
            self.persist_profile_data()

        if self.should_save_plots:
            self.save_plots()

    ####################################################################################
    # MISC / HELPER
    ####################################################################################
    def add_experiment_metadata_to_data_array(self, data_array: xr.DataArray):
        metadata_keys = [
            "strategy",
            "reference_wavelength",
            "r_min",
            "r_max",
            "ratio_numerator",
            "ratio_denominator",
            "instrument_factor",
            "midpoint_potential",
            "z",
            "temperature",
            "measurement_order",
            "measure_thickness",
            "channel_register",
            "population_register",
            "n_deriv",
            "warp_n_basis",
            "warp_order",
            "warp_lambda",
            "smooth_lambda",
            "smooth_n_breaks",
            "smooth_order",
            "rough_lambda",
            "rough_n_breaks",
            "rough_order",
        ]

        metadata_dict = {k: getattr(self, k) for k in metadata_keys}

        return data_array.assign_attrs(**metadata_dict)


def run_analysis(experiment_dir, log_level):
    """
    Analyze a stack of ratiometric pharynx images
    """
    log_map = {1: logging.INFO, 2: logging.DEBUG}
    if log_level > 0:
        # TODO: fix logging so debug messages don't come from matplotlib
        logging.basicConfig(
            format="%(asctime)s %(levelname)s:%(message)s",
            level=log_map[log_level],
            datefmt="%I:%M:%S",
        )
    Experiment(Path(experiment_dir)).full_pipeline()


if __name__ == "__main__":
    run_analysis(sys.argv[1], 1)
