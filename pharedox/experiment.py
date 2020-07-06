"""
This module contains the Experiment class. This class is the object that orchestrates
the analysis pipeline for redox imaging experiments.
"""

import datetime
import logging
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
from strictyaml import (
    Bool,
    Float,
    Int,
    Map,
    MapPattern,
    CommaSeparated,
    Str,
    load,
)
from tqdm import tqdm
from yaml import YAMLError

from pharedox import image_processing as ip
from pharedox import io as pio
from pharedox import plots, profile_processing, utils


class Experiment:
    """
    This class orchestrates the analysis pipeline for our redox imaging experiments.
    """

    experiment_schema = Map(
        {
            "pipeline": Map(
                {
                    "strategy": Str(),
                    "channel_order": CommaSeparated(Str()),
                    "trimmed_profile_length": Int(),
                    "untrimmed_profile_length": Int(),
                    "seg_threshold": Int(),
                    "measurement_order": Int(),
                    "measure_thickness": Float(),
                    "reference_wavelength": Str(),
                    "image_register": Int(),
                    "channel_register": Int(),
                    "population_register": Int(),
                    "trimmed_regions": MapPattern(Str(), CommaSeparated(Float())),
                    "untrimmed_regions": MapPattern(Str(), CommaSeparated(Float())),
                }
            ),
            "redox": Map(
                {
                    "ratio_numerator": Str(),
                    "ratio_denominator": Str(),
                    "r_min": Float(),
                    "r_max": Float(),
                    "instrument_factor": Float(),
                    "midpoint_potential": Float(),
                    "z": Int(),
                    "temperature": Float(),
                }
            ),
            "registration": Map(
                {
                    "n_deriv": Float(),
                    "warp_n_basis": Float(),
                    "warp_order": Float(),
                    "warp_lambda": Float(),
                    "smooth_lambda": Float(),
                    "smooth_n_breaks": Float(),
                    "smooth_order": Float(),
                    "rough_lambda": Float(),
                    "rough_n_breaks": Float(),
                    "rough_order": Float(),
                }
            ),
            "output": Map(
                {
                    "should_save_plots": Bool(),
                    "should_save_profile_data": Bool(),
                    "should_save_summary_data": Bool(),
                }
            ),
        }
    )

    seg_images: xr.DataArray = None
    rot_fl: xr.DataArray = None
    rot_seg: xr.DataArray = None

    midlines: xr.DataArray = None

    untrimmed_profiles: xr.DataArray = None
    trimmed_profiles: xr.DataArray = None

    channel_warps: xr.DataArray = None
    std_warps: xr.DataArray = None

    def __init__(self, exp_dir):
        self.experiment_dir = Path(exp_dir)
        self.settings_path = self.experiment_dir.joinpath("settings.yaml")
        try:
            with open(self.settings_path, "r") as f:
                self.config = load(f.read(), self.experiment_schema).data
        except YAMLError:
            raise ValueError("Incorrectly specified config file.")

        self.experiment_id = self.experiment_dir.stem

        # compute the filenames/paths for this experiment
        self.raw_img_stack_path = self.experiment_dir.joinpath(
            self.experiment_id + ".tif"
        )
        self.movement_path = self.experiment_dir.joinpath(
            self.experiment_id + "-mvmt.csv"
        )
        self.indexer_path = self.experiment_dir.joinpath(
            self.experiment_id + "-indexer.csv"
        )
        self.processed_images_dir = self.experiment_dir.joinpath("processed_images")
        self.rot_seg_dir = self.processed_images_dir.joinpath("rot_seg")
        self.rot_fl_dir = self.processed_images_dir.joinpath("rot_fl")
        self.fl_imgs_dir = self.processed_images_dir.joinpath("fluorescent_images")
        self.orig_images_path = self.processed_images_dir.joinpath("images.nc")
        self.seg_images_path = self.processed_images_dir.joinpath("seg_images.nc")
        self.aligned_images_path = self.processed_images_dir.joinpath(
            "aligned_images.nc"
        )
        self.aligned_seg_images_path = self.processed_images_dir.joinpath(
            "aligned_seg_images.nc"
        )

        self.strains = pio.load_strain_map_from_disk(self.indexer_path)

        # load images
        self.images = self._load_raw_images()

        # try to load masks
        try:
            self.load_masks()
        except IOError:
            logging.info("No masks found in experiment directory")
            pass

    # Computed Filepaths
    @property
    def fig_dir(self):
        return self.analysis_dir.joinpath("figs")

    @property
    def untrimmed_profile_data_path(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-untrimmed_profile_data.nc"
        )

    @property
    def trimmed_profile_data_path(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-trimmed_profile_data.nc"
        )

    @property
    def channel_warp_data_path(self):
        return self.analysis_dir.joinpath(self.experiment_id + "-channel_warps.nc")

    @property
    def std_warp_data_path(self):
        return self.analysis_dir.joinpath(self.experiment_id + "-std_warps.nc")

    @property
    def untrimmed_profile_data_csv_path(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-untrimmed_profile_data.csv"
        )

    @property
    def trimmed_profile_data_csv_path(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-trimmed_profile_data.csv"
        )

    @property
    def untrimmed_region_data_path(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-untrimmed_region_data.csv"
        )

    @property
    def trimmed_region_data_path(self):
        return self.analysis_dir.joinpath(
            self.experiment_id + "-trimmed_region_data.csv"
        )

    @property
    def analysis_dir(self) -> Path:
        date_str = datetime.datetime.now().strftime("%Y-%m-%d")
        strategy = self.config["pipeline"]["strategy"]
        if len(strategy) > 0:
            suffix = f"_{strategy}"
        else:
            suffix = ""
        analysis_dir_ = self.experiment_dir.joinpath(
            "analyses", utils.get_valid_filename(f"{date_str}{suffix}"),
        )
        # analysis_dir_.mkdir(parents=True, exist_ok=True)
        return analysis_dir_

    def _load_raw_images(self):
        """
        This returns the raw (non-median-subtracted) images
        """
        logging.info(f"Loading image data from {self.raw_img_stack_path}")
        raw_image_data = pio.load_tiff_as_hyperstack(
            img_stack_path=self.raw_img_stack_path,
            strain_map=self.strains,
            movement_path=self.movement_path,
            channel_order=self.config["pipeline"]["channel_order"],
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
        movement_path = self.experiment_dir.joinpath(self.experiment_id + "-mvmt.csv")
        try:
            df = pd.read_csv(movement_path)
            df = df.pivot_table(
                index="animal", columns=["region", "pair"], values="movement"
            )
            df = df.stack("pair")
            return df
        except FileNotFoundError:
            logging.warning(f"Tried to access {movement_path}; file was not found")
            return None

    def make_analysis_dir(self) -> None:
        logging.info(f"Making analysis directory at {self.analysis_dir}")
        self.analysis_dir.mkdir(parents=True, exist_ok=True)

    @property
    def trimmed_summary_table(self):
        df = profile_processing.summarize_over_regions(
            self.trimmed_profiles,
            regions=self.config["pipeline"]["trimmed_regions"],
            **self.config["redox"],
        )
        return df

    @property
    def untrimmed_summary_table(self):
        df = profile_processing.summarize_over_regions(
            self.untrimmed_profiles,
            regions=self.config["pipeline"]["untrimmed_regions"],
            **self.config["redox"],
        )
        return df

    ####################################################################################
    # PIPELINE
    ####################################################################################

    def full_pipeline(self):
        logging.info(f"Starting full pipeline run for {self.experiment_dir}")

        self.make_analysis_dir()

        logging.info(f"Saving fluorescent images to {self.fl_imgs_dir}")
        pio.save_images_xarray_to_tiffs(
            self.images, self.fl_imgs_dir, prefix=self.experiment_id
        )

        self.segment_pharynxes()
        self.register_images()
        self.align_and_center()
        self.calculate_midlines()
        self.measure_under_midlines()
        self.register_profiles()
        self.trim_data()
        self.calculate_redox()
        self.do_manual_ap_flips()
        self.persist_to_disk()

        logging.info(f"Finished full pipeline run for {self.experiment_dir}")

        return self

    def run_neuron_pipeline(self):
        logging.info(
            f"Starting full neuron analysis pipeline run for {self.experiment_dir}"
        )
        self.make_analysis_dir()
        df = ip.measure_under_labels(self.images, self.seg_images).reset_index()

        df.to_csv(self.analysis_dir / (self.experiment_id + "-neuron_analysis.csv"))

    def segment_pharynxes(self):
        if self.seg_images is not None:
            logging.info("masks have been specified. skipping mask generation")
            self.save_masks()
            return
        else:
            logging.info("Generating masks")
            self.seg_images = ip.segment_pharynxes(
                self.images, wvl=self.config["pipeline"]["reference_wavelength"]
            )
            self.save_masks()

    def register_images(self):
        if self.config["pipeline"]["image_register"]:
            logging.info("Registering Images")
            self.images = ip.register_all_images(self.images, self.seg_images)

    def align_and_center(self):
        logging.info("Centering and rotating pharynxes")
        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(
            self.images, self.seg_images,
        )

        logging.info(f"Saving rotated FL images to {self.aligned_images_path}")
        pio.save_profile_data(self.rot_fl, self.aligned_images_path)

        logging.info(f"Saving rotated masks to {self.aligned_seg_images_path}")
        pio.save_profile_data(self.rot_seg, self.aligned_seg_images_path)

    def calculate_midlines(self):
        logging.info("Calculating midlines")
        self.midlines = ip.calculate_midlines(self.rot_seg, degree=4)

    def measure_under_midlines(self):
        logging.info("Measuring under midlines")
        self.untrimmed_profiles = ip.measure_under_midlines(
            self.rot_fl,
            self.midlines,
            n_points=self.config["pipeline"]["untrimmed_profile_length"],
            order=self.config["pipeline"]["measurement_order"],
            thickness=float(self.config["pipeline"]["measure_thickness"]),
        )
        self.untrimmed_profiles = profile_processing.align_pa(self.untrimmed_profiles)
        self.untrimmed_profiles = self.add_experiment_metadata_to_data_array(
            self.untrimmed_profiles
        )

        # subtract the image medians from the profile data
        logging.info("Subtracting image medians from profile data")
        self.untrimmed_profiles = ip.subtract_medians(
            self.untrimmed_profiles, self.images
        )

    def register_profiles(self):

        if self.config["pipeline"]["population_register"]:
            logging.info("Standardizing profiles")
            (
                self.untrimmed_profiles,
                self.std_warps,
            ) = profile_processing.register_profiles_pop(
                self.untrimmed_profiles,
                self.config["redox"],
                **self.config["registration"],
            )

        if self.config["pipeline"]["channel_register"]:
            logging.info("Channel-Registering profiles")
            (
                self.untrimmed_profiles,
                self.channel_warps,
            ) = profile_processing.channel_register(
                self.untrimmed_profiles,
                redox_params=self.config["redox"],
                reg_params=self.config["registration"],
            )

    def trim_data(self):
        logging.info("Trimming intensity data")
        self.trimmed_profiles = self.add_experiment_metadata_to_data_array(
            profile_processing.trim_profiles(
                self.untrimmed_profiles,
                self.config["pipeline"]["seg_threshold"],
                ref_wvl=self.config["pipeline"]["reference_wavelength"],
            )
        )

    def calculate_redox(self):
        logging.info("Calculating redox measurements")

        redox_params = self.config["redox"]

        # Images
        self.images = utils.add_derived_wavelengths(self.images, **redox_params)
        self.rot_fl = utils.add_derived_wavelengths(self.rot_fl, **redox_params)

        # profiles
        self.trimmed_profiles = utils.add_derived_wavelengths(
            self.trimmed_profiles, **redox_params
        )

        self.untrimmed_profiles = utils.add_derived_wavelengths(
            self.untrimmed_profiles, **redox_params
        )

    def do_manual_ap_flips(self):
        # TODO finish implementation
        logging.info("skipping manual AP flips - not implemented")

    def flip_at(self, idx):
        # TODO use get_axis_num here and flip with xarray along correct dimension
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
            (self.images, self.orig_images_path),
            (self.rot_fl, self.aligned_images_path),
            (self.seg_images, self.seg_images_path),
            (self.rot_seg, self.aligned_seg_images_path),
        ]
        for img, path in imgs_paths:
            if img is not None:
                logging.info(f"Saving images to {path}")
                img.to_netcdf(path)

    # def load_tiff_as_hyperstack(self):
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
                for title, fig in plots.generate_wvl_pair_timepoint_profile_plots(data):
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
            imgs = utils.add_derived_wavelengths(self.images, **self.config["redox"])
            with PdfPages(mvmt_annotation_img_path) as pdf:
                for i in tqdm(range(self.images.animal.size)):
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
                                    wavelength=self.config["redox"]["ratio_numerator"],
                                    pair=pair,
                                    timepoint=tp,
                                )
                                / self.rot_fl.sel(
                                    wavelength=self.config["redox"][
                                        "ratio_denominator"
                                    ],
                                    pair=pair,
                                    timepoint=tp,
                                )
                            )[i]
                            I = self.rot_fl.sel(
                                wavelength=self.config["redox"]["ratio_numerator"],
                                pair=pair,
                                timepoint=tp,
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
                                *self.midlines.sel(pair=pair, timepoint=tp,)[i]
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
            f"Saving untrimmed profile data to {self.untrimmed_profile_data_path}"
        )
        pio.save_profile_data(self.untrimmed_profiles, self.untrimmed_profile_data_path)

        logging.info(f"Saving trimmed profile data to {self.trimmed_profile_data_path}")
        pio.save_profile_data(self.trimmed_profiles, self.trimmed_profile_data_path)

        # Warps, if necessary
        if self.config["pipeline"]["channel_register"]:
            logging.info(f"Saving channel warp data to {self.channel_warp_data_path}")
            self.channel_warps.to_netcdf(self.channel_warp_data_path)

        if self.config["pipeline"]["population_register"]:
            logging.info(f"Saving channel warp data to {self.std_warp_data_path}")
            self.std_warps.to_netcdf(self.std_warp_data_path)

        # Now, save the profile data in "tidy" format as CSV
        profile_processing.to_dataframe(self.trimmed_profiles, "value").to_csv(
            self.trimmed_profile_data_csv_path
        )
        profile_processing.to_dataframe(self.untrimmed_profiles, "value").to_csv(
            self.untrimmed_profile_data_csv_path
        )

    def save_summary_data(self):
        # Persist the region means
        logging.info(
            f"Saving untrimmed region means to {self.untrimmed_region_data_path}"
        )
        self.untrimmed_summary_table.to_csv(self.untrimmed_region_data_path)
        logging.info(f"Saving trimmed region means to {self.trimmed_region_data_path}")
        self.trimmed_summary_table.to_csv(self.trimmed_region_data_path)

    def save_masks(self):
        logging.info(f"writing masks to {self.seg_images_path}")
        pio.save_profile_data(self.seg_images, self.seg_images_path)
        logging.info(f"Saved masks to {self.seg_images_path}")

    def load_masks(self):
        self.seg_images = pio.load_profile_data(self.seg_images_path)
        logging.info(f"Loaded masks from {self.seg_images_path}")

    def persist_to_disk(self):
        logging.info(f"Saving {self.experiment_id} inside {self.experiment_dir}")

        if self.config["output"]["should_save_summary_data"]:
            self.save_summary_data()

        if self.config["output"]["should_save_profile_data"]:
            self.persist_profile_data()

        if self.config["output"]["should_save_plots"]:
            self.save_plots()

    ####################################################################################
    # MISC / HELPER
    ####################################################################################
    def add_experiment_metadata_to_data_array(self, data_array: xr.DataArray):
        params = {}
        params.update(self.config["pipeline"])
        params.update(self.config["redox"])
        params.update(self.config["registration"])

        to_remove = ["channel_order", "trimmed_regions", "untrimmed_regions"]
        for k in to_remove:
            del params[k]

        return data_array.assign_attrs(**params)


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
    import sys

    run_analysis(sys.argv[1], 1)
