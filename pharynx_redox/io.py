import re
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Union
import os

import numpy as np
import pandas as pd
import xarray as xr
from skimage.external import tifffile
from sklearn.cluster import KMeans

from pharynx_redox import utils


def load_profile_data(path: Union[Path, str]) -> xr.DataArray:
    logging.info("Loading data from %s" % path)

    # data = xr.load_dataarray(path).set_index(animal=["experiment_id", "animal"])
    data = xr.load_dataarray(path)

    try:
        data = data.rename({"animal_": "animal"})
    except ValueError:
        pass

    return data


def save_profile_data(
    profile_data: xr.DataArray, path: Union[Path, str]
) -> xr.DataArray:
    logging.info("Saving data to %s" % path)
    try:
        profile_data.reset_index("animal").to_netcdf(path)
    except KeyError:
        profile_data.to_netcdf(path)


def _parse_illum_setting(ilum_setting: str) -> str:
    """
    Extract wavelength from `Illumination Setting` string from MetaMorph.
    
    Parameters
    ----------
    ilum_setting : str
        the `Illumination Setting` string found as the value in the property tag with
        the "_IllumSetting_" id.
    
    Returns
    -------
    str
        If there is a number in the string, that number is returned. If `transmitted
        light` is in the string, "TL" is returned. Otherwise, the string itself is 
        returned unchanged.
    """
    try:
        return re.search(r"(\d+)", ilum_setting).group(0)
    except AttributeError:
        if "transmitted light" in ilum_setting.lower():
            return "TL"
        else:
            return ilum_setting


def get_metadata_from_tiff(image_path: Path) -> List[Dict]:
    """
    Get metadata from each plane of a TIFF image stack.
    
    Parameters
    ----------
    image_path : Path
        the path to the TIFF image stack
    
    Returns
    -------
    List[Dict]
        a list where each entry is a dictionary containing the metadata for that image
    """

    # Each image has a "description" string (formatted as XML) associated with it

    # The keys (index 0) here refer to the id property in the XML, each of which has a
    # corresponding value property

    # The labels (index 1) refer to the column name that the values will be stored under

    # the functions are how we should process the string stored in the "value" property
    # for example, ``lambda x: x`` simply returns the string itself.
    metadata_keyfuncs = [
        ("_IllumSetting_", "wavelength", _parse_illum_setting),
        (
            "Exposure Time",
            "exposure",
            lambda x: int(re.search(r"(\d+) ms", x).group(1)),
        ),
        ("Prior Stage X", "stage_x", int),
        ("Prior Stage Y", "stage_y", int),
        ("Prior Z", "stage_z", int),
        (
            "acquisition-time-local",
            "time",
            lambda x: datetime.strptime(x, "%Y%m%d %H:%M:%S.%f"),
        ),
        ("camera-binning-x", "bin_x", int),
        ("camera-binning-y", "bin_y", int),
        # TODO: gain
    ]

    with tifffile.TiffFile(str(image_path)) as tif:
        all_metadata = []
        for page in tif.pages:
            descr = str(page.tags["image_description"].value, "utf-8")
            metadata = {}
            for key, label, fn in metadata_keyfuncs:
                try:
                    val = fn(re.search(rf'id="{key}".*value="(.*)"', descr).group(1))
                except AttributeError:
                    raise AttributeError(
                        f"Could not find key[{key}] in plane description\n{descr}"
                    )
                metadata[label] = val
            all_metadata.append(metadata)

    md = {}
    md["data"] = all_metadata
    md["types"] = {
        "stage_x": np.int16,
        "stage_y": np.int16,
        "stage_z": np.int16,
        "time": "datetime64[us]",
        "bin_x": np.uint8,
        "bin_y": np.uint8,
        "exposure": np.uint8,
    }
    return md


def load_tiff_from_disk(image_path: Path) -> np.ndarray:
    """
    Load a tiff file from disk

    Parameters
    ----------
    image_path: Path
        the path to the image file

    Returns
    -------
    img_data: np.ndarray of shape
        the image stack as a numpy array with the following dimensions::

            (n_images, height, width)
    metadata: List[Dict]
        the metadata associated with each image. Only returned if return_metadata=True

    """
    img_data = tifffile.imread(str(image_path))

    try:
        metadata = get_metadata_from_tiff(str(image_path))
    except AttributeError:
        metadata = None

    return img_data, metadata


def save_images_xarray_to_disk(
    imgs: xr.DataArray, dir_path: str, prefix: str = None, suffix: str = None
) -> None:
    """
    Save the given image stack to disk inside the given directory, separated by
    wavelength and pair

    Parameters
    ----------
    imgs
        the images to save to disk
    dir_path
        the directory to save the images in;
    prefix: optional
        a string with which to prepend the filenames
    suffix: optional
        a string with which to append the filenames

    Returns
    -------
    """
    # TODO: test
    dir_path = Path(dir_path)
    dir_path.mkdir(parents=True, exist_ok=True)

    prefix = f"{prefix}-" if prefix else ""
    suffix = f"-{suffix}" if suffix else ""

    for pair in imgs.pair.data:
        for wvl in imgs.wavelength.data:
            final_path = dir_path.joinpath(f"{prefix}wvl={wvl}_pair={pair}{suffix}.tif")
            if imgs.data.dtype == np.bool:
                data = np.uint8(imgs.sel(wavelength=wvl, pair=pair).data * 255)
            else:
                data = imgs.sel(wavelength=wvl, pair=pair).data

            tifffile.imsave(str(final_path), data)


def load_and_restack_img_set(
    dir_path: Union[str, Path], template_img_stack: xr.DataArray
):
    """
    Given a directory containing multiple .TIFF files, each labeled with wavelength
    and pair, load them into an xarray in memory.
    
    Parameters
    ----------
    dir_path : Union[str, Path]
        the directory containing the images. Must contain *only* those images to load
        into the xarray.
    """
    # TODO: test
    new_array = template_img_stack.copy()

    for fn in os.listdir(dir_path):
        try:
            fn = dir_path.joinpath(fn)
            wvl = re.search("wvl=(.+)_pair", str(fn)).group(1)
            pair = int(re.search("pair=(\d+)_?", str(fn)).group(1))

            img = load_tiff_from_disk(fn)
            new_array.loc[dict(wavelength=wvl, pair=pair)] = img
        except AttributeError:
            logging.warn(
                f"Encountered non-standard file in segmentation directory: {fn}"
            )
            continue

    return new_array


def process_imaging_scheme_str(imaging_scheme: str, delimiter="/") -> [(str, int)]:
    """
    Split the imaging scheme string by the given delimiter, and return
    [(wavelength, nth_occurrence), ...]

    Examples
    --------
    >>> process_imaging_scheme_str("TL/470/410/470/410", delimiter='/')
    [("TL", 0), ("470", 0), ("410", 0), ("470", 1), ("410", 1)]

    Parameters
    ----------
    imaging_scheme
        A string of wavelengths which indicate the order in which images were taken, separated by the delimiter
    delimiter
        a string which separates the wavelengths in the imaging scheme string

    Returns
    -------
    list
        [(wavelength, nth_occurrence), ...]
    """
    return utils.create_occurrence_count_tuples(imaging_scheme.split(delimiter))


def load_images(
    img_stack_path: Path,
    channel_order: [str],
    strain_map: [str] = None,
    indexer_path: Path = None,
    movement_path: Path = None,
    dtype=np.uint16,
) -> xr.DataArray:
    # Check Arguments
    if (indexer_path is None) and (strain_map is None):
        raise ValueError(
            "either `indexer_path` or `strain_map` must be supplied. Neither was."
        )
    if (indexer_path is not None) and (strain_map is not None):
        raise ValueError(
            "Only one of `indexer_path` or `strain_map` may be given. Both were."
        )

    if indexer_path:
        strain_map = load_strain_map_from_disk(indexer_path)

    # Load Data
    imgs, metadata = load_tiff_from_disk(img_stack_path)

    # Calculate dimensions of hyperstack
    wvls, pairs = zip(*utils.create_occurrence_count_tuples(channel_order))
    n_frames = imgs.shape[0]
    n_frames_per_animal = len(channel_order)
    n_wvls = len(np.unique(channel_order))
    n_animals = len(strain_map)
    n_pairs = np.max(pairs) + 1
    n_timepoints = int(n_frames / (n_frames_per_animal * n_animals))

    # Set up empty hyperstack
    da = xr.DataArray(
        np.full(
            (n_animals, n_timepoints, n_pairs, n_wvls, imgs.shape[-2], imgs.shape[-1]),
            np.nan,
            dtype=dtype,
        ),
        dims=["animal", "timepoint", "pair", "wavelength", "y", "x"],
        coords={
            "wavelength": np.unique(channel_order),
            "strain": ("animal", strain_map),
            "animal": ("animal", np.arange(n_animals, dtype=np.uint16)),
        },
    )

    # Set up empty metadata coords if necessary
    if metadata:
        metadata_df = pd.DataFrame(metadata["data"])
        all_coords = {
            metadata_key: np.zeros(
                (n_animals, n_timepoints, n_pairs, n_wvls), dtype=metadata_type
            )
            for metadata_key, metadata_type in metadata["types"].items()
        }

    # fill in the hyperstack
    frame = 0
    for timepoint in range(n_timepoints):
        for animal in range(n_animals):
            for wvl, pair in utils.create_occurrence_count_tuples(channel_order):
                # Assign image data from current frame to correct index
                da.loc[
                    dict(animal=animal, timepoint=timepoint, pair=pair, wavelength=wvl)
                ] = imgs[frame]

                if metadata:
                    # build the metadata coords for current frame using metadata DF
                    wvl_idx = np.where(np.unique(channel_order) == wvl)[0][0]
                    for key in metadata["types"].keys():
                        all_coords[key][
                            animal, timepoint, pair, wvl_idx
                        ] = metadata_df.loc[frame][key]

                frame += 1

    if metadata:
        da = da.assign_coords(
            {
                k: (("animal", "timepoint", "pair", "wavelength"), v)
                for k, v in all_coords.items()
            }
        )
        logging.info(f"Assigned coordinates to image data ({list(all_coords.keys())})")

    # Now, assign movement annotations to coordinates if possible
    default_mvmt_regions = ["posterior", "anterior", "sides_of_tip", "tip"]
    default_mvmt_metadata = {
        r: np.zeros((da.animal.size, da.pair.size, da.timepoint.size), dtype=np.uint8)
        for r in default_mvmt_regions
    }

    try:
        mvmt = pd.read_csv(movement_path)

        regions = list(mvmt.columns.drop(["experiment", "animal", "pair", "timepoint"]))
        try:
            regions.remove("notes")
        except ValueError:
            logging.info("no notes in movement file")

        mvmt_metadata = {
            r: np.zeros(
                (da.animal.size, da.pair.size, da.timepoint.size), dtype=np.uint8
            )
            for r in regions
        }

        logging.info("Adding movement annotations to image data")
        for animal in mvmt.animal.unique():
            for pair in mvmt.pair.unique():
                for timepoint in mvmt.timepoint.unique():
                    for region in regions:
                        idx = mvmt.index[
                            (mvmt["animal"] == animal)
                            & (mvmt["pair"] == pair)
                            & (mvmt["timepoint"] == timepoint)
                        ]
                        mvmt_metadata[region][animal, pair, timepoint] = mvmt.loc[idx][
                            region
                        ].values[0]
        logging.info(f"Loaded movement file from {movement_path}")

    except IOError:
        logging.warn(
            f"A movement file path ({movement_path}) was supplied but could not be found. All movement ({default_mvmt_regions}) assigned to 0."
        )
        mvmt_metadata = default_mvmt_metadata

    except ValueError:
        logging.info(
            f"No movement file supplied. All movement ({default_mvmt_regions}) assigned to 0."
        )
        mvmt_metadata = default_mvmt_metadata

    mvmt_coords = {
        f"mvmt-{r}": (("animal", "pair", "timepoint"), mvmt_labels)
        for r, mvmt_labels in mvmt_metadata.items()
    }

    da = da.assign_coords(mvmt_coords)

    return da


def save_split_images_to_disk(images: xr.DataArray, prefix: str, dir_path: str) -> None:
    """
    Save the given image container to disk, splitting the images by wavelength and pair.

    Parameters
    ----------
    images
        The DataArray containing the images
    prefix
        The name that will be used to save the images
    dir_path
        The directory to save the images in

    Returns
    -------
    None

    """
    dir_path = Path(dir_path)
    dir_path.mkdir(parents=True, exist_ok=True)
    for wvl in images.wavelength.data:
        for pair in images.pair.data:
            tifffile.imsave(
                str(dir_path.joinpath(f"{prefix}-{wvl}-{pair}.tif").absolute()),
                images.sel(wavelength=wvl, pair=pair).data,
            )


def load_strain_map_from_disk(strain_map_path: Path) -> np.ndarray:
    """
    Load strain map from disk, generate a 1D array where the index corresponds to the
    strain of the worm at that index

    Parameters
    ----------
    strain_map_path
        the path to strain map

        The strain map is a CSV file which tells the system which maps animal number to
        strain. It has three columns: ``Strain, Start_Animal, and End_Animal``.

        A valid strain map might look like this::

            Strain,Start_Animal,End_Animal
            HD233,1,60
            SAY47,61,123

    Returns
    -------
    numpy.ndarray
        an array of strings where the index corresponds to the strain of the worm at
        that index
    """

    strain_map_df = pd.read_csv(strain_map_path)
    strain_map_df.rename(columns=lambda x: x.strip(), inplace=True)

    # Check that all numbers are unique
    idx = strain_map_df[["start_animal", "end_animal"]].values.ravel()
    if len(idx) != len(np.unique(idx)):
        raise AttributeError("All indices in the indexer must be unique.")

    return np.concatenate(
        [
            np.repeat(x.strain, x.end_animal - x.start_animal + 1)
            for x in strain_map_df.itertuples()
        ]
    ).flatten()


def load_all_rot_fl(meta_dir: Path) -> xr.DataArray:
    """
    TODO: Documentation

    Returns
    -------

    """
    return xr.concat(
        [load_profile_data(p) for p in meta_dir.glob("**/*/*_rot_fl.nc")], dim="spec"
    )


def load_all_rot_seg() -> xr.DataArray:
    """
    TODO: Documentation

    Returns
    -------

    """
    return xr.load_dataarray("../data/paired_ratio/all_rot_seg.nc")


if __name__ == "__main__":
    load_strain_map_from_disk(
        "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47-indexer.csv"
    )
