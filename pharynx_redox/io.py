import re
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Union
import os
import itertools

import numpy as np
import pandas as pd
import xarray as xr
from skimage.external import tifffile
from sklearn.cluster import KMeans

from pharynx_redox import utils


def load_profile_data(path: Union[Path, str]) -> xr.DataArray:
    """Load the data at the given path into an xr.DataArray.

    Right now, this is just a wrapper around `xr.load_dataarray` but I made it just so we'd be using a standarzied IO
    interface for the data.

    Parameters
    ----------
    path : Union[Path, str]
        path to the profile data

    Returns
    -------
    xr.DataArray
        the data
    """
    data = xr.load_dataarray(path)
    logging.info("Loaded data from %s" % path)
    return data


def save_profile_data(profile_data: xr.DataArray, path: Union[Path, str]) -> None:
    """Save the data to disk.

    Right now, this is just a wrapper around `xr.to_netcdf` but I made it just so we'd be using a standarzied IO
    interface for the data.

    Parameters
    ----------
    profile_data: xr.DataArray
        the data to save
    path : Union[Path, str]
        the path to save the data at
    """
    logging.info("Saving data to %s" % path)
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
    
    Raises
    ------
    AttributeError
        The given image file does not contain metadata in the expected format
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
        "time": "datetime64[ns]",
        "bin_x": np.uint8,
        "bin_y": np.uint8,
        "exposure": np.uint8,
    }
    return md


def load_tiff_from_disk(image_path: Path) -> (np.ndarray, dict):
    """
    Load a tiff file from disk

    Parameters
    ----------
    image_path: Path
        the path to the image file

    Returns
    -------
    img_data: np.ndarray
        the image stack as a numpy array with the following dimensions::

            (n_images, height, width)

    metadata: dict
        the metadata associated with each image. Only returned if return_metadata=True

    """
    img_data = tifffile.imread(str(image_path))

    try:
        metadata = get_metadata_from_tiff(str(image_path))
    except AttributeError:
        metadata = None

    return img_data, metadata


def save_images_xarray_to_disk(
    imgs: xr.DataArray,
    dir_path: str,
    prefix: str = None,
    suffix: str = None,
    z_axis: str = "animal",
) -> None:
    """
    Save the given image stack to disk inside the given directory, separated by each axis except for
    the given z-axis.

    For example, say you've imaged 60 animals, each with 3 timepoints and 2 pairs. In this case,
    you will most likely want to scroll through the animals in FIJI/etc. So you would use `z_axis='animal'`
    to output stacks with 60 frames - one frame per animal. You would have each combination of (timepoint, pair, wavelength)
    in its own stack.

    However, say you've imaged 6 animals, each with 100 timepoints and 2 pairs. In this case, you would
    use `z_axis='timepoint'`. This would result in each each stack getting 100 frames - one per timepoint. Each
    (animal, pair, wavelength) combination would get its own file.

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
    z_axis:


    Returns
    -------
    """
    # TODO: test
    dir_path = Path(dir_path)
    dir_path.mkdir(parents=True, exist_ok=True)

    prefix = f"{prefix}-" if prefix else ""
    suffix = f"-{suffix}" if suffix else ""

    # remove the z-axis from the dimensons over which we will iterate
    dims = list(imgs.dims)
    for d in [z_axis, "x", "y"]:
        dims.remove(d)

    # Build a list of dictionaries with each combination of unique values for the remaining dimensions
    selector_values = list(
        itertools.product(*[np.unique(imgs[dim]).tolist() for dim in dims])
    )
    selectors = [dict(zip(dims, vals)) for vals in selector_values]

    for selector in selectors:
        selector_str = "-".join([f"{k}={v}" for k, v in selector.items()])
        path = dir_path.joinpath(f"{prefix}{selector_str}{suffix}.tiff")
        data = imgs.sel(**selector).values

        if data.dtype == np.bool:
            data = data.astype(np.uint8) * 255

        tifffile.imsave(str(path), data)


def load_images(
    img_stack_path: Union[Path, str],
    channel_order: [str],
    strain_map: [str] = None,
    indexer_path: Path = None,
    movement_path: Path = None,
    dtype=np.uint16,
) -> xr.DataArray:

    img_stack_path = Path(img_stack_path)

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
            "exp_animal": ("animal", np.arange(n_animals, dtype=np.uint16)),
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
        if movement_path is None:
            logging.info(
                f"No movement path specified. All movement ({default_mvmt_regions}) assigned to 0."
            )
        else:
            logging.warn(
                f"Movement file incorrectly specified. All movement ({default_mvmt_regions}) assigned to 0."
            )
        mvmt_metadata = default_mvmt_metadata

    mvmt_coords = {
        f"mvmt-{r}": (("animal", "pair", "timepoint"), mvmt_labels)
        for r, mvmt_labels in mvmt_metadata.items()
    }

    da = da.assign_coords(mvmt_coords)
    da.name = img_stack_path.stem

    return da


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
