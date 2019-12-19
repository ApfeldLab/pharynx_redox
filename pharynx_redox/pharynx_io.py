import re
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import xarray as xr
from skimage.external import tifffile

from . import utils


def load_profile_data(path: Union[Path, str]) -> xr.DataArray:
    logging.info("Loading data from %s" % path)
    return xr.load_dataarray(path)


def save_profile_data(
    profile_data: xr.DataArray, path: Union[Path, str]
) -> xr.DataArray:
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
    """

    # Each image has a "description" string (formatted as XML) associated with it

    # The keys here refer to the id property in the XML, each of which has a
    # corresponding value property

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
        return all_metadata


def load_tiff_from_disk(image_path: Path, return_metadata=False) -> np.ndarray:
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

    if return_metadata:
        return img_data, get_metadata_from_tiff(str(image_path))
    else:
        return img_data


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
    dir_path = Path(dir_path)
    dir_path.mkdir(parents=True, exist_ok=True)

    if prefix is not None:
        prefix = f"{prefix}-"
    if suffix is not None:
        suffix = f"-{suffix}"

    for pair in imgs.pair.data:
        for wvl in imgs.wavelength.data:
            final_path = dir_path.joinpath(
                f"{prefix}_wvl={wvl}_pair={pair}_{suffix}.tif"
            )
            if imgs.data.dtype == np.bool:
                data = np.uint8(imgs.sel(wavelength=wvl, pair=pair).data * 255)
            else:
                data = imgs.sel(wavelength=wvl, pair=pair).data

            tifffile.imsave(str(final_path), data)


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
    intercalated_image_stack_path: Path, strain_map: [str] = None, indexer_path=None
) -> xr.DataArray:
    """
    Loads the images specified by the path into an `xarray.DataArray <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.html#xarray-dataarray/>`_,
    organized by strain, wavelength, and pair.
    
    Parameters
    ----------
    intercalated_image_stack_path
        the path to the raw image stack
    imaging_scheme
        a string denoting the order of the wavelengths used during imaging.

        Transmitted light should be represented as ``TL``. The wavelengths should be
        separated by ``/``.

        E.g.::

            TL/470/410/470/410
    strain_map
        a list of strain names, corresponding to the strain of each animal. The length
        must therefore be the same as the number of animals imaged. Overridden by indexer_path. If None, indexer_path must be given.
    indexer_path
        if given, use the indexer to load the strain map instead of passing it explicity. Overrides the strain_map parameter. If None, strain_map must be given.

    Returns
    -------
    img_stack: xr.DataArray
        A multi-dimensional array of the form::

            (frame, wavelength, pair, height, width)

        The data is returned in the form of an `xarray.DataArray` object. This is very
        similar to (and indeed uses for its implementation) a numpy ndarray. The major
        difference (exploited by this code-base) is that dimensions may be accessed by
        labels in addition to the traditional index access.

        For example, all transmitted-light images for the strain HD233 may be accessed
        as follows::

            >> all_images = load_images(intercalated_image_stack_path, imaging_scheme, strains)
            >> all_images.data.shape
            (123, 3, 2, 130, 174)
            >> only_tl = all_images.sel(strain='HD233', wavelength='TL', pair=0)
            >> only_tl.data.shape
            (60, 130, 174)
    """
    imgdata, metadata = load_tiff_from_disk(
        intercalated_image_stack_path, return_metadata=True
    )

    if indexer_path:
        strain_map = load_strain_map_from_disk(indexer_path)

    # Organize the metadata
    df = pd.DataFrame(metadata)
    df["animal"] = df.groupby(["stage_x", "stage_y"], sort=False).ngroup(ascending=True)
    df["pair"] = df.groupby(["animal", "wavelength"]).cumcount()

    n_animals = len(df["animal"].unique())
    n_pairs = len(df["pair"].unique())
    n_wvls = len(df["wavelength"].unique())

    imgdata_reshaped = np.full(
        (n_animals, n_pairs, n_wvls, imgdata.shape[-2], imgdata.shape[-1]),
        np.nan,
        dtype=np.int16,
    )

    metadata_keys = [
        ("stage_x", np.int16),
        ("stage_y", np.int16),
        ("stage_z", np.int16),
        ("time", "datetime64[us]"),
        ("bin_x", np.uint8),
        ("bin_y", np.uint8),
        ("exposure", np.uint8),
    ]

    all_coords = {
        k: np.empty((n_animals, n_pairs, n_wvls), dtype=dtype)
        for k, dtype in metadata_keys
    }

    for animal in df["animal"].unique():
        for pair in df["pair"].unique():
            for wvl_idx, wvl in enumerate(df["wavelength"].unique()):
                try:
                    idx = df.index[
                        (df["animal"] == animal)
                        & (df["wavelength"] == wvl)
                        & (df["pair"] == pair)
                    ][0]
                    imgdata_reshaped[animal, pair, wvl_idx] = imgdata[idx]
                    for key, _ in metadata_keys:
                        all_coords[key][animal, pair, wvl_idx] = df.loc[idx][key]
                except IndexError:
                    # This happens when, for example, TL is only imaged once per pair
                    continue

    da = xr.DataArray(
        imgdata_reshaped,
        dims=["animal", "pair", "wavelength", "y", "x"],
        coords={
            "wavelength": df["wavelength"].unique(),
            "strain": ("animal", strain_map),
            "experiment_id": "2019-12-10_PD4793_ts",
            "time": (("animal", "pair", "wavelength"), all_coords["time"]),
            "stage_x": (("animal", "pair", "wavelength"), all_coords["stage_x"]),
            "stage_y": (("animal", "pair", "wavelength"), all_coords["stage_y"]),
            "stage_z": (("animal", "pair", "wavelength"), all_coords["stage_z"]),
            # "exposure": (
            #     ("animal", "pair", "wavelength"),
            #     all_coords["exposure"],
            #     dict(units="ms"),
            # ),
        },
    )

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
