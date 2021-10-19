import itertools
import logging
import pickle
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import pandera as pa
import tifffile
import xarray as xr
from skimage import io as skio


def get_image_metadata(
    img_path: Union[Path, str], acquisition_method: str
) -> pd.DataFrame:
    acquisition_method_map = {
        "acquire": get_image_metadata_metamorph_acquire,
        "mda": get_image_metadata_metamorph_mda,
    }

    return acquisition_method_map[acquisition_method](img_path)


def get_image_metadata_metamorph_mda(img_path: Union[Path, str]):
    raise NotImplementedError


def get_image_metadata_metamorph_acquire(img_path: Union[Path, str]) -> pd.DataFrame:
    """
    Return a DataFrame of the image stack's metdata where each row is indexed by the
    frame number. The columns come from MetaMorph's metadata XML as formatted when using
    the `Acquire` window in MetaMorph.

    Parameters
    ----------
    img_path
        the location of the image

    Returns
    -------
    metadata
        a DataFrame containing the image metadata, indexed by frame

    """
    type_lookup = {
        "string": str,
        "int": int,
        "bool": bool,
        "float": float,
        "guid": str,
        "colorref": str,
        "float-array": np.ndarray,
        "time": lambda t: datetime.strptime(t, "%Y%m%d %H:%M:%S.%f"),
    }
    metadata = defaultdict(list)
    with tifffile.TiffFile(str(img_path)) as tif:
        for page in tif.pages:
            xml_string = page.tags["ImageDescription"].value
            try:
                root = ET.fromstring(xml_string)
            except ET.ParseError:
                raise AttributeError(
                    "Invalid metadata in file. Might have been processed by ImageJ?"
                )

            for prop in root.findall("prop"):
                if prop.get("id") == "Description":
                    descr_strings = prop.get("value").splitlines()
                    for s in descr_strings:
                        try:
                            k, v = s.split(": ")
                            metadata[k].append(v)
                        except ValueError:
                            pass
                else:
                    metadata[prop.get("id")].append(
                        type_lookup[prop.get("type")](prop.get("value"))
                    )

            for x in root.find("PlaneInfo"):
                metadata[x.get("id")].append(type_lookup[x.get("type")](x.get("value")))

    metadata = pd.DataFrame(metadata)
    metadata.index.name = "frame"

    return metadata


def save_midlines(path: Union[Path, str], midlines: xr.DataArray):
    """
    Save the midlines to the given path

    Parameters
    ----------
    path
        the path to save the midlines
    midlines
        the midlines to save

    Returns
    -------
    None

    """
    logging.info(f"saving midlines to {path}")
    with open(path, "wb") as f:
        pickle.dump(midlines, f, pickle.HIGHEST_PROTOCOL)


def load_midlines(path: Union[Path, str]) -> xr.DataArray:
    with open(path, "rb") as f:
        midlines = pickle.load(f)
    logging.info(f"loaded midlines from {path}")
    return midlines


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


def save_images_xarray_to_tiffs(
    imgs: xr.DataArray,
    dir_path: str,
    prefix: str = None,
    suffix: str = None,
    z_axis: str = "animal",
) -> None:
    """
    Save the given image stack to disk inside the given directory, separated by each
    axis except for the given z-axis.

    For example, say you've imaged 60 animals, each with 3 timepoints and 2 pairs. In
    this case, you will most likely want to scroll through the animals in FIJI/etc. So
    you would use `z_axis='animal'` to output stacks with 60 frames - one frame per
    animal. You would have each combination of (timepoint, pair, wavelength) in its own
    stack.

    However, say you've imaged 6 animals, each with 100 timepoints and 2 pairs. In this
    case, you would use `z_axis='timepoint'`. This would result in each each stack
    getting 100 frames - one per timepoint. Each (animal, pair, wavelength) combination
     would get its own file.

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


def read_image_raw(
    image_path: Union[Path, str],
) -> np.ndarray:
    """
    Read an image from disk, returning a `np.ndarray`. Right now, just a wrapper around
    `scikit-image`'s `imread` function.

    Parameters
    ----------
    image_path
        the path to the image

    Returns
    -------
    image
        an nD image
    """
    return skio.imread(str(image_path))


def validate_frame_metadata(metadata: pd.DataFrame):

    schema = pa.DataFrameSchema(
        {
            "frame": pa.Column(pa.Int, pa.Check.greater_than_or_equal_to(0)),
            "animal": pa.Column(pa.Int, pa.Check.greater_than_or_equal_to(0)),
            "wavelength": pa.Column(pa.String, coerce=True),
            "timepoint": pa.Column(pa.Int),
            "pair": pa.Column(pa.Int, pa.Check.greater_than_or_equal_to(0)),
            "strain": pa.Column(pa.String),
            "notes": pa.Column(pa.String, nullable=True, coerce=True),
        }
    )

    metadata = schema.validate(metadata)
    groupby_dims = ["animal", "timepoint", "pair"]
    assert (
        metadata.groupby(groupby_dims)["strain"].nunique().max() == 1
    ), "more than one strain was given to a wavelength set"

    return metadata


def validate_movement_annotations(mvmt_tbl: pd.DataFrame) -> pd.DataFrame:
    mvmt_tbl_schema = pa.DataFrameSchema(
        {
            "animal": pa.Column(pa.Int, pa.Check.greater_than_or_equal_to(0)),
            "timepoint": pa.Column(pa.Int, pa.Check.greater_than_or_equal_to(0)),
            "pair": pa.Column(pa.Int, pa.Check.greater_than_or_equal_to(0)),
            "mvmt-*": pa.Column(pa.Int, regex=True),
        },
        strict=True,
    )

    return mvmt_tbl_schema.validate(mvmt_tbl)


def load_tiff_as_hyperstack(
    img_stack_path: Union[Path, str],
    manual_metadata: Union[Path, str, pd.DataFrame],
    mvmt_metadata: Union[Path, str] = None,
    z_dim: int = 0,
    y_dim: int = 1,
    x_dim: int = 2,
) -> xr.DataArray:
    """
    Convert a 3D image stack to a 6D "hyperstack" with the dimensions:

        animal, timepoint, pair, wavelength, y, x

    If movement metadata is supplied, it will be added to the resultant `xr.DataArray`
    as coordinates.

    Parameters
    ----------
    img_stack_path
        the path to a 3D array on disk. If a path is given, attempts to open the given
        path into a 3D numpy array. Can open the following formats:

            .tiff, .tif, .stk

    manual_metadata
        the path to the metadata table (CSV). The table should have the following
        columns:

            frame, animal, wavelength, pair, timepoint
    mvmt_metadata
        the path to a movement metadata table (CSV). The table should have the following
        columns:

            animal, pair, timepoint, mvmt-*

        there may be an arbitrary number of `mvmt-*` columns, where `*` is replaced with
        a signifier (e.g. `posterior`).

    z_dim
        the dimension of the z-axis in the given 3D stack
    y_dim
        the dimension of the y-axis in the given 3D stack
    x_dim
        the dimension of the x-axis in the given 3D stack

    Returns
    -------
    hyperstack
        a 6D xarray, with the dimensions:

            animal, timepoint, pair, wavelength, y, x

    """
    img_data_3d = read_image_raw(img_stack_path)

    if isinstance(manual_metadata, (str, Path)):
        manual_metadata = pd.read_csv(manual_metadata, dtype={"notes": str})

    manual_metadata = validate_frame_metadata(manual_metadata)

    # metadata_cols = ["strain"]
    # try:
    #     image_metadata = get_image_metadata_metamorph_acquire(img_stack_path)
    #
    #     frame_metadata = (
    #         manual_metadata.set_index("frame")
    #         .join(image_metadata, rsuffix="-auto_detected")
    #         .reset_index()
    #     )
    # except ValueError:
    #     frame_metadata = manual_metadata

    # indices are 0-indexed
    n_animals = manual_metadata["animal"].max() + 1
    n_timepoints = manual_metadata["timepoint"].max() + 1
    n_pairs = manual_metadata["pair"].max() + 1
    n_wvls = len(manual_metadata["wavelength"].unique())

    # Set up empty hyperstack
    dim_order = ["animal", "timepoint", "pair", "wavelength", "y", "x"]
    da = xr.DataArray(
        np.zeros(
            (
                n_animals,
                n_timepoints,
                n_pairs,
                n_wvls,
                img_data_3d.shape[y_dim],
                img_data_3d.shape[x_dim],
            ),
            dtype=np.float64,
        ),
        dims=dim_order,
        coords={
            "wavelength": manual_metadata.wavelength.unique(),
            "strain": ("animal", manual_metadata.groupby("animal").strain.first()),
            "exp_animal": ("animal", manual_metadata.animal.unique()),
        },
    )
    da = da.assign_coords(
        {"frame": manual_metadata.set_index(dim_order[:-2])["frame"].to_xarray()}
    )

    # METADATA
    if mvmt_metadata is not None:
        mvmt_metadata = Path(mvmt_metadata)
        if mvmt_metadata.exists():
            mvmt_metadata = validate_movement_annotations(pd.read_csv(mvmt_metadata))
            mvmt_xarray = mvmt_metadata.set_index(
                ["animal", "timepoint", "pair"]
            ).to_xarray()
            da = da.assign_coords(mvmt_xarray)

    # add image data to hyperstack
    # this is the slowest piece of the function, but can't think of a better way to do it
    for frame in range(img_data_3d.shape[z_dim]):
        indexer = dict(
            manual_metadata.iloc[frame][["animal", "wavelength", "pair", "timepoint"]]
        )
        da.loc[indexer] = img_data_3d[frame]

    return da
