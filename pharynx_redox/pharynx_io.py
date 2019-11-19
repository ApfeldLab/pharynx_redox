from pathlib import Path
from typing import Union
import numpy as np
import pandas as pd
import xarray as xr
from skimage.external import tifffile

from pharynx_redox import utils


def load_profile_data(profile_data_path: Union[Path, str]) -> xr.DataArray:
    try:
        return xr.load_dataarray(profile_data_path).set_index(
            spec=["experiment", "animal"]
        )
    except ValueError:
        return xr.load_dataarray(profile_data_path)


def save_profile_data(
    profile_data: xr.DataArray, profile_data_path: Union[Path, str]
) -> xr.DataArray:
    profile_data.reset_index("spec").to_netcdf(profile_data_path)


def load_tiff_from_disk(image_path: Path) -> np.ndarray:
    """
    Load a tiff file from disk

    Parameters
    ----------
    image_path: Path
        the path to the image file

    Returns
    -------
    img: np.ndarray of shape
        the image stack as a numpy array with the following dimensions::

            (n_images, height, width)

    """
    return tifffile.imread(str(image_path))


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
    intercalated_image_stack_path: Path,
    imaging_scheme: str,
    strain_map: [str] = None,
    indexer_path=None,
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
    if indexer_path is not None:
        strain_map = load_strain_map_from_disk(indexer_path)

    intercalated_image_stack = load_tiff_from_disk(intercalated_image_stack_path)
    lambdas_with_counts = process_imaging_scheme_str(imaging_scheme, "/")
    lambdas = np.array([l[0] for l in lambdas_with_counts])
    pairs = [l[1] for l in lambdas_with_counts]
    unique_lambdas = [
        lambdas[idx] for idx in sorted(np.unique(lambdas, return_index=True)[1])
    ]
    n_animals = intercalated_image_stack.shape[0] // len(lambdas)
    img_height = intercalated_image_stack.shape[1]
    img_width = intercalated_image_stack.shape[2]

    tmp_img_stack = np.reshape(
        intercalated_image_stack, (n_animals, len(lambdas), img_height, img_width)
    )

    reshaped_img_stack = np.empty(
        (
            n_animals,
            len(np.unique(lambdas)),
            np.max(pairs) + 1,
            intercalated_image_stack.shape[1],
            intercalated_image_stack.shape[2],
        ),
        dtype=intercalated_image_stack.dtype,
    )

    for i, (wvl, pair) in enumerate(lambdas_with_counts):
        j = np.where(lambdas == wvl)[0][0]
        reshaped_img_stack[:, j, pair, :, :] = tmp_img_stack[:, i, :, :]

    exp_id = intercalated_image_stack_path.stem

    spec = pd.MultiIndex.from_arrays(
        [np.repeat(exp_id, n_animals), np.arange(n_animals)],
        names=("experiment", "animal"),
    )
    return xr.DataArray(
        reshaped_img_stack,
        dims=["spec", "wavelength", "pair", "y", "x"],
        coords={
            "wavelength": unique_lambdas,
            "spec": spec,
            "strain": ("spec", strain_map),
        },
    )


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


def load_all_rot_fl() -> xr.DataArray:
    """
    TODO: Documentation

    Returns
    -------

    """
    return xr.load_dataarray("../data/paired_ratio/all_rot_fl.nc")


def load_all_rot_seg() -> xr.DataArray:
    """
    TODO: Documentation

    Returns
    -------

    """
    return xr.load_dataarray("../data/paired_ratio/all_rot_seg.nc")
