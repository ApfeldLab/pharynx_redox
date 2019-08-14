from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from skimage import io as sk_io
from skimage.external import tifffile

from pharynx_analysis import utils


def load_tiff_from_disk(image_path: Path) -> np.ndarray:
    """ Load a tiff stack from disk.

    :param image_path: path to the image stack
    :return: a numpy array with dimensions (frame, height, width)
    """
    return sk_io.imread(str(image_path))


def save_images_xarray_to_disk(
    imgs: xr.DataArray, dir_path: str, stem: str, suffix: str
):
    dir_path = Path(dir_path)
    dir_path.mkdir(parents=True, exist_ok=True)
    for pair in imgs.pair.data:
        for wvl in imgs.wavelength.data:
            final_path = dir_path.joinpath(f"{stem}-{wvl}-{pair}-{suffix}.tif")
            if imgs.data.dtype == np.bool:
                data = np.uint8(imgs.sel(wavelength=wvl, pair=pair).data * 255)
            else:
                data = imgs.sel(wavelength=wvl, pair=pair).data

            tifffile.imsave(str(final_path), data)


def process_imaging_scheme_str(imaging_scheme_str: str, delimiter="/") -> [(str, int)]:
    """Split the imaging scheme string by the given delimiter, and return [(wavelength, nth_occurrence), ...]

    :param imaging_scheme_str: A string of wavelengths which indicate the order in which images were taken, separated by the delimiter
    :param delimiter: a string which separates the wavelengths in the imaging scheme string
    :return: a list [(wavelength, nth_occurrence), ...]
    """
    return utils.create_occurrence_count_tuples(imaging_scheme_str.split(delimiter))


def load_images(
    intercalated_image_stack_path: str, imaging_scheme: str, strain_map: [str]
) -> xr.DataArray:
    """ Loads the images specified by the path into an xarray.DataArray, organized by strain and wavelength.
    
    The wavelengths are split up according to the `imaging_scheme`, which specifies the order of the wavelengths taken
    during imaging. Transmitted light should be represented as "TL". The wavelengths should be separated by
    a single "/".

    Example `imaging_scheme`: "TL/470/410/470".

    The resultant stack is split up by strain. This information is provided by the `strains` parameter, which should be
    a list containing strings corresponding to the strain of each animal. Thus, the length of the `strains` list must
    be the same as the number of animals imaged.


    HOW TO ACCESS THE RESULTANT DATA:

    The data is returned in the form of an `xarray.DataArray` object. This is very similar to (and indeed uses for its
    implementation) a numpy ndarray. The major difference (exploited by this code-base) is that dimensions may be
    accessed by labels in addition to the traditional index access.

    The data is organized as a 5-d matrix according to the following structure:
        (frame, wavelength, pair, height, width)

    For example, all transmitted-light images for the strain HD233 may be accessed as follows:

        >> all_images = load_images(intercalated_image_stack_path, imaging_scheme, strains)
        >> all_images.data.shape
        (123, 3, 2, 130, 174)
        >> only_tl = all_images.sel(strain='HD233', wavelength='TL', pair=0)
        >> only_tl.data.shape
        (60, 130, 174)


    :param intercalated_image_stack_path: the path to the raw image stack
    :param imaging_scheme: a string denoting the order of the wavelengths used during imaging
    :param strain_map: a list of strain names, corresponding to the strain of each animal
    :return: An `xarray.DataArray` with the form (frame, wavelength, y, x)

    """
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

    return xr.DataArray(
        reshaped_img_stack,
        dims=["strain", "wavelength", "pair", "y", "x"],
        coords={"wavelength": unique_lambdas, "strain": strain_map},
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
    """ Load strain map from disk, generate a 1D array where the index corresponds to the strain of the worm at that index

    The strain map is a CSV file which tells the system which maps animal number to strain. It has three columns:
    Strain, Start_Animal, and End_Animal.

    A valid strain map might look like this:
        Strain,Start_Animal,End_Animal
        HD233,1,60
        SAY47,61,123

    :param strain_map_path: path to strain map
    :return: a numpy array of strings where the index corresponds to the strain of the worm at that index
    """
    strain_map_df = pd.read_csv(strain_map_path)
    return np.concatenate(
        [
            np.repeat(x.strain, x.end_animal - x.start_animal + 1)
            for x in strain_map_df.itertuples()
        ]
    ).flatten()
