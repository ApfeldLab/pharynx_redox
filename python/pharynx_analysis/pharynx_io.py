import numpy as np
import pandas as pd
from skimage import io as sk_io
import xarray as xr


def load_tiff_from_disk(image_path: str) -> np.ndarray:
    """ Load a tiff stack from disk.

    :param image_path: path to the image stack
    :return: a numpy array with dimensions (frame, height, width)
    """
    return sk_io.imread(image_path)


def load_images(intercalated_image_stack_path: str, imaging_scheme: str, strain_map: [str]) -> xr.DataArray:
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
    implementation) a numpy ndarray. The major difference (exploited by this code-base; in fact there are many detailed
    in the xarray documentation) is that dimensions may be accessed by labels in addition to the traditional index
    access.

    The data is organized as a 4-d matrix according to the following structure:
        (frame, wavelength, height, width)

    For example, all transmitted-light images for the strain HD233 may be accessed as follows:

        >> all_images = load_images(intercalated_image_stack_path, imaging_scheme, strains)
        >> all_images.data.shape
        (123, 5, 130, 174)
        >> only_tl = all_images.sel(strain='HD233', wavelength='TL')
        >> only_tl.data.shape
        (60, 130, 174)


    :param intercalated_image_stack_path: the path to the raw image stack
    :param imaging_scheme: a string denoting the order of the wavelengths used during imaging
    :param strain_map: a list of strain names, corresponding to the strain of each animal
    :return: An `xarray.DataArray` with the form (frame, wavelength, y, x)

    """
    intercalated_image_stack = load_tiff_from_disk(intercalated_image_stack_path)
    lambdas = imaging_scheme.split("/")
    n_animals = intercalated_image_stack.shape[0] // len(lambdas)

    reshaped_img_stack = np.reshape(
        intercalated_image_stack,
        (n_animals, len(lambdas), intercalated_image_stack.shape[1], intercalated_image_stack.shape[2]))

    return xr.DataArray(reshaped_img_stack,
                        dims=['strain', 'wavelength', 'y', 'x'],
                        coords={'wavelength': lambdas, 'strain': strain_map})


def load_strain_map(strain_map_path: str) -> np.ndarray:
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
        [np.repeat(x.Strain, x.End_Animal - x.Start_Animal + 1) for x in strain_map_df.itertuples()]).flatten()
