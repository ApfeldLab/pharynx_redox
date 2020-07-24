import logging
import os
import re
import subprocess
import sys
import typing
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage as ndi
from skimage.measure import label, regionprops

from pharedox import profile_processing as pp


def midlines_xarray_to_napari(midlines: xr.DataArray, n: int) -> None:
    """
    Convert the xarray midline storage to the format that napari expects

    Parameters
    ----------
    midlines
    n

    Returns
    -------

    """

    def get_midline_pts(midline, n=10):
        y, x = midline.linspace(n)
        return np.stack([x, y]).T

    lines = xr.apply_ufunc(
        get_midline_pts,
        midlines,
        n,
        vectorize=True,
        output_core_dims=[["vertex", "xy"]],
    ).transpose("vertex", "animal", "timepoint", "pair", "xy")
    lines_list = []
    for a in lines.animal:
        for t in lines.timepoint:
            for p in lines.pair:
                line_i = lines.sel(animal=a, timepoint=t, pair=p).values
                line_i = np.hstack((np.tile([a, t, p], (n, 1)), line_i))
                lines_list.append(line_i)
    return lines_list


def stack_to_hyperstack(data, channel_order, strains, metadata=None):
    """Create a 6D 'Hyperstack' from a 3D stack of images

    channel_order: list
        the order of the channels in the stack
    strains: list
        the strain each animal
    """
    wvls, pairs = zip(*create_occurrence_count_tuples(channel_order))

    # calculate dimensions of hyperstack
    n_frames = data.shape[0]
    n_frames_per_animal = len(channel_order)
    n_wvls = len(np.unique(channel_order))
    n_animals = len(strains)
    n_pairs = np.max(pairs) + 1
    n_timepoints = int(n_frames / (n_frames_per_animal * n_animals))

    # create empty hyperstack with correct dimensions
    da = xr.DataArray(
        np.full(
            (n_animals, n_timepoints, n_pairs, n_wvls, data.shape[-2], data.shape[-1]),
            np.nan,
            dtype=data.dtype,
        ),
        dims=["animal", "timepoint", "pair", "wavelength", "y", "x"],
        coords={"wavelength": np.unique(channel_order), "strain": ("animal", strains),},
    )

    # fill hyperstack with data
    frame = 0
    for timepoint in range(n_timepoints):
        for animal in range(n_animals):
            for wvl, pair in zip(*[wvls, pairs]):
                # Assign image data from current frame to correct index
                da.loc[
                    dict(animal=animal, timepoint=timepoint, pair=pair, wavelength=wvl)
                ] = data[frame]

                frame += 1

    return da


def open_folder(path):
    if sys.platform == "darwin":
        subprocess.check_call(["open", "--", path])
    elif sys.platform == "linux2":
        subprocess.check_call(["xdg-open", "--", path])
    elif sys.platform == "win32":
        subprocess.check_call(["explorer", path])


def requires_matlab(func):
    pass


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


def mm2inch(*tupl):
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


def custom_round(x, base=5):
    return int(base * round(float(x) / base))


def round_down(num, divisor):
    return num - (num % divisor)


def parse_reg_param_filename(s: Path):
    """
    Extract the parameters used in the registration parameter sweep into a dictionary

    Example
    -------
    >>> s = Path("/Users/sean/code/pharedox/data/registration_param_sweep/n_deriv=2.0_rough_lambda=0.01_rough_n_breaks=300.0_rough_order=6.0_smooth_lambda=8.0_smooth_n_breaks=100.0_smooth_order=6.0_warp_lambda=10000.0_warp_n_basis=10.0_warp_order=4.0.nc")
    >>> parse_reg_param_filename(s)
    >>> {'n_deriv': '2.0',
    ...  'rough_lambda': '0.01',
    ...  'rough_n_breaks': '300.0',
    ...  'rough_order': '6.0',
    ...  'smooth_lambda': '8.0',
    ...  'smooth_n_breaks': '100.0',
    ...  'smooth_order': '6.0',
    ...  'warp_lambda': '10000.0',
    ...  'warp_n_basis': '10.0',
    ...  'warp_order': '4.0'}

    Parameters
    ----------
    s : Path
        the filename for the registration data

    Returns
    -------
    dict
        a dictionary mapping parameter keys to values
    """
    s = s.stem
    param_vals = re.findall(r"=(\d+\.\d+)_?", s)
    param_keys = re.split(r"=(\d+\.\d+)_?", s)[::2]
    return dict(zip(param_keys, param_vals))


def validate_pharynx_mask(mask: xr.DataArray):
    """
    Validate that the given pharyngeal mask image satisfies all requirements for further
    analysis.

    Parameters
    ----------
    mask : xr.DataArray
        the pharyngeal mask to validate
    """
    raise NotImplementedError


def send_data_to_matlab(data: typing.Union[xr.DataArray, np.ndarray], var_name: str):
    """
    Send data to MATLAB session currently running. Useful for debugging MATLAB code.

    Parameters
    ----------
    data : typing.Union[xr.DataArray, np.ndarray]
        The data to send to MATLAB
    var_name : str
        The name of the variable that will appear in MATLAB workspace
    """

    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Data not sent to MATLAB")

    engines = matlab.engine.find_matlab()
    if len(engines) == 0:
        raise EnvironmentError(
            """
            No MATLAB instances currently running OR the MATLAB instance is not shared...
            (1) Ensure MATLAB is running.
            (2) Make sure that the instance is shared. To do this, run (in MATLAB): matlab.engine.shareEngine
            (3) Try this function again 
            """
        )

    eng = matlab.engine.connect_matlab(engines[0])
    eng.workspace[var_name] = matlab.double(data.values.tolist())


def scale_region_boundaries(regions: dict, profile_length: int):
    """
    Scale the region boundaries

    the length of the profiles may be changed as a parameter of the experiment. As such,
    region boundaries are expressed through percentages of profile length, then scaled
    to match the length of the profile vectors.

    Examples
    --------

        >>> regions = {'a': [0, 0.1], 'b': [.2, .4]}
        >>> scale_region_boundaries(regions, 10)
        {'a': [0, 1], 'b': [2, 4]}


    Parameters
    ----------
    regions
        a dictionary mapping region names to [left_bound, right_bound] expressed as
        percentages of the profile
    profile_length
        the length of the profile that the resultant scaled regions will be used on

    Returns
    -------
    dict
        the scaled region boundary map

    """
    return {
        region: np.int_(profile_length * np.asarray(regions[region]))
        for region in regions.keys()
    }


def rmse(u, v, axis=0):
    """
    Calculate the root mean squared error (RMSE) between the two given vectors

    Parameters
    ----------
    u
        the first vector
    v
        the second vector
    axis
        the axis over which the mean should be taken

    Returns
    -------
    float
        RMSE between the two vectors

    """
    return np.sqrt(np.mean((u - v) ** 2, axis=axis))


def jaccard(im1, im2):
    """
    Computes the Jaccard metric, a measure of set similarity.

    Parameters
    ----------
    im1 : array-like, bool
        Any array of arbitrary size. If not boolean, will be converted.
    im2 : array-like, bool
        Any other array of identical size. If not boolean, will be converted.

    Returns
    -------
    jaccard : float
        Jaccard metric returned is a float on range [0,1].
        Maximum similarity = 1
        No similarity = 0
    """
    im1 = np.asarray(im1).astype(np.bool)
    im2 = np.asarray(im2).astype(np.bool)

    intersection = np.logical_and(im1, im2)
    union = np.logical_or(im1, im2)

    return intersection.sum() / float(union.sum())


def create_occurrence_count_tuples(l: typing.Iterable) -> [(typing.Any, int)]:
    """
    Given a list of things, return a list of tuples ``(item, nth_occurrence)``

    .. todo::
        test

    Parameters
    ----------
    l
        the list of things

    Returns
    -------
    occurrence_count_list
        a list of tuples ``(item, nth_occurrence)``

    """
    count_tuples = []
    c = Counter()
    for item in l:
        c.update([item])
        count_tuples.append((item, c[item] - 1))
    return count_tuples


def figure_to_np_array(fig):
    """Given a figure, return a 3D numpy array, where the dimensions are (height, width, RGB)"""
    fig.tight_layout(pad=0)
    fig.canvas.draw()
    data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return data


def calc_max_bbox(
    rot_seg_stack: xr.DataArray, ref_pair: int = 0, ref_wvl: str = "410"
) -> (int, int, int, int):
    """
    Calculate the smallest bounding-box that encapsulates all of the the pharynxes in
    the given set of rotated, segmented pharynxes.


    Parameters
    ----------
    rot_seg_stack
        the rotated and segmented stack of pharynx images
    ref_pair
        the pair to consider when calculating bounding boxes
    ref_wvl
        the wavelength to consider when calculating bounding boxes

    Returns
    -------
    a tuple of ``(min_row, min_col, max_row, max_col)``

    """
    b_boxes = []

    for i in range(rot_seg_stack.animal.size):
        props = regionprops(
            label(rot_seg_stack.sel(pair=ref_pair, wavelength=ref_wvl).isel(animal=i))
        )[0]
        b_boxes.append(props.bbox)

    b_boxes = np.vstack(b_boxes)
    min_row = np.min(b_boxes[:, 0])
    min_col = np.min(b_boxes[:, 1])
    max_row = np.max(b_boxes[:, 2])
    max_col = np.max(b_boxes[:, 3])

    return min_row, min_col, max_row, max_col


def get_mvmt_pair_i(mvmt, pair):
    return mvmt.loc[pd.IndexSlice[:, pair], :]


def get_valid_filename(s):
    """
    Return the given string converted to a string that can be used for a clean
    filename. Remove leading and trailing spaces; convert other spaces to
    underscores; and remove anything that is not an alphanumeric, dash,
    underscore, or dot.

    From: https://github.com/django/django/blob/master/django/utils/text.py

    Examples
    --------

        >>> get_valid_filename("john's portrait in 2004.jpg")
        'johns_portrait_in_2004.jpg'

    Parameters
    ----------
    s
        a string to convert to valid filename

    Returns
    -------
    str
        a valid filename

    """
    s = str(s).strip().replace(" ", "_")
    return re.sub(r"(?u)[^-\w.]", "", s)


def z_transform(fd):
    """
    Z-standardize a functional data grid object

    computed by taking (x-mean(x))/std(x)

    Parameters
    ----------
    fd

    Returns
    -------

    """
    try:
        data = np.squeeze(fd.data_matrix)
    except AttributeError:
        data = fd
    means = np.mean(data, axis=1)
    stds = np.std(data, axis=1)
    try:
        fd.data_matrix = ((data.T - means) / stds).T
    except AttributeError:
        return ((data.T - means) / stds).T
    return fd


def measure_shifted_midlines(
    experiment,
    shift_range: typing.Iterable[float],
    shift_steps: int,
    n_points: int = 200,
) -> (xr.DataArray, typing.List[float], typing.List[int]):
    """
    Measure under shifted midlines for use in synthetic movement analysis

    First, **only non-moving animals** are selected for analysis (this guarantees that
    for 0-shifts, there will be no errors).

    shifted midlines are then generated according to the following pseudocode::

        for each pair:
            for each animal:
                for each shift (dx):
                    measure under 410 and 470 using the 410 midline, storing the resultant
                    measurements as pair 0

                    shift the 470 midline by dx

                    measure under 410 and 470 using 410-dx midline, storing the resultant
                    measurements as pair 1

    Parameters
    ----------
    experiment
        The experiment to analyze
    shift_range
        a tuple of ``(min, max)`` indicating the bounds of the shifts
    shift_steps
        the number of shifts to use within the bounds of ``shift_range``
    n_points
        the number of points to measure under each midline

    Returns
    -------
    (xr.DataArray, typing.List[float], typing.List[int])
        a tuple containing:

            - the measurements
            - a list of shifts (where the index corresponds to the first dimension of
              the measurements (the animal))
            - a list of original indices (i.e. the experimental indices) for each
              synthetic index

    """
    shifts = np.linspace(*shift_range, shift_steps)

    ex_meas = experiment.trimmed_raw_profiles

    df = experiment.movement
    df = pd.DataFrame(df.to_records())

    # indexed by pair
    stationary_animals = np.array(
        [
            df[
                (df.pair == pair) & (df.anterior == 0) & (df.posterior == 0)
            ].animal.values
            for pair in df.pair.unique()
        ]
    )

    measurements = xr.DataArray(
        np.zeros(
            (
                len(np.concatenate(stationary_animals)) * len(shifts),
                ex_meas.wavelength.size,
                ex_meas.pair.size,
                n_points,
            )
        ),
        dims=["animal", "wavelength", "pair", "position"],
        coords={"wavelength": ex_meas.wavelength, "pair": ex_meas.pair},
    )

    all_shifts = []
    orig_idx = []
    new_animal_idx = 0
    for pair in range(stationary_animals.shape[0]):
        for orig_animal_idx in stationary_animals[pair]:
            midline = experiment.midlines[orig_animal_idx]["410"][pair]

            # TODO: generalize over wavelengths
            i410 = experiment.rot_fl.sel(wavelength="410", pair=pair)[
                orig_animal_idx
            ].values
            i470 = experiment.rot_fl.sel(wavelength="470", pair=pair)[
                orig_animal_idx
            ].values

            for shift in shifts:
                all_shifts.append(shift)
                orig_idx.append(orig_animal_idx)

                xs, ys = midline.linspace(n=n_points)

                # First, measure under 410 and 470 like normal, (will be pair 0)
                unshifted_410 = ndi.map_coordinates(i410, np.stack([ys, xs]), order=1)
                unshifted_470 = ndi.map_coordinates(i470, np.stack([ys, xs]), order=1)

                # now, shift the coordinates of the midline and measure under 470
                # (will be pair 1)
                shifted_470 = ndi.map_coordinates(
                    i470, np.stack([ys, xs + shift]), order=1
                )

                measurements[new_animal_idx].loc[
                    {"pair": 0, "wavelength": "410"}
                ] = unshifted_410
                measurements[new_animal_idx].loc[
                    {"pair": 0, "wavelength": "470"}
                ] = unshifted_470

                measurements[new_animal_idx].loc[
                    {"pair": 1, "wavelength": "410"}
                ] = unshifted_410
                measurements[new_animal_idx].loc[
                    {"pair": 1, "wavelength": "470"}
                ] = shifted_470

                new_animal_idx += 1

    return measurements, all_shifts, orig_idx


def expand_dimension(
    data: xr.DataArray,
    dim: str,
    new_coords: typing.Dict[str, typing.Union[xr.DataArray, np.ndarray]],
):
    all_coords = np.append(data[dim].data, list(new_coords.keys()))
    data = data.reindex(**{dim: all_coords})
    for coord, new_data in new_coords.items():
        try:
            # if new_data is an xr.DataArray
            data.loc[{dim: coord}] = new_data.values
        except AttributeError:
            # if new_data is an np.ndarray
            data.loc[{dim: coord}] = new_data
    return data


def add_derived_wavelengths(
    data,
    r_min=0.852,
    r_max=6.65,
    instrument_factor=0.171,
    midpoint_potential=-265.0,
    z=2,
    temperature=22.0,
    ratio_numerator="410",
    ratio_denominator="470",
):
    """ 
    Add "derived" wavelengths to the given DataArray. These derived wavelengths are
    the ratio (`r`), the fraction oxidized (`oxd`), and the reduction potential (`e`).

    Parameters
    ----------
    data : [type]
        [description]
    r_min : [type]
        [description]
    r_max : [type]
        [description]
    instrument_factor : [type]
        [description]
    midpoint_potential : [type]
        [description]
    z : [type]
        [description]
    temperature : [type]
        [description]
    numerator : str, optional
        [description], by default "410"
    denominator : str, optional
        [description], by default "470"

    Returns
    -------
    [type]
        [description]
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = data.sel(wavelength=ratio_numerator) / data.sel(
            wavelength=ratio_denominator
        )
        oxd = pp.r_to_oxd(
            r, r_min=r_min, r_max=r_max, instrument_factor=instrument_factor
        )
        e = pp.oxd_to_redox_potential(
            oxd, midpoint_potential=midpoint_potential, z=z, temperature=temperature
        )

    # Need to add time coordinates to these dimensions
    # time comes from avg(time(410), time(470))
    try:
        t1 = data.sel(wavelength=ratio_numerator).time
        t2 = data.sel(wavelength=ratio_denominator).time
        time = t1 + ((t1 - t2) / 2)
    except AttributeError:
        time = None

    for wvl, wvl_data in zip(["r", "oxd", "e"], [r, oxd, e]):
        if wvl in data.wavelength.values:
            data.loc[dict(wavelength=wvl)] = wvl_data
        else:
            data = expand_dimension(data, "wavelength", {wvl: wvl_data})

        if time is not None:
            data.time.loc[dict(wavelength=wvl)] = time

    return data


def git_version() -> str:
    """
    Return the current git revision.
    
    Stolen from Numpy's internal code at https://github.com/numpy/numpy/blob/578f4e7dca4701637284c782d8c74c0d5b688341/setup.py#L65
    
    Returns
    -------
    str
        the current git revision
    """

    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ["SYSTEMROOT", "PATH"]:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env["LANGUAGE"] = "C"
        env["LANG"] = "C"
        env["LC_ALL"] = "C"
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(["git", "rev-parse", "HEAD"])
        GIT_REVISION = out.strip().decode("ascii")
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def setup_logging(loglevel="info"):
    import logging

    logmap = {"info": logging.INFO, "debug": logging.DEBUG}

    logging.basicConfig(
        format="%(asctime)s %(levelname)s:%(message)s",
        level=logmap[loglevel],
        datefmt="%I:%M:%S",
    )
