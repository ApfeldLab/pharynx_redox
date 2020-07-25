import logging
import os
import re
import subprocess
import sys
import typing
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage as ndi

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


def open_folder(path):
    if sys.platform == "darwin":
        subprocess.check_call(["open", "--", path])
    elif sys.platform == "linux2":
        subprocess.check_call(["xdg-open", "--", path])
    elif sys.platform == "win32":
        subprocess.check_call(["explorer", path])


def requires_matlab(func):
    # TODO use this as a decorator
    pass


def validate_pharynx_mask(mask: xr.DataArray):
    """
    Validate that the given pharyngeal mask image satisfies all requirements for further
    analysis.

    Parameters
    ----------
    mask : xr.DataArray
        the pharyngeal mask to validate
    """
    # TODO implement this
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
    data: xr.DataArray,
    r_min: float = 0.852,
    r_max: float = 6.65,
    instrument_factor: float = 0.171,
    midpoint_potential: float = -265.0,
    z: int = 2,
    temperature: float = 22.0,
    ratio_numerator: str = "410",
    ratio_denominator: str = "470",
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
    ratio_numerator
    ratio_denominator

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
