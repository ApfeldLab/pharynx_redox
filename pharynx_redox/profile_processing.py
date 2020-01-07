import collections
from dataclasses import dataclass
from typing import Union
import logging

import scipy
import numpy as np
import xarray as xr
import pandas as pd
from sklearn.preprocessing import scale

from . import utils
from . import constants

import matlab.engine


def summarize_over_regions(
    data: Union[np.ndarray, xr.DataArray],
    regions: dict,
    value_name: str = None,
    rescale: bool = False,
):
    """
    Summarize the profile data over region boundaries, storing the resultant data in a pandas DataFrame.

    The functional observations stored as profile data is difficult to analyze with
    traditional statistical methods. To make it easier, we can average the values within
    regions of the profile which correspond to salient biological regions (different
    muscle cells, etc).

    This image demonstrates the idea, where averages would be taken inside the shaded
    areas, and the value assigned to the label of that shaded region.

    .. image:: _static/ex_i410.png



    Parameters
    ----------
    data
        the profile data to be summarized. maximally 2-dimensional where the shape is
        (observation, position)
    regions
        a dictionary mapping region name to [left_boundary, right_boundary]. This can be
        unscaled (where boundaries are expressed as percentages) or scaled (where
        boundaries are expressed as integer indices of the profile vector). If unscaled,
        the ``rescale`` parameter should be set to ``True``.
    value_name
        the name of the value being measured (e.g. "410", "E", "OxD"). If left as None,
        defaults to "value". This will be the heading of the column storing the
        measurements in the resultant DataFrame
    rescale
        whether the regions should be rescaled, see the ``regions`` parameter for more
        info

    Returns
    -------
    pd.DataFrame
        a pandas DataFrame containing the mean value of each profile within the given
        region boundaries

        ====  =======  ========  ========
          ..  value    region      animal
        ====  =======  ========  ========
           0  7856.16  pm3              0
           1  7598.52  pm3              1
           2  9302.02  pm3              2
        ...   ...      ...       ...
        ====  =======  ========  ========

    """
    if rescale:
        regions = utils.scale_region_boundaries(regions, data.shape[-1])
    if value_name is None:
        value_name = "value"

    if type(data) == np.ndarray:
        data = xr.DataArray(data, dims=("observations", "position"))

    dfs = []
    for region, bounds in regions.items():
        reg_df = (
            data[dict(position=range(bounds[0], bounds[1]))]
            .mean(dim="position", skipna=True)
            .to_pandas()
            .to_frame()
        )
        # reg_df = reg_df.reset_index()
        reg_df["region"] = region
        # reg_df["animal"] = range(len(reg_df))
        reg_df.rename({0: value_name}, inplace=True, axis="columns")
        # reg_df.drop(reg_df.columns[0], axis=1, inplace=True)
        dfs.append(reg_df)
    return pd.concat(dfs)


def smooth_profile_data(
    profile_data: Union[np.ndarray, xr.DataArray],
    l: float = 1.12,
    order: float = 4.0,
    n_basis: float = 200,
    n_eval_pts: int = 1000,
):
    """
    Smooth profile data by fitting smoothing B-splines

    Implemented in MATLAB as smooth_profiles
    """

    smooth_profile_data = xr.DataArray(
        0,
        dims=profile_data.dims,
        coords={
            "animal": profile_data.animal,
            "wavelength": profile_data.wavelength,
            "pair": profile_data.pair,
            "position": np.arange(n_eval_pts),
        },
    )
    for wvl in profile_data.wavelength.values:
        for pair in profile_data.pair.values:
            raw_data = matlab.double(
                profile_data.sel(wavelength=wvl, pair=pair)
                .values.astype(np.float)
                .T.tolist()
            )

            eng = matlab.engine.start_matlab()
            sm_data = np.array(
                eng.smooth_profiles(raw_data, l, order, n_basis, n_eval_pts)
            ).T
            smooth_profile_data.loc[dict(pair=pair, wavelength=wvl)] = sm_data
    return smooth_profile_data


def register_profiles_pairs(
    profile_data: xr.DataArray, eng: matlab.engine.MatlabEngine = None, **reg_params
) -> xr.DataArray:

    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()

    for pair in profile_data.pair:
        reg_pair = register_profiles(profile_data.sel(pair=pair), eng=eng, **reg_params)
        reg_profile_data.loc[dict(pair=pair)] = reg_pair

    return reg_profile_data


def register_profiles(
    profile_data: xr.DataArray,
    eng: matlab.engine.MatlabEngine = None,
    n_deriv: float = 2.0,
    rough_lambda: float = 10.0 ** 0.05,
    rough_n_breaks: float = 300.0,
    rough_order: float = 4.0,
    smooth_lambda: float = 10.0 ** 2,
    smooth_n_breaks: float = 100.0,
    smooth_order: float = 4.0,
    warp_lambda: float = 5.0e3,
    warp_n_basis: float = 30.0,
    warp_order: float = 4.0,
) -> xr.DataArray:
    """
    Register the 470nm channel into the 410nm channel profile data
    
    Parameters
    ----------
    profile_data : xr.DataArray
        The data to register. Must be maximally 3-dimensional. The dimensions are (animal, wavelength, position along profile).
    eng: matlab.engine.MatlabEngine
        The MATLAB engine to use for registration. If None, a new engine is created.
    """
    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()

    i410 = matlab.double(profile_data.sel(wavelength="410").values.tolist())
    i470 = matlab.double(profile_data.sel(wavelength="470").values.tolist())

    resample_resolution = float(profile_data.position.size)

    # Call the MATLAB subroutine
    r410, r470, _ = eng.channel_register(
        i410,
        i470,
        resample_resolution,
        warp_n_basis,
        warp_order,
        warp_lambda,
        smooth_lambda,
        smooth_n_breaks,
        smooth_order,
        rough_lambda,
        rough_n_breaks,
        rough_order,
        n_deriv,
        nargout=3,
    )
    r410, r470 = np.array(r410).T, np.array(r470).T

    reg_profile_data.loc[dict(wavelength="410")] = r410
    reg_profile_data.loc[dict(wavelength="470")] = r470
    reg_profile_data = utils.add_derived_wavelengths(reg_profile_data)

    return reg_profile_data


def trim_profile(
    profile: Union[np.ndarray, xr.DataArray], threshold: float, new_length: int
):
    """
    Trim the given profile data by finding the first/last values where the profile
    crosses the specified threshold, then interpolating to fit the given new length.

    .. note::
        Uses linear interpolation

    Parameters
    ----------
    profile
        the data to trim
    threshold
        the threshold
    new_length
        the length of the resultant interpolated profiles

    Returns
    -------

    """
    first = np.argmax(profile > threshold)
    last = len(profile) - np.argmax(np.flip(profile > threshold))

    trimmed = profile[first : last + 1]
    new_xs = np.linspace(0, len(trimmed), new_length)
    old_xs = np.arange(0, len(trimmed))

    return np.interp(new_xs, old_xs, trimmed)


def get_trim_boundaries(
    data: xr.DataArray, ref_wvl: str = "410", thresh: float = 2000.0
) -> (np.ndarray, np.ndarray):
    """
    Find the "left" and "right" indices to use to trim intensity profiles given a
    threshold.

    Essentially, we find the first index where the intensity profile crosses the given
    threshold and call that the "left", then do the same on the reversed profile and
    call that the "right".


    Parameters
    ----------
    data
        the intensity profile data (potentially containing multiple wavelengths)
    ref_wvl
        the wavelength to use to calculate boundaries
    thresh
        the threshold

    Returns
    -------
    (np.ndarray, np.ndarray)
        the (left, right) bounds for each profile, where the index in the array
        corresponds to the index of the animal in ``data``.

    """
    prof_len = data.position.size
    l_bound = np.argmax(data.sel(wavelength=ref_wvl) >= thresh, axis=2).data - 1
    r_bound = (
        prof_len
        - np.argmax(
            np.flip(data.sel(wavelength=ref_wvl), axis=2) >= thresh, axis=2
        ).data
    )
    return l_bound, r_bound


def trim_profiles(
    intensity_data: xr.DataArray, threshold: float, ref_wvl: str = "410"
) -> xr.DataArray:
    """
    Trim the background away from the profiles.
    
    Parameters
    ----------
    intensity_data : xr.DataArray
        the profile data to trim
    threshold : float
        the threshold under which data will be thrown away
    ref_wvl : str, optional
        the wavelength to be used to calculate trim boundaries. Other wavelengths will
        be trimmed using these boundaries. By default "410"
    
    Returns
    -------
    xr.DataArray
        the trimmed profiles
    """
    trimmed_intensity_data = intensity_data.copy()

    l, r = get_trim_boundaries(intensity_data, ref_wvl=ref_wvl, thresh=threshold)

    # trimmed_intensity_data.to_netcdf("/Users/sean/Desktop/test_trimming-untrimmed.nc")

    for img_idx in intensity_data.animal:
        for wvl_idx in range(intensity_data.wavelength.size):
            wvl = intensity_data.wavelength.data[wvl_idx]
            if "tl" not in wvl.lower():
                for pair in range(intensity_data.pair.size):
                    selector = dict(wavelength=wvl, pair=pair, animal=img_idx)
                    data = intensity_data.sel(selector).data

                    trimmed = data[l[img_idx, pair] : r[img_idx, pair]]
                    new_xs = np.linspace(0, len(trimmed), intensity_data.position.size)
                    old_xs = np.arange(0, len(trimmed))
                    resized = np.interp(new_xs, old_xs, trimmed)

                    trimmed_intensity_data.loc[selector] = resized

    return trimmed_intensity_data


def r_to_oxd(
    r: Union[np.ndarray, xr.DataArray],
    r_min: float = 0.852,
    r_max: float = 6.65,
    instrument_factor: float = 0.171,
):
    """
    Convert ratios to OxD

    Parameters
    ----------
    r
    r_min
    r_max
    instrument_factor

    Returns
    -------

    """
    return (r - r_min) / ((r - r_min) + instrument_factor * (r_max - r))


def oxd_to_redox_potential(
    oxd: Union[np.ndarray, xr.DataArray],
    midpoint_potential: float = -265.0,
    z: float = 2.0,
    temperature: float = 22.0,
):
    """
    Convert OxD to redox potential

    .. warning::
        May contain ``NaN`` values

    Parameters
    ----------
    oxd
    midpoint_potential
    z
    temperature

    Returns
    -------

    """
    return midpoint_potential - (
        8314.462 * (273.15 + temperature) / (z * 96485.3415)
    ) * np.log((1 - oxd) / oxd)


def find_optimal_regions(initial_regions, err_data, min_width=10, rescale_regions=True):
    try:
        err_data = err_data.values
    except AttributeError:
        pass

    def avg_err_for_bounds(region_bounds):
        l, r = int(region_bounds[0]), int(region_bounds[1])
        if np.abs(r - l) < min_width:
            return np.inf
        if r < l:
            l, r = r, l
        return np.nanmean(err_data[:, int(l) : int(r)])

    if rescale_regions:
        initial_regions = utils.scale_region_boundaries(
            initial_regions, err_data.shape[-1]
        )

    opt_regions = {}
    for region, bounds in initial_regions.items():
        opt_bounds = scipy.optimize.brute(
            avg_err_for_bounds,
            ranges=(slice(bounds[0], bounds[1], 1), slice(bounds[0], bounds[1], 1)),
            finish=None,
        )
        opt_regions[region] = opt_bounds.astype(np.int)
    return opt_regions
