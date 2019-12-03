import collections
from dataclasses import dataclass
from typing import Union

import scipy
import numpy as np
import xarray as xr
import pandas as pd
import skfda
import matlab.engine
from skfda.representation.basis import BSpline
from skfda.preprocessing.smoothing import BasisSmoother
from sklearn.preprocessing import scale

from pharynx_redox import utils, constants


@dataclass
class RegData:
    """
    A `dataclass <https://docs.python.org/3/library/dataclasses.html>`_ which acts as a
    container for registered data.

    Attributes
    ----------
    r410: skfda.representation.FDataBasis
        A functional object containing the registered 410nm profile data
    r470: skfda.representation.FDataBasis
        A functional object containing the registered 470nm profile data
    warp: skfda.representation.FDataGrid
        A functional object containing the warp function that were used to warp the 470nm profile data into the 410nm profile data
    reg_data: xr.DataArray
        The quantized registered profile data

    """

    r410: [skfda.representation.FDataBasis]
    r470: [skfda.representation.FDataBasis]
    warps: [skfda.representation.FDataGrid]
    reg_data: xr.DataArray


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
            .to_pandas().to_frame()
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
            "spec": profile_data.spec,
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
            sm_data = np.array(
                constants.matlab_engine.smooth_profiles(
                    raw_data, l, order, n_basis, n_eval_pts
                )
            ).T
            smooth_profile_data.loc[dict(pair=pair, wavelength=wvl)] = sm_data
    return smooth_profile_data


def register_profiles_matlab(
    raw_profile_data,
    n_deriv: int = 2,
    rough_lambda: float = 10e-3,
    smooth_lambda: float = 10e-1,
    warp_lambda: float = 10e-1,
    smooth_n_breaks: float = 100.0,
    rough_n_breaks: float = 300.0,
    warp_n_basis: float = 300.0,
    warp_order: float = 4.0,
    smooth_order: float = 4.0,
    rough_order: float = 4.0,
    resample_resolution: int = None,
) -> (xr.DataArray, np.ndarray):
    eng = matlab.engine.start_matlab()
    reg_profile_data = raw_profile_data.copy()

    # Registration Parameters
    if resample_resolution is None:
        resample_resolution = reg_profile_data.shape[-1]

    for pair in raw_profile_data.pair.values:
        print(f"Registering Pair {pair}")
        i410 = matlab.double(
            raw_profile_data.sel(wavelength="410", pair=pair).values.T.tolist()
        )
        i470 = matlab.double(
            raw_profile_data.sel(wavelength="470", pair=pair).values.T.tolist()
        )

        # Call the MATLAB subroutine
        r410, r470, warpfd, wfd, d410, reg_d470 = eng.channel_register_discrete(
            i410,
            i470,
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
            resample_resolution,
            nargout=6,
        )
        r410, r470 = np.array(r410).T, np.array(r470).T
        warpfd = np.array(warpfd)
        wfd = np.array(wfd)
        d410 = np.array(d410)
        reg_d470 = np.array(reg_d470)

        reg_profile_data.loc[dict(pair=pair, wavelength="410")] = r410
        reg_profile_data.loc[dict(pair=pair, wavelength="470")] = r470

    return reg_profile_data


def scale_by_wvl(data: xr.DataArray) -> xr.DataArray:
    """
    Apply Z-transform scaling on a per-wavelength basis

    Parameters
    ----------
    data
        the data to scale

    Returns
    -------
    xr.DataArray
        the z-scaled data

    """
    scaled = data.copy()
    for wvl in scaled.wavelength:
        for pair in scaled.pair:
            unscaled = scaled.sel(wavelength=wvl, pair=pair).values
            scaled.loc[dict(wavelength=wvl, pair=pair)] = scale(unscaled, axis=1)
    return scaled


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


# def create_empty_profile_xarray(template=None, position_size=None):
#     return xr.DataArray(
#         0,
#         dims=["spec", "wavelength", "pair", "position"],
#         coords={
#             "spec": intensity_data.spec,
#             "wavelength": intensity_data.wavelength,
#             "pair": intensity_data.pair,
#             "position": np.arange(new_length),
#             "strain": ("spec", intensity_data.strain),
#         },
#     )


def trim_profiles(
    intensity_data: xr.DataArray,
    threshold: float,
    new_length: int,
    ref_wvl: str = "410",
):
    """
    Parameters
    ----------
    ref_wvl
    intensity_data
    threshold
    new_length

    Returns
    -------

    """
    trimmed_intensity_data = xr.DataArray(
        0,
        dims=["spec", "wavelength", "pair", "position"],
        coords={
            "spec": intensity_data.spec,
            "wavelength": intensity_data.wavelength,
            "pair": intensity_data.pair,
            "position": np.arange(new_length),
            "strain": ("spec", intensity_data.strain),
        },
    )

    l, r = get_trim_boundaries(intensity_data, ref_wvl=ref_wvl, thresh=threshold)

    for img_idx in range(intensity_data.spec.size):
        for wvl_idx in range(intensity_data.wavelength.size):
            wvl = intensity_data.wavelength.data[wvl_idx]
            if "tl" not in wvl.lower():
                for pair in range(intensity_data.pair.size):
                    data = (
                        intensity_data.sel(wavelength=wvl, pair=pair)
                        .isel(spec=img_idx)
                        .data
                    )

                    trimmed = data[l[img_idx, pair] : r[img_idx, pair]]
                    new_xs = np.linspace(0, len(trimmed), new_length)
                    old_xs = np.arange(0, len(trimmed))
                    resized = np.interp(new_xs, old_xs, trimmed)

                    trimmed_intensity_data[img_idx, wvl_idx, pair, :] = resized

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


if __name__ == "__main__":
    profile_data = xr.load_dataarray(
        "/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/analyses/2019-08-26_single_unreg/2017_02_22-HD233_SAY47-profile_data.nc"
    )
    register_profiles_matlab(profile_data)
