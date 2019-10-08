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

from pharynx_redox import utils


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
            .to_pandas()
        )
        reg_df = reg_df.reset_index()
        reg_df["region"] = region
        reg_df["animal"] = range(len(reg_df))
        reg_df.rename({0: value_name}, inplace=True, axis="columns")
        reg_df.drop(reg_df.columns[0], axis=1, inplace=True)
        dfs.append(reg_df)
    return pd.concat(dfs)


def smooth_profile_data(
    profile_data: Union[np.ndarray, xr.DataArray],
    order: int = 5,
    nbasis: int = 200,
    smoothing_parameter: float = 1e-8,
    ret_basis: bool = True,
):
    """
    Smooth profile data by fitting smoothing B-splines

    Parameters
    ----------
    profile_data
        the data to smooth, maximally 2-dimensional of shape (animals, position)
    order
        the order of the basis polynomials
    nbasis
        the number of basis polynomials
    smoothing_parameter
        the smoothing parameter. Larger values increase smoothness. If 0, no smoothing
        is performed. Try values logarithmically.

    Returns
    -------
    skfda.representation.grid.FDataBasis
        the smoothed data represented in Basis form
    skfda.representation.basis.FDataBasis
        the basis used for smoothed data

    """
    fd = skfda.FDataGrid(profile_data)
    basis = BSpline(order=order, nbasis=nbasis)
    smoother = BasisSmoother(
        basis, smoothing_parameter=smoothing_parameter, method="qr", return_basis=True
    )
    fd_smooth = smoother.fit_transform(fd)

    return fd_smooth, basis


def register_pair(
    data1,
    data2,
    smooth_lambda: float = 1e-5,
    rough_lambda: float = 1e-7,
    warp_lam: float = 1e-1,
    smooth_nbasis: int = 64,
    rough_nbasis: int = 200,
    grid_dim: int = 7,
):
    """
    TODO: documentation

    Parameters
    ----------
    data1
    data2
    smooth_lambda
    rough_lambda
    warp_lam
    smooth_nbasis
    rough_nbasis
    grid_dim

    Returns
    -------

    """
    pos_size = data1.shape[-1]

    f1_sm, _ = smooth_profile_data(
        data1, smoothing_parameter=smooth_lambda, nbasis=smooth_nbasis
    )
    f2_sm, _ = smooth_profile_data(
        data2, smoothing_parameter=smooth_lambda, nbasis=smooth_nbasis
    )

    f1_rough, _ = smooth_profile_data(
        data1, smoothing_parameter=rough_lambda, nbasis=rough_nbasis
    )
    f2_rough, _ = smooth_profile_data(
        data2, smoothing_parameter=rough_lambda, nbasis=rough_nbasis
    )

    # Resample
    f1_sm = f1_sm.to_grid(np.linspace(*f1_sm.domain_range[0], pos_size))
    f2_sm = f2_sm.to_grid(np.linspace(*f2_sm.domain_range[0], pos_size))

    # We want to register with the derivatives
    d1_sm = f1_sm.derivative()
    d2_sm = f2_sm.derivative()

    # Calculate the warp
    warp = skfda.preprocessing.registration.elastic_registration_warping(
        d1_sm, d2_sm, lam=warp_lam, grid_dim=grid_dim
    )

    # Apply inverse warp to f2
    warp_inv = skfda.preprocessing.registration.invert_warping(warp)
    f2_rough = f2_rough.compose(warp_inv)

    return f1_rough, f2_rough, warp


# noinspection PyUnresolvedReferences
def register_profiles(
    raw_profile_data: xr.DataArray,
    smooth_lambda: float = 1e-5,
    rough_lambda: float = 1e-7,
    warp_lam: float = 1e-1,
    smooth_nbasis: int = 64,
    rough_nbasis: int = 200,
    warp_to_mean: bool = False,
) -> RegData:
    """
    Register 470 channels into 410 using elastic registration.

    First, the data is smoothed using a B-spline basis with the specified number of
    basis functions. The *derivatives* of these functions are then computed for use in
    registration.

    The 470nm profile data is registered into the 410nm data by calculating warping
    functions (which map old x-coordinates to new ones) using these smooth derivatives.

    Finally, the warping functions are applied to a "rougher" fit of the data (again,
    this rough fit is a B-spline fit using the specified number of basis functions and
    smoothing parameter).

    TODO: generalize over wavelengths

    Parameters
    ----------
    smooth_lambda
        The smoothing constraint for the "smooth" data, which will be used for
        registration
    rough_lambda
        The smoothing constraint for the "rough" data, which is what will be returned
        as the measurement data
    warp_lam
        The smoothing constraint for the warp function
    raw_profile_data
        The profile data to register
    smooth_nbasis
        the number of basis functions for the "smooth" data, which will be used for
        registration
    rough_nbasis
        the number of basis function for the "rough" data, which is what will be
        returned after warping as the measurement data
    warp_to_mean
        if True, warps all profiles to the mean of the 410nm (pair 0) profile data


    Returns
    -------
    RegData
        A named tuple containing::

            (r410, r470, warp, reg_profile_data)

        r410, r470, and warp are lists indexed by pair, reg_profile_data is an
        xr.DataArray with the same structure as raw_profile_data


    """
    reg_profile_data = raw_profile_data.copy()

    pos_size = raw_profile_data.position.size

    # For a two-pair image set, these lists will be two items long. The 0th item will
    # contain a functional data object containing ALL observations for pair 0, etc.
    f410 = []
    f470 = []
    warps = []

    # Register each pair separately
    for pair in raw_profile_data.pair.data:
        i410 = raw_profile_data.sel(wavelength="410", pair=pair).values
        i470 = raw_profile_data.sel(wavelength="470", pair=pair).values

        f410_sm, _ = smooth_profile_data(
            i410, smoothing_parameter=smooth_lambda, nbasis=smooth_nbasis
        )
        f470_sm, _ = smooth_profile_data(
            i470, smoothing_parameter=smooth_lambda, nbasis=smooth_nbasis
        )

        f410_rough, _ = smooth_profile_data(
            i410, smoothing_parameter=rough_lambda, nbasis=rough_nbasis
        )
        f470_rough, _ = smooth_profile_data(
            i470, smoothing_parameter=rough_lambda, nbasis=rough_nbasis
        )

        # We will register with the derivatives
        d410_sm = f410_sm.derivative()
        d470_sm = f470_sm.derivative()

        # First, register all to mean if specified
        if warp_to_mean:
            template = f410_rough.mean().to_grid()

            # Calculate the warp using the derivatives of the smooth data
            warp_d410_sm = skfda.preprocessing.registration.elastic_registration_warping(
                d410_sm.to_grid(), template, lam=warp_lam
            )
            warp_d470_sm = skfda.preprocessing.registration.elastic_registration_warping(
                d470_sm.to_grid(), template, lam=warp_lam
            )

            # Apply the calculated warps to the rough fits
            warp_d410_sm_inv = skfda.preprocessing.registration.invert_warping(
                warp_d410_sm
            )
            warp_d470_sm_inv = skfda.preprocessing.registration.invert_warping(
                warp_d470_sm
            )

            f410_rough = f410_rough.compose(warp_d410_sm_inv)
            f470_rough = f470_rough.compose(warp_d470_sm_inv)

            # Re-smooth so we can do pair-wise registration later
            f410_sm, _ = smooth_profile_data(
                np.squeeze(f410_rough.to_grid().data_matrix),
                smoothing_parameter=smooth_lambda,
                nbasis=smooth_nbasis,
            )
            f470_sm, _ = smooth_profile_data(
                np.squeeze(f470_rough.to_grid().data_matrix),
                smoothing_parameter=smooth_lambda,
                nbasis=smooth_nbasis,
            )

            d410_sm = f410_sm.derivative()
            d470_sm = f470_sm.derivative()

        # Calculate the pair-wise warp
        warp = skfda.preprocessing.registration.elastic_registration_warping(
            d410_sm.to_grid(), d470_sm.to_grid(), lam=warp_lam
        )

        # Apply inverse warp to 470
        warp_inv = skfda.preprocessing.registration.invert_warping(warp)
        f470_rough = f470_rough.compose(warp_inv)

        # Keep track of the functional objects to return
        f410.append(f410_rough)
        f470.append(f470_rough)
        warps.append(warp)

        # Save the quantized data
        reg_profile_data.loc[dict(pair=pair, wavelength="410")] = np.squeeze(
            f410_rough.to_grid(
                np.linspace(*f470_sm.domain_range[0], pos_size)
            ).data_matrix
        )
        reg_profile_data.loc[dict(pair=pair, wavelength="470")] = np.squeeze(
            f470_rough.to_grid(
                np.linspace(*f470_sm.domain_range[0], pos_size)
            ).data_matrix
        )

    return RegData(f410, f470, warps, reg_profile_data)


def register_profiles_matlab(raw_profile_data) -> (xr.DataArray, np.ndarray):
    eng = matlab.engine.start_matlab()
    reg_profile_data = raw_profile_data.copy()

    # Warps indexed by pair
    all_warps = []

    # Registration Parameters
    resample_resolution = reg_profile_data.shape[-1]

    warp_n_basis = 300.0
    warp_order = 4.0
    warp_lambda = 10 ** 5

    smooth_lambda = 10 ** 2
    smooth_n_breaks = 100.0
    smooth_order = 4.0

    rough_lambda = 10 ** 0.05
    rough_n_breaks = 300.0
    rough_order = 4.0

    n_deriv = 2.0

    for pair in raw_profile_data.pair.values:
        print(f"Registering Pair {pair}")
        i410 = matlab.double(
            raw_profile_data.sel(wavelength="410", pair=pair).values.T.tolist()
        )
        i470 = matlab.double(
            raw_profile_data.sel(wavelength="470", pair=pair).values.T.tolist()
        )

        # _, _, r410, r470, warps = eng.ChannelRegister(i410, i470, nargout=5)
        r410, r470, warps = eng.channel_register(
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

        reg_profile_data.loc[dict(pair=pair, wavelength="410")] = r410
        reg_profile_data.loc[dict(pair=pair, wavelength="470")] = r470

        all_warps.append(np.array(warps))

    return reg_profile_data, np.array(all_warps)


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

