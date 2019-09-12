import collections
from typing import Union

import numpy as np
import xarray as xr
import pandas as pd
import skfda
from skfda.representation.basis import BSpline
from skfda.preprocessing.smoothing import BasisSmoother

from sklearn.preprocessing import scale

from pharynx_redox import utils

RegData = collections.namedtuple("RegData", "r410 r470 warp reg_data")
RegData.__doc__ = """
    A collection class for storing the results of data registration
"""


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
            .mean(dim="position")
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
    order=5,
    nbasis=200,
    smoothing_parameter=1e-8,
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
    skfda.representation.grid.FDataGrid
        the smoothed data
    skfda.representation.basis.FDataBasis
        the basis representation of the smoothed data

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
    smooth_lambda=1e-5,
    rough_lambda=1e-7,
    warp_lam=1e-1,
    smooth_nbasis=64,
    rough_nbasis=200,
    grid_dim=7,
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
    Register 470 channels into 410 using elastic registration

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
        if True, warps all profiles to the mean


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

    # lists indexed by pair
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

        # Resample
        # f410_sm = f410_sm.to_grid(np.linspace(*f410_sm.domain_range[0], pos_size))
        # f470_sm = f470_sm.to_grid(np.linspace(*f470_sm.domain_range[0], pos_size))

        # # Rescale with Z-transform
        # z410_sm = utils.z_transform(f410_sm)
        # z470_sm = utils.z_transform(f470_sm)

        # We will register with the derivatives of the z-transformed data
        d410_sm = f410_sm.derivative()
        d470_sm = f470_sm.derivative()

        # First, register all to mean if specified
        if warp_to_mean:
            template = f410_rough.mean().to_grid()

            warp_d410_sm = skfda.preprocessing.registration.elastic_registration_warping(
                d410_sm.to_grid(), template, lam=warp_lam
            )
            warp_d470_sm = skfda.preprocessing.registration.elastic_registration_warping(
                d470_sm.to_grid(), template, lam=warp_lam
            )

            warp_d410_sm_inv = skfda.preprocessing.registration.invert_warping(
                warp_d410_sm
            )
            warp_d470_sm_inv = skfda.preprocessing.registration.invert_warping(
                warp_d470_sm
            )

            f410_rough = f410_rough.compose(warp_d410_sm_inv)
            f470_rough = f470_rough.compose(warp_d470_sm_inv)

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


def scale_by_wvl(data: xr.DataArray) -> xr.DataArray:
    scaled = data.copy()
    for wvl in scaled.wavelength:
        for pair in scaled.pair:
            unscaled = scaled.sel(wavelength=wvl, pair=pair).values
            scaled.loc[dict(wavelength=wvl, pair=pair)] = scale(unscaled, axis=1)
    return scaled


def trim_profile(profile, threshold, new_length):
    """
    Parameters
    ----------
    profile
    threshold
    new_length

    Returns
    -------

    """
    first = np.argmax(profile > threshold)
    last = len(profile) - np.argmax(np.flip(profile > threshold))

    trimmed = profile[first : last + 1]
    new_xs = np.linspace(0, len(trimmed), new_length)
    old_xs = np.arange(0, len(trimmed))

    return np.interp(new_xs, old_xs, trimmed)


def get_trim_boundaries(data, ref_wvl="410", thresh=2000):
    prof_len = data.position.size
    l_bound = np.argmax(data.sel(wavelength=ref_wvl) >= thresh, axis=2).data - 1
    r_bound = (
        prof_len
        - np.argmax(
            np.flip(data.sel(wavelength=ref_wvl), axis=2) >= thresh, axis=2
        ).data
    )
    return l_bound, r_bound


def trim_profiles(intensity_data, threshold, new_length, ref_wvl="410"):
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
        np.zeros(
            (
                intensity_data.strain.size,
                intensity_data.wavelength.size,
                intensity_data.pair.size,
                new_length,
            )
        ),
        dims=["strain", "wavelength", "pair", "position"],
        coords={
            "strain": intensity_data.strain,
            "wavelength": intensity_data.wavelength,
            "pair": intensity_data.pair,
        },
    )

    l, r = get_trim_boundaries(intensity_data, ref_wvl=ref_wvl, thresh=threshold)

    for img_idx in range(intensity_data.strain.size):
        for wvl_idx in range(intensity_data.wavelength.size):
            wvl = intensity_data.wavelength.data[wvl_idx]
            if "tl" not in wvl.lower():
                for pair in range(intensity_data.pair.size):
                    data = (
                        intensity_data.sel(wavelength=wvl, pair=pair)
                        .isel(strain=img_idx)
                        .data
                    )

                    trimmed = data[l[img_idx, pair] : r[img_idx, pair]]
                    new_xs = np.linspace(0, len(trimmed), new_length)
                    old_xs = np.arange(0, len(trimmed))
                    resized = np.interp(new_xs, old_xs, trimmed)

                    trimmed_intensity_data[img_idx, wvl_idx, pair, :] = resized

    return trimmed_intensity_data


def r_to_oxd(r, r_min=0.852, r_max=6.65, instrument_factor=0.171):
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


def oxd_to_redox_potential(oxd, midpoint_potential=-265, z=2, temperature=22):
    """
    Convert OxD to redox potential

    NOTE: may return NaN

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
