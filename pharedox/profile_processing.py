import logging
import warnings
from typing import Dict, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal, spatial

from pharedox import utils


def to_dataframe(data: xr.DataArray, *args, **kwargs) -> pd.DataFrame:
    """
    Replacement for `xr.DataArray.to_dataframe` that adds the attrs for the given
    DataArray into the resultant DataFrame.
    
    Parameters
    ----------
    data : xr.DataArray
        the data to convert to DataFrame
    
    Returns
    -------
    pd.DataFrame
        a pandas DataFrame containing the data in the given DataArray, including the
        global attributes
    """
    df = data.to_dataframe(*args, **kwargs)
    for k, v in data.attrs.items():
        df[k] = v
    return df


def align_pa(
    intensity_data: xr.DataArray,
    reference_wavelength: str = "410",
    reference_pair: int = 0,
    reference_timepoint: int = 0,
) -> xr.DataArray:
    """
    Given intensity profile data, flip each animal along their anterior-posterior axis
    if necessary, so that all face the same direction

    Parameters
    ----------
    intensity_data
        the data to align
    reference_wavelength: optional
        the wavelength to calculate the alignment for
    reference_pair: optional
        the pair to calculate the alignment for
    reference_timepoint
        the timepoint to calculate the alignment for

    Returns
    -------
    aligned_intensity_data
        the PA-aligned intensity data

    Notes
    -----
    The alignments are calculated for a single wavelength and pair for each animal, then
    applied to all wavelengths and pairs for that animal.

    The algorithm works as follows:

        - take the derivative of the (trimmed) intensity profiles (this accounts for
          differences in absolute intensity between animals)
        - use the first animal in the stack as the reference profile
        - for all animals:

           - compare a forward and reverse profile to the reference profile (using the
             cosine-similarity metric)
           - keep either the forward or reverse profile accordingly

        - finally, determine the location of the peaks in the *average* profile

            - reverse all profiles if necessary (this will be necessary if the first
              animal happens to be reversed)

    """
    data = intensity_data

    ref_data = data.sel(
        wavelength=reference_wavelength,
        pair=reference_pair,
        timepoint=reference_timepoint,
    )
    ref_profile = ref_data.isel(animal=0).data

    ref_vecs = np.tile(ref_profile, (data.animal.size, 1))
    unflipped = data.sel(
        wavelength=reference_wavelength,
        pair=reference_pair,
        timepoint=reference_timepoint,
    ).data
    flipped = np.fliplr(unflipped)

    # cosine-similarity measurements
    should_flip = (
        spatial.distance.cdist(ref_vecs, unflipped, "cosine")[0, :]
        > spatial.distance.cdist(ref_vecs, flipped, "cosine")[0, :]
    )

    # Do the actual flip
    # position needs to be reindexed, otherwise xarray freaks out
    intensity_data[should_flip] = np.flip(
        intensity_data[should_flip].values, axis=intensity_data.get_axis_num("position")
    )
    intensity_data = intensity_data.reindex(
        position=np.linspace(0, 1, intensity_data.position.size)
    )

    mean_intensity = trim_profile(
        np.mean(
            intensity_data.sel(
                wavelength=reference_wavelength,
                pair=reference_pair,
                timepoint=reference_timepoint,
            ),
            axis=0,
        ).data,
        threshold=2000,
        new_length=100,
    )

    # parameters found experimentally
    # TODO these could use some tweaking
    peaks, _ = signal.find_peaks(
        mean_intensity, distance=0.2 * len(mean_intensity), prominence=200, wlen=10
    )

    if len(peaks) < 2:
        return intensity_data

    if peaks[0] < len(mean_intensity) - peaks[1]:
        intensity_data = np.flip(
            intensity_data, axis=intensity_data.get_axis_num("position")
        )

    return intensity_data


def summarize_over_regions(
    data: xr.DataArray,
    regions: Dict,
    rescale: bool = True,
    value_name: str = "value",
    pointwise: Union[bool, str] = False,
    **redox_params,
):
    if pointwise == "both":
        # recursively call this function for pointwise=T/F and concat the results
        return pd.concat(
            [
                summarize_over_regions(
                    data, regions, rescale, value_name, pointwise=False
                ),
                summarize_over_regions(
                    data, regions, rescale, value_name, pointwise=True
                ),
            ]
        )

    if rescale:
        regions = utils.scale_region_boundaries(regions, data.shape[-1])

    # Ensure that derived wavelengths are present
    data = utils.add_derived_wavelengths(data, **redox_params)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        all_region_data = []

        for _, bounds in regions.items():
            if isinstance(bounds, (int, float)):
                all_region_data.append(data.interp(position=bounds))
            else:
                all_region_data.append(
                    data.sel(position=slice(bounds[0], bounds[1])).mean(
                        dim="position", skipna=True
                    )
                )

    region_data = xr.concat(all_region_data, pd.Index(regions.keys(), name="region"))
    region_data = region_data.assign_attrs(**data.attrs)

    region_data.loc[dict(wavelength="r")] = region_data.sel(
        wavelength=redox_params["ratio_numerator"]
    ) / region_data.sel(wavelength=redox_params["ratio_denominator"])
    region_data.loc[dict(wavelength="oxd")] = r_to_oxd(
        region_data.sel(wavelength="r"),
        r_min=redox_params["r_min"],
        r_max=redox_params["r_max"],
        instrument_factor=redox_params["instrument_factor"],
    )
    region_data.loc[dict(wavelength="e")] = oxd_to_redox_potential(
        region_data.sel(wavelength="oxd"),
        midpoint_potential=redox_params["midpoint_potential"],
        z=redox_params["z"],
        temperature=redox_params["temperature"],
    )

    df = to_dataframe(region_data, value_name)
    df["pointwise"] = pointwise

    try:
        df.set_index(["experiment_id"], append=True, inplace=True)
    except ValueError:
        pass

    return df


def smooth_profile_data(
    profile_data: Union[np.ndarray, xr.DataArray],
    lambda_: float = 100.0,
    order: float = 4.0,
    n_basis: float = 100.0,
    n_deriv=0.0,
    eng=None,
):
    """
    Smooth profile data by fitting smoothing B-splines

    Implemented in MATLAB as smooth_profiles
    """

    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed. Skipping smoothing.")
        return profile_data

    if eng is None:
        eng = matlab.engine.start_matlab()

    resample_resolution = profile_data.position.size

    return xr.apply_ufunc(
        lambda x: np.array(
            eng.smooth_profiles(
                matlab.double(x.tolist()),
                resample_resolution,
                n_basis,
                order,
                lambda_,
                n_deriv,
            )
        ).T,
        profile_data,
        input_core_dims=[["position"]],
        output_core_dims=[["position"]],
        vectorize=True,
    )


def standardize_profiles(
    profile_data: xr.DataArray, redox_params, eng=None, **reg_kwargs,
) -> (xr.DataArray, xr.DataArray):

    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Skipping registration.")
        return profile_data

    if eng is None:
        eng = matlab.engine.start_matlab()

    std_profile_data = profile_data.copy()
    std_warp_data = profile_data.copy().isel(wavelength=0)

    for tp in profile_data.timepoint:
        for pair in profile_data.pair:
            selector = dict(timepoint=tp, pair=pair)

            data = std_profile_data.sel(**selector)
            std_data_, warp_ = _standardize_profiles(
                profile_data=data, eng=eng, **reg_kwargs
            )
            std_profile_data.loc[selector] = std_data_
            std_warp_data.loc[selector] = np.array(warp_).T

    std_profile_data = std_profile_data.assign_attrs(**reg_kwargs)
    std_profile_data = utils.add_derived_wavelengths(std_profile_data, **redox_params)

    return std_profile_data, std_warp_data


def _standardize_profiles(
    profile_data: xr.DataArray,
    eng=None,
    n_deriv: float = 2.0,
    rough_lambda: float = 0.01,
    rough_n_breaks: float = 300.0,
    rough_order: float = 4.0,
    smooth_lambda: float = 10.0 ** 2,
    smooth_n_breaks: float = 100.0,
    smooth_order: float = 4.0,
    warp_lambda: float = 5.0e3,
    warp_n_basis: float = 30.0,
    warp_order: float = 4.0,
    ratio_numerator: str = "410",
    ratio_denominator: str = "470",
) -> (xr.DataArray, xr.DataArray):
    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Skipping registration.")
        return profile_data

    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()

    i410 = matlab.double(profile_data.sel(wavelength=ratio_numerator).values.tolist())
    i470 = matlab.double(profile_data.sel(wavelength=ratio_denominator).values.tolist())

    resample_resolution = float(profile_data.position.size)

    r410, r470, warp_data = eng.pop_register(
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

    reg_profile_data.loc[dict(wavelength=ratio_numerator)] = r410
    reg_profile_data.loc[dict(wavelength=ratio_denominator)] = r470

    return reg_profile_data, warp_data


def channel_register(
    profile_data: xr.DataArray, redox_params: dict, reg_params: dict, eng=None,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Perform channel-registration on the given profile data

    Parameters
    ----------
    profile_data
        the data to register
    redox_params
        the redox parameters
    reg_params
        the registration parameters
    eng
        the MATLAB engine (optional)

    Returns
    -------
    reg_data: xr.DataArray
        the registered data
    warp_data: xr.DataArray
        the warp functions used to register the data

    """

    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Skipping registration.")
        return profile_data, None

    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()
    # get rid of the wavelength dimension for warp data
    warp_data = profile_data.copy().isel(wavelength=0)

    for pair in profile_data.pair:
        for timepoint in profile_data.timepoint:
            reg, warp_ = _channel_register(
                profile_data.sel(pair=pair, timepoint=timepoint), eng=eng, **reg_params
            )
            reg_profile_data.loc[dict(pair=pair, timepoint=timepoint)] = reg
            warp_data.loc[dict(pair=pair, timepoint=timepoint)] = np.array(warp_).T

    reg_profile_data = utils.add_derived_wavelengths(reg_profile_data, **redox_params)

    return reg_profile_data, warp_data


def _channel_register(
    profile_data: xr.DataArray,
    eng=None,
    n_deriv: float = 2.0,
    rough_lambda: float = 0.01,
    rough_n_breaks: float = 300.0,
    rough_order: float = 4.0,
    smooth_lambda: float = 100.0,
    smooth_n_breaks: float = 100.0,
    smooth_order: float = 4.0,
    warp_lambda: float = 5.0e3,
    warp_n_basis: float = 30.0,
    warp_order: float = 4.0,
) -> (np.ndarray, np.ndarray):
    """
    Register the 470nm channel into the 410nm channel profile data
    
    Parameters
    ----------
    profile_data : xr.DataArray
        The data to register. Must be maximally 3-dimensional. The dimensions are (animal, wavelength, position along profile).
    eng: matlab.engine.MatlabEngine
        The MATLAB engine to use for registration. If None, a new engine is created.
    """

    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Skipping registration.")
        return profile_data

    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()

    i410 = matlab.double(profile_data.sel(wavelength="410").values.tolist())
    i470 = matlab.double(profile_data.sel(wavelength="470").values.tolist())

    resample_resolution = float(profile_data.position.size)

    # Call the MATLAB subroutine
    r410, r470, warp_data = eng.channel_register(
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

    return reg_profile_data, warp_data


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
    data_reversed = data.reindex(position=list(reversed(data.position)))
    l_bound = (data.sel(wavelength=ref_wvl) >= thresh).argmax(dim="position").data - 1
    r_bound = (
        prof_len
        - (data_reversed.sel(wavelength=ref_wvl) >= thresh).argmax(dim="position").data
    ) - 1
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

    for i, img_idx in enumerate(intensity_data.animal):
        for wvl_idx in range(intensity_data.wavelength.size):
            wvl = intensity_data.wavelength.data[wvl_idx]
            if "tl" not in wvl.lower():
                for pair in range(intensity_data.pair.size):
                    for tp in intensity_data.timepoint.values:
                        selector = dict(
                            wavelength=wvl, pair=pair, animal=img_idx, timepoint=tp
                        )
                        data = intensity_data.sel(selector).data
                        l_i, r_i = l[i, tp, pair], r[i, tp, pair]
                        try:
                            trimmed = data[l_i:r_i]
                            new_xs = np.linspace(
                                0, len(trimmed), intensity_data.position.size
                            )
                            old_xs = np.arange(0, len(trimmed))
                            resized = np.interp(new_xs, old_xs, trimmed)

                            trimmed_intensity_data.loc[selector] = resized
                        except ValueError:
                            logging.warning(
                                f"trim boundaries close ({np.abs(r_i - l_i)}) for (animal: {i}, wvl: {wvl}, pair: {pair}) - skipping trimming this animal"
                            )

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
    # We can get NaN ratios because of background subtraction, this is expected
    # so we suppress the warnings here
    with np.errstate(invalid="ignore"):
        return midpoint_potential - (
            8314.462 * (273.15 + temperature) / (z * 96485.3415)
        ) * np.log((1 - oxd) / oxd)
