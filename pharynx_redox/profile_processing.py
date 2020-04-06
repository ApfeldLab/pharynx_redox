import collections
from dataclasses import dataclass
from typing import Union, Dict
import logging
import warnings

import scipy
import numpy as np
import xarray as xr
import pandas as pd
from scipy import spatial, signal
from sklearn.preprocessing import scale
from tqdm.auto import tqdm

from pharynx_redox import utils
from pharynx_redox import constants

from numba import vectorize, int64


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
    reference_pair: optional
        the pair to calculate the alignment for
    reference_wavelength: optional
        the wavelength to calculate the alignment for

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

    # do the actual cosine-similarity measurements
    should_flip = (
        spatial.distance.cdist(ref_vecs, unflipped, "cosine")[0, :]
        > spatial.distance.cdist(ref_vecs, flipped, "cosine")[0, :]
    )

    # position needs to be reindexed, otherwise xarray freaks out
    # Axis=4 because that is the index of `position`
    intensity_data[should_flip] = np.flip(intensity_data[should_flip], axis=4).reindex(
        position=np.arange(intensity_data.position.size)
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
    peaks, _ = signal.find_peaks(
        mean_intensity, distance=0.2 * len(mean_intensity), prominence=200, wlen=10
    )

    if len(peaks) < 2:
        return intensity_data

    if peaks[0] < len(mean_intensity) - peaks[1]:
        intensity_data = np.flip(intensity_data, axis=3)

    return intensity_data


@vectorize
def get_mvmt(a, b):
    if a == False:
        return b
    if b == False:
        return a
    else:
        return np.nan


def summarize_over_regions_old(
    data: xr.DataArray,
    regions: dict,
    value_name: str = None,
    rescale: bool = False,
    add_attrs: bool = True,
):
    """
    Summarize the profile data over region boundaries, storing the resultant data in a 
    pandas DataFrame.

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
    add_attrs
        add attributes from each observation to the table

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

    dfs = []
    for region, bounds in regions.items():
        reg_df = data[dict(position=range(bounds[0], bounds[1]))].mean(
            dim="position", skipna=True
        )
        try:
            reg_df = reg_df.to_pandas().to_frame()
        except:
            reg_df = pd.DataFrame({value_name: [reg_df.values[()]]})

        reg_df["region"] = region
        reg_df.rename({0: value_name}, inplace=True, axis="columns")
        reg_df = reg_df.reset_index()

        try:
            reg_df["strain"] = data.strain.values
        except AttributeError:
            logging.warning("No strain info in data, not adding to summarization table")

        try:
            reg_df["time"] = data.time.values
        except AttributeError:
            logging.warning("No timestamp in data, not adding to summarization table")

        if add_attrs:
            for attr, val in data.attrs.items():
                reg_df[attr] = val

        for region in ["posterior", "anterior", "sides_of_tip", "tip"]:
            try:
                reg_df[f"mvmt-{region}"] = data[f"mvmt-{region}"]
            except KeyError:
                pass

        dfs.append(reg_df)
    return pd.concat(dfs)


def summarize_over_regions(
    data: xr.DataArray,
    regions: Dict,
    rescale: bool = True,
    value_name: str = "value",
    pointwise: Union[bool, str] = False,
    ratio_numerator: str = "410",
    ratio_denominator: str = "470",
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
    try:
        data = utils.add_derived_wavelengths(
            data, numerator=ratio_numerator, denominator=ratio_denominator
        )
    except ValueError:
        pass

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        region_data = xr.concat(
            [
                data[dict(position=range(bounds[0], bounds[1]))].mean(
                    dim="position", skipna=True
                )
                for _, bounds in regions.items()
            ],
            pd.Index(regions.keys(), name="region"),
        )
    region_data = region_data.assign_attrs(**data.attrs)

    if not pointwise:
        try:
            region_data.loc[dict(wavelength="r")] = region_data.sel(
                wavelength=ratio_numerator
            ) / region_data.sel(wavelength=ratio_denominator)
            region_data.loc[dict(wavelength="oxd")] = r_to_oxd(
                region_data.sel(wavelength="r"),
                r_min=data.r_min,
                r_max=data.r_max,
                instrument_factor=data.instrument_factor,
            )
            region_data.loc[dict(wavelength="e")] = oxd_to_redox_potential(
                region_data.sel(wavelength="oxd"),
                midpoint_potential=data.midpoint_potential,
                z=data.z,
                temperature=data.temperature,
            )
        except ValueError:
            pass

    df = to_dataframe(region_data, value_name)
    df["pointwise"] = pointwise

    try:
        df.set_index(["experiment_id",], append=True, inplace=True)
    except ValueError:
        pass

    return df


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

    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed. Skipping smoothing.")
        return profile_data

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


def _register_profiles_pop(
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
) -> xr.DataArray:
    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Skipping registration.")
        return profile_data

    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()
    warp_data = []

    i410 = matlab.double(profile_data.sel(wavelength="410").values.tolist())
    i470 = matlab.double(profile_data.sel(wavelength="470").values.tolist())

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

    reg_profile_data.loc[dict(wavelength="410")] = r410
    reg_profile_data.loc[dict(wavelength="470")] = r470

    return reg_profile_data, warp_data


def register_profiles_pop(
    profile_data: xr.DataArray,
    eng=None,
    progress_bar=False,
    ratio_numerator="410",
    ratio_denominator="470",
    **reg_kwargs,
) -> xr.DataArray:
    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Skipping registration.")
        return profile_data

    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()
    warp_data = []
    disable_progress = not progress_bar
    for tp in tqdm(profile_data.timepoint, disable=disable_progress, desc="timepoint"):
        for pair in tqdm(profile_data.pair, disable=disable_progress, desc="pair"):
            selector = dict(timepoint=tp, pair=pair)

            data = reg_profile_data.sel(**selector)
            reg_data, warp_data = _register_profiles_pop(
                profile_data=data, eng=eng, **reg_kwargs
            )
            reg_profile_data.loc[selector] = reg_data
    reg_profile_data = reg_profile_data.assign_attrs(**reg_kwargs)
    reg_profile_data = utils.add_derived_wavelengths(
        reg_profile_data, numerator=ratio_numerator, denominator=ratio_denominator
    )
    return reg_profile_data, warp_data


def channel_register(
    profile_data: xr.DataArray,
    eng=None,
    ratio_numerator="410",
    ratio_denominator="470",
    **reg_params,
) -> xr.DataArray:

    try:
        import matlab.engine
    except ImportError:
        logging.warn("MATLAB engine not installed! Skipping registration.")
        return profile_data

    if eng is None:
        eng = matlab.engine.start_matlab()

    reg_profile_data = profile_data.copy()
    warp_data = []

    for pair in profile_data.pair:
        for timepoint in profile_data.timepoint:
            reg, warp_ = _channel_register(
                profile_data.sel(pair=pair, timepoint=timepoint), eng=eng, **reg_params
            )
            reg_profile_data.loc[dict(pair=pair, timepoint=timepoint)] = reg
            # TODO keep track of warp data

    reg_profile_data = utils.add_derived_wavelengths(
        reg_profile_data, numerator=ratio_numerator, denominator=ratio_denominator
    )
    return reg_profile_data, warp_data


def _channel_register(
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
    # axis=3 b/c that is where `position` is after selecting wavelength
    l_bound = np.argmax(data.sel(wavelength=ref_wvl) >= thresh, axis=3).data - 1
    r_bound = (
        prof_len
        - np.argmax(
            np.flip(data.sel(wavelength=ref_wvl), axis=2) >= thresh, axis=3
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
                            logging.warn(
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
