import matlab
import matlab.engine
import numpy as np
import xarray as xr
from scipy.interpolate import UnivariateSpline


# noinspection PyUnresolvedReferences
def register_profiles(raw_profile_data, f_n_basis=64, f_order=4, f_deriv_penalty=2, f_lambda=1e-6,
                      warp_n_basis=6, warp_order=4, warp_deriv_penalty=2, warp_lambda=5e3):
    reg_profile_data = raw_profile_data.copy()
    eng = matlab.engine.start_matlab()
    for pair in raw_profile_data.pair.data:
        data410 = matlab.double(raw_profile_data.sel(pair=pair, wavelength='410').values.tolist())
        data470 = matlab.double(raw_profile_data.sel(pair=pair, wavelength='470').values.tolist())

        r410, r470 = map(np.asarray, eng.register_pairs(
            data410, data470,
            np.float(f_n_basis), np.float(f_order), np.float(f_deriv_penalty), np.float(f_lambda),
            np.float(warp_n_basis), np.float(warp_order), np.float(warp_deriv_penalty), np.float(warp_lambda),
            nargout=2))

        reg_profile_data.loc[dict(pair=pair, wavelength='410')] = r410.T
        reg_profile_data.loc[dict(pair=pair, wavelength='470')] = r470.T

    return reg_profile_data


def smooth_profile(profile_data: xr.DataArray, s: int = 1e6):
    """
    Smooth the given profile data using B-splines constrained by the given smoothing factor

    Parameters
    ----------
    profile_data
        the data to smooth
    s
        the smoothing factor (see scipy.UnivariateSpline for more details). Default chosen by eye.

    Returns
    -------
    sm
        The smoothed data

    """
    sm = profile_data.copy()
    xs = np.arange(profile_data.position.size)

    for strain in range(profile_data.strain.size):
        for wvl in profile_data.wavelength:
            for pair in profile_data.pair:
                rough_data = profile_data.sel(wavelength=wvl, pair=pair).isel(strain=strain)
                sm.loc[dict(wavelength=wvl, pair=pair)][strain] = UnivariateSpline(xs, rough_data, s=s)(xs)

    return sm


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

    trimmed = profile[first:last + 1]
    new_xs = np.linspace(0, len(trimmed), new_length)
    old_xs = np.arange(0, len(trimmed))

    return np.interp(new_xs, old_xs, trimmed)


def get_trim_boundaries(data, ref_wvl='410', thresh=2000):
    prof_len = data.position.size
    l_bound = np.argmax(data.sel(wavelength=ref_wvl) >= thresh, axis=2).data - 1
    r_bound = prof_len - np.argmax(np.flip(data.sel(wavelength=ref_wvl), axis=2) >= thresh, axis=2).data
    return l_bound, r_bound


def trim_profiles(intensity_data, threshold, new_length, ref_wvl='410'):
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
            (intensity_data.strain.size, intensity_data.wavelength.size, intensity_data.pair.size, new_length)
        ),
        dims=['strain', 'wavelength', 'pair', 'position'],
        coords={'strain': intensity_data.strain, 'wavelength': intensity_data.wavelength, 'pair': intensity_data.pair}
    )

    l, r = get_trim_boundaries(intensity_data, ref_wvl=ref_wvl, thresh=threshold)

    for img_idx in range(intensity_data.strain.size):
        for wvl_idx in range(intensity_data.wavelength.size):
            wvl = intensity_data.wavelength.data[wvl_idx]
            if 'tl' not in wvl.lower():
                for pair in range(intensity_data.pair.size):
                    data = intensity_data.sel(wavelength=wvl, pair=pair).isel(strain=img_idx).data

                    trimmed = data[l[img_idx, pair]:r[img_idx, pair]]
                    new_xs = np.linspace(0, len(trimmed), new_length)
                    old_xs = np.arange(0, len(trimmed))
                    resized = np.interp(new_xs, old_xs, trimmed)

                    trimmed_intensity_data[img_idx, wvl_idx, pair, :] = resized

    return trimmed_intensity_data


def r_to_oxd(r, r_min=0.852, r_max=6.65, instrument_factor=0.171):
    """

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
    """Convert OxD to redox potential

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
    return midpoint_potential - (8314.462 * (273.15 + temperature) / (z * 96485.3415)) * np.log((1 - oxd) / oxd)
