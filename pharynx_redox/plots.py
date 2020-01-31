"""
Miscellaneous and sundry plotting functions for to please your visual cortex 
"""

from pprint import pformat
from typing import Dict, Tuple, Union

from tqdm.auto import tqdm

from . import profile_processing as pp
from . import data_analysis as da
from . import utils

import matplotlib.colors
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from matplotlib import cm, gridspec, colors, image
from matplotlib.gridspec import GridSpec
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.weightstats import DescrStatsW
from matplotlib.backends.backend_pdf import PdfPages


def plot_stage_layout(
    image_data: xr.DataArray, pair: int = 0
) -> sns.axisgrid.FacetGrid:
    """
    Shows a scatter plot where each point is an animal located on the imaging stage and
    the points are colored by strain.

    A useful visualization to make sure that the strain map is accurate.

    .. image:: _static/plot_stage_layout.png
    
    Parameters
    ----------
    image_data : xr.DataArray
        The image data acquired by metamorph.
    pair : int
        The image pair to display
    
    Returns
    -------
    seaborn.axisgrid.FacetGrid
        The grid object returned by seaborns's lmplot 
    
    See Also
    --------
    pharynx_io.load_images
    seaborn.lmplot
    """
    df = pd.DataFrame(
        dict(
            stage_x=image_data.sel(wavelength="410", pair=1).stage_x,
            stage_y=image_data.sel(wavelength="410", pair=1).stage_y,
            strain=image_data.sel(wavelength="410", pair=1).strain,
        )
    )

    return sns.lmplot(x="stage_x", y="stage_y", data=df, hue="strain", fit_reg=False)


def cdf_plot(data, *args, **kwargs):
    """
    Plot a CDF, compatible with Seaborn's FacetGrid

    data
        1-D vector of numbers to plot the CDF of
    *args
        ignored
    **kwargs
        keyword arguments passed onto ``plt.step``
    """
    ecdf = sm.distributions.ECDF(data)
    x = np.linspace(min(data), max(data))
    y = ecdf(x)
    plt.step(x, y, **kwargs)


def add_regions_to_axis(
    ax,
    regions: dict,
    label_dist_bottom_percent: float = 0.03,
    label_x_offset_percent: float = 0.005,
    alpha: float = 0.03,
    hide_labels: bool = False,
    xs=None,
    color="black",
    **kwargs,
):
    """
    TODO: Documentation

    Parameters
    ----------
    ax
        the axis to add the regions to
    regions
        the region dictionary, formatted as such::

                {
                    'pm3': [1, 10],
                    'pm4': [12, 30],
                    ...
                }

    label_dist_bottom_percent
        the distance from the bottom of the axis that the region labels should be placed, expressed as a percentage of the axis height
    label_x_offset_percent
        the distance from the left of the region annotation, expressed as a percentage of the axis length
    alpha
        the opacity of the region annotations (0 = transparent, 1=opaque)
    hide_labels
        if True, does not add labels to regions
    kwargs
        these will be passed onto ``ax.axvspan``

    """
    min_y, max_y = ax.get_ylim()
    min_x, max_x = ax.get_xlim()

    text_y = ((max_y - min_y) * label_dist_bottom_percent) + min_y

    text_x_offset = (max_x - min_x) * label_x_offset_percent

    for region, bounds in regions.items():
        ax.axvspan(bounds[0], bounds[1], alpha=alpha, color=color, linewidth=0, **kwargs)
        if not hide_labels:
            ax.annotate(region, xy=(bounds[0] + text_x_offset, text_y))


def plot_profile_avg_with_bounds(
    data, ax=None, confint_alpha=0.05, label=None, xs=None, **kwargs
):
    """
    TODO: Documentation

    Parameters
    ----------
    data
    ax
    confint_alpha
    label
    kwargs

    Returns
    -------

    """
    if ax is None:
        _, ax = plt.subplots()

    if xs is not None:
        ax.plot(xs, np.nanmean(data, axis=0), label=label, **kwargs)
    else:
        ax.plot(np.nanmean(data, axis=0), label=label, **kwargs)
    lower, upper = DescrStatsW(data).tconfint_mean(alpha=confint_alpha)
    if xs is None:
        xs = np.arange(len(lower))

    kwargs.pop("linestyle", None)
    kwargs.pop("linewidth", None)
    kwargs.pop("lw", None)
    ax.fill_between(xs, lower, upper, alpha=0.3, lw=0, **kwargs)

    return ax

def plot_profile_avg_with_sem_bounds(
    data, ax=None, label=None, xs=None, **kwargs
):
    """
    TODO: Documentation

    Parameters
    ----------
    data
    ax
    label
    kwargs

    Returns
    -------

    """
    if ax is None:
        _, ax = plt.subplots()

    mean = np.nanmean(data, axis=0)
    if xs is not None:
        ax.plot(xs, mean, label=label, **kwargs)
    else:
        ax.plot(mean, label=label, **kwargs)
    
    sem = stats.sem(data)

    lower, upper = mean-sem, mean+sem
    if xs is None:
        xs = np.arange(len(lower))

    kwargs.pop("linestyle", None)
    kwargs.pop("linewidth", None)
    kwargs.pop("lw", None)
    ax.fill_between(xs, lower, upper, alpha=0.3, lw=0, **kwargs)

    return ax


def plot_profile_avg(data, ax=None, label=None, xs=None, **kwargs):
    if ax is None:
        _, ax = plt.subplots()

    if xs is not None:
        ax.plot(xs, np.nanmean(data, axis=0), label=label, **kwargs)
    else:
        ax.plot(np.nanmean(data, axis=0), label=label, **kwargs)

    return ax


def imshow_ratio_normed(
    ratio_img,
    fl_img,
    profile_data=None,
    prob=0.999,
    cmap="coolwarm",
    r_min=None,
    r_max=None,
    i_min=0,
    i_max=None,
    clip=True,
    ax=None,
    colorbar=False,
    colorbar_kwargs_dict={},
    **imshow_kwargs,
):
    """
    Show the given ratio image, first converting to HSV and setting the "V" (value) channel to be the given (normalized) intensity image

    Parameters
    ----------
    ratio_img
        the ratio image to display
    fl_img
        the fluorescent intensity image with which to "value-correct" the ratio image. A good choice here is the max value of both intensity channels used in the ratio.
    profile_data
        the midline profile data corresponding to the ratio image. This is used to center and to choose min/max values for the ratio colormap.
    prob
        The "confidence interval" around the center of the ratio values to include in the colormap. For example, 0.95 translates to a min/max of mean(ratio) +/- (1.96*std(ratio))
    cmap
        The colormap used to display the ratio image. Diverging colormaps are a good choice here (default is RdBu_r).
    r_min
        The minimum value for the ratio colormap. If None, uses the `prob` parameter (see its description), and requires `profile_data`.
    r_max
        The maximum value for the ratio colormap. If None, uses the `prob` parameter (see its description), and requires `profile_data`.
    i_min
        The intensity to map to 0 in the value channel
    i_max
        The intensity to map to 1 in the value channel
    clip
        Whether or not the value channel should be clipped to [0, 1] before converting back to RGB. Leaving this as True is a sane default.
    ax
        If given, the image is plotted on this axis. If ``None``, this function uses the pyplot interface.
    imshow_args
        keyword arguments that will be passed along to the ``imshow`` function
    """

    # Convert ratio to RGB
    if profile_data is None:
        if (r_min is None) or (r_max is None):
            raise ValueError("r_min and r_max must be set if profile_data is not given")
    else:
        r_mean = np.mean(profile_data)
        r_std = np.std(profile_data)
        Z = stats.norm.ppf(
            prob
        )  # this converts probability -> "Z" value (e.g. 0.95 -> 1.96)
        window = r_std * Z
        r_min_ = r_mean - window
        r_max_ = r_mean + window
    if r_min is None:
        r_min = r_min_
    if r_max is None:
        r_max = r_max_

    norm_ratio = colors.Normalize(vmin=r_min, vmax=r_max)

    cmap = cm.get_cmap(cmap)

    img_rgba = cmap(norm_ratio(ratio_img))

    # Now convert RGB to HSV, using intensity image as the "V" (value)
    if i_max is None:
        i_max = np.max(fl_img)

    norm_fl = colors.Normalize(vmin=i_min, vmax=i_max, clip=clip)
    hsv_img = colors.rgb_to_hsv(img_rgba[:, :, :3])  # ignore the "alpha" channel
    hsv_img[:, :, -1] = norm_fl(fl_img)

    if ax is None:
        ax = plt.gca()

    # Convert HSV back to RGB and plot
    img_rgba = colors.hsv_to_rgb(hsv_img)

    im = ax.imshow(img_rgba, **imshow_kwargs)
    im.cmap = cmap
    im.norm = norm_ratio

    if colorbar:
        add_img_colorbar(ax, **colorbar_kwargs_dict)

    return im


def add_img_colorbar(ax, position="right", size="5%", pad=0.05, **colorbar_kwargs):
    try:
        axes_img = [
            obj for obj in ax.get_children() if isinstance(obj, image.AxesImage)
        ][0]
    except IndexError:
        raise ValueError(
            "No image found in axis children. This method only works for axes with images."
        )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes(position, size=size, pad=pad)
    return plt.colorbar(
        cm.ScalarMappable(norm=axes_img.norm, cmap=axes_img.cmap),
        cax=cax,
        **colorbar_kwargs,
    )


def registration_diagnostic_plot_stack(fl, raw_prof, reg_prof, filepath, **reg_params):
    with PdfPages(filepath) as pdf:
        for i in range(fl.spec.size):
            f = registration_diagnostic_plot(fl, raw_prof, reg_prof, i, **reg_params)
            pdf.savefig(f)
            plt.close(f)


def registration_diagnostic_plot(fl, raw_prof, reg_prof, idx, **params) -> plt.Figure:
    if "pair" in fl.dims:
        n_pairs = fl.pair.size
    else:
        n_pairs = 1

    fig_width, fig_height = 12, 12
    fig_scale = 1

    fig, axes = plt.subplots(
        4, n_pairs, figsize=(fig_scale * fig_width, fig_scale * fig_height)
    )

    colors = {"410": "tab:blue", "470": "tab:orange"}

    ###############
    # RATIO IMAGE #
    ###############
    for pair in range(n_pairs):
        ax = axes[0, pair]
        colormap_profile_buffer = 10
        imshow_ratio_normed(
            fl.sel(wavelength="r", pair=pair)[idx],
            fl.sel(wavelength="410", pair=pair)[idx],
            profile_data=raw_prof.sel(wavelength="r", pair=pair)[idx][
                colormap_profile_buffer:-colormap_profile_buffer
            ],
            prob=0.999,
            i_max=2000,
            colorbar=True,
            ax=ax,
        )

        ax.set_xlim(45, 120)
        ax.set_ylim(50, 80)

    #####################
    # INTENSITY PROFILE #
    #####################
    for pair in range(n_pairs):
        # Intensity Profile
        ax = axes[1, pair]
        xs = np.linspace(0, 1, raw_prof.position.size)
        # REG 410 ("unregistered, but smooth")
        ax.scatter(
            xs,
            raw_prof.sel(wavelength="410", pair=pair)[idx],
            color="k",
            label="raw 410",
            s=1,
        )
        # RAW 410 (to test hyper-fine smoothing)
        ax.plot(
            xs,
            reg_prof.sel(wavelength="410", pair=pair)[idx],
            lw=1,
            color=colors["410"],
            label="r410",
            linestyle="-",
        )
        # RAW 470
        ax.plot(
            xs,
            raw_prof.sel(wavelength="470", pair=pair)[idx],
            linestyle="-",
            lw=1,
            color=colors["470"],
            label="470",
        )
        # REG 470
        ax.plot(
            xs,
            reg_prof.sel(wavelength="470", pair=pair)[idx],
            linestyle="--",
            lw=1,
            color=colors["470"],
            label="r470",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 2.5e4)

        # UNREGISTERED d(I) PROFILE
        ax = ax.twinx()
        # RAW 410 (to test hyper-fine smoothing)
        ax.plot(
            xs,
            reg_prof.sel(wavelength="410", pair=pair)[idx].differentiate(
                coord="position"
            ),
            lw=1,
            color=colors["410"],
            label="r410",
            linestyle="-",
        )
        # RAW 470
        ax.plot(
            xs,
            raw_prof.sel(wavelength="470", pair=pair)[idx].differentiate(
                coord="position"
            ),
            linestyle="-",
            lw=1,
            color=colors["470"],
            label="470",
        )
        # REG 470
        ax.plot(
            xs,
            reg_prof.sel(wavelength="470", pair=pair)[idx].differentiate(
                coord="position"
            ),
            linestyle="--",
            lw=1,
            color=colors["470"],
            label="r470",
        )
        ax.axhline(0, linestyle="--", color="lightgray", lw=1)

        ax.set_xlim(0, 1)
        ax.set_ylim(-5e3, 1e3)
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

    #################
    # RATIO PROFILE #
    #################

    # get appropriate y-limits for ratios
    buffer = int(raw_prof.position.size * 0.30)
    r_min = raw_prof.sel(wavelength="r")[..., buffer:-buffer].min()
    r_max = raw_prof.sel(wavelength="r")[..., buffer:-buffer].max()

    for pair in range(n_pairs):
        # Intensity Profile
        ax = axes[2, pair]
        xs = np.linspace(0, 1, raw_prof.position.size)

        ax.plot(
            xs, raw_prof.sel(wavelength="r", pair=pair)[idx], label="raw", color="k"
        )
        ax.plot(
            xs,
            reg_prof.sel(wavelength="r", pair=pair)[idx],
            label="raw",
            color="tab:red",
            linestyle="-",
        )
        ax.set_xlim(0, 1)
        autoscale_percent_buffer = 0.0
        ax.set_ylim(
            r_min - (r_min * autoscale_percent_buffer),
            r_max + (r_max * autoscale_percent_buffer),
        )
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

    ax = axes[3, 0]
    ax.set_axis_off()

    #################
    # Parameter box #
    #################
    textstr = "\n".join(f"{k}={v}" for k, v in params.items())
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    ax.text(
        0.01,
        0.97,
        textstr,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=props,
    )

    axes[3, 1].set_axis_off()

    plt.tight_layout()

    return fig
