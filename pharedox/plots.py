"""
Miscellaneous and sundry plotting functions for to please your visual cortex 
"""

import typing
import warnings
from pathlib import Path
from typing import Dict, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from matplotlib import cm, colors, gridspec, image, transforms
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
from skimage.measure import label, regionprops
from statsmodels.stats.weightstats import DescrStatsW
from tqdm.auto import tqdm

from pharedox import constants
from pharedox import data_analysis as da


def imshow_r_stack(
    imgs: xr.DataArray,
    profile_data: xr.DataArray,
    output_dir: Union[str, Path],
    per_animal_cmap: bool = True,
    fl_wvl: str = "410",
    cmap: str = "coolwarm",
    width: int = 80,
    height: int = 30,
    colorbar=True,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    center = (np.array(imgs.shape[-2:]) / 2).astype(np.int)
    wpad = int(width / 2)
    hpad = int(height / 2)
    for tp in tqdm(imgs.timepoint.values, leave=False, desc="timepoint"):
        for pair in tqdm(imgs.pair.values, leave=False, desc="pair"):
            filepath = output_dir.joinpath(f"timepoint={tp}_pair={pair}.pdf")
            with PdfPages(filepath) as pdf:
                i = 0
                for animal in tqdm(imgs.animal.values, desc="animal", leave=False):
                    fig, ax = plt.subplots()
                    selector = dict(animal=animal, timepoint=tp, pair=pair)
                    im, cbar = imshow_ratio_normed(
                        imgs.sel(wavelength="r", **selector),
                        imgs.sel(wavelength=fl_wvl, **selector),
                        profile_data=profile_data.sel(wavelength="r", **selector),
                        prob=0.999,
                        colorbar=colorbar,
                        i_min=0,
                        i_max=3000,
                        cmap=cmap,
                        ax=ax,
                    )
                    ax.set_xlim(center[1] - wpad, center[1] + wpad)
                    ax.set_ylim(center[0] - hpad, center[0] + hpad)
                    ax.set_title(str(selector))
                    pdf.savefig()
                    if (i % 20) == 0:
                        plt.close("all")
                    i += 1


def generate_wvl_pair_timepoint_profile_plots(data: xr.DataArray, ignored_wvls=None):
    """
    For each wavelength and pair in the given data, this function plots a line plot with
    each color representing a unique strain. The line is the mean value across animals
    for that strain, and the shaded regions are the 95% confidence intervals
    
    Parameters
    ----------
    data
    ignored_wvls
    """
    if ignored_wvls is None:
        ignored_wvls = ["TL"]
    strains = np.unique(data.strain.values)
    cmap = plt.get_cmap("Set2")
    colormap = dict(zip(strains, cmap.colors))

    wvls = list(map(lambda x: x.lower(), data.wavelength.values))
    for wvl in ignored_wvls:
        try:
            wvls.remove(wvl.lower())
        except ValueError:
            continue

    for wvl in wvls:
        for pair in data.pair.values:
            for tp in data.timepoint.values:
                fig, ax = plt.subplots()
                for strain in strains:
                    strain_data = data.where(data["strain"] == strain, drop=True)
                    ax.plot(
                        strain_data.sel(wavelength=wvl, pair=pair, timepoint=tp).T,
                        color=colormap[strain],
                        alpha=0.5,
                    )

                title = f"wavelength = {wvl} ; pair = {pair} ; timepoint = {tp}"
                ax.set_title(title)
                ax.legend(
                    [
                        plt.Line2D([0], [0], color=color, lw=4)
                        for color in cmap.colors[: len(strains)]
                    ],
                    strains,
                )
                yield title, fig


def generate_avg_wvl_pair_profile_plots(
    data: xr.DataArray, ignored_wvls: typing.List[str] = None
):
    """
    For each wavelength and pair in the given data, this function plots a line plot with
    each color representing a unique strain. The line is the mean value across animals
    for that strain, and the shaded regions are the 95% confidence intervals
    
    Parameters
    ----------
    ignored_wvls
    data : [type]
        [description]
    """
    if ignored_wvls is None:
        ignored_wvls = ["TL"]
    strains = np.unique(data.strain.values)
    cmap = plt.get_cmap("Set2")
    colormap = dict(zip(strains, cmap.colors))
    wvls = list(map(lambda x: x.lower(), data.wavelength.values))
    for wvl in ignored_wvls:
        try:
            wvls.remove(wvl.lower())
        except ValueError:
            continue
    for wvl in wvls:
        for pair in data.pair.values:
            for tp in data.timepoint.values:
                fig, ax = plt.subplots()
                for strain in np.unique(data.strain.values):
                    strain_data = data.where(data["strain"] == strain, drop=True)
                    plot_profile_avg_with_bounds(
                        strain_data.sel(wavelength=wvl, pair=pair, timepoint=tp),
                        label=strain,
                        ax=ax,
                        color=colormap[strain],
                    )
                title = f"wavelength = {wvl} ; pair = {pair} ; timepoint = {tp}"
                ax.set_title(title)
                ax.legend()
                yield title, fig


def plot_err_with_region_summaries(
    data: xr.DataArray,
    measure_regions: Dict,
    display_regions=None,
    ax=None,
    profile_color="black",
    label=None,
):
    st_color = "k"
    mv_color = "tab:red"

    if ax is None:
        _, ax = plt.subplots()

    if display_regions is None:
        display_regions = measure_regions

    df = da.fold_v_point_table(data, measure_regions)
    df_avgs = df.reset_index().groupby("region").agg(["mean", "sem"]).reset_index()

    xs = np.linspace(0, 1, data.position.size)

    plot_profile_avg_with_sem_bounds(
        100 * da.fold_error(data), xs=xs, ax=ax, color=profile_color, label=label
    )

    for region, region_err_mean, region_err_sem in zip(
        df_avgs["region"],
        df_avgs["fold_error_region"][1]["mean"],
        df_avgs["fold_error_region"][1]["sem"],
    ):
        try:
            ax.axhline(
                100 * region_err_mean,
                *display_regions[region],
                color=profile_color,
                alpha=1,
                lw=2,
                solid_capstyle="butt",
            )
            ax.errorbar(
                x=np.mean(display_regions[region]),
                y=100 * region_err_mean,
                yerr=100 * region_err_sem,
                color=profile_color,
                elinewidth=0.5,
                capsize=1,
                capthick=0.5,
            )
        except:
            continue

    ax.set_xlim(0, 1)

    add_regions_to_axis(
        ax, display_regions, alpha=0.3, hide_labels=True, skip=["medial_axis"]
    )
    ax.set_xlabel("position along midline")


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
    io.load_tiff_as_hyperstack
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


def ecdf_(data):
    """ Compute ECDF """
    x = np.sort(data)
    n = x.size
    y = np.arange(1, n + 1) / n
    return (x, y)


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
    # ecdf = sm.distributions.ECDF(data)
    x, y = ecdf_(data)
    # x = np.linspace(min(data), max(data), len(data))
    # logging.debug(x)
    # y = ecdf(x)
    # plt.step(x, y, **kwargs)
    # sns.kdeplot(data, cumulative=True)
    plt.step(x, y, **kwargs)


def add_regions_to_axis(
    ax,
    regions: dict,
    skip=None,
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
    skip
        the regions to skip plotting
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
    if skip is None:
        skip = []
    min_y, max_y = ax.get_ylim()
    min_x, max_x = ax.get_xlim()

    text_y = ((max_y - min_y) * label_dist_bottom_percent) + min_y

    text_x_offset = (max_x - min_x) * label_x_offset_percent

    for region, bounds in regions.items():
        if region in skip:
            continue
        ax.axvspan(
            bounds[0], bounds[1], alpha=alpha, color=color, linewidth=0, **kwargs
        )
        if not hide_labels:
            ax.annotate(region, xy=(bounds[0] + text_x_offset, text_y))


def add_region_bars_to_axis(
    ax, regions, skip=None, bar_height=8, bar_width=1, fontsize=3
):
    if skip is None:
        skip = []

    for region, region_bounds in regions.items():
        if region in skip:
            continue

        yy = -0.01

        ax.annotate(
            "",
            xy=(region_bounds[0], yy),
            xycoords=("data", "axes fraction"),
            xytext=(region_bounds[1], yy),
            textcoords=("data", "axes fraction"),
            arrowprops=dict(
                arrowstyle="-",
                connectionstyle=f"bar,armA=-{bar_height},armB=-{bar_height},fraction=0.0",
                capstyle="butt",
                joinstyle="miter",
                lw=bar_width,
            ),
            annotation_clip=False,
        )
        ax.annotate(
            region,
            xy=((region_bounds[0] + region_bounds[1]) / 2, yy - 0.08),
            xycoords=("data", "axes fraction"),
            ha="center",
            fontsize=fontsize,
        )

    ax.xaxis.labelpad = 25


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
        ax = plt.gca()

    if xs is None:
        try:
            # if the data is an xr.DataArray
            xs = data.position
        except ValueError:
            # if it's a numpy array
            xs = np.arange(len(data))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax.plot(xs, np.nanmean(data, axis=0), label=label, **kwargs)

        lower, upper = DescrStatsW(data).tconfint_mean(alpha=confint_alpha)

    kwargs.pop("linestyle", None)
    kwargs.pop("linewidth", None)
    kwargs.pop("lw", None)
    ax.fill_between(xs, lower, upper, alpha=0.3, lw=0, **kwargs)

    return ax


def plot_profile_avg_with_sem_bounds(data, ax=None, label=None, xs=None, **kwargs):
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

    with np.errstate(invalid="ignore"):
        sem = stats.sem(data)

        lower, upper = mean - sem, mean + sem
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


def imgs_to_rgb(
    imgs,
    r_min,
    r_max,
    cmap="coolwarm",
    i_min=0,
    i_max=None,
    i_wvls=["410", "470"],
    ratio_numerator="410",
    ratio_denominator="470",
):
    if i_max is None:
        i_max = np.max(imgs.sel(wavelength=["410", "470"]))

    try:
        R = imgs.sel(wavelength="R")
    except KeyError:
        R = imgs.sel(wavelength=ratio_numerator) / imgs.sel(
            wavelength=ratio_denominator
        )

    norm_ratio = colors.Normalize(vmin=r_min, vmax=r_max)
    cmap = cm.get_cmap(cmap)
    img_rgba = cmap(norm_ratio(R))
    norm_fl = colors.Normalize(vmin=i_min, vmax=i_max, clip=True)
    hsv_img = colors.rgb_to_hsv(img_rgba[..., :3])  # ignore the "alpha" channel
    hsv_img[..., -1] = norm_fl(imgs.sel(wavelength=i_wvls).max(dim="wavelength"))
    img_rgba = colors.hsv_to_rgb(hsv_img)
    return img_rgba


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
    Show the given ratio image, first converting to HSV and setting the "V" (value) 
    channel to be the given (normalized) intensity image

    Parameters
    ----------
    ratio_img
        the ratio image to display
    fl_img
        the fluorescent intensity image with which to "value-correct" the ratio image. 
        A good choice here is the max value of both intensity channels used in the 
        ratio.
    profile_data
        the midline profile data corresponding to the ratio image. This is used to 
        center and to choose min/max values for the ratio colormap.
    prob
        The "confidence interval" around the center of the ratio values to include in 
        the colormap. For example, 0.95 translates to a min/max of 
        mean(ratio) +/- (1.96*std(ratio))
    cmap
        The colormap used to display the ratio image. Diverging colormaps are a good 
        choice here (default is RdBu_r).
    r_min
        The minimum value for the ratio colormap. If None, uses the `prob` parameter 
        (see its description), and requires `profile_data`.
    r_max
        The maximum value for the ratio colormap. If None, uses the `prob` parameter 
        (see its description), and requires `profile_data`.
    i_min
        The intensity to map to 0 in the value channel
    i_max
        The intensity to map to 1 in the value channel
    clip
        Whether or not the value channel should be clipped to [0, 1] before converting 
        back to RGB. Leaving this as True is a sane default.
    ax
        If given, the image is plotted on this axis. If ``None``, this function uses the
        pyplot interface.
    colorbar
        show the colorbar or not
    imshow_args
        keyword arguments that will be passed along to the ``imshow`` function
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Convert ratio to RGB
        if profile_data is None:
            if (r_min is None) or (r_max is None):
                raise ValueError(
                    "r_min and r_max must be set if profile_data is not given"
                )
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
            cbar = add_img_colorbar(ax, **colorbar_kwargs_dict)
            return im, cbar
        else:
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


def plot_pharynx_R_imgs(
    img: xr.DataArray,
    mask: xr.DataArray,
    crop: bool = True,
    crop_pad: int = 10,
    cmap_normalization: str = "frame",
    cmap: str = "coolwarm",
    fig_kwargs=None,
):
    """
    Generate a figure which has ratio images broken up by timepoint and pair
    
    Parameters
    ----------
    img : xr.DataArray
        The image to display. Should contain a single animal, and the `r` and `410`
        wavelengths.
    mask : xr.DataArray
        The mask with which a ROI will be used for calculated the average R value of
        the pharynx
    crop : bool, optional
        Whether the image should be cropped, by default True
    crop_pad : int, optional
        The padding for the crop, as number of pixels on each side of
        the bounding box surrounding the pharynx, by default 10
    cmap_normalization : str, optional
        How the colormap should be normalized, by default "frame". "frame" means each
        timepoint and pair will be normalized separately
    cmap : str, optional
        The colormap to use, by default "coolwarm"
    fig_kwargs : [type], optional
        Keyword arguments to be passed to `matplotlib.pyplot.figure`, by default None
    
    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    NotImplementedError
        [description]
    """
    if "animal" in img.dims:
        raise ValueError(
            f"Image must contain single animal. Given stack contains {img.animal.size} animals"
        )

    valid_cmap_normalizations = ["frame", "animal"]
    if cmap_normalization not in valid_cmap_normalizations:
        raise ValueError(
            f"`cmap_normaliztion` must be one of {valid_cmap_normalizations} (given <{cmap_normalization}>)"
        )

    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(ncols=img.pair.size, nrows=img.timepoint.size, figure=fig)

    for i, tp in enumerate(img.timepoint.values):
        for j, pair in enumerate(img.pair.values):
            ax = fig.add_subplot(gs[i, j])
            ax.set_title(f"timepoint = {tp} | Pair={pair}")
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

            R = img.sel(wavelength="r", timepoint=tp, pair=pair).values
            I = img.sel(wavelength="410", timepoint=tp, pair=pair).values
            M = mask.sel(timepoint=tp, pair=pair).astype(bool).values

            P_r = R[M]

            if cmap_normalization == "frame":
                r_min = np.mean(P_r) - 1.96 * np.std(P_r)
                r_max = np.mean(P_r) + 1.96 * np.std(P_r)

                i_max = 0.30 * regionprops(label(M), intensity_image=I)[0].max_intensity
            if cmap_normalization == "animal":
                raise NotImplementedError

            if crop:
                rp = regionprops(label(M), intensity_image=R)[0]
                (min_row, min_col, max_row, max_col) = rp.bbox
                ax.set_ylim(max_row + crop_pad, min_row - crop_pad)
                ax.set_xlim(min_col - crop_pad, max_col + crop_pad)

            imshow_ratio_normed(
                R, I, cmap=cmap, i_min=0, i_max=i_max, r_min=r_min, r_max=r_max, ax=ax
            )
    return fig


def plot_multiple_pop_errors(
    data_and_labels,
    ylim=None,
    xlim=None,
    add_regions=True,
    figsize=(20, 10),
    dpi=100,
    regions=constants.untrimmed_regions,
):
    """Plot multiple error profiles and their corresponding labels

    Parameters
    ----------
    data_and_labels : List[Tuple(xr.DataArray, str)]
        A list of (data, label), one for each data set to plot
    ylim : Tuple(float, float), optional
        The y limits of the plot, by default None
    xlim : Tuple(float, float), optional
        The x-limits of the plot, by default None
    add_regions : bool, optional
        Whether or not to add regions to the plot, by default True
    figsize : Tuple(float, floa), optional
        The figure size, by default (20, 10)
    dpi : int, optional
        The DPI of the plot, by default 100
    regions : dict, optional
        The regions to plot, if enabled, by default constants.untrimmed_regions

    Returns
    -------
    fig, ax
        The figure and axis object
    """

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    for data, label in data_and_labels:
        xs = np.linspace(0, 1, data.position.size)
        plot_profile_avg_with_bounds(
            da.fold_error(data.sel(timepoint=0)), xs=xs, ax=ax, label=label
        )

    if ylim:
        ax.set_ylim(*ylim)
    if xlim:
        ax.set_xlim(*xlim)

    if add_regions:
        add_regions_to_axis(ax, regions, alpha=0.3, hide_labels=True)

    ax.set_ylabel("Absolute Error (%)")
    ax.set_xlabel("position along midline")

    ax.legend(loc="lower left")

    return fig, ax


def plot_multiple_pop_wvl(
    data_and_labels,
    wvl="r",
    ylim=None,
    xlim=None,
    add_regions=True,
    figsize=(20, 10),
    dpi=100,
    regions=constants.untrimmed_regions,
):
    """Plot multiple error profiles and their corresponding labels

    Parameters
    ----------
    data_and_labels : List[Tuple(xr.DataArray, str)]
        A list of (data, label), one for each data set to plot
    ylim : Tuple(float, float), optional
        The y limits of the plot, by default None
    xlim : Tuple(float, float), optional
        The x-limits of the plot, by default None
    add_regions : bool, optional
        Whether or not to add regions to the plot, by default True
    figsize : Tuple(float, floa), optional
        The figure size, by default (20, 10)
    dpi : int, optional
        The DPI of the plot, by default 100
    regions : dict, optional
        The regions to plot, if enabled, by default constants.untrimmed_regions

    Returns
    -------
    fig, ax
        The figure and axis object
    """

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    for data, label in data_and_labels:
        # xs = np.linspace(0, 1, data.position.size)
        xs = data.position
        plot_profile_avg_with_bounds(
            data.sel(timepoint=0, wavelength=wvl).mean(dim="pair"),
            xs=xs,
            ax=ax,
            label=label,
        )

    if ylim:
        ax.set_ylim(*ylim)
    if xlim:
        ax.set_xlim(*xlim)

    if add_regions:
        add_regions_to_axis(ax, regions, alpha=0.3, hide_labels=True)

    ax.set_ylabel(wvl)
    ax.set_xlabel("position along midline")

    ax.legend(loc="lower left")

    return fig, ax
