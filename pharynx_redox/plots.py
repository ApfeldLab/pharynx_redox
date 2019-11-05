from pprint import pformat
from typing import Dict, Tuple, Union

from tqdm.auto import tqdm

from pharynx_redox.profile_processing import scale_by_wvl
from pharynx_redox import data_analysis as da
from pharynx_redox import utils

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


def plot_paired_experiment_summary(experiment):
    """
    TODO: Documentation

    Parameters
    ----------
    experiment

    Returns
    -------

    """
    fig, axes = plt.subplots(5, 2, figsize=(20, 20))
    plot_average_by_strain_and_pair(
        experiment.trimmed_intensity_data.sel(wavelength="410"),
        regions=experiment.scaled_regions,
        axes=[axes[0, 0], axes[0, 1]],
        title="410",
        legend=True,
    )
    plot_average_by_strain_and_pair(
        experiment.trimmed_intensity_data.sel(wavelength="470"),
        regions=experiment.scaled_regions,
        axes=[axes[1, 0], axes[1, 1]],
        title="470",
    )
    plot_average_by_strain_and_pair(
        experiment.trimmed_intensity_data.sel(wavelength="r"),
        regions=experiment.scaled_regions,
        axes=[axes[2, 0], axes[2, 1]],
        title="410/470",
    )
    plot_average_by_strain_and_pair(
        experiment.trimmed_intensity_data.sel(wavelength="oxd"),
        regions=experiment.scaled_regions,
        axes=[axes[3, 0], axes[3, 1]],
        title="OxD",
    )
    plot_average_by_strain_and_pair(
        experiment.trimmed_intensity_data.sel(wavelength="e"),
        regions=experiment.scaled_regions,
        axes=[axes[4, 0], axes[4, 1]],
        title="E",
    )
    plt.tight_layout()
    return fig, axes


def plot_individual_profile_data_by_strain_and_pair(
    profile_data,
    cmin,
    cmax,
    cmap_name,
    linewidth=1,
    ylim=None,
    alpha=1,
    figsize=None,
    cmap_boundary_trim=0,
):
    """
    TODO: Documentation

    Parameters
    ----------
    profile_data
    cmin
    cmax
    cmap_name
    linewidth
    ylim
    alpha
    figsize : object
    cmap_boundary_trim

    Returns
    -------

    """
    cmap = cm.get_cmap(name=cmap_name)
    norm = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax)

    strains = np.unique(profile_data.strain.data)

    if ylim is None:
        ylim = [profile_data.min(), profile_data.max()]

    n_strains = len(strains)
    if "pair" in profile_data.dims:
        _, axes = plt.subplots(n_strains, 2, figsize=figsize)

        for strain, ax in zip(strains, axes):
            for i in range(2):
                data = profile_data.isel(pair=i).sel(strain=strain)
                color_data = data[
                    :,
                    cmap_boundary_trim : profile_data.position.size
                    - cmap_boundary_trim,
                ]
                means_normed = norm(color_data.mean(dim="position").data)
                colors = cmap(means_normed)
                colors[:, 3] = alpha
                for animal_idx in range(data.shape[0]):
                    ax[i].plot(
                        data[animal_idx], c=colors[animal_idx], linewidth=linewidth
                    )
                    ax[i].set_ylim(ylim)
                ax[i].set_title(f"{strain} (Pair {i})")
                if cmap_boundary_trim > 0:
                    boundary_line_args = {
                        "color": "k",
                        "linewidth": 1,
                        "alpha": 0.1,
                        "linestyle": "--",
                    }
                    ax[i].axvline(cmap_boundary_trim, **boundary_line_args)
                    ax[i].axvline(
                        profile_data.position.size - cmap_boundary_trim,
                        **boundary_line_args,
                    )
    else:
        fig, ax = plt.subplots(figsize=figsize)
        for strain in strains:
            data = profile_data.sel(strain=strain)
            color_data = data[
                :, cmap_boundary_trim : profile_data.position.size - cmap_boundary_trim
            ]
            means_normed = norm(color_data.mean(dim="position").data)
            colors = cmap(means_normed)
            colors[:, 3] = alpha
            for animal_idx in range(data.shape[0]):
                ax.plot(data[animal_idx], c=colors[animal_idx], linewidth=linewidth)
                ax.set_ylim(ylim)
            if cmap_boundary_trim > 0:
                boundary_line_args = {
                    "color": "k",
                    "linewidth": 1,
                    "alpha": 0.1,
                    "linestyle": "--",
                }
                ax.axvline(cmap_boundary_trim, **boundary_line_args)
                ax.axvline(
                    profile_data.position.size - cmap_boundary_trim,
                    **boundary_line_args,
                )

    plt.tight_layout()


def plot_average_by_strain_and_pair(
    profile_data, ylim=None, regions=None, axes=None, title=None, legend=False
):
    """
    TODO: Documentation

    Parameters
    ----------
    profile_data
    ylim
    regions
    axes
    title
    legend

    Returns
    -------

    """
    strains = np.unique(profile_data.strain.data)

    if "pair" in profile_data.dims:
        n_pairs = profile_data.pair.size

        if axes is None:
            _, axes = plt.subplots(n_pairs, 1, figsize=(10, 10))

        for i, ax in zip(range(n_pairs), axes):

            for strain in strains:
                data = profile_data.isel(pair=i).sel(strain=strain)
                ax.plot(
                    np.mean(data, axis=0),
                    label=f"{strain} (n={profile_data.strain.sel(strain=strain).size})",
                )
                lower, upper = DescrStatsW(data).tconfint_mean()
                xs = np.arange(len(lower))
                ax.fill_between(xs, lower, upper, alpha=0.4)
                if legend:
                    ax.legend()
                if title:
                    ax.set_title(f"{title}-{i}")
                else:
                    ax.set_title(f"Pair {i}")
                if ylim:
                    ax.set_ylim(ylim)
            if regions:
                add_regions_to_axis(ax, regions)

        return axes

    else:
        if axes is None:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        else:
            ax = axes
        for strain in strains:
            data = profile_data.sel(strain=strain)
            ax.plot(np.mean(data, axis=0), label=strain)
            lower, upper = DescrStatsW(data).tconfint_mean()
            xs = np.arange(len(lower))
            ax.fill_between(xs, lower, upper, alpha=0.4)
            if legend:
                ax.legend()
            ax.set_title(title)
            if ylim:
                ax.set_ylim(ylim)
            if regions:
                add_regions_to_axis(ax, regions)

        return ax


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
        ax.axvspan(bounds[0], bounds[1], alpha=alpha, color=color, **kwargs)
        if not hide_labels:
            ax.annotate(region, xy=(bounds[0] + text_x_offset, text_y))


def plot_multi_profile_by_wvl_and_pair(data, color="royalblue", alpha=0.3):
    """
    TODO: Documentation

    Parameters
    ----------
    data
    color
    alpha

    Returns
    -------

    """
    fl_wvls = list(filter(lambda x: x != "TL", data.wavelength))

    fig, axes = plt.subplots(
        len(fl_wvls), data.pair.size, sharex="col", sharey="row", figsize=(20, 10)
    )
    fig.subplots_adjust(wspace=0, hspace=0.1)
    for row, wvl in enumerate(fl_wvls):
        if wvl != "TL":
            for col, pair in enumerate(data.pair.data):
                axes[row][col].plot(
                    data.sel(wavelength=wvl, pair=pair).T, color=color, alpha=alpha
                )
                axes[row][col].set_title(f"{wvl.data}-{pair}")

    return fig, axes


def plot_profile_avg_by_strain(profile_data, ax_title=None, alpha=0.05, ax=None):
    """
    TODO: Documentation

    Parameters
    ----------
    profile_data
    ax_title
    alpha
    ax

    Returns
    -------

    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))

    ax.set_title(ax_title)
    for strain in np.unique(profile_data.strain.data):
        d = profile_data.sel(strain=strain)
        ax.plot(np.mean(d, axis=0), label=strain)
        lower, upper = DescrStatsW(d).tconfint_mean(alpha=alpha)
        xs = np.arange(len(lower))
        ax.fill_between(xs, lower, upper, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
    plt.subplots_adjust(right=0.7)

    return plt.gcf(), ax


def boxplot_by_strain(y, summary_table):
    """
    TODO: Documentation

    Parameters
    ----------
    y
    summary_table

    Returns
    -------

    """
    return sns.catplot(y=y, x="strain", data=summary_table, kind="box")


def single_animal_diagnostic_plot(
    animal_idx,
    rot_fl,
    midlines,
    profile_data: Union[Tuple[xr.DataArray, xr.DataArray], xr.DataArray],
    wvl_limits: Dict[str, Tuple[float, float]] = None,
):
    """
    TODO: Documentation

    Parameters
    ----------
    animal_idx
    rot_fl
    midlines
    profile_data
    wvl_limits

    Returns
    -------

    """
    fig = plt.figure(constrained_layout=True, figsize=(15, 15))
    grid_spec = GridSpec(5, 4, figure=fig)
    midline_xs = np.arange(40, 120)

    ax_idx = 0
    if type(profile_data) == tuple:
        pairs = profile_data[0].pair.data
        strains = profile_data[0].strain.data
    else:
        pairs = profile_data.pair.data
        strains = profile_data.strain.data

    # Plot images in first [n_pairs] rows
    for pair in np.arange(pairs.size):
        i410 = rot_fl.sel(wavelength="410", pair=pair).isel(spec=animal_idx)
        i470 = rot_fl.sel(wavelength="470", pair=pair).isel(spec=animal_idx)

        # Plot 410 image
        ax = fig.add_subplot(grid_spec[pair, 0])
        ax.imshow(i410)
        ax.plot(
            midline_xs, midlines[animal_idx]["410"][pair](midline_xs), color="orange"
        )
        ax.set_title(f"410-{pair}")

        # Plot 470 image
        ax = fig.add_subplot(grid_spec[pair, 1])
        ax.imshow(i470)
        ax.plot(midline_xs, midlines[animal_idx]["470"][pair](midline_xs), color="r")
        ax.set_title(f"470-{pair}")

        # Plot Ratio image
        ax = fig.add_subplot(grid_spec[pair, 2])
        ax.imshow(i410 / i470)
        ax.plot(
            midline_xs,
            midlines[animal_idx]["410"][pair](midline_xs),
            color="orange",
            label="410",
            alpha=0.5,
        )
        ax.plot(
            midline_xs,
            midlines[animal_idx]["470"][pair](midline_xs),
            color="r",
            label="470",
            alpha=0.5,
        )
        ax.set_title(f"(410/470)-{pair}")
        ax.legend()
        ax_idx += 1

    if type(profile_data) == tuple:
        starts = [0, 2]
        ends = [2, 4]
        times = [0, 1]

    else:
        starts = [0]
        ends = [4]
        times = [0]

    ax_idx_save = ax_idx
    if type(profile_data) != tuple:
        profile_data = (profile_data,)
    for i in times:
        # Plot intensity profiles
        prof_data = profile_data[i]
        ax = fig.add_subplot(grid_spec[ax_idx, starts[times[i]] : ends[times[i]]])
        for pair in np.arange(pairs.size):
            ax.plot(
                prof_data.sel(wavelength="410", pair=pair).isel(spec=animal_idx),
                label=f"410-{pair}",
            )
            ax.plot(
                prof_data.sel(wavelength="470", pair=pair).isel(spec=animal_idx),
                label=f"470-{pair}",
            )
            if wvl_limits:
                ax.set_ylim(wvl_limits.get("intensity", (0, 1.6e4)))
        ax.legend()
        ax.set_title("Intensity")
        ax_idx += 1

        # Plot R, OxD, and E
        for wvl in ["r", "oxd", "e"]:
            ax = fig.add_subplot(grid_spec[ax_idx, starts[times[i]] : ends[times[i]]])
            for pair in np.arange(pairs.size):
                ax.plot(
                    prof_data.sel(wavelength=wvl, pair=pair).isel(spec=animal_idx),
                    label=f"{pair}",
                )
                if wvl_limits:
                    ax.set_ylim(wvl_limits.get(wvl))
                else:
                    ax.set_ylim(
                        [
                            np.nanquantile(prof_data.sel(wavelength=wvl).data, 0.01),
                            np.nanquantile(prof_data.sel(wavelength=wvl).data, 0.99),
                        ]
                    )
            ax.legend()
            ax.set_title(wvl)
            ax_idx += 1

        ax_idx = ax_idx_save

    plt.suptitle(f"Animal {animal_idx} ({strains[animal_idx]})")

    return fig


def plot_reg_diagnostic(ex_raw, ex_reg, i, col_w=6, row_h=3, **params):
    """
    TODO: Documentation

    Parameters
    ----------
    ex_raw
    ex_reg
    i
    col_w
    row_h
    params

    Returns
    -------

    """
    n_cols = 3
    n_rows = 4
    f = plt.figure(figsize=(col_w * n_cols, row_h * n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols, figure=f)

    def gridify_axis(ax):
        ax.grid(b=True, which="major", axis="x")

    r_clim = [1, 2]
    rr_clim = [0, 0.2]

    sc_raw = scale_by_wvl(ex_raw.untrimmed_profiles)
    sc_reg = scale_by_wvl(ex_reg.untrimmed_profiles)

    for pair in [0, 1]:
        col = pair
        # Raw Intensity
        ax = f.add_subplot(gs[0, col])
        ax.plot(sc_raw.sel(wavelength="410", pair=pair).isel(spec=i), label="410")
        ax.plot(sc_raw.sel(wavelength="470", pair=pair).isel(spec=i), label="470")
        ax.plot(sc_reg.sel(wavelength="470", pair=pair).isel(spec=i), label="r470")
        ax.legend(loc="upper left")
        ax.set_title(fr"I_{pair} Raw")
        ax.set_ylim([-2, 2])
        gridify_axis(ax)

        # Diff Raw Intensity
        ax = f.add_subplot(gs[1, col])
        ax.plot(
            np.diff(sc_raw.sel(wavelength="410", pair=pair).isel(spec=i)), label="410"
        )
        ax.plot(
            np.diff(sc_raw.sel(wavelength="470", pair=pair).isel(spec=i)), label="470"
        )
        ax.plot(
            np.diff(sc_reg.sel(wavelength="470", pair=pair).isel(spec=i)), label="r470"
        )
        ax.set_title(fr"diff(I_{pair})")
        ax.set_ylim([-0.3, 0.3])
        gridify_axis(ax)

        # (Raw, Reg) Ratios
        ax = f.add_subplot(gs[2, col])
        ax.plot(
            ex_raw.untrimmed_profiles.sel(wavelength="410", pair=pair).isel(spec=i)
            / ex_raw.untrimmed_profiles.sel(wavelength="470", pair=pair).isel(spec=i),
            color="k",
            label="raw",
        )

        ax.plot(
            ex_reg.untrimmed_profiles.sel(wavelength="410", pair=pair).isel(spec=i)
            / ex_reg.untrimmed_profiles.sel(wavelength="470", pair=pair).isel(spec=i),
            color="r",
            label="reg",
        )
        ax.legend(loc="upper right")
        ax.set_ylim([1, 2])
        ax.set_title(fr"410/470 ({pair})")
        gridify_axis(ax)

    # Error
    ax = f.add_subplot(gs[3, 1])
    ax.plot(
        np.abs(
            1
            - (
                ex_raw.untrimmed_profiles.sel(wavelength="410", pair=0).isel(spec=i)
                / ex_raw.untrimmed_profiles.sel(wavelength="470", pair=0).isel(spec=i)
            )
            / (
                ex_raw.untrimmed_profiles.sel(wavelength="410", pair=1).isel(spec=i)
                / ex_raw.untrimmed_profiles.sel(wavelength="470", pair=1).isel(spec=i)
            )
        ),
        label="Raw0/Raw1",
        color="black",
    )
    ax.plot(
        np.abs(
            1
            - (
                ex_reg.untrimmed_profiles.sel(wavelength="410", pair=1).isel(spec=i)
                / ex_reg.untrimmed_profiles.sel(wavelength="470", pair=1).isel(spec=i)
            )
            / (
                ex_raw.untrimmed_profiles.sel(wavelength="410", pair=0).isel(spec=i)
                / ex_raw.untrimmed_profiles.sel(wavelength="470", pair=0).isel(spec=i)
            )
        ),
        label="Reg1/Raw0",
        color="purple",
    )
    gridify_axis(ax)
    ax.set_ylim(rr_clim)
    ax.set_title("|1-(r0/r1)|")
    ax.legend(loc="upper right")

    ax = f.add_subplot(gs[3, 0])
    ax.plot(
        np.abs(
            1
            - (
                ex_raw.untrimmed_profiles.sel(wavelength="410", pair=0).isel(spec=i)
                / ex_raw.untrimmed_profiles.sel(wavelength="470", pair=0).isel(spec=i)
            )
            / (
                ex_raw.untrimmed_profiles.sel(wavelength="410", pair=1).isel(spec=i)
                / ex_raw.untrimmed_profiles.sel(wavelength="470", pair=1).isel(spec=i)
            )
        ),
        label="Raw0/Raw1",
        color="black",
    )
    ax.plot(
        np.abs(
            1
            - (
                ex_reg.untrimmed_profiles.sel(wavelength="410", pair=0).isel(spec=i)
                / ex_reg.untrimmed_profiles.sel(wavelength="470", pair=0).isel(spec=i)
            )
            / (
                ex_raw.untrimmed_profiles.sel(wavelength="410", pair=1).isel(spec=i)
                / ex_raw.untrimmed_profiles.sel(wavelength="470", pair=1).isel(spec=i)
            )
        ),
        label="Reg0/Raw1",
        color="red",
    )
    gridify_axis(ax)
    ax.set_ylim(rr_clim)
    ax.set_title("|1-(r_REG/r_RAW)|")
    ax.legend(loc="upper right")

    # Images
    def clip_im(ax):
        ax.set_xlim(45, 125)
        ax.set_ylim(50, 80)

    ax = f.add_subplot(gs[0, 2])
    R0 = ex_raw.rot_fl.sel(wavelength="410", pair=0).isel(spec=i) / ex_raw.rot_fl.sel(
        wavelength="470", pair=0
    ).isel(spec=i)
    im = ax.imshow(R0)
    im.set_clim(r_clim)
    clip_im(ax)
    ax.set_title("R0")
    ax.plot(*ex_raw.midlines[i]["410"][0].linspace(), color="r", linewidth=3)

    ax = f.add_subplot(gs[1, 2])
    R1 = ex_raw.rot_fl.sel(wavelength="410", pair=1).isel(spec=i) / ex_raw.rot_fl.sel(
        wavelength="470", pair=1
    ).isel(spec=i)
    im = ax.imshow(R1)
    im.set_clim(r_clim)
    clip_im(ax)
    ax.set_title("R1")
    ax.plot(*ex_raw.midlines[i]["410"][1].linspace(), color="r", linewidth=3)

    ax = f.add_subplot(gs[2, 2])
    RR = np.abs(1 - (R0 / R1))
    im = ax.imshow(RR)
    im.set_clim(rr_clim)
    clip_im(ax)
    ax.set_title("|1 - R0/R1|")

    try:
        m0 = ex_raw.movement.loc[pd.IndexSlice[i, 0], :]
        m1 = ex_raw.movement.loc[pd.IndexSlice[i, 1], :]
        plt.suptitle(
            f"{ex_raw.experiment_id} --(Animal {i})-- mvmt(T:{m0['tip']}/{m1['tip']} | SoT:{m0['sides_of_tip']}/{m1['sides_of_tip']} | A:{m0['anterior']}/{m1['anterior']} | P:{m0['posterior']}/{m1['posterior']})"
        )
    except:
        pass

    # Parameter box
    textstr = "\n".join([f"{k}={v}" for k, v in params.items()])
    props = dict(boxstyle="round", facecolor="wheat")
    ax = f.add_subplot(gs[3, 2])
    ax.text(
        0.05,
        0.95,
        textstr,
        transform=ax.transAxes,
        fontsize=14,
        verticalalignment="top",
        bbox=props,
    )
    ax.axis("off")

    f.tight_layout(rect=[0, 0.03, 1, 0.95])

    return f


def save_reg_diagnostics(ex_raw, ex_reg, output_filepath, **params):
    """
    TODO: Documentation

    Parameters
    ----------
    ex_raw
    ex_reg
    output_filepath
    params

    Returns
    -------

    """
    with PdfPages(output_filepath) as pdf:
        for i in tqdm(range(ex_raw.strains.size)):
            f = plot_reg_diagnostic(ex_raw, ex_reg, i, **params)
            pdf.savefig()
            plt.close(f)


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
        fig, ax = plt.subplots()

    if xs is not None:
        ax.plot(xs, np.nanmean(data, axis=0), label=label, **kwargs)
    else:
        ax.plot(np.nanmean(data, axis=0), label=label, **kwargs)
    lower, upper = DescrStatsW(data).tconfint_mean(alpha=confint_alpha)
    if xs is None:
        xs = np.arange(len(lower))

    kwargs.pop("linestyle", None)
    kwargs.pop("linewidth", None)
    ax.fill_between(xs, lower, upper, alpha=0.3, **kwargs)

    return ax


def loose_movement_stratified_plots(
    mvmt: pd.DataFrame,
    raw_prof: xr.DataArray,
    reg_prof: xr.DataArray,
    fname: str = None,
    param_dict: dict = None,
) -> None:
    """
    Plot Movement error according to "loose" movement stratifications

    Parameters
    ----------
    mvmt
        a DataFrame containing the movement calls, expected in the following format:

        ====  ========================  ====================  =============================  ==============================  =================================  ========================  =============================  ==============================  =================================  ========================
          ..  ('experiment', '', '')      ('animal', '', '')    ('movement', 0, 'anterior')    ('movement', 0, 'posterior')    ('movement', 0, 'sides_of_tip')    ('movement', 0, 'tip')    ('movement', 1, 'anterior')    ('movement', 1, 'posterior')    ('movement', 1, 'sides_of_tip')    ('movement', 1, 'tip')
        ====  ========================  ====================  =============================  ==============================  =================================  ========================  =============================  ==============================  =================================  ========================
           0  2017_02_22-HD233_SAY47                       0                              0                               0                                  0                         1                              0                               0                                  0                         0
           1  2017_02_22-HD233_SAY47                       1                              0                               0                                  1                         0                              1                               1                                  1                         1
           2  2017_02_22-HD233_SAY47                       2                              0                               0                                  0                         0                              0                               0                                  0                         0
        ====  ========================  ====================  =============================  ==============================  =================================  ========================  =============================  ==============================  =================================  ========================

        The tuples in the headers indicate the use of a pandas Multi-Index.
        TODO: write function that will munge the original movement file into this format

    raw_prof
        The raw measurement profiles
    reg_prof
        the registered measurement profiles
    fname
        the file name to save the resultant plot. If none, does not save plot
    param_dict
        the relevant parameters used in the analysis

    Returns
    -------

    """
    fig, axes = plt.subplots(1, 2, figsize=(18, 5), sharey="all")

    ylims = (0.007, 0.05)

    # Posterior
    roi = "posterior"

    for ax, roi in zip(axes, ["anterior", "posterior"]):
        ax.set_title(f"{roi} (+ maybe others)")
        ax.set_ylim(*ylims)

        m_0_0, m_0_1, m_1_0, m_1_1 = da.split_by_movement_types(mvmt, roi)
        moving = pd.concat([m_0_1, m_1_0]).drop_duplicates().reset_index(drop=True)

        # stationary
        plot_profile_avg_with_bounds(
            da.get_resid_rr(raw_prof)[m_0_0.index.values],
            ax=ax,
            label=f"Raw, Stationary (n={len(m_0_0.index.values)})",
        )
        plot_profile_avg_with_bounds(
            da.get_resid_rr(reg_prof)[m_0_0.index.values],
            ax=ax,
            label=f"Reg., Stationary (n={len(m_0_0.index.values)})",
        )
        # moving
        plot_profile_avg_with_bounds(
            da.get_resid_rr(raw_prof)[moving.index.values],
            ax=ax,
            label=f"Raw, Moving (n={len(moving.index.values)})",
        )
        plot_profile_avg_with_bounds(
            da.get_resid_rr(reg_prof)[moving.index.values],
            ax=ax,
            label=f"Reg., Moving (n={len(moving.index.values)})",
        )

        ax.legend(loc="upper left")
        ax.set_xlabel("A-P Axis")

    axes[0].set_ylabel("error")
    plt.tight_layout()

    if fname:
        plt.savefig(fname)

    if param_dict:
        param_str = pformat(param_dict)
        plt.suptitle(param_str)


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
            raise ValueError('r_min and r_max must be set if profile_data is not given')
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

    # Convert HSV back to RGB and plot
    img_rgba = colors.hsv_to_rgb(hsv_img)
    if ax is None:
        im = plt.imshow(img_rgba, **imshow_kwargs)
    else:
        im = ax.imshow(img_rgba, **imshow_kwargs)

    im.cmap = cmap
    im.norm = norm_ratio

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

