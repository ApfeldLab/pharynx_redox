from typing import Dict, Tuple, Union

import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import xarray as xr
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from statsmodels.stats.weightstats import DescrStatsW


def plot_paired_experiment_summary(experiment):
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
    figsize
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
        fig, axes = plt.subplots(n_strains, 2, figsize=figsize)

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
    strains = np.unique(profile_data.strain.data)

    if "pair" in profile_data.dims:
        n_pairs = profile_data.pair.size

        if axes is None:
            fig, axes = plt.subplots(n_pairs, 1, figsize=(10, 10))

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
    ax, regions, label_dist_bottom_percent=0.03, label_x_offset_percent=0.005, alpha=0.1
):
    min_y, max_y = ax.get_ylim()
    min_x, max_x = ax.get_xlim()
    text_y = ((max_y - min_y) * label_dist_bottom_percent) + min_y

    text_x_offset = (max_x - min_x) * label_x_offset_percent

    for region, bounds in regions.items():
        ax.axvspan(bounds[0], bounds[1], alpha=alpha)
        ax.annotate(region, xy=(bounds[0] + text_x_offset, text_y))


def plot_multi_profile_by_wvl_and_pair(data, color="royalblue", alpha=0.3):
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
    return sns.catplot(y=y, x="strain", data=summary_table, kind="box")


def single_animal_diagnostic_plot(
    animal_idx,
    rot_fl,
    midlines,
    profile_data: Union[Tuple[xr.DataArray, xr.DataArray], xr.DataArray],
    wvl_limits: Dict[str, Tuple[float, float]] = None,
):
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
        i410 = rot_fl.sel(wavelength="410", pair=pair).isel(strain=animal_idx)
        i470 = rot_fl.sel(wavelength="470", pair=pair).isel(strain=animal_idx)

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
                prof_data.sel(wavelength="410", pair=pair).isel(strain=animal_idx),
                label=f"410-{pair}",
            )
            ax.plot(
                prof_data.sel(wavelength="470", pair=pair).isel(strain=animal_idx),
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
                    prof_data.sel(wavelength=wvl, pair=pair).isel(strain=animal_idx),
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


def plot_profile_avg_with_bounds(
    data, ax=None, confint_alpha=0.05, label=None, **kwargs
):
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(np.mean(data, axis=0), label=label, **kwargs)
    lower, upper = DescrStatsW(data).tconfint_mean(alpha=confint_alpha)
    xs = np.arange(len(lower))
    ax.fill_between(xs, lower, upper, alpha=0.3, **kwargs)

    return ax
