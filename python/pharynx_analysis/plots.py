import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from statsmodels.stats.weightstats import DescrStatsW


def plot_individual_profile_data_by_strain_and_pair(profile_data, cmin, cmax, cmap_name, linewidth=1, ylim=None,
                                                    alpha=1, figsize=None, cmap_boundary_trim=0):
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
    fig, axes = plt.subplots(n_strains, 2, figsize=figsize)

    for strain, ax in zip(strains, axes):
        for i in range(2):
            data = profile_data.isel(pair=i).sel(strain=strain)
            color_data = data[:, cmap_boundary_trim:profile_data.position.size - cmap_boundary_trim]
            means_normed = norm(color_data.mean(dim='position').data)
            colors = cmap(means_normed)
            colors[:, 3] = alpha
            for animal_idx in range(data.shape[0]):
                ax[i].plot(data[animal_idx], c=colors[animal_idx], linewidth=linewidth)
                ax[i].set_ylim(ylim)
            ax[i].set_title(f'{strain} (Pair {i})')
            if cmap_boundary_trim > 0:
                boundary_line_args = {
                    'color': 'k', 'linewidth': 1, 'alpha': 0.1,
                    'linestyle': '--'

                }
                ax[i].axvline(cmap_boundary_trim, **boundary_line_args)
                ax[i].axvline(profile_data.position.size - cmap_boundary_trim, **boundary_line_args)

    plt.tight_layout()


def plot_average_by_strain_and_pair(profile_data, ylim=None, sup_title=None):
    strains = np.unique(profile_data.strain.data)

    if 'pair' in profile_data.dims:
        n_pairs = profile_data.pair.size

        fig, axes = plt.subplots(n_pairs, 1, figsize=(10, 10))

        for i, ax in zip(range(n_pairs), axes):
            for strain in strains:
                data = profile_data.isel(pair=i).sel(strain=strain)
                ax.plot(np.mean(data, axis=0), label=strain)
                lower, upper = DescrStatsW(data).tconfint_mean()
                xs = np.arange(len(lower))
                ax.fill_between(xs, lower, upper, alpha=0.4)
                ax.legend()
                ax.set_title(f'Pair: {i}')
                if ylim:
                    ax.set_ylim(ylim)
    else:
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        for strain in strains:
            data = profile_data.sel(strain=strain)
            ax.plot(np.mean(data, axis=0), label=strain)
            lower, upper = DescrStatsW(data).tconfint_mean()
            xs = np.arange(len(lower))
            ax.fill_between(xs, lower, upper, alpha=0.4)
            ax.legend()
            if ylim:
                ax.set_ylim(ylim)

    if sup_title:
        fig.suptitle(sup_title)
