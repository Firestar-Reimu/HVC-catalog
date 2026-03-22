import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils.aperture import SkyEllipticalAperture


# Plot a 2-panel figure for one cluster: moment 0 image + spectra
def plot_cluster_and_spectrum(
    i,
    cluster_params,
    crafts_velocity_axis_list,
    crafts_spectra,
    v_fit_list,
    y_model_list,
    hi4pi_velocity_axis_list,
    hi4pi_spectra,
    VMIN,
    VMAX,
    VPEAK,
    ID=None,
    *,
    figsize=(10, 4),
    cmap="viridis_r",
    contour_level=0.5,
    contour_color="white",
    ellipse_edgecolor="red",
    grid_kwargs=None,
    legend_kwargs=None,
    colorbar_label="Moment 0 (K km/s)",
    cluster_title="Cluster",
    spectrum_title="Spectrum",
    suptitle_fmt="{i}: {id}",
    save_dir=None,
    save_name_fmt="{i}_{id}.png",
    dpi=300,
    show=True,
):
    """
    Plot a 2-panel figure for one cluster:
      Left: moment 0 image + contour + fitted ellipse
      Right: spectra (CRAFTS + Gaussian Fit + HI4PI) + velocity markers

    Parameters
    ----------
    i : int
        Index of the cluster to plot.
    cluster_params : sequence
        cluster_params[i] should be a tuple:
        (moment_0_crop, extent, xc, yc, a, b, theta)
    crafts_velocity_axis_list, crafts_spectra, v_fit_list, y_model_list,
    hi4pi_velocity_axis_list, hi4pi_spectra : list/array-like
        Spectral data containers indexed by i.
    VMIN, VMAX, VPEAK : array-like
        Velocity marker arrays indexed by i.
    ID : array-like or None
        Optional ID list/array indexed by i (used in suptitle and filename).
    save_dir : str or None
        If provided, save the figure into this directory.
    show : bool
        If True, call plt.show(). If False, do not show (useful in scripts).

    Returns
    -------
    fig, ax : matplotlib Figure and Axes array (shape (2,))
    """
    moment_0_crop, extent, xc, yc, a, b, theta = cluster_params[i]

    if grid_kwargs is None:
        grid_kwargs = dict(color="white", ls="--", alpha=0.5)
    if legend_kwargs is None:
        legend_kwargs = dict(fontsize=8)

    fig, ax = plt.subplots(1, 2, figsize=figsize)

    # --- Left panel: moment 0 map ---
    im = ax[0].imshow(
        moment_0_crop,
        extent=extent,
        origin="lower",
        cmap=cmap,
        interpolation="none",
    )

    # contour at contour_level * max
    level_val = contour_level * np.max(moment_0_crop)
    ax[0].contour(
        moment_0_crop,
        levels=[level_val],
        extent=extent,
        colors=contour_color,
        linewidths=1,
        linestyles="solid",
    )

    ellipse = Ellipse(
        (xc, yc),
        width=2 * a,
        height=2 * b,
        angle=theta,
        edgecolor=ellipse_edgecolor,
        facecolor="none",
        lw=2,
    )
    ax[0].add_patch(ellipse)

    ax[0].set_xlabel("Right Ascension (FK5)")
    ax[0].set_ylabel("Declination (FK5)")
    ax[0].tick_params(axis="both", which="both", labelsize=10)
    ax[0].set_title(cluster_title)
    ax[0].grid(**grid_kwargs)

    # colorbar for left panel
    fig.colorbar(im, ax=ax[0], label=colorbar_label)

    # --- Right panel: spectra ---
    ax[1].plot(
        crafts_velocity_axis_list[i],
        crafts_spectra[i],
        lw=0.5,
        alpha=0.5,
        label="CRAFTS",
    )
    ax[1].plot(v_fit_list[i], y_model_list[i], lw=1.5, label="Gaussian Fit")
    ax[1].plot(hi4pi_velocity_axis_list[i], hi4pi_spectra[i], lw=1, label="HI4PI")

    ax[1].axvline(VMIN[i], color="cyan", linestyle="--")
    ax[1].axvline(VMAX[i], color="cyan", linestyle="--")
    ax[1].axvline(VPEAK[i], color="magenta", linestyle="-.")

    ax[1].set_xlabel("Velocity LSR (km/s)")
    ax[1].set_ylabel("Flux Density (K)")
    ax[1].yaxis.set_label_position("right")
    ax[1].yaxis.tick_right()
    ax[1].set_title(spectrum_title)
    ax[1].legend(**legend_kwargs)

    # --- Titles / layout ---
    fig.tight_layout()

    the_id = None
    if ID is not None:
        the_id = ID[i]

    if the_id is not None:
        fig.suptitle(suptitle_fmt.format(i=i, id=the_id), fontsize=16)
    else:
        fig.suptitle(str(i), fontsize=16)

    # Need to re-adjust after suptitle sometimes
    fig.subplots_adjust(top=0.85)

    # --- Save if needed ---
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        safe_id = the_id if the_id is not None else "unknown"
        save_name = save_name_fmt.format(i=i, id=safe_id)
        out_path = os.path.join(save_dir, save_name)
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, ax


# Plot moment maps for multiple clusters by indices
def plot_clusters_and_spectra(
    indices,
    *,
    show=True,
    **kwargs,
):
    """
    Plot moment maps for multiple clusters by indices.
    Returns list of (fig, ax).
    """
    figs = []
    for i in indices:
        fig, ax = plot_cluster_and_spectrum(i, show=show, **kwargs)
        figs.append((fig, ax))
    return figs


def plot_moments(
    moments,
    titles,
    ellipse_params,
    wcs,
    mapping,
    eps,
    minPts,
    MIN_SNR,
    MIN_SIZE,
    MIN_FWHM,
    figsize=(10, 10),
    index=True,
    colorbar=None,
    suptitle=None,
    save_fig=False,
    annotate_offset=(0, 0),
    annotate_text_kwargs=None,
):
    n = len(moments)

    fig, axes = plt.subplots(
        n,
        1,
        figsize=figsize,
        sharex=True,
        sharey=True,
        layout="compressed",
        subplot_kw={"projection": wcs.celestial},
    )

    axes = np.atleast_1d(axes)

    base_text_kwargs = dict(
        color="black",
        fontsize=8,
        ha="center",
        va="center",
        zorder=1000,
        bbox=None,
    )
    if annotate_text_kwargs:
        base_text_kwargs.update(annotate_text_kwargs)

    for ax, moment, title in zip(axes, moments, titles):
        data = getattr(moment, "value", moment)

        norm = None if np.all(np.isnan(data)) else simple_norm(data, percent=95)
        im = ax.imshow(data, norm=norm, cmap="viridis_r", origin="lower")

        for i in range(1, len(mapping) + 1):
            idx = mapping[i]
            xc, yc, a, b, theta = ellipse_params[idx]

            aper = SkyEllipticalAperture(
                SkyCoord(ra=xc * u.deg, dec=yc * u.deg),
                a * u.deg,
                b * u.deg,
                theta=(90 - theta) * u.deg,
            )
            aper_to_pixels = aper.to_pixel(wcs.celestial)
            aper_to_pixels.plot(ax=ax, color="red", lw=1)

            pos = aper_to_pixels.positions

            if index is True:
                ax.annotate(
                    str(i),
                    xy=(pos[0], pos[1]),
                    xycoords="data",
                    xytext=annotate_offset,
                    textcoords="offset points",
                    **base_text_kwargs,
                )

        ax.grid(linestyle="--")
        ax.set_title(title)

        lon = ax.coords[0]
        lat = ax.coords[1]
        lon.set_major_formatter("dd")
        lat.set_major_formatter("dd")
        lon.set_axislabel("RA")
        lat.set_axislabel("Dec")

        if colorbar is not None:
            fig.colorbar(im, ax=ax, label=colorbar[titles.index(title)])
        else:
            fig.colorbar(im, ax=ax)

    if suptitle is not None:
        fig.suptitle(suptitle)
    else:
        fig.suptitle(
            f"eps = {np.round(eps, 3)}, minPts = {minPts}, MIN_SNR = {MIN_SNR}, "
            f"MIN_SIZE = {MIN_SIZE}, MIN_FWHM = {MIN_FWHM}"
        )

    if save_fig is True:
        fig.savefig(f"./Figures/{suptitle}.png", dpi=300)

    rect = Rectangle(
        (0, 0),
        1,
        1,
        transform=fig.transFigure,
        facecolor="none",
        edgecolor="black",
        linewidth=1.5,
        linestyle="-",
    )
    fig.add_artist(rect)

    plt.show()
    return fig, axes
