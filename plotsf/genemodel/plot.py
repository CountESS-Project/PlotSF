from typing import Mapping, Sequence, Optional, Any, Union
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.figure
import pandas as pd
import numpy as np


def plot_figure(
    df: pd.DataFrame, genes: Union[Mapping, Sequence[Mapping]]
) -> matplotlib.figure.Figure:
    """Wrapper for single- or multi-gene figures.

    If multiple genes are provided, each will be rendered as its own axis in a vertical stack and each will be labeled
    with a single uppercase letter, as would be suitable for a multi-panel figure in a manuscript.

    Parameters
    ----------
    df: pd.DataFrame
    genes: Union[Mapping, Sequence[Mapping]]

    Returns
    -------
    matplotlib.figure.Figure
        The created figure object.

    """
    if isinstance(genes, Mapping):  # convert single genes into a list of length 1
        genes = [genes]
    fig = plt.figure(figsize=(10, 2 * len(genes)), tight_layout="true")
    spec = gridspec.GridSpec(ncols=1, nrows=len(genes), figure=fig)
    for i, g in enumerate(genes):
        ax = fig.add_subplot(spec[i, 0])
        if len(genes) > 1:
            ax.text(0, 1, chr(ord("A") + i), transform=ax.transAxes, size="xx-large")
        plot_variants(ax, df, g)
        if i == 0:
            fig.legend()
        if i == len(genes) - 1:
            ax.set_xlabel("amino acid position")
    return fig


def plot_gene_body(
    ax: plt.axes,
    gene: Mapping[str, Any],
    gene_model_height: float,
    vertical_offset: float = 0,
) -> None:
    """

    Parameters
    ----------
    ax
    gene
    gene_model_height
    vertical_offset

    Returns
    -------

    """
    gene_model_y = -gene_model_height + vertical_offset
    gene_rect = patches.Rectangle(
        xy=(1, gene_model_y),
        width=gene["length"],
        height=gene_model_height,
        facecolor="white",
        edgecolor="black",
        linewidth=1,
    )
    ax.add_artist(gene_rect)

    if "domains" in gene:
        for domain in gene["domains"]:
            domain_width = domain["end"] - domain["start"] + 1
            rect = patches.Rectangle(
                xy=(domain["start"], gene_model_y),
                width=domain_width,
                height=gene_model_height,
                facecolor=domain["color"],
                edgecolor="black",
                linewidth=1,
            )
            ax.add_artist(rect)
            ax.annotate(
                s=domain["name"],
                xy=(
                    domain["start"] + domain_width / 2,
                    gene_model_y + gene_model_height / 2,
                ),
                ha="center",
                va="center",
                color=domain["textcolor"],
            )

    if "endlabels" in gene:
        for endlabel in gene["endlabels"]:
            # spacing = gene["length"] * 0.01     # TODO: fix so that the labels match up across multiple gene models
            spacing = 10
            if endlabel["side"] == "left":
                endlabel_x = -spacing
                endlabel_ha = "right"
            elif endlabel["side"] == "right":
                endlabel_x = gene["length"] + spacing
                endlabel_ha = "left"
            else:
                raise ValueError("endlabel side must be one of 'left' or 'right'")

            ax.annotate(
                s=endlabel["text"],
                xy=(endlabel_x, gene_model_y + gene_model_height / 2),
                ha=endlabel_ha,
                va="center",
                color=endlabel["textcolor"],
            )


def format_axes_for_gene(
    ax: plt.axes, gene: Mapping[str, Any], gene_model_height: float
) -> None:
    """

    Parameters
    ----------
    ax
    gene
    gene_model_height

    Returns
    -------

    """
    ax.set_xlim((1, gene["length"] + 1))
    ax.set_ylim((-gene_model_height, 1))
    ax.set_title(gene["name"], style="italic")

    ax.get_yaxis().set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)


def plot_lollipops(
    ax: plt.axes,
    positions: Sequence[float],
    heights: Union[float, Sequence[float]],
    **kwargs: Optional,
) -> None:
    """

    Parameters
    ----------
    ax
    positions
    label
    heights
    kwargs: Optional
        Keyword arguments passed to :py:func:`matplotlib.pyplot.scatter`. Useful options include label
        (used for figure legends), symbol, color, edgecolor, and zorder.

    Returns
    -------

    """
    if isinstance(heights, Sequence):
        if len(heights) != len(positions):
            raise ValueError(
                f"unequal number of lollipop heights ({len(heights)}) and positions ({len(positions)})"
            )
    else:
        heights = [heights] * len(positions)

    for p, h in zip(positions, heights):
        ax.plot([p, p], [0.0, h], color="black", linestyle="-", linewidth=1, zorder=1)
    ax.scatter(x=positions, y=heights, **kwargs)


def plot_variants(ax: plt.axes, df: pd.DataFrame, gene: Mapping[str, Any]):
    """

    Parameters
    ----------
    ax: plt.axes
        A matplotlib axes object used to generate the plot.
    df: pd.DataFrame
        A pandas DataFrame containing variant information.
    gene: Mapping[str, Any]
        A mapping (e.g. dictionary) containing gene information, including the name, length, and any domains.

    Returns
    -------
    None

    """
    gene_model_height = 0.3
    vus_positions = df.loc[
        np.logical_and(df["clinsig"].isin(["C", "VUS"]), df["gene"] == gene["name"])
    ]["aa_position"].values
    path_positions = df.loc[
        np.logical_and(df["clinsig"].isin(["P", "LP"]), df["gene"] == gene["name"])
    ]["aa_position"].values
    ben_positions = df.loc[
        np.logical_and(df["clinsig"].isin(["B", "LB"]), df["gene"] == gene["name"])
    ]["aa_position"].values

    # plot the gene
    format_axes_for_gene(ax, gene, gene_model_height)
    plot_gene_body(ax, gene, gene_model_height)

    # plot VUS and conflicting positions
    plot_lollipops(
        ax,
        positions=vus_positions,
        heights=0.5,
        # label=f"VUS/Conflicting (n={len(vus_positions)})",
        label=f"VUS/Conflicting",
        marker="o",
        color="#1b9e77",
        edgecolors="black",
        zorder=2,
    )
    # plot pathogenic variants
    plot_lollipops(
        ax,
        positions=path_positions,
        heights=0.7,
        # label=f"Pathogenic/Likely pathogenic (n={len(path_positions)})",
        label=f"Pathogenic/Likely pathogenic",
        marker="o",
        color="#d95f02",
        edgecolors="black",
        zorder=2,
    )
    plot_lollipops(
        ax,
        positions=ben_positions,
        heights=0.3,
        # label=f"Benign/Likely benign (n={len(ben_positions)})",
        label=f"Benign/Likely benign",
        marker="o",
        color="#7570b3",
        edgecolors="black",
        zorder=2,
    )
