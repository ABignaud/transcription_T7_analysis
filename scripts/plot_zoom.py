#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.hic as bch
import bacchus.io as bcio
import bacchus.plot as bcp
import cooler
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
from os.path import join
import pandas as pd
from scipy import ndimage
import seaborn as sns

# Import snakemake parameters.
annotation_file = snakemake.input.annotation
mat_wt_file = str(snakemake.input.cool_wt)
mat_rf_file = str(snakemake.input.cool_rf)
rna_wt_file = snakemake.input.rna_wt
rna_rf_file = snakemake.input.rna_rf
cov_wt_file = snakemake.input.cov_wt
cov_rf_file = snakemake.input.cov_rf
gc_file = snakemake.input.gc
epod_file = snakemake.input.EPODs
res = int(snakemake.params.res)
cmap = snakemake.params.cmap
outdir = str(snakemake.params.outdir)
width = snakemake.params.width
pos = snakemake.params.positions

# Create outdir if necessary.
os.makedirs(outdir, exist_ok=True)


def import_annotation_gff(annotation_file):
    """Function to create a table of the gene positions from the gff file."""
    annotation = pd.DataFrame(
        columns=["type", "start", "end", "strand", "name"]
    )
    with open(annotation_file, "r") as file:
        for line in file:
            # Header.
            if line.startswith("#"):
                continue
            # Stop at the fasta sequences.
            elif line.startswith(">"):
                break
            else:
                line = line.split("\t")
                if line[2] in ["gene", "tRNA", "rRNA"]:
                    if line[2] == "gene":
                        name = line[8].split("Name=")[-1].split(";")[0]
                        # Extract gene position.
                        annot = {
                            "type": line[2],
                            "start": int(line[3]),
                            "end": int(line[4]),
                            "strand": line[6],
                            "gene_name": name,
                        }
                        annotation = annotation.append(annot, ignore_index=True)
    return annotation


# Import files.
annotation = import_annotation_gff(annotation_file)
mat_wt = cooler.Cooler(f"{mat_wt_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_wt[np.isnan(mat_wt)] = 0
mat_rf = cooler.Cooler(f"{mat_rf_file}::/resolutions/{res}").matrix(
    balance=True, sparse=False
)[:]
mat_rf[np.isnan(mat_rf)] = 0
rna_wt, _ = bcio.extract_big_wig(rna_wt_file, binning=100)
rna_rf, _ = bcio.extract_big_wig(rna_rf_file, binning=100)
cov_wt, _ = bcio.extract_big_wig(cov_wt_file, binning=1000)
cov_rf, _ = bcio.extract_big_wig(cov_wt_file, binning=1000)
gc_content = pd.read_csv(gc_file, sep="\t", header=None).iloc[:, 1]
epod = pd.read_csv(epod_file, sep="\t", header=None)

# Make the rotation of the matrix to plot them
mat_wt_rot = copy.copy(mat_wt)
mat_wt_rot = bch.interpolate_white_lines(mat_wt_rot)
mat_wt_rot[np.isnan(mat_wt_rot)] = 0
mat_wt_rot = ndimage.rotate(mat_wt_rot, 45, reshape=True)
mat_rf_rot = copy.copy(mat_rf)
mat_rf_rot = bch.interpolate_white_lines(mat_rf_rot)
mat_rf_rot[np.isnan(mat_rf_rot)] = 0
mat_rf_rot = ndimage.rotate(mat_rf_rot, 45, reshape=True)


def plot_region(
    M1,
    M1_rot,
    M2,
    M2_rot,
    rna1,
    rna2,
    cov1,
    cov2,
    annotation,
    binning,
    width,
    zoom_ini,
    outfile,
    split=False,
):
    rna_binning = 100
    pal = sns.color_palette("Paired")

    # Defined values to specify the borders of the matrices and tables
    zoom = [zoom // binning for zoom in zoom_ini]
    width = width / 1000
    row1 = len(M1_rot) // 2 - int(np.sqrt(2) * width)
    row2 = len(M1_rot) // 2 + int(np.sqrt(2) * width)
    col1 = int(zoom[0] * np.sqrt(2))
    col2 = int(zoom[1] * np.sqrt(2))
    annotation_zoom = annotation[
        np.logical_and(
            annotation.start > zoom_ini[0], annotation.end < zoom_ini[1]
        )
    ]
    ymax = np.nanmax(
        np.concatenate(
            (
                rna1[zoom_ini[0] // rna_binning : zoom_ini[1] // rna_binning],
                rna2[zoom_ini[0] // rna_binning : zoom_ini[1] // rna_binning],
            )
        )
    )
    max_cov = np.nanmax(cov1) * 1.02
    min_cov = np.nanmin(cov1) * 0.98

    # Define panels depending on split or not.
    if split and ymax > 1500:
        a = 4
        fig, ax = plt.subplots(
            7,
            2,
            figsize=(16, 12),
            gridspec_kw={"height_ratios": [1, 30, 1, 6, 12, 3, 3]},
            sharex=False,
        )

        # Parameters to merge panel 2 and 3.
        d = 0.02
        D = 0.04
        for j in range(2):
            ax[a - 1, j].set_xlim(zoom_ini[0] // 1000, zoom_ini[1] // 1000)
            ax[a - 1, j].tick_params(axis="both", labelsize=15)
            ax[a - 1, j].spines["bottom"].set_visible(False)
            ax[a - 1, j].spines["top"].set_visible(False)
            ax[a - 1, j].spines["right"].set_visible(False)
            ax[a - 1, j].set_ylim(ymax - 525, ymax)
            ax[a - 1, j].get_xaxis().set_visible(False)

            # Add the small cut between the two panels
            kwargs = dict(
                transform=ax[a - 1, j].transAxes, color="k", clip_on=False
            )
            ax[a - 1, j].plot((-d, +d), (-D, +D), **kwargs)  # top-left diagonal
            kwargs.update(
                transform=ax[a, j].transAxes
            )  # switch to the bottom axes
            ax[a, j].plot(
                (-d, +d), (1 - d, 1 + d), **kwargs
            )  # bottom-left diagonal

        # Plot the RNA on the split panel
        ax[a - 1, 0].fill_between(
            np.arange(
                zoom_ini[0] // 1000, zoom_ini[1] // 1000, rna_binning / 1000
            ),
            rna1[zoom_ini[0] // rna_binning : zoom_ini[1] // rna_binning],
            color="black",
        )
        ax[a - 1, 1].fill_between(
            np.arange(
                zoom_ini[0] // 1000, zoom_ini[1] // 1000, rna_binning / 1000
            ),
            rna2[zoom_ini[0] // rna_binning : zoom_ini[1] // rna_binning],
            color="black",
        )

    else:
        a = 3
        fig, ax = plt.subplots(
            6,
            2,
            figsize=(16, 12),
            gridspec_kw={"height_ratios": [1, 30, 1, 15, 3, 3]},
            sharex=False,
        )

    # Plot expressed genes - blue are forward - red are reversed
    for j in range(2):
        ax[0, j].set_xlim(zoom_ini[0], zoom_ini[1])
        ax[0, j].get_xaxis().set_visible(False)
        ax[0, j].get_yaxis().set_visible(False)
        ax[2, j].set_xlim(zoom_ini[0], zoom_ini[1])
        ax[2, j].get_xaxis().set_visible(False)
        ax[2, j].get_yaxis().set_visible(False)
    pos = 1
    for i in range(len(annotation_zoom)):
        # Extract annotation information
        annot = annotation_zoom.iloc[
            i,
        ]
        strand = annot.strand
        start = annot.start // rna_binning
        end = annot.end // rna_binning
        name = annot.gene_name
        # Defined color depending on the strand
        if strand == "+":
            color = pal.as_hex()[1]
        else:
            color = pal.as_hex()[5]
        # Print it only if it transcribed (10% most transcribed genes)
        if np.mean(rna1[start:end]) >= 120.7628998139441:
            pos = pos * -1
            for j in range(2):
                ax[0, j].add_patch(
                    patches.Rectangle(
                        (start * 100, 0),
                        (end - start) * 100,
                        1,
                        edgecolor=color,
                        facecolor=color,
                        fill=True,
                    )
                )
                ax[0, j].text(
                    x=(start + ((end - start) / 2)) * 100,
                    y=pos + 0.5,
                    s=name,
                    rotation=90 * pos,
                    wrap=True,
                )
    color = pal.as_hex()[3]
    for i in epod.index:
        start = epod.loc[i, 3]
        end = epod.loc[i, 4]
        if (
            (start > zoom_ini[0])
            and (start < zoom_ini[1])
            or (end > zoom_ini[0])
            and (end < zoom_ini[1])
        ):
            start = max(start, zoom_ini[0]) // rna_binning
            end = min(end, zoom_ini[1]) // rna_binning
            for j in range(2):
                ax[2, j].add_patch(
                    patches.Rectangle(
                        (start * 100, 0),
                        (end - start) * 100,
                        1,
                        edgecolor=color,
                        facecolor=color,
                        fill=True,
                    )
                )
    #                 ax[2, j].text(
    #                     x=(start + ((end - start) / 2)) * 100,
    #                     y=pos + 0.5,
    #                     s="EPOD",
    #                     rotation=90 * pos,
    #                     wrap=True,
    #                 )

    # Plot the matrices
    ax[1, 0].set_ylabel("Genomic distance (kb)", fontsize=15)
    for j in range(2):
        ax[1, j].get_xaxis().set_visible(False)
        ax[1, j].tick_params(axis="both", labelsize=15)

    # Plot the matrice 1
    im = ax[1, 0].imshow(
        M1_rot[row1:row2, col1:col2],
        cmap=cmap,
        interpolation="none",
        vmin=0,
        vmax=np.nanpercentile(M1[zoom[0] : zoom[1], zoom[0] : zoom[1]], 97),
        extent=(col1 / np.sqrt(2), col2 / np.sqrt(2), width, -width),
    )

    # Plot the matrice 2
    ax[1, 1].imshow(
        M2_rot[row1:row2, col1:col2],
        cmap=cmap,
        interpolation="none",
        vmin=0,
        vmax=np.nanpercentile(M1[zoom[0] : zoom[1], zoom[0] : zoom[1]], 97),
        extent=(col1 / np.sqrt(2), col2 / np.sqrt(2), width, -width),
    )

    # RNAseq legends and plot at the bottom panel.
    ax[a, 0].set_ylabel("RNA count (CPM)", fontsize=15)
    for j in range(2):
        #         ax[a, j].set_xlabel("Genomic coordinates (kb)", size=15)
        ax[a, j].get_xaxis().set_visible(False)
        ax[a, j].tick_params(axis="both", labelsize=15)
        ax[a, j].set_xlim(zoom_ini[0] // 1000, zoom_ini[1] // 1000)
        ax[a, j].spines["top"].set_visible(False)
        ax[a, j].spines["right"].set_visible(False)
        if split:
            if ymax > 1500:
                ax[a, j].set_ylim(0, 1050)
            else:
                ax[a, j].set_ylim(0, 1500)
        else:
            ax[a, j].set_ylim(0, ymax)

    # Plot RNA at the bottom panel
    ax[a, 0].fill_between(
        np.arange(zoom_ini[0] // 1000, zoom_ini[1] // 1000, rna_binning / 1000),
        rna1[zoom_ini[0] // rna_binning : zoom_ini[1] // rna_binning],
        color="black",
    )
    ax[a, 1].fill_between(
        np.arange(zoom_ini[0] // 1000, zoom_ini[1] // 1000, rna_binning / 1000),
        rna2[zoom_ini[0] // rna_binning : zoom_ini[1] // rna_binning],
        color="black",
    )

    # GC content
    gc_binning = 100
    for j in range(2):
        ax[a + 1, j].plot(
            np.arange(
                zoom_ini[0] // 1000, zoom_ini[1] // 1000, gc_binning / 1000
            ),
            gc_content[
                (zoom_ini[0] + 250)
                // gc_binning : (zoom_ini[1] + 250)
                // gc_binning
            ],
            color="black",
        )
        ax[a + 1, j].set_xlim(zoom_ini[0] // 1000, zoom_ini[1] // 1000)
        ax[a + 1, j].tick_params(axis="both", labelsize=15)
        ax[a + 1, j].get_xaxis().set_visible(False)
    ax[a + 1, 0].set_ylabel("GC", fontsize=15)

    # Coverage
    cov_binning = 1000
    ax[a + 2, 0].plot(
        np.arange(zoom_ini[0] // 1000, zoom_ini[1] // 1000, cov_binning / 1000),
        cov1[
            (zoom_ini[0] + 250)
            // cov_binning : (zoom_ini[1] + 250)
            // cov_binning
        ],
        color="black",
    )
    ax[a + 2, 1].plot(
        np.arange(zoom_ini[0] // 1000, zoom_ini[1] // 1000, cov_binning / 1000),
        cov2[
            (zoom_ini[0] + 250)
            // cov_binning : (zoom_ini[1] + 250)
            // cov_binning
        ],
        color="black",
    )
    for j in range(2):
        ax[a + 2, j].set_xlim(zoom_ini[0] // 1000, zoom_ini[1] // 1000)
        ax[a + 2, j].set_ylim(min_cov, max_cov)
        ax[a + 2, j].set_xlabel("Genomic coordinates (kb)", size=15)
        ax[a + 2, j].tick_params(axis="both", labelsize=15)
    ax[a + 2, 0].set_ylabel("HiC\ncoverage\n(CPM)", fontsize=15)

    # Colorbar
    cbar = fig.colorbar(
        im, ax=ax.ravel().tolist(), shrink=0.3, anchor=(1.2, 0.75)
    )
    cbar.ax.tick_params(labelsize=15)

    # Save the fig and adjust margins
    if split and ymax > 1500:
        plt.subplots_adjust(wspace=0.2, hspace=0.1)
    else:
        plt.subplots_adjust(wspace=0.2, hspace=0.1)
    plt.savefig(outfile, bbox_inches="tight", dpi=200)


for split in [True, False]:
    if split:
        outfile = join(
            outdir,
            f"region_{pos[0]}_{pos[1]}_split.pdf",
        )
    else:
        outfile = join(
            outdir,
            f"region_{pos[0]}_{pos[1]}.pdf",
        )
    position = [int(p) * 1000 for p in pos]
    plot_region(
        mat_wt,
        mat_wt_rot,
        mat_rf,
        mat_rf_rot,
        rna_wt,
        rna_rf,
        cov_wt,
        cov_rf,
        annotation,
        res,
        width,
        position,
        outfile,
        split=split,
    )
