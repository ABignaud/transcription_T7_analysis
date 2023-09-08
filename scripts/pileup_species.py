#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.io as bcio
import bacchus.hic as bch
import bacchus.plot as bcp
import bacchus.transcription as bct
import cooler
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import scipy.stats as sst
from os.path import dirname, join

# Make sure that the figures directory exists.
os.makedirs(dirname(snakemake.output.pileup), exist_ok=True)
os.makedirs(dirname(snakemake.output.pileup_zoom), exist_ok=True)

# Import parameters
binning = int(snakemake.params.binning)
label = str(snakemake.params.species)
circular = str(snakemake.params.circular)
window = int(snakemake.params.window)
threshold = int(snakemake.params.threshold)
unit_length = int(snakemake.params.unit_length)

# Import RNA tracks.
rna, chrom_start = bcio.extract_big_wig(
    file=str(snakemake.input.rna),
    binning=binning,
    circular=circular,
    sigma=None,
    ztransform=None,
)
# Ugly loop to have the chrom length...
chrom_start_size = {}
for i, name in enumerate(chrom_start):
    if i != 0:
        chrom_start_size[prev_name] = {
            "start": start,
            "length": chrom_start[name] - start,
        }
    prev_name = name
    start = chrom_start[name]
chrom_start_size[prev_name] = {
    "start": start,
    "length": len(rna) - start,
}

# RNA at 1bp.
rna_1bp, chrom_start_1bp = bcio.extract_big_wig(
    file=str(snakemake.input.rna),
    binning=1,
    circular=circular,
    sigma=None,
    ztransform=None,
)

# Extract annotation from GFF.
annotation = pd.DataFrame(
    columns=["type", "chr", "start", "end", "strand", "name", "tss", "tts", "rpkm"]
)
with open(snakemake.input.annotation, "r") as file:
    for line in file:
        if line.startswith("#"):
            continue
        # Stop if fasta sequence at the end or empty line.
        elif line.startswith(">") or (line.startswith("\n")):
            break
        else:
            line = line.split("\t")
            if line[2] == "CDS":
                name = line[8].split("Name=")[-1].split(";")[0]
                annot = {
                    "type": line[2],
                    "chr": line[0],
                    "start": int(line[3]),
                    "end": int(line[4]),
                    "strand": line[6],
                    "name": name,
                }
                annotation = annotation.append(annot, ignore_index=True)

print(annotation)

# Find TSS/TTS positions.
for i in annotation.index:
    if annotation.loc[i, "strand"] == "+":
        annotation.loc[i, "tss"] = annotation.loc[i, "start"]
        annotation.loc[i, "tts"] = annotation.loc[i, "end"]
    else:
        annotation.loc[i, "tts"] = annotation.loc[i, "start"]
        annotation.loc[i, "tss"] = annotation.loc[i, "end"]

# Compute the RPKM values from the genes.
for i in annotation.index:
    chr_start = chrom_start_1bp[annotation.loc[i, "chr"]]
    annotation.loc[i, "rpkm"] = np.nanmean(
        rna_1bp[
            chr_start
            + annotation.loc[i, "start"] : chr_start
            + annotation.loc[i, "end"]
        ]
    )

# Compute coding density.
n = 0
for i in annotation.index:
    n += annotation.loc[i, "end"] - annotation.loc[i, "start"]
coding_density = n / len(rna_1bp)

# Plot RPKM distribution.
plt.hist(annotation.rpkm, bins=25)
plt.text(
    x=np.nanpercentile(annotation.rpkm, 80),
    y=len(annotation.rpkm) // 25,
    s=f"corr={coding_density:.2f}",
)
plt.savefig(snakemake.output.rpkm)

# Compute pileup
for tu_length in [0, 3000]:
    for threshold2 in np.arange(5, 30, 5):
        rna_pileup_pos, _, pileup_pos, _ = bct.pileup_genes(
            clr=cooler.Cooler(f"{snakemake.input.hic}::/resolutions/{binning}"),
            annotation=annotation,
            rna=rna,
            chrom_start_size=chrom_start_size,
            window_size=25000,
            binning=binning,
            threshold=threshold2,
            neg="detrend",
            tu_length=tu_length,
            operation="mean",
            circular=circular,
        )

        # Plot pileup
        # Parameters
        ax_kb = 1000
        window_plot = 25000 // ax_kb

        # Plot
        fig, ax = plt.subplots(
            2, 1, figsize=(8, 13), gridspec_kw={"height_ratios": [7, 3]}
        )
        # RNA
        ax[1].axvline(
            0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4
        )
        ax[1].tick_params(axis="both", labelsize=14)
        ax[1].plot(
            np.arange(
                -window_plot,
                window_plot + (1 * binning / ax_kb),
                binning / ax_kb,
            ),
            rna_pileup_pos,
        )
        ax[1].set_ylabel("Transcription (CPM)", fontsize=15)
        ax[1].set_xlabel("Genomic distance (kb)", fontsize=15)
        # HiC
        if tu_length == 0:
            vmax = 0.012
        else:
            vmax = 0.008
        ax[0].get_xaxis().set_visible(False)
        im = ax[0].imshow(
            pileup_pos**0.8,
            cmap="Reds",
            vmin=0,
            vmax=vmax,
            extent=[-window_plot, window_plot, window_plot, -window_plot],
        )
        # ax[0].axvline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
        # ax[0].axhline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
        ax[0].tick_params(axis="both", labelsize=14)
        ax[0].set_ylabel("Genomic distance (kb)", fontsize=15)

        # Add colorbar, title and savefig
        fig.colorbar(
            im, ax=ax.ravel().tolist(), shrink=0.33, anchor=(1.3, 0.75)
        )
        plt.subplots_adjust(hspace=0.1)
        ax[0].set_title(label, fontsize=20)
        plt.savefig(
            join(
                snakemake.params.out_dir,
                f"pileup_pos_{threshold2}_TU{tu_length}.pdf",
            ),
            dpi=200,
        )

        hic_signal = bch.compute_hic_signal(
            pileup_pos, binning=500, start=0, stop=5000
        )
        print(
            f"{label} (Gene) - Pearson correlation: {sst.pearsonr(hic_signal[10:40], rna_pileup_pos[10:40])[0]:.2f}"
        )
        print(
            f"{label} (Gene) - p-value: {sst.pearsonr(hic_signal[10:40], rna_pileup_pos[10:40])[1]:.5f}"
        )

# Log ratio pileup
rna_pileup_pos, rna_pileup_neg, pileup_pos, pileup_neg = bct.pileup_genes(
    clr=cooler.Cooler(f"{snakemake.input.hic}::/resolutions/{binning}"),
    annotation=annotation,
    rna=rna,
    chrom_start_size=chrom_start_size,
    window_size=window,
    binning=binning,
    threshold=threshold,
    neg="random-neighbor",
    tu_length=unit_length,
    operation="mean",
    circular=circular,
)

# Print correlation:
hic_signal = bch.compute_hic_signal(pileup_pos, binning=500, start=0, stop=5000)
print(
    f"{label} (TU) - Pearson correlation: {sst.pearsonr(hic_signal[10:window//500 - 10], rna_pileup_pos[10:window//500 - 10])[0]:.2f}"
)
print(
    f"{label} (TU) - p-value: {sst.pearsonr(hic_signal[10:window//500 - 10], rna_pileup_pos[10:window//500 - 10])[1]:.5f}"
)

# Plot pileup pos.
# Plot pileup
# Parameters
ax_kb = 1000
window_plot = 25000 // ax_kb

# Plot pos
fig, ax = plt.subplots(
    2, 1, figsize=(8, 13), gridspec_kw={"height_ratios": [7, 3]}
)
# RNA
ax[1].axvline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
ax[1].tick_params(axis="both", labelsize=14)
ax[1].plot(
    np.arange(
        -window_plot,
        window_plot + (1 * binning / ax_kb),
        binning / ax_kb,
    ),
    rna_pileup_pos[window // 500 - 50 : window // 500 + 51],
)
ax[1].set_ylabel("Transcription (CPM)", fontsize=15)
ax[1].set_xlabel("Genomic distance (kb)", fontsize=15)
# HiC
ax[0].get_xaxis().set_visible(False)
im = ax[0].imshow(
    pileup_pos[
        window // 500 - 50 : window // 500 + 51,
        window // 500 - 50 : window // 500 + 51,
    ]
    ** 0.8,
    cmap="Reds",
    vmin=0.003,
    vmax=0.007,
    extent=[-window_plot, window_plot, window_plot, -window_plot],
)
# ax[0].axvline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
# ax[0].axhline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
ax[0].tick_params(axis="both", labelsize=14)
ax[0].set_ylabel("Genomic distance (kb)", fontsize=15)

# Add colorbar, title and savefig
fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.33, anchor=(1.3, 0.75))
plt.subplots_adjust(hspace=0.1)
ax[0].set_title(label, fontsize=20)
plt.savefig(snakemake.output.pileup_pos_TU, dpi=200)

# Plot neg
fig, ax = plt.subplots(
    2, 1, figsize=(8, 13), gridspec_kw={"height_ratios": [7, 3]}
)
# RNA
ax[1].axvline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
ax[1].tick_params(axis="both", labelsize=14)
ax[1].plot(
    np.arange(
        -window_plot,
        window_plot + (1 * binning / ax_kb),
        binning / ax_kb,
    ),
    rna_pileup_neg[window // 500 - 50 : window // 500 + 51],
)
ax[1].set_ylabel("Transcription (CPM)", fontsize=15)
ax[1].set_xlabel("Genomic distance (kb)", fontsize=15)
# HiC
ax[0].get_xaxis().set_visible(False)
im = ax[0].imshow(
    pileup_neg[
        window // 500 - 50 : window // 500 + 51,
        window // 500 - 50 : window // 500 + 51,
    ]
    ** 0.8,
    cmap="Reds",
    vmin=0,
    vmax=np.nanpercentile(
        pileup_pos[
            window // 500 - 50 : window // 500 + 51,
            window // 500 - 50 : window // 500 + 51,
        ]
        ** 0.8,
        99,
    ),
    extent=[-window_plot, window_plot, window_plot, -window_plot],
)
# ax[0].axvline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
# ax[0].axhline(0, color="black", linestyle="dashed", linewidth=1.5, alpha=0.4)
ax[0].tick_params(axis="both", labelsize=14)
ax[0].set_ylabel("Genomic distance (kb)", fontsize=15)

# Add colorbar, title and savefig
fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.33, anchor=(1.3, 0.75))
plt.subplots_adjust(hspace=0.1)
ax[0].set_title(label, fontsize=20)
plt.savefig(snakemake.output.pileup_neg_TU, dpi=200)

# Plot pileup genes

# Plot the pileup ratio.
bcp.pileup_plot(
    pileup=pileup_pos,
    pileup_control=pileup_neg,
    gen_tracks=[rna_pileup_pos],
    gen_tracks_control=[rna_pileup_neg],
    binning=binning,
    window=50000,
    ratio="log",
    out_file=snakemake.output.pileup,
    title=label,
    dpi=200,
    vmax=0.2,
)
bcp.pileup_plot(
    pileup=pileup_pos,
    pileup_control=pileup_neg,
    gen_tracks=[rna_pileup_pos],
    gen_tracks_control=[rna_pileup_neg],
    binning=binning,
    window=25000,
    ratio="log",
    out_file=snakemake.output.pileup_zoom,
    title=label,
    dpi=200,
    vmax=0.2,
)
bcp.pileup_plot(
    pileup=pileup_pos,
    pileup_control=pileup_neg,
    gen_tracks=[rna_pileup_pos],
    gen_tracks_control=[rna_pileup_neg],
    binning=binning,
    window=10000,
    ratio="log",
    out_file=snakemake.output.pileup_zoom2,
    title=label,
    dpi=200,
    vmax=0.2,
)
