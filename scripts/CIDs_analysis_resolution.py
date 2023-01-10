#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bacchus.directional as bcd
import bacchus.insulation as bci
import bacchus.io as bcio
import cooler
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as mc
import numpy as np
import os
import pandas as pd
import scipy.sparse as sp
import scipy.stats as st
import seaborn as sns
from random import choices

mat_file = snakemake.input.mat
rna_file = snakemake.input.rna
annotation_file = snakemake.input.annotation
cmap = snakemake.params.cmap
dpi = snakemake.params.dpi
text_file = str(snakemake.output.text_file)
full_mat_file = str(snakemake.output.full_mat)
zoom_1_file = str(snakemake.output.zoom1)
zoom_2_file = str(snakemake.output.zoom2)
full_mat_di_file = str(snakemake.output.full_mat_di)
zoom_1_di_file = str(snakemake.output.zoom1_di)
zoom_2_di_file = str(snakemake.output.zoom2_di)
out_files_mat = [full_mat_file, zoom_1_file, zoom_2_file]
out_files_mat_di = [full_mat_di_file, zoom_1_di_file, zoom_2_di_file]
bor_trans = str(snakemake.output.bor_trans)
trans_bor = str(snakemake.output.trans_bor)

# Create outdir if necessary.
os.makedirs(str(snakemake.params.outdir), exist_ok=True)


def write_stats(text_file, text, writing_type="a"):
    with open(text_file, writing_type) as out:
        out.write(f"{text}\n")


# Compute CIDs at 5kb resolution with 100kb window size.
M5000 = cooler.Cooler(f"{mat_file}::/resolutions/{5000}").matrix(
    balance=True, sparse=True
)[:]
di5000 = bcd.directional_index(M5000, 20)
l5000 = bcd.di_borders(di5000)
write_stats(
    text_file,
    f"Numbers of CIDs at 5kb with 100kb window size: {len(l5000)}",
    "w",
)

# Compute CIDs at 2kb resolution with 50kb window size.
M2000 = cooler.Cooler(f"{mat_file}::/resolutions/{2000}").matrix(
    balance=True, sparse=True
)[:]
di2000 = bcd.directional_index(M2000, 25)
l2000 = bcd.di_borders(di2000)
write_stats(
    text_file, f"Numbers of CIDs at 2kb with 50kb window size: {len(l2000)}"
)

# Compute CIDs at 1kb resolution using insulation score.
M1000 = cooler.Cooler(f"{mat_file}::/resolutions/{1000}").matrix(
    balance=True, sparse=True
)[:]
di1000 = bcd.directional_index(M1000, 20)
l1000 = bcd.di_borders(di1000)
write_stats(
    text_file, f"Numbers of CIDs at 1kb with 50kb window size: {len(l1000)}"
)

mat = cooler.Cooler(f"{mat_file}::/resolutions/{1000}").matrix(
    balance=True, sparse=False
)[:]
mat[np.isnan(mat)] = 0
final_borders_1000, lri_1000 = bci.get_insulation_score(
    mat, [10, 15, 20, 25, 30]
)
write_stats(
    text_file,
    f"{len(final_borders_1000)} CIDs have been detected with insulation score at 1kb resolution.",
)


rna, _ = bcio.extract_big_wig(rna_file, ztransform=False)

annotation = pd.DataFrame(
    columns=["type", "start", "end", "strand", "gene_name", "rpkm"]
)

n = 0
with open(annotation_file, "r") as file:
    for line in file:
        if line.startswith("#"):
            continue
        elif line.startswith(">"):
            break
        else:
            line = line.split("\t")
            if line[2] in ["gene", "tRNA", "rRNA"]:
                if line[2] == "gene":
                    name = line[8].split("Name=")[-1].split(";")[0]
                    annot = {
                        "type": line[2],
                        "start": int(line[3]),
                        "end": int(line[4]),
                        "strand": line[6],
                        "gene_name": name,
                        "rpkm": np.mean(rna[int(line[3]) : int(line[4])]),
                    }
                    annotation = annotation.append(annot, ignore_index=True)


###
# Print matrices with CIDs borders (Fig 1e).
###
starts = [0, 2500, 3000]
ends = [len(mat), 3000, 3500]
col = ["blue", "lime", "cyan"]

for i in range(3):
    out_file, start, end = out_files_mat[i], starts[i], ends[i]

    vmax = 99
    title = None
    fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=dpi)

    # Axis values
    scaling_factor = 1.0
    axis = "kb"

    # No end values given.
    if end == 0:
        end = len(mat)

    # Display plots
    im = ax.imshow(
        mat[start:end, start:end],
        cmap=cmap,
        vmin=0,
        vmax=np.nanpercentile(mat, vmax),
        extent=(
            start * scaling_factor,
            end * scaling_factor,
            end * scaling_factor,
            start * scaling_factor,
        ),
    )

    # Legend
    ax.set_xlabel(f"Genomic coordinates ({axis:s})", fontsize=16)
    ax.set_ylabel(f"Genomic coordinates ({axis:s})", fontsize=16)
    ax.tick_params(axis="both", labelsize=16)

    # Title
    if title is not None:
        ax.set_title(title, size=18)

    # Colorbar
    cbar = plt.colorbar(im, shrink=0.33, anchor=(0, 0.5))
    cbar.ax.tick_params(labelsize=16)

    first = True
    for i in [x * 5 for x in l5000]:
        if i in np.arange(start, end):
            #         plt.axvline(i, lw=1, ls='dashed', c='blue')
            if first:
                s = i
                first = False
                ax.axvline(
                    s,
                    ymin=(end - s) / (end - start),
                    ymax=1,
                    c=col[0],
                    lw=1,
                    zorder=0.5,
                )
                ax.axhline(
                    s,
                    xmin=0,
                    xmax=1 - ((end - s) / (end - start)),
                    c=col[0],
                    lw=1,
                    zorder=0.5,
                )
            else:
                e = i
                ax.add_patch(
                    patches.Rectangle(
                        (s, s), e - s, e - s, edgecolor=col[0], fill=False
                    )
                )
                s = i
    ax.axvline(
        s, ymin=(end - s) / (end - start), ymax=0, c=col[0], lw=1, zorder=0.5
    )
    ax.axhline(
        s,
        xmin=1,
        xmax=1 - ((end - s) / (end - start)),
        c=col[0],
        lw=1,
        zorder=0.5,
    )

    first = True
    for i in [x * 5 for x in l2000]:
        if i in np.arange(start, end):
            #         plt.axvline(i, lw=1, ls='dashed', c='cyan')
            if first:
                s = i
                first = False
                ax.axvline(
                    s,
                    ymin=(end - s) / (end - start),
                    ymax=1,
                    c=col[1],
                    lw=1,
                    zorder=0.5,
                )
                ax.axhline(
                    s,
                    xmin=0,
                    xmax=1 - ((end - s) / (end - start)),
                    c=col[1],
                    lw=1,
                    zorder=0.5,
                )
            else:
                e = i
                ax.add_patch(
                    patches.Rectangle(
                        (s, s), e - s, e - s, edgecolor=col[1], fill=False
                    )
                )
                s = i
    ax.axvline(
        s, ymin=(end - s) / (end - start), ymax=0, c=col[1], lw=1, zorder=0.5
    )
    ax.axhline(
        s,
        xmin=1,
        xmax=1 - ((end - s) / (end - start)),
        c=col[1],
        lw=1,
        zorder=0.5,
    )

    first = True
    for i in final_borders_1000:
        if i in np.arange(start, end):
            #         plt.axvline(i, lw=1, ls='dashed', c='magenta')
            if first:
                s = i
                first = False
                ax.axvline(
                    s,
                    ymin=(end - s) / (end - start),
                    ymax=1,
                    c=col[2],
                    lw=1,
                    zorder=0.5,
                )
                ax.axhline(
                    s,
                    xmin=0,
                    xmax=1 - ((end - s) / (end - start)),
                    c=col[2],
                    lw=1,
                    zorder=0.5,
                )
            else:
                e = i
                ax.add_patch(
                    patches.Rectangle(
                        (s, s), e - s, e - s, edgecolor=col[2], fill=False
                    )
                )
                s = i
    ax.axvline(
        s, ymin=(end - s) / (end - start), ymax=0, c=col[2], lw=1, zorder=0.5
    )
    ax.axhline(
        s,
        xmin=1,
        xmax=1 - ((end - s) / (end - start)),
        c=col[2],
        lw=1,
        zorder=0.5,
    )

    # Savefig
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
    plt.close()


###
# Print full matrices with CIDs borders using DI only (Supp Fig 1f).
###

for i in range(3):
    out_file, start, end = out_files_mat_di[i], starts[i], ends[i]

    vmax = 99
    title = None
    fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=dpi)

    # Axis values
    scaling_factor = 1.0
    axis = "kb"

    # No end values given.
    if end == 0:
        end = len(mat)

    # Display plots
    im = ax.imshow(
        mat[start:end, start:end],
        cmap=cmap,
        vmin=0,
        vmax=np.percentile(mat, vmax),
        extent=(
            start * scaling_factor,
            end * scaling_factor,
            end * scaling_factor,
            start * scaling_factor,
        ),
    )

    # Legend
    ax.set_xlabel(f"Genomic coordinates ({axis:s})", fontsize=16)
    ax.set_ylabel(f"Genomic coordinates ({axis:s})", fontsize=16)
    ax.tick_params(axis="both", labelsize=16)

    # Title
    if title is not None:
        ax.set_title(title, size=18)

    # Colorbar
    cbar = plt.colorbar(im, shrink=0.33, anchor=(0, 0.5))
    cbar.ax.tick_params(labelsize=16)

    first = True
    for i in [x * 5 for x in l5000]:
        if i in np.arange(start, end):
            #         plt.axvline(i, lw=1, ls='dashed', c='blue')
            if first:
                s = i
                first = False
                ax.axvline(
                    s,
                    ymin=(end - s) / (end - start),
                    ymax=1,
                    c=col[0],
                    lw=1,
                    zorder=0.5,
                )
                ax.axhline(
                    s,
                    xmin=0,
                    xmax=1 - ((end - s) / (end - start)),
                    c=col[0],
                    lw=1,
                    zorder=0.5,
                )
            else:
                e = i
                ax.add_patch(
                    patches.Rectangle(
                        (s, s), e - s, e - s, edgecolor=col[0], fill=False
                    )
                )
                s = i
    ax.axvline(
        s, ymin=(end - s) / (end - start), ymax=0, c=col[0], lw=1, zorder=0.5
    )
    ax.axhline(
        s,
        xmin=1,
        xmax=1 - ((end - s) / (end - start)),
        c=col[0],
        lw=1,
        zorder=0.5,
    )

    first = True
    for i in [x * 5 for x in l2000]:
        if i in np.arange(start, end):
            # plt.axvline(i, lw=1, ls='dashed', c='cyan')
            if first:
                s = i
                first = False
                ax.axvline(
                    s,
                    ymin=(end - s) / (end - start),
                    ymax=1,
                    c=col[1],
                    lw=1,
                    zorder=0.5,
                )
                ax.axhline(
                    s,
                    xmin=0,
                    xmax=1 - ((end - s) / (end - start)),
                    c=col[1],
                    lw=1,
                    zorder=0.5,
                )
            else:
                e = i
                ax.add_patch(
                    patches.Rectangle(
                        (s, s), e - s, e - s, edgecolor=col[1], fill=False
                    )
                )
                s = i
    ax.axvline(
        s, ymin=(end - s) / (end - start), ymax=0, c=col[1], lw=1, zorder=0.5
    )
    ax.axhline(
        s,
        xmin=1,
        xmax=1 - ((end - s) / (end - start)),
        c=col[1],
        lw=1,
        zorder=0.5,
    )

    first = True
    for i in [x * 5 for x in l1000]:
        if i in np.arange(start, end):
            # plt.axvline(i, lw=1, ls='dashed', c='cyan')
            if first:
                s = i
                first = False
                ax.axvline(
                    s,
                    ymin=(end - s) / (end - start),
                    ymax=1,
                    c=col[2],
                    lw=1,
                    zorder=0.5,
                )
                ax.axhline(
                    s,
                    xmin=0,
                    xmax=1 - ((end - s) / (end - start)),
                    c=col[2],
                    lw=1,
                    zorder=0.5,
                )
            else:
                e = i
                ax.add_patch(
                    patches.Rectangle(
                        (s, s), e - s, e - s, edgecolor=col[2], fill=False
                    )
                )
                s = i
    ax.axvline(
        s, ymin=(end - s) / (end - start), ymax=0, c=col[2], lw=1, zorder=0.5
    )
    ax.axhline(
        s,
        xmin=1,
        xmax=1 - ((end - s) / (end - start)),
        c=col[2],
        lw=1,
        zorder=0.5,
    )

    # Savefig
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
    plt.close()

###
# Transcription at borders.
###


def get_borders_transcriptions(annotation, borders, size, step):
    borders_transcriptions = []
    for i in annotation.index:
        if annotation.loc[i, "strand"] == "-":
            pos = annotation.loc[i, "end"]
        else:
            pos = annotation.loc[i, "start"]
        for j in borders:
            a = j * size + step
            if abs(pos - a) < 2500:
                borders_transcriptions.append(annotation.loc[i, "rpkm"])
    return borders_transcriptions


borders_trans_5000 = get_borders_transcriptions(annotation, l5000, 5000, 2500)
borders_trans_2000 = get_borders_transcriptions(annotation, l2000, 5000, 1000)
borders_trans_1000 = get_borders_transcriptions(
    annotation, final_borders_1000, 1000, 500
)
data = {
    "RPKM (log)": np.log(
        list(annotation.rpkm)
        + list(np.log(borders_trans_5000))
        + list(np.log(borders_trans_2000))
        + list(np.log(borders_trans_1000))
    ),
    "Genes": list(np.repeat(f"All genes\nn={len(annotation)}", len(annotation)))
    + list(
        np.repeat(f"5kb\nn={len(borders_trans_5000)}", len(borders_trans_5000))
    )
    + list(
        np.repeat(f"2kb\nn={len(borders_trans_2000)}", len(borders_trans_2000))
    )
    + list(
        np.repeat(f"1kb\nn={len(borders_trans_1000)}", len(borders_trans_1000))
    ),
}
data = pd.DataFrame(data)
data.replace([np.inf, -np.inf, np.nan], 0, inplace=True)

sns.violinplot(x="Genes", y="RPKM (log)", data=data, palette="tab10")
plt.xlabel("Resolution of the matrix to detect the borders", size=14)
plt.ylabel("RPKM (log)", size=14)
plt.title(
    "Genes transcription of genes at less\nthan 5kb of a detected border",
    size=16,
)
plt.xticks(size=12)
plt.yticks(size=12)
plt.savefig(bor_trans, bbox_inches="tight")
plt.close()

write_stats(
    text_file,
    f"p-value between 5kb borders and whole genome: {st.mannwhitneyu(np.log(list(annotation.rpkm)), np.log(borders_trans_5000))[1]}",
)
write_stats(
    text_file,
    f"p-value between 2kb borders and whole genome: {st.mannwhitneyu(np.log(list(annotation.rpkm)), np.log(borders_trans_2000))[1]}",
)
write_stats(
    text_file,
    f"p-value between 1kb borders and whole genome: {st.mannwhitneyu(np.log(list(annotation.rpkm)), np.log(borders_trans_1000))[1]}",
)

##
# Transcription depending on distance borders
###

annotation["dist_5kb"] = 0
annotation["dist_2kb"] = 0
annotation["dist_1kb"] = 0
annotation["log1p_rpkm"] = np.log1p(annotation["rpkm"])
for i in annotation.index:
    value = np.inf
    if annotation.loc[i, "strand"] == "+":
        for x in l5000:
            pos = annotation.loc[i, "start"]
            if pos < x * 5000:
                a = x * 5000 - pos
            elif pos > (x + 1) * 5000:
                a = pos - (x + 1) * 5000
            else:
                a = -1
            value = min(a, value)
    else:
        for x in l5000:
            pos = annotation.loc[i, "end"]
            if pos < x * 5000:
                a = x * 5000 - pos
            elif pos > (x + 1) * 5000:
                a = pos - (x + 1) * 5000
            else:
                a = -1
            value = min(a, value)
    annotation.loc[i, "dist_5kb"] = value

for i in annotation.index:
    value = np.inf
    if annotation.loc[i, "strand"] == "+":
        for x in l2000:
            pos = annotation.loc[i, "start"]
            if pos < x * 5000:
                a = x * 5000 - pos
            elif pos > (x + 0.4) * 5000:
                a = pos - (x + 0.4) * 5000
            else:
                a = -1
            value = min(a, value)
    else:
        for x in l2000:
            pos = annotation.loc[i, "end"]
            if pos < x * 5000:
                a = x * 5000 - pos
            elif pos > (x + 0.4) * 5000:
                a = pos - (x + 0.4) * 5000
            else:
                a = -1
            value = min(a, value)
    annotation.loc[i, "dist_2kb"] = value

for i in annotation.index:
    value = np.inf
    if annotation.loc[i, "strand"] == "+":
        for x in final_borders_1000:
            pos = annotation.loc[i, "start"]
            if pos < x * 1000:
                a = x * 1000 - pos
            elif pos > (x + 1) * 1000:
                a = pos - (x + 1) * 1000
            else:
                a = -1
            value = min(a, value)
    else:
        for x in final_borders_1000:
            pos = annotation.loc[i, "end"]
            if pos < x * 1000:
                a = x * 1000 - pos
            elif pos > (x + 1) * 1000:
                a = pos - (x + 1) * 1000
            else:
                a = -1
            value = min(a, value)
    annotation.loc[i, "dist_1kb"] = value

bin_means_l3, bin_edges_l3, binnumber_l3 = st.binned_statistic(
    annotation.dist_1kb,
    annotation.rpkm,
    statistic="mean",
    bins=41,
    range=(-5000, 200000),
)
bin_means_l2, bin_edges_l2, binnumber_l2 = st.binned_statistic(
    annotation.dist_2kb,
    annotation.rpkm,
    statistic="mean",
    bins=41,
    range=(-5000, 200000),
)
bin_means_l1, bin_edges_l1, binnumber_l1 = st.binned_statistic(
    annotation.dist_5kb,
    annotation.rpkm,
    statistic="mean",
    bins=41,
    range=(-5000, 200000),
)


def bootstrap_down(data):
    l = np.zeros((10000))
    for i in range(10000):
        l[i] = np.mean(choices(data, k=len(data)))
    return np.percentile(l, 5)


def bootstrap_up(data):
    l = np.zeros((10000))
    for i in range(10000):
        l[i] = np.mean(choices(data, k=len(data)))
    return np.percentile(l, 95)


bin_perc5_l1, bin_edges_l1, binnumber_l1 = st.binned_statistic(
    annotation.dist_5kb,
    annotation.rpkm,
    statistic=bootstrap_down,
    bins=41,
    range=(-5000, 200000),
)
bin_perc95_l1, bin_edges_l1, binnumber_l1 = st.binned_statistic(
    annotation.dist_5kb,
    annotation.rpkm,
    statistic=bootstrap_up,
    bins=41,
    range=(-5000, 200000),
)
bin_perc5_l2, bin_edges_l2, binnumber_l2 = st.binned_statistic(
    annotation.dist_2kb,
    annotation.rpkm,
    statistic=bootstrap_down,
    bins=41,
    range=(-5000, 200000),
)
bin_perc95_l2, bin_edges_l2, binnumber_l2 = st.binned_statistic(
    annotation.dist_2kb,
    annotation.rpkm,
    statistic=bootstrap_up,
    bins=41,
    range=(-5000, 200000),
)
bin_perc5_l3, bin_edges_l3, binnumber_l3 = st.binned_statistic(
    annotation.dist_1kb,
    annotation.rpkm,
    statistic=bootstrap_down,
    bins=41,
    range=(-5000, 200000),
)
bin_perc95_l3, bin_edges_l3, binnumber_l3 = st.binned_statistic(
    annotation.dist_1kb,
    annotation.rpkm,
    statistic=bootstrap_up,
    bins=41,
    range=(-5000, 200000),
)

plt.plot(bin_edges_l1[1:] / 1000, bin_means_l1, label="5kb")
plt.fill_between(
    bin_edges_l1[1:] / 1000, bin_perc5_l1, bin_perc95_l1, alpha=0.5
)
plt.plot(bin_edges_l2[1:] / 1000, bin_means_l2, label="2kb")
plt.fill_between(
    bin_edges_l2[1:] / 1000, bin_perc5_l2, bin_perc95_l2, alpha=0.5
)
plt.plot(bin_edges_l3[1:] / 1000, bin_means_l3, label="1kb")
plt.fill_between(
    bin_edges_l3[1:] / 1000, bin_perc5_l3, bin_perc95_l3, alpha=0.5
)
plt.xlabel("Distance from the closest border (kb)", size=14)
plt.ylabel("Mean RPKM", size=14)
plt.xlim(0, 50)
plt.title(
    "Genes transcription depending on\nthe distance to the closest border",
    size=16,
)
plt.xticks(size=12)
plt.yticks(size=12)
plt.legend()
plt.savefig(trans_bor, bbox_inches="tight")
plt.close()
