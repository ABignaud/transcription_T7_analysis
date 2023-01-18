#!/bin/env snakemake -s

# Rules to plot matrices and ratio of the contact map of T7 system contact map.

# Figure 3
rule multi_prom_gapr:
    input:
        mat_t7 = join(OUT_DIR, 'HiC', 'T7.mcool'),
        mat_t7_2p100 = join(OUT_DIR, 'HiC', 'T7_2P100.mcool'),
        mat_t7_2p60 = join(OUT_DIR, 'HiC', 'T7_2P60.mcool'),
        mat_t7_2p100_dv = join(OUT_DIR, 'HiC', 'T7_2P100_DV.mcool'),
        mat_t7_2p60_dv = join(OUT_DIR, 'HiC', 'T7_2P60_DV.mcool'),
        mat_t7_2p100_cv = join(OUT_DIR, 'HiC', 'T7_2P100_CV.mcool'),
        rna_t7 = join(OUT_DIR, 'RNA_tracks', 'T7_rna_unstranded.bw'),
        rna_t7_2p100 = join(OUT_DIR, 'RNA_tracks', 'T7_2P100_rna_unstranded.bw'),
        rna_t7_2p60 = join(OUT_DIR, 'RNA_tracks', 'T7_2P60_rna_unstranded.bw'),
        rna_t7_2p100_dv = join(OUT_DIR, 'RNA_tracks', 'T7_2P100_DV_rna_unstranded.bw'),
        rna_t7_2p60_dv = join(OUT_DIR, 'RNA_tracks', 'T7_2P60_DV_rna_unstranded.bw'),
        rna_t7_2p100_cv = join(OUT_DIR, 'RNA_tracks', 'T7_2P100_CV_rna_unstranded.bw'),
        chip_RNA_t7 = join(OUT_DIR, 'ChIP_tracks', 'CCchip08.bw'),
        chip_RNA_t7_2p100 = join(OUT_DIR, 'ChIP_tracks', 'CCchip09.bw'),
        chip_RNA_t7_2p60 = join(OUT_DIR, 'ChIP_tracks', 'CCchip11.bw'),
        chip_RNA_t7_2p100_dv = join(OUT_DIR, 'ChIP_tracks', 'CCchip10.bw'),
        chip_RNA_t7_2p60_dv = join(OUT_DIR, 'ChIP_tracks', 'CCchip12.bw'),
        chip_RNA_t7_2p100_cv = join(OUT_DIR, 'ChIP_tracks', 'CCchip06.bw'),
        chip_gapr_t7 = join(OUT_DIR, 'ChIP_tracks', 'CCchip16.bw'),
        chip_gapr_t7_2p100 = join(OUT_DIR, 'ChIP_tracks', 'CCchip18.bw'),
        chip_gapr_t7_2p60 = join(OUT_DIR, 'ChIP_tracks', 'CCchip20.bw'),
        chip_gapr_t7_2p100_dv = join(OUT_DIR, 'ChIP_tracks', 'CCchip19.bw'),
        chip_gapr_t7_2p60_dv = join(OUT_DIR, 'ChIP_tracks', 'CCchip21.bw'),
        chip_gapr_t7_2p100_cv = join(OUT_DIR, 'ChIP_tracks', 'CCchip17.bw'),
        chip_gapr_control = join(OUT_DIR, 'ChIP_tracks', 'CCchip14.bw'),
    params:
        cmap = config['cmap'],
        res = 1000,
        tracks_res = 100,
        outdir = join(OUT_DIR, 'figures', 'T7_system', 'multi_prom'),
    output:
        plot = join(OUT_DIR, 'figures', 'T7_system', 'multi_prom', 'multi_prom_gapr.pdf'),
        txt = join(OUT_DIR, 'figures', 'T7_system', 'multi_prom', 'multi_prom_corr.txt'),
    conda: '../envs/bacchus.yaml'
    threads: 1
    script: '../scripts/multi_prom_gapr.py'
