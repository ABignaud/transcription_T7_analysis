#!/bin/env snakemake -s

# Rules to plot matrices and ratio of the contact map of T7 system contact map.

# Figure 2
rule rifampicine_fw:
    input:
        mat_t7_glu = join(OUT_DIR, 'HiC', 'T7_glu.mcool'),
        mat_t7_ara = join(OUT_DIR, 'HiC', 'T7_ara.mcool'),
        mat_t7_glu_rif = join(OUT_DIR, 'HiC', 'T7_glu_rif.mcool'),
        mat_t7_ara_rif = join(OUT_DIR, 'HiC', 'T7_ara_rif.mcool'),
        rna_t7_glu = join(OUT_DIR, 'RNA_tracks', 'T7_glu_rna_unstranded.bw'),
        rna_t7_ara = join(OUT_DIR, 'RNA_tracks', 'T7_ara_rna_unstranded.bw'),
        rna_t7_glu_rif = join(OUT_DIR, 'RNA_tracks', 'T7_glu_rif_rna_unstranded.bw'),
        rna_t7_ara_rif = join(OUT_DIR, 'RNA_tracks', 'T7_ara_rif_rna_unstranded.bw'),
    params:
        cmap = config['cmap'],
        res = 1000,
        res_ratio = 5000,
        width = 64000,
        outdir = join(OUT_DIR, 'figures', 'T7_system', 'T7_fw'),
    output:
        matrices = join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'matrices.pdf'),
        zoom = join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'zoom_matrices.pdf'),
        zoom_ratio = join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'zoom_ratio.pdf'),
        rna_only = join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'rna_only.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: config['threads_serp']
    script: '../scripts/rifampicine_for.py'

# Supplementary Figure 2a
rule rifampicine_rv:
    input:
        mat_t7_glu = join(OUT_DIR, 'HiC', 'T7_rv_glu.mcool'),
        mat_t7_ara = join(OUT_DIR, 'HiC', 'T7_rv_ara.mcool'),
        mat_t7_glu_rif = join(OUT_DIR, 'HiC', 'T7_rv_glu_rif.mcool'),
        mat_t7_ara_rif = join(OUT_DIR, 'HiC', 'T7_rv_ara_rif.mcool'),
        rna_t7_glu = join(OUT_DIR, 'RNA_tracks', 'T7_rv_glu_rna_unstranded.bw'),
        rna_t7_ara = join(OUT_DIR, 'RNA_tracks', 'T7_rv_ara_rna_unstranded.bw'),
        rna_t7_glu_rif = join(OUT_DIR, 'RNA_tracks', 'T7_rv_glu_rif_rna_unstranded.bw'),
        rna_t7_ara_rif = join(OUT_DIR, 'RNA_tracks', 'T7_rv_ara_rif_rna_unstranded.bw'),
        chip_t7_glu = join(OUT_DIR, 'ChIP_tracks', 'CCC01.bw'),
        chip_t7_ara = join(OUT_DIR, 'ChIP_tracks', 'CCC02.bw'),
        chip_t7_glu_rif = join(OUT_DIR, 'ChIP_tracks', 'CCC03.bw'),
        chip_t7_ara_rif_R1 = join(OUT_DIR, 'ChIP_tracks', 'CCC04.bw'),
        chip_t7_ara_rif_R2 = join(OUT_DIR, 'ChIP_tracks', 'CCchip05.bw'),
    params:
        cmap = config['cmap'],
        res = 1000,
        res_ratio = 5000,
        width = 64000,
        outdir = join(OUT_DIR, 'figures', 'T7_system', 'T7_rev'),
    output:
        matrices = join(OUT_DIR, 'figures', 'T7_system', 'T7_rev', 'matrices.pdf'),
        zoom = join(OUT_DIR, 'figures', 'T7_system', 'T7_rev', 'zoom_matrices.pdf'),
        zoom_ratio = join(OUT_DIR, 'figures', 'T7_system', 'T7_rev', 'zoom_ratio.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: config['threads_serp']
    script: '../scripts/rifampicine_rev.py'

rule all_fig2:
    input:
        join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'matrices.pdf'),
        join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'zoom_matrices.pdf'),
        join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'zoom_ratio.pdf'),
        join(OUT_DIR, 'figures', 'T7_system', 'T7_fw', 'rna_only.pdf'),
        join(OUT_DIR, 'figures', 'T7_system', 'T7_rev', 'matrices.pdf'),
        join(OUT_DIR, 'figures', 'T7_system', 'T7_rev', 'zoom_matrices.pdf'),
        join(OUT_DIR, 'figures', 'T7_system', 'T7_rev', 'zoom_ratio.pdf'),
    output:
        touch(join(TMP, 'all_fig2.done'))


