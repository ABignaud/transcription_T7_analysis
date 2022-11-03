#!/bin/env snakemake -s

# Rules to plot matrices and ratio of the contact map of mutants or drug treated
# bacteria.

rule rifampicine:
    input:
        mat_t7_glu = join(OUT_DIR, 'HiC', 'T7_fw_glu.mcool'),
        mat_t7_ara = join(OUT_DIR, 'HiC', 'T7_fw_ara.mcool'),
        mat_t7_glu_rif = join(OUT_DIR, 'HiC', 'T7_fw_glu_rif.mcool'),
        mat_t7_ara_rif = join(OUT_DIR, 'HiC', 'T7_fw_ara_rif.mcool'),
        rna_t7_glu = join(OUT_DIR, 'RNA_tracks', 'T7_fw_glu_unstranded.bw'),
        rna_t7_ara = join(OUT_DIR, 'RNA_tracks', 'T7_fw_ara_unstranded.bw'),
        rna_t7_glu_rif = join(OUT_DIR, 'RNA_tracks', 'T7_fw_glu_rif_unstranded.bw'),
        rna_t7_ara_rif = join(OUT_DIR, 'RNA_tracks', 'T7_fw_ara_rif_unstranded.bw'),
    params:
        cmap = config['cmap'],
        res = 1000,
    conda: '../envs/bacchus.yaml'
    output:
        matrices = join(OUT_DIR, 'figures', 'T7_fw', 'matrices.pdf'),
        zoom = join(OUT_DIR, 'figures', 'T7_fw', 'zoom.pdf'),
        zoom_ratio = join(OUT_DIR, 'figures', 'T7_fw', 'zoom_ratio.pdf'),
    script: '../scripts/rifampicine.py'

rule rifampicine_reverse:

rule chlorophenicol:
rule stop_codon:
rule matP:
rule novobiocin:

