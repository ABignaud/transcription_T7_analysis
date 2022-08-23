#!/bin/env snakemake -s

# Generates pileups plot from the output of previous tracks.

# rule corr:
#     input:
#         hic = join(OUT_DIR, '{species}', 'cool', '{species}.mcool'),
#         rna = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_unstranded.bw'),
#     output:
#         mat_fig = join(OUT_DIR, '{species}', 'figures', '{species}_map.pdf'),
#         corr_fig = join(OUT_DIR, '{species}', 'figures', '{species}_corr.pdf'),
#     params:
#         binning = lambda w: config[w.species]['res'][0],
#         species = '{species}',
#         circular = lambda w: config[w.species]['circular'],
#     conda: '../envs/bacchus.yaml'
#     script: '../scripts/corr.py'


rule pileup:
    input:
        hic = join(OUT_DIR, '{species}', 'cool', '{species}.mcool'),
        rna = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_unstranded.bw'),
        annotation = lambda w: join(REF_DIR, config[w.species]['annotation']),
    output:
        pileup_pos_TU = join(OUT_DIR, '{species}', 'figures', '{species}_pileup_pos_TU.pdf'),
        pileup_neg_TU = join(OUT_DIR, '{species}', 'figures', '{species}_pileup_neg_TU.pdf'),
        pileup = join(OUT_DIR, '{species}', 'figures', '{species}_pileup.pdf'),
        pileup_zoom = join(OUT_DIR, '{species}', 'figures', '{species}_pileup_zoom.pdf'),
        pileup_zoom2 = join(OUT_DIR, '{species}', 'figures', '{species}_pileup_zoom2.pdf'),
        rpkm = join(OUT_DIR, '{species}', 'figures', '{species}_rpkm.pdf'),
    params:
        binning = lambda w: config[w.species]['res'][0],
        species = '{species}',
        out_dir = join(OUT_DIR, '{species}', 'figures'),
        circular = lambda w: config[w.species]['circular'],
        window = lambda w: config[w.species]['window'],
        threshold = lambda w: config[w.species]['threshold'],
        unit_length = lambda w: config[w.species]['unit_length'],
    conda: '../envs/bacchus.yaml'
    script: '../scripts/pileup.py'
