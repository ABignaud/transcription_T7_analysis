#!/bin/env snakemake -s

# Rules to plot the pileup from the genomics tracks. 

# Do the correlations between hic and rna tracks.
rule corr:
    input:
        hic = join(OUT_DIR, 'HiC', '{species}.mcool'),
        rna = join(OUT_DIR,  'RNA_tracks', '{species}_rna_unstranded.bw'),
    output:
        mat_fig = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_map.pdf'),
        corr_fig = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_corr.pdf'),
    params:
        binning = lambda w: config[w.species]['res'][0],
        species = '{species}',
        circular = lambda w: config[w.species]['circular'],
    conda: '../envs/bacchus.yaml'
    script: '../scripts/corr_hic_rna.py'

# Generates pileups plot of HEGs from the output of previous tracks for each 
# species.
rule pileup:
    input:
        hic = join(OUT_DIR, 'HiC', '{species}.mcool'),
        rna = join(OUT_DIR,  'RNA_tracks', '{species}_rna_unstranded.bw'),
        annotation = lambda w: join(REF_DIR, config[w.species]['annotation']),
    output:
        pileup_pos_TU = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_pos_TU.pdf'),
        pileup_neg_TU = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_neg_TU.pdf'),
        pileup = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup.pdf'),
        pileup_zoom = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_zoom.pdf'),
        pileup_zoom2 = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_zoom2.pdf'),
        rpkm = join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_rpkm.pdf'),
    params:
        binning = lambda w: config[w.species]['res'][0],
        species = '{species}',
        out_dir = join(OUT_DIR, '{species}', 'figures'),
        circular = lambda w: config[w.species]['circular'],
        window = lambda w: config[w.species]['window'],
        threshold = lambda w: config[w.species]['threshold'],
        unit_length = lambda w: config[w.species]['unit_length'],
    conda: '../envs/bacchus.yaml'
    script: '../scripts/pileup_species.py'


rule all_pileup:
    input:
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_map.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_corr.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_pos_TU.pdf'), species=species),
        expand(join(OUT_DIR, 'figures', 'pileup', '{species}', '{species}_pileup_neg_TU.pdf'), species=species),
    output: join(TMP, 'pileup.done')
