#!/bin/env snakemake -s

# Rules to generates contact map plot.

# Generate GC content using dnaglider
rule generate_gc:
    input:
        fasta = join(REF_DIR, config['E_coli']['ref'])
    output:
        gc = join(OUT_DIR, 'E_coli', 'ref', 'gc_content.tsv')
    threads: config['threads']
    conda: "../envs/bacchus.yaml"
    shell: 
        """
        go get -u github.com/cmdoret/dnaglider/dnaglider
        dnaglider \
            -fasta {input.fasta} \
            -fields GC \
            -out {output.gc} \
            -stride 100 \
            -threads {threads} \
            -window 500
        """

# Rule plot the matrix.
rule plot_matrix:
    input:
        mat = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
    params:
        res = lambda w: int(w.res) * 1000,
        cmap = config['cmap'],
        title = 'E. coli WT - binning {res}kb'
    output:
        join(OUT_DIR, 'figures', 'E_coli_contact_map', 'WT_{res}kb.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: config["threads"]
    script: '../scripts/plot_matrix.py'

# EPODs positions download from Freddolino et al., PLOS, 2021.
# Rule plot the zoom of matrices with different metrics.
rule plot_zoom:
    input:
        annotation = join(REF_DIR, config['E_coli']['annotation']),
        cool_wt = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
        cool_rf = join(OUT_DIR, 'HiC', 'E_coli_rif.mcool'),
        rna_wt = join(OUT_DIR, 'RNA_tracks', 'E_coli_rna_unstranded.bw'),
        rna_rf = join(OUT_DIR, 'RNA_tracks', 'E_coli_rif_rna_unstranded.bw'),
        gc = join(OUT_DIR, 'E_coli', 'ref', 'gc_content.tsv'),
        EPODs = join(REF_DIR, config['epods']),
    params:
        res = config['E_coli']['annotation'][0],
        cmap = config['cmap'],
        width = config['zoom_width'],
        positions = lambda w: [w.pos1, w.pos2],
        out_dir = join(OUT_DIR, 'figures', 'E_coli_contact_map'),
    output:
        join(OUT_DIR, 'figures', 'E_coli_contact_map', 'mat_zoom', 'region_{pos1}_{pos2}.pdf'),
    conda: '../envs/bacchus.yaml'
    threads: config["threads"]
    script: '../scripts/plot_matrix.py'

# Run Macrodomain and CIDs detection using directional index and comapre it to 
# Lioy dataset.
rule CIDs_comparison_Lioy:
    input:
        mat = join(OUT_DIR, 'HiC', 'E_coli.mcool'),
    params:
        res = 5000,
        cmap = config['cmap'],
    output:
        macrodomain = join(OUT_DIR, 'figures', 'E_coli_contact_map', 'CIDs', 'macrodomaim.pdf'),
        cid = join(OUT_DIR, 'figures', 'E_coli_contact_map', 'CIDs', 'CIDs.pdf'),
        compare = join(OUT_DIR, 'figures', 'E_coli_contact_map', 'CIDs', 'comaprison_Lioy.pdf'),
    conda: '../envs/bacchus.yaml'
    script: '../scripts/CIDs_analysis_Lioy.py'

# Run high resolution CIDs analysis based on Directional index and insulation 
# score.
rule CIDs_analysis:

# Detect TIDs and compute some metrics on them 
rule TIDs_analysis: