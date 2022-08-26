#!/bin/env snakemake -s

# Rule to generates the RNA tracks

# Build bowtie2 index of the reference genome.
# rule bt2_index:
#   input: lambda w: config['ref'][w.species]
#   output: touch(join(OUT_DIR, '{species}', 'ref', 'bt2_index.done'))
#   params:
#     idx = join(OUT_DIR, '{species}', 'ref', 'genome')
#   # singularity: "docker://koszullab/hicstuff:v3.1.0"
#   # conda: "../envs/hic_processing.yaml"
#   shell: "bowtie2-build {input} {params.idx}"


# # Make splits from Hi-C fastq files to speed up mapping. [between 1 and 999]
# N_SPLITS = 4
# split_names = [f'part_{s:03}' for s in range(1, N_SPLITS + 1)] #base names of split files


# Map the reads and sort them.
rule align_rna:
  input:
    index_flag = lambda w: join(OUT_DIR, samples.species[w.rna_library], 'ref', 'bt2_index.done'),
    R1 = lambda w: join(samples.folder[w.rna_library], '{rna_library}_R1.fq.gz'),
    R2 = lambda w: join(samples.folder[w.rna_library], '{rna_library}_R2.fq.gz'),
  output: 
    join(OUT_DIR, 'bam', '{rna_library}.bam'),
  params:
    index = lambda w: join(OUT_DIR, samples.species[w.rna_library], 'ref' , 'genome'),
    bt2_presets = config['bowtie2_args'],
  threads: config['threads']
#   singularity: "docker://koszullab/hicstuff:v3.1.0"
#   conda: "../envs/hic_processing.yaml"
  log: "logs/mapping/{rna_library}.log"
  shell:
    """
    bowtie2 {params.bt2_presets} \
            -p {threads} \
            -x {params.index} \
            --maxins 1000 \
            -1 {input.R1} \
            -2 {input.R2} 2> {log} | \
      samtools sort -@ {threads} -O BAM - | \
      samtools fixmate -@ {threads} --output-fmt bam -m - - | \
      samtools markdup -@ {threads} --output-fmt bam -r - - | \
      samtools view -@ 16 --output-fmt bam -f 2 -q 10 -1 -b - | \
      samtools sort -@ 16 --output-fmt bam -l 9 -o {output} 
    """


# Build RNA tracks
rule rna_coverage:
    input:
        bam = join(OUT_DIR, 'bam', '{rna_library}.bam'),
    output:
        unstranded = join(OUT_DIR, 'tracks', '{rna_library}_unstranded.bw'),
        forward = join(OUT_DIR, 'tracks', '{rna_library}_forward.bw'),
        reverse = join(OUT_DIR, 'tracks', '{rna_library}_reverse.bw'),
    threads: config['threads']
    shell:
        """
        bamCoverage --bam {input} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
            --ignoreDuplicates
        bamCoverage --bam {input} \
            --outFileName {output.forward} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
            --ignoreDuplicates\
            --filterRNAstrand forward
        bamCoverage --bam {input} \
            --outFileName {output.reverse} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
            --ignoreDuplicates\
            --filterRNAstrand reverse
        """


# Merge multiple bam files into one.
rule merge_rna_bam:
    input:
        lambda w: [join(
                OUT_DIR, 'bam', f'{i}.bam'
            ) for i in samples[
                (samples.species == w.species) & (samples.type == 'rna')
            ].index]
    output: join(OUT_DIR, '{species}', 'bam', 'rna_merge.bam')
    threads: config["threads"]
    params: 
        lambda w: " ".join([join(
                OUT_DIR, 'bam', f'{i}.bam'
            ) for i in samples[
                (samples.species == w.species) & (samples.type == 'rna')
            ].index])
    shell:'samtools merge -l 9 -O BAM -@ {threads} -f {output} {params}'
        

# RNA coverage for species merge.
rule rna_coverage:
    input:
        join(OUT_DIR, '{species}', 'bam', 'rna_merge.bam'),
    output:
        unstranded = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_unstranded.bw'),
        forward = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_forward.bw'),
        reverse = join(OUT_DIR,  '{species}', 'tracks', '{species}_rna_reverse.bw'),
    threads: config['threads']
    shell:
        """
        bamCoverage --bam {input} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
            --ignoreDuplicates
        bamCoverage --bam {input} \
            --outFileName {output.forward} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
            --ignoreDuplicates\
            --filterRNAstrand forward
        bamCoverage --bam {input} \
            --outFileName {output.reverse} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads \
            --ignoreDuplicates\
            --filterRNAstrand reverse
        """