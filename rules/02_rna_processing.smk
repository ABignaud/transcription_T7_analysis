#!/bin/env snakemake -s

# Rule to generates the RNA tracks


# Map the reads and sort them.
rule align_rna:
  input:
    index_flag = lambda w: join(TMP, 'ref', f'{samples.species[w.rna_library]}_bt2_index.done'),
    R1 = join(TMP, 'split_reads', '{rna_library}_R1', '{rna_library}_R1.{split}.fq.gz'),
    R2 = join(TMP, 'split_reads', '{rna_library}_R2', '{rna_library}_R2.{split}.fq.gz'),
  output: 
    join(TMP, 'align', '{rna_library}', '{rna_library}.{split}.bam'),
  params:
    index = lambda w: join(TMP, 'ref', f'{samples.species[w.rna_library]}_genome'),
    bt2_presets = config['bowtie2_args'],
  threads: config['threads']
  conda: "../envs/gen_tracks.yaml"
  log: join(OUT_DIR, 'logs', 'rna', '{rna_library}_{split}.log')
  shell:
    """
    bowtie2 {params.bt2_presets} \
            -p {threads} \
            -x {params.index} \
            --maxins 1000 \
            -1 {input.R1} \
            -2 {input.R2} 2> {log} | \
        samtools sort -@ {threads} -n -O BAM - | \
        samtools fixmate -@ {threads} --output-fmt bam -m - - | \
        samtools sort -@ {threads} -O BAM - | \
        samtools markdup -@ {threads} --output-fmt bam -r - - | \
        samtools view -@ {threads} --output-fmt bam -f 2 -q 10 -1 -b - | \
        samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
    """

# rule align_SE_rna:
#   input:
#     index_flag = lambda w: join(TMP, 'ref', f'{samples.species[w.rna_se_library]}_bt2_index.done'),
#     fq = join(TMP, 'split_reads', '{rna_se_library}', '{rna_se_library}.{split}.fq.gz'),
#   output: 
#     join(TMP, 'align', '{rna_se_library}',  '{rna_se_library}.{split}.bam'),
#   params:
#     index = lambda w: join(TMP, 'ref', f'{samples.species[w.rna_se_library]}_genome'),
#     bt2_presets = config['bowtie2_args'],
#   threads: config['threads']
#   conda: "../envs/gen_tracks.yaml"
#   log: join(OUT_DIR, 'logs', 'mapping', '{rna_se_library}_{split}.log')
#   shell:
#     """
#     bowtie2 {params.bt2_presets} \
#             -p {threads} \
#             -x {params.index} \
#             --maxins 1000 \
#             -U {input.fq} 2> {log} | \
#         samtools sort -@ {threads} -O BAM - | \
#         samtools markdup -@ {threads} --output-fmt bam -r - - | \
#         samtools view -@ {threads} --output-fmt bam -q 10 -1 -b - | \
#         samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
#     """

rule merge_split_rna_alignments:
  input:
    expand(
      join(TMP, 'align', '{{rna_library}}',  '{{rna_library}}.{split}.bam'),
      split=split_names
    )
  output: join(TMP, 'bam', '{rna_library}.bam')
  threads: config['threads_small']
  conda: "../envs/gen_tracks.yaml"
  shell:
    """
    samtools merge -n -O BAM -@ {threads} - {input} | 
        samtools sort -@ {threads} --output-fmt bam -l 9 -o {output}
    """

# Build RNA tracks
rule rna_coverage:
    input:
        bam = join(TMP, 'bam', '{rna_library}.bam'),
    output:
        unstranded = join(OUT_DIR, 'RNA_tracks', '{rna_library}_unstranded.bw'),
        fw = join(OUT_DIR, 'RNA_tracks', '{rna_library}_forward.bw'),
        rv = join(OUT_DIR, 'RNA_tracks', '{rna_library}_reverse.bw'),
    threads: config['threads_small']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools index {input.bam} -@ {threads}
        bamCoverage --bam {input.bam} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200
        bamCoverage --bam {input.bam} \
            --outFileName {output.fw} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --filterRNAstrand forward
        bamCoverage --bam {input.bam} \
            --outFileName {output.rv} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --filterRNAstrand reverse
        """


# Merge multiple bam files into one.
rule merge_rna_bam:
    input:
        lambda w: [join(
              TMP, 'bam',  f'{i}.bam'
            ) for i in samples[
                (samples.condition == w.condition) & (samples.type == 'rna')
            ].index]
    output: join(TMP, 'bam', '{condition}.bam'),
    threads: config["threads_small"]
    conda: "../envs/gen_tracks.yaml"
    params: 
        lambda w: " ".join([join(
                TMP, 'bam',  f'{i}.bam'
            ) for i in samples[
                (samples.condition == w.condition) & (samples.type == 'rna')
            ].index])
    shell:'samtools merge -l 9 -O BAM -@ {threads} -f {output} {params}'
        

# RNA coverage for species merge.
rule rna_coverage_merge:
    input:
        join(TMP, 'bam', '{condition}.bam'),
    output:
        unstranded = join(OUT_DIR,  'RNA_tracks', '{condition}_rna_unstranded.bw'),
        fw = join(OUT_DIR,  'RNA_tracks', '{condition}_rna_forward.bw'),
        rv = join(OUT_DIR,  'RNA_tracks', '{condition}_rna_reverse.bw'),
    threads: config['threads_small']
    conda: "../envs/gen_tracks.yaml"
    shell:
        """
        samtools index {input} -@ {threads}
        bamCoverage --bam {input} \
            --outFileName {output.unstranded} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200
        bamCoverage --bam {input} \
            --outFileName {output.fw} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --filterRNAstrand forward
        bamCoverage --bam {input} \
            --outFileName {output.rv} \
            --binSize 1 \
            --numberOfProcessors {threads} \
            --normalizeUsing CPM \
            --extendReads 200 \
            --filterRNAstrand reverse
        """
