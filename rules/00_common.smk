#!/bin/env snakemake -s

# Rules to generates the build bowtie2 index and split fastq for alignement.


# Build bowtie2 index of the reference genome.
rule bt2_index:
    input: lambda w: join(REF_DIR, config[w.species]['ref']),
    output: touch(join(TMP, 'ref', '{species}_bt2_index.done')),
    params:
        idx = join(TMP, 'ref', '{species}_genome'),
    threads: config['threads']
    conda: "../envs/hic_processing.yaml"
    shell: "bowtie2-build --threads {threads} {input} {params.idx}"


# Make splits from Hi-C fastq files to speed up mapping. [between 1 and 999].
N_SPLITS = config['n_splits']
split_names = [f'part_{s:03}' for s in range(1, N_SPLITS + 1)] #base names of split files.

rule split_fastq:
    input: join(FASTQ_DIR, '{library}_R{end}.fq.gz'),
    output: expand(join(TMP, 'split_reads', '{{library}}_R{{end}}', '{{library}}_R{{end}}.{split}.fq.gz'), split=split_names),
    params:
        n_splits = N_SPLITS,
        split_dir = join(TMP, 'split_reads', "{library}_R{end}"),
    message: "Splitting {wildcards.library}_{wildcards.end} into {params.n_splits} split fastq."
    conda: "../envs/split_fastq.yaml"
    threads: config['threads_split']
    shell:
      """
      mkdir -p {params.split_dir}
      # 100 split fastqs will be created with name pattern 001.fq - 100.fq
      seqkit split2 -p {params.n_splits} \
          -w 0 \
          -f \
          -j {threads} \
          -1 {input} \
          -O {params.split_dir}
      """


# Make split om single ends libraries.
# rule split_fastq_se:
#     input: join(FASTQ_DIR, '{rna_se_library}.fq.gz'),
#     output: expand(join(TMP, 'split_reads', '{{rna_se_library}}', '{{rna_se_library}}.{split}.fq.gz'), split=split_names),
#     params:
#         n_splits = N_SPLITS,
#         split_dir = join(TMP, 'split_reads', "{rna_se_library}"),
#     message: "Splitting {wildcards.rna_se_library} into {params.n_splits} split fastq."
#     conda: "../envs/split_fastq.yaml"
#     threads: config['threads_split']
#     shell:
#       """
#       mkdir -p {params.split_dir}
#       # 100 split fastqs will be created with name pattern 001.fq - 100.fq
#       seqkit split2 -p {params.n_splits} \
#           -w 0 \
#           -f \
#           -j {threads} \
#           -1 {input} \
#           -O {params.split_dir}
#       """
