# T7 promoter transcription analysis

In this repository, we present the code for the analysis of study of the transcription's impact on Escherichia coli chromosome described in this [paper](https://www.biorxiv.org/content/biorxiv/early/2022/09/16/2022.09.16.507559.full.pdf).

## Requirements

- [Conda](https://docs.conda.io/en/latest/), Package, dependency and environment management for any language.
- [Snakemake](https://snakemake.readthedocs.io/en/stable/), Workflow management system.
- [Go](https://go.dev/doc/install), programming language.
- [Dnaglider](https://github.com/cmdoret/dnaglider), Command line utility to compute sliding window genome statistics from a fasta file.

## Usage

To launch the analysis clone this repository.

```bash
git clone https://github.com/ABignaud/transcription_T7_analysis.git
cd transcription_T7_analysis
snakemake -j 32 --use-conda 
```

## Citation

[Transcriptional units form the elementary constraining building blocks of the bacterial chromosome](https://www.biorxiv.org/content/biorxiv/early/2022/09/16/2022.09.16.507559.full.pdf), A. Bignaud, C. Cockram, E. Allemand, J. Mozziconnacci, O. Espeli, R. Koszul, *BioArXiv*, 2022.
