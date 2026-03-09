# dyna-clust-predict-am

`dyna-clust-predict-am` is a **special-use adaptation** of
`dyna-clust-predict` for predicting similarity cut-offs in
**arbuscular mycorrhizal (AM) fungi** using the
**SSU V4 WANDA-AML2 fragment**.

It predicts optimal sequence similarity cut-offs for classification and
clustering using **vsearch global alignment** and
**F-measure optimisation**, following Vu *et al.* (2022).

## Overview

This repository predicts optimal sequence similarity thresholds for an
AM-fungi-focused V4 workflow based on EUKARYOME SSU V4 input data.

Main pipelines in this repository:

-   `scripts/02_predict_glom_cutoffs_V4.sh`
  Combined preparation + global cutoff prediction for Glomeromycota at
  family, genus, and species ranks.
-   `scripts/03_predict_euk_cutoffs_V4.sh`
  Combined preparation + global cutoff prediction for Eukaryome at
  kingdom, phylum, class, and order ranks.

## Attribution

This project builds directly on ideas and implementation patterns from:

- [dnabarcoder](https://github.com/vuthuyduong/dnabarcoder)
- [dyna-clust-predict](https://github.com/LukeLikesDirt/dyna-clust-predict)

Relative to those upstream projects, this repository is specialised for
AM fungi and the V4 WANDA-AML2 marker fragment.

Core characteristics retained/adapted here include:

-   Uses **vsearch** for global alignment
-   Supports scalable parallel computation
-   Integrates similarity prediction within an R-based workflow

## Pipeline steps

  `scripts/01_prepare_euk_V4.sh`            Download/process EUKARYOME V4
                                            and write `eukaryome_V4.*`

  `scripts/02_predict_glom_cutoffs_V4.sh`   Build Glomeromycota subsets,
                                            compute family similarity matrix,
                                            predict global family/genus/species cutoffs

  `scripts/03_predict_euk_cutoffs_V4.sh`    Build required global Eukaryome
                                            ID sets, subset on the fly,
                                            compute per-rank similarity,
                                            predict global kingdom/phylum/class/order cutoffs

All scripts must be run from the **project root directory**.

### R modules

  `R/utils.R`               Shared utility functions (e.g., `is_identified()`)

  `R/reformat.R`            FASTA header parsing and taxonomy extraction

  `R/check_annotations.R`   Infraspecific annotation standardisation

  `R/subset.R`              Balanced taxon subset generation

  `R/compute_sim.R`         Pairwise similarity computation using vsearch

  `R/predict.R`             Cut-off prediction (parallel or sequential)
  
## Directory structure

  data/
  ├── eukaryome_V4.fasta
  ├── eukaryome_V4.classification
  ├── glomeromycota_*_V4.fasta
  ├── glomeromycota_*_V4.classification
  ├── *.cutoffs.json
  └── *.predicted

## Environment setup

A conda environment specification is provided in `environment.yml`.

[mamba](https://mamba.readthedocs.io/) is strongly recommended over
`conda` for faster dependency resolution.

### Local installation

``` bash
mamba env create -f environment.yml
conda activate dyna_clust_predict
```

### HPC (SLURM)

Run from the project root:

``` bash
sbatch scripts/00_create_dyna_clust_env.sh
conda activate dyna_clust_predict
```

## Quick start

From the project root:

```bash
bash scripts/01_prepare_euk_V4.sh
bash scripts/02_predict_glom_cutoffs_V4.sh
bash scripts/03_predict_euk_cutoffs_V4.sh
```

Both prediction scripts write timestamped logs under `logs/`.

# Key parameters

## Sequence selection

Used in: `subset.R` (via `04_prepare_subsets.sh`)

These parameters control taxonomic balance and sampling constraints
prior to similarity prediction.

    --min_subgroups   INT    Minimum unique child taxa per parent taxon (default: 10)
    --min_sequences   INT    Minimum sequences per parent taxon after proportion cap (default: 30)
    --max_sequences   INT    Maximum sequences per parent taxon; excess are balanced round-robin downsampled across child taxa (default: 25000)
    --max_proportion  FLOAT  Maximum fraction a child taxon may represent (default: 0.5)

Example:

``` bash
Rscript R/subset.R \
  --fasta_in input.fasta \
  --classification_in input.classification \
  --min_subgroups 10 \
  --min_sequences 30 \
  --max_sequences 25000 \
  --max_proportion 0.5 \
  --output_dir output
```

## Similarity prediction

Used in: `predict.R` (via `06a_predict_cutoffs.sh` / `06b_predict_cutoffs_region.sh`)

These parameters select the rank combination to predict:

    --rank          STR    Target rank(s), comma-separated (e.g. species,genus) [required]
    --higher_rank   STR    Parent rank(s) for local prediction, comma-separated (e.g. genus,family).
                           Omit to run a single global prediction across all sequences. [default: ""]

These parameters control the similarity sweep:

    --start_threshold   FLOAT   Starting similarity threshold (default: 0.0)
    --end_threshold     FLOAT   Ending similarity threshold (default: 1.0)
    --step              FLOAT   Threshold step size (default: 0.001)
    --min_cutoff        FLOAT   Min cutoff value to report in output (default: 0.0)

The same filtering thresholds used during sequence selection (`subset.R`) also apply here, but as
strict dataset filters rather than balancing rules — datasets that fail a threshold are
skipped entirely rather than downsampled:

    --min_group_no    INT    Min unique child taxa required to report a cutoff (default: 10)
    --min_seq_no      INT    Min sequences required to report a cutoff (default: 30)
    --max_seq_no      INT    Max sequences per dataset; excess is randomly sampled (default: 25000)
    --max_proportion  FLOAT  Skip datasets where the dominant group exceeds this fraction (default: 1.0)

Execution controls:

    --run_parallel   yes/no   Parallel dataset processing via furrr/future (default: yes)
    --n_cpus         INT      Workers (parallel) or vsearch threads (sequential) (default: all−1)
    --tmp_dir        DIR      Directory for temporary vsearch output (default: ./tmp)

Example:

``` bash
Rscript R/predict.R \
  --input data/full_ITS/eukaryome_ITS.fasta \
  --classification data/full_ITS/eukaryome_ITS.classification \
  --rank species \
  --higher_rank genus \
  --start_threshold 0.9 \
  --end_threshold 1.0 \
  --step 0.001 \
  --run_parallel yes \
  --n_cpus 80 \
  --out data/full_ITS \
  --prefix eukaryome_ITS
```

## Citation

Vu, D., Nilsson, R. H., & Verkley, G. J. (2022). Dnabarcoder: An open‐source software package for analysing and predicting DNA sequence similarity cutoffs for fungal sequence identification. Molecular Ecology Resources, 22(7), 2793-2809 https://doi.org/10.1111/1755-0998.13651
