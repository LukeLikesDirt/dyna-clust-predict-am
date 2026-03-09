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
AM-fungi-focused V4 workflow based on EUKARYOME v.20 (Tedersoo *et al.* 2024) SSU input data trimmed to V4 using the WANDA-AML2 primer pair.

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

## Pipeline steps

  `scripts/01_prepare_euk_V4.sh`            Download EUKARYOME and extract full-legth V4
                                            fragments covering WANDA-AML2

  `scripts/02_predict_glom_cutoffs_V4.sh`   Build Glomeromycota subsets,
                                            compute similarity matrix,
                                            predict global family/genus/species cutoffs

  `scripts/03_predict_euk_cutoffs_V4.sh`    Build Eukaryome subsets
                                            Icompute similarity matrix,
                                            predict global kingdom/phylum/class/order cutoffs

All scripts must be run from the **project root directory**.

### R modules

  `R/utils.R`               Shared utility functions (e.g., `is_identified()`)

  `R/reformat.R`            FASTA header parsing and taxonomy extraction

  `R/check_annotations.R`   Infraspecific annotation standardisation

  `R/subset.R`              Balanced taxon subset generation

  `R/compute_sim.R`         Pairwise similarity computation using vsearch

  `R/predict.R`             Cut-off prediction (parallel or sequential)
  
```

## Citation

Tedersoo, L., Hosseyni Moghaddam, M. S., Mikryukov, V., Hakimzadeh, A., Bahram, M., Nilsson, R. H., ... & Anslan, S. (2024). EUKARYOME: the rRNA gene reference database for identification of all eukaryotes. Database, 2024, baae043. https://doi.org/10.1093/database/baae043

Vu, D., Nilsson, R. H., & Verkley, G. J. (2022). Dnabarcoder: An open‐source software package for analysing and predicting DNA sequence similarity cutoffs for fungal sequence identification. Molecular Ecology Resources, 22(7), 2793-2809 https://doi.org/10.1111/1755-0998.13651
