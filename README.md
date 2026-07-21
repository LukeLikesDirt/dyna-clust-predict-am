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

Predicted similarity cut-offs from this pipeline are provided in `output/cutoffs_glom_V4.txt` and can be used directly for classification and clustering without re-running the pipelines. Reference sequences and taxonomy (FASTA and classification files) will be made available on Figshare; the README will be updated with links once they are published.

This project builds directly on ideas and implementation patterns from [dnabarcoder](https://github.com/vuthuyduong/dnabarcoder) and [dyna-clust-predict](https://github.com/LukeLikesDirt/dyna-clust-predict)

Main pipelines in this repository:

- `scripts/02_predict_glom_cutoffs_V4.sh`
  Combined preparation + global cutoff prediction for Glomeromycota at
  order, family, genus, and species ranks.
- `scripts/03_predict_endo_cutoffs_V4.sh`
  Combined preparation + global cutoff prediction for Endogonomycetes
  (class within phylum Mucoromycota) at order, family, genus, and species
  ranks, following Tedersoo *et al.* (2024).
- `scripts/04_predict_fun_cutoffs_V4.sh`
  Combined preparation + global cutoff prediction for Fungi at
  phylum, class, order, family, genus, and species ranks.
- `scripts/05_predict_euk_cutoffs_V4.sh`
  Combined preparation + global cutoff prediction for Eukaryome at
  kingdom, phylum, class, and order ranks.

## Pipeline steps

  `scripts/01_prepare_euk_V4.sh`            Download EUKARYOME and extract full-legth V4
                                            fragments covering WANDA-AML2

  `scripts/02_predict_glom_cutoffs_V4.sh`   Build Glomeromycota subsets,
                                            compute similarity matrix,
                                            predict global order/family/genus/species cutoffs

  `scripts/03_predict_endo_cutoffs_V4.sh`   Build Endogonomycetes subsets,
                                            compute similarity matrix,
                                            predict global order/family/genus/species cutoffs

  `scripts/04_predict_fun_cutoffs_V4.sh`    Build Fungi subsets,
                                            compute similarity matrix,
                                            predict global phylum/class/order/family/genus/species cutoffs

  `scripts/05_predict_euk_cutoffs_V4.sh`    Build Eukaryome subsets,
                                            compute similarity matrix,
                                            predict global kingdom/phylum/class/order cutoffs

All scripts must be run from the **project root directory**.

### Alphanumeric taxon codes (Endogonomycetes and Glomeromycota)

Tedersoo *et al.* (2024, *MycoKeys* 107: 273–325) propose an alphanumeric
coding system (e.g. `Densosporales.fam02.gen01`) to communicate undescribed
family- and genus-level groups in Glomeromycota and Endogonomycetes ahead of
formal description. Since the paper's authors also curate EUKARYOME, these
codes appear directly in EUKARYOME's placeholder taxonomy fields.
`R/reformat.R` preserves these codes (converting `.` to `_`, e.g.
`Densosporales_fam02_gen01`) rather than collapsing them to "unidentified",
so they are retained as usable taxonomic groups throughout the pipeline.
Genuine unplaced/incertae sedis placeholders (e.g. `.reg` codes) are still
treated as unidentified.

### R modules

  `R/utils.R`               Shared utility functions (e.g., `is_identified()`)

  `R/reformat.R`            FASTA header parsing and taxonomy extraction

  `R/check_annotations.R`   Infraspecific annotation standardisation

  `R/subset.R`              Balanced taxon subset generation

  `R/compute_sim.R`         Pairwise similarity computation using vsearch

  `R/predict.R`             Cut-off prediction (parallel or sequential)

## Citations

Tedersoo, L., Hosseyni Moghaddam, M. S., Mikryukov, V., Hakimzadeh, A., Bahram, M., Nilsson, R. H., ... & Anslan, S. (2024). EUKARYOME: the rRNA gene reference database for identification of all eukaryotes. Database, 2024, baae043. https://doi.org/10.1093/database/baae043

Vu, D., Nilsson, R. H., & Verkley, G. J. (2022). Dnabarcoder: An open‐source software package for analysing and predicting DNA sequence similarity cutoffs for fungal sequence identification. Molecular Ecology Resources, 22(7), 2793-2809. https://doi.org/10.1111/1755-0998.13651
