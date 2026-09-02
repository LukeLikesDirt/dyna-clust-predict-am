# dyna-clust-predict-am

Predicts sequence similarity cut-offs for classifying/clustering
**arbuscular mycorrhizal (AM) fungi**, using **vsearch global alignment**
and **F-measure optimisation** (following Vu *et al.* 2022). A special-use
adaptation of `dyna-clust-predict`, parameterised for AM fungi and specific V4 primer pairs.

Reference: EUKARYOME SSU v2.0 (Tedersoo *et al.* 2024), trimmed to each primer-pair amplicon.

## Primer pairs

| Name | Primers | Amplicon | Length window |
|---|---|---|---|
| `wanda_aml2` | WANDA (5'-CAGCCGCGGTAATTCCAGCT-3'); AML2 (5'-GAACCCAAACACTTTGGTTTCC-3') | SSU V4 | 440-540bp |
| `amv45nf_amdgr` | AMV4.5NF (5'-AAGCTCGTAGTTGAATTTCG-3'); AMDGR (5'-CCCAACTATCCCTATTAATCAT-3') | SSU V4 (nested) | 150-280bp |

Scales to new primer pairs -- see "Adding a primer
pair" below.

## Output

`output/cutoffs_am.txt` -- master table, every primer pair, every dataset
(Glomeromycota, Endogonomycetes, Eukaryome), family/genus/species.
`output/cutoffs_am_<primer_set>.txt` -- same, one primer pair only.

Columns: `rank`, `higher_rank`, `dataset`, `cut-off`,
`confidence`, `sequence number`, `group number`, `max proportion`,
`primers`, `primer_string`.

Glomeromycota and Endogonomycetes cut-offs give lineage-specific thresholds for classifying and clustering AM fungi. Eukaryome cut-offs (global, across all eukaryotes) are a rough fallback threshold for everything else in a metabarcoding dataset — non-target taxa that still need clustering/classification before they're filtered out.

## Pipeline

Run from the project root.

| Script | Does |
|---|---|
| `01_prepare_reference.sh <primer_set>` | Download EUKARYOME (once, shared), trim to the primer pair -> `data/<primer_set>/` |
| `02_predict_cutoffs.sh <primer_set> <taxon>` | Subset by taxon, compute similarity, predict cut-offs per rank |
| `03_predict_euk_cutoffs.sh <primer_set>` | Global cut-offs across all Eukaryota (balanced hierarchical sampling, `R/subset.R`) |
| `05_combine_cutoffs.sh` | Combine every primer pair/taxon found into `output/cutoffs_am*.txt` |
| `run_all.sh <primer_set> [taxon ...]` | Runs 01 -> 02 (each taxon) -> 03 for one primer pair |

Taxon groups run through 02: glomeromycota, endogonomycetes. Eukaryome uses a per-rank similarity matrix instead of one shared matrix, so it has its own script and config (03_predict_euk_cutoffs.sh, config/taxa/eukaryome.conf) rather than going through 02.

### Adding a primer pair

1. `config/primers/<name>.conf` (copy an existing one). Calibrate the
   length window empirically -- see the comment in
   `config/primers/amv45nf_amdgr.conf`. Set `PRIMER_LABEL` to
   `"<forward name>–<reverse name>"` (en dash).
2. `scripts/run_all.sh <name>`
3. `scripts/05_combine_cutoffs.sh` -- auto-discovers the new primer set.

### Adding a taxon group

`config/taxa/<name>.conf` (copy an existing one). `CASCADE_BASE_RANK`
empty = filter directly from the taxon match (Glomeromycota,
Endogonomycetes); set to e.g. `"phylum"` to require identification at a
coarser rank first before going finer (Fungi).

## Alphanumeric taxon codes

Tedersoo *et al.* (2024) assign alphanumeric codes to undescribed
Glomeromycota/Endogonomycetes family/genus groups (e.g.
`Densosporales.fam02.gen01`). `R/reformat.R` keeps these as usable taxa. Incertae sedis placeholders collapse into "unidentified", since one incertae sedis label can lump together several unrelated groups under the same parent.

## R modules

`R/utils.R` -- `is_identified()` and shared helpers
`R/reformat.R` -- FASTA header parsing, taxonomy extraction
`R/check_annotations.R` -- infraspecific annotation standardisation
`R/dereplicate_lca.R` -- collapses identical sequences, resolves taxonomy to the LCA
`R/combine_lca.R` -- merges the full and partial reference sets, resolving overlaps by LCA
`R/subset.R` -- balanced taxon subsampling (used by Eukaryome only)
`R/compute_sim.R` -- pairwise similarity via vsearch
`R/predict.R` -- cut-off prediction
`R/combine_cutoffs.R` -- combines per-taxon predictions (invoke via `05_combine_cutoffs.sh`)
`R/consolidate_cutoffs.R` -- gap-filling and monotonicity repair for nested cut-off tables (not currently used by the pipeline)

Superseded pre-refactor scripts are in `scripts/backups/` (untracked).

## Citations

Tedersoo, L., Hosseyni Moghaddam, M. S., Mikryukov, V., Hakimzadeh, A., Bahram, M., Nilsson, R. H., ... & Anslan, S. (2024). EUKARYOME: the rRNA gene reference database for identification of all eukaryotes. Database, 2024, baae043. https://doi.org/10.1093/database/baae043

Vu, D., Nilsson, R. H., & Verkley, G. J. (2022). Dnabarcoder: An open‐source software package for analysing and predicting DNA sequence similarity cutoffs for fungal sequence identification. Molecular Ecology Resources, 22(7), 2793-2809. https://doi.org/10.1111/1755-0998.13651
