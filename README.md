# Genetic control and tissue specificity of transcriptomic rhythmicity: analysis code

Code for every analysis, figure and supplementary table in the manuscript.
Scripts are numbered in run order; each opens with a header giving what it does,
how it is run, and which figure or table it produces. Everything is in R.

## Setup

- **R ≥ 4.4** with: tidyverse, data.table, genio, broom, lmtest, future,
  future.apply, lsa, dryR, MetBrewer, corrplot, ggExtra, ggrepel, patchwork
  (≥ 1.2), openxlsx, httr, jsonlite.
- **`00_config.R`** holds every path and analysis constant (TSS window, MAF,
  cd-eQTL threshold, minimum genotype group size, classification FDR, the CCG
  list, the two example loci). Each script sources it, and it sources
  **`00_helperFunctions.R`**, the shared model-fitting functions.
- Three environment variables relocate the tree without editing code:
  `CIRCADIAN_ROOT` (input data), `CIRCADIAN_RESULTS` (intermediate results) and
  `CIRCADIAN_SUBMISSION` (figures and tables).
- `*.sh` files are SLURM wrappers. Submit them from this folder (`sbatch 01_rhythmicGenes.sh`);
  logs go to `logs/`. Scripts without a wrapper run in seconds to minutes
  (`Rscript 05_supTable4_coefficients.R`), but still submit them to a compute node.
- Any script that loads the genotypes (`loadGenotypes()`) reads the full 8.25M-variant
  plink set and needs about 100 GB of memory.

## Inputs

| Input | Source |
|---|---|
| Normalised expression and covariates (5 genotype PCs, 10 PEER factors, sex, WGS platform, PCR protocol) | GTEx v10, `GTEx_Analysis_v10_eQTL_expression_matrices.tar`, `GTEx_Analysis_v10_eQTL_covariates.tar` |
| Genotypes and variant lookup table | GTEx WGS (dbGaP controlled access) |
| Donor inferred phase (`Data/CHIRAL/DIP.RData`) | CHIRAL, Talamanca et al. 2023 |
| Sample size per tissue (`Data/PhenotypeFile/GTEx_v10_sample_size.csv`) | GTEx v10 portal |
| Pan-UK Biobank associations of the cd-eQTL lead variants (`Results/phewas/inferred/`: `phewasClean/*.tsv`, `Pan-uk-phenotype.csv`, `Neale_rsID_key.csv`) | Pan-UK Biobank, European ancestry, 1,085 phenotypes |
| Associations named by Chen et al. (`Data/published/associationListPublished.csv`) | Chen et al. 2025 |

## Run order

| Step | Script | Wrapper | What it does |
|---|---|---|---|
| 01 | `01_rhythmicGenes.R` | `.sh`, array 1-50 | Cosinor vs linear LRT for every gene in each tissue; variants within 50 kb of rhythmic gene TSSs |
| 02 | `02_cdeQTLscan.R` | `.sh`, array 1-50 | Three nested models per gene-variant pair (MAF ≥ 10%); the interaction LRT is the cd-eQTL test |
| 03 | `03_geneSymbolMap.R` | none (needs internet) | Ensembl ID to gene symbol map (Ensembl REST) |
| 04 | `04_cdeQTLtables.R` | `.sh` | Rhythmic genes (Bonferroni < 0.05 within tissue); cd-eQTL (p < 1e-4); the lead variant per gene x tissue |
| 05 | `05_supTable4_coefficients.R` | none | Interaction-model coefficients of every cd-eQTL |
| 06 | `06_figure2_heatmap.R` | `.sh` | Interaction p of each cd-eQTL gene in every tissue; heatmap |
| 07 | `07_binghamClassification.R` | `.sh` | Bingham et al. amplitude / acrophase F-tests on residualised per-genotype cosinors; zero-amplitude gating; classification |
| 08 | `08_supTable6_classification.R` | none | Classification table |
| 09 | `09_methodComparison.R` | `.sh` | Per-genotype amplitude and acrophase: interaction model vs residualised cosinor |
| 10 | `10_suppFigure2_methodComparison.R` | none | Their agreement (Pearson r) |
| 11 | `11_dryRvalidation.R` | `.sh` | dryR model selection at the 54 loci |
| 12 | `12_suppFigure1_dryRmodelGrid.R` | `.sh` (R 4.4.2) | BIC weight of every dryR model |
| 13 | `13_supTable5_dryR.R` | none | dryR validation table |
| 14 | `14_multiTissueRegression.R` | `.sh` | Interaction model of each cd-eQTL in every tissue; tissue x tissue cosine similarity (49 rhythmic tissues) |
| 15 | `15_tissueSimilarity.R` | none | CCG vs CRG Wilcoxon on mean similarity (all four coefficients, and interaction terms only); example-locus matrices |
| 16 | `16_figure4_suppFigure3.R` | none | Similarity heatmaps for BMAL1 and TMED10 |
| 17 | `17_figure3.R` | `.sh` | Example loci: residualised expression, confidence clocks, amplitude and acrophase across tissues |
| 18 | `18_supTables7to10.R` | `.sh` | Example loci: per-genotype amplitude and acrophase in the 49 tissues |
| 19 | `19_figure5_phewas.R` | none | PheWAS (p < 0.05/1,085) |
| 20 | `20_chenComparison.R` | `.sh` | Interaction model applied to the Chen et al. associations; AIF1 example |
| 21 | `21_chenComparison.R` | `.sh` | Chen et al. Supplementary Data 2 tests within the TSS window, merged per tissue with the cd-eQTL scan |
| 21 | `21_supplementaryWorkbook.R` | none | All supplementary tables in one workbook with legends |
| 22 | `22_chenPipeline.R` | `.sh` | Chen et al. rhyQTL pipeline (genotype-group, cosinor, dryR, G-test, HANOVA filters) applied to the 54 independent cd-eQTL |

Dependencies: 02 needs 01; 03 needs 01; 04 needs 02 and 03; all later steps need 04.
Beyond that: 08 needs 07, 10 needs 09, 12 and 13 need 11, 15 needs 14, 16 needs 15,
17 needs 14, and 21 needs every table.

## Figures and tables

| Item | Script | Output (under `Writing/Submission/`) |
|---|---|---|
| Figure 1a | drawn by hand (schematic) | |
| Figure 1b | 04 | `Figures/Figure_1b_RhythmicGenes.jpeg` |
| Figure 2 | 06 | `Figures/Figure_2_cdEQTL_heatmap.jpeg` |
| Figure 3a-b | 17 | `Figures/Figure_3{a,b}_*_expression`, `_confidenceClock` |
| Figure 3c-f | 17 (data from 14) | `Figures/Figure_3{c,d}_*_amplitude`, `Figure_3{e,f}_*_phase` |
| Figure 4 | 16 (data from 14) | `Figures/Figure_4_cosineSimilarity` |
| Figure 5 | 19 | `Figures/Figure_5_PheWAS` |
| Supplementary Figure 1 | 12 (data from 11) | `Supplementary_Figures/SuppFigure_1_dryRModelSelection` |
| Supplementary Figure 2 | 10 (data from 09) | `Supplementary_Figures/SuppFigure_2{a,b}_EffectSizeComparison_*` |
| Supplementary Figure 3 | 16 (data from 15) | `Supplementary_Figures/SuppFigure_3_cosineSimilarityInteraction` |
| Supplementary Figure 4 | 20 | `Supplementary_Figures/SuppFigure_4_AIF1_ChenComparison` |
| Supplementary Tables 1-3 | 04 | `Supplementary_Tables/SupTable_{1,2,3}_*.csv` |
| Supplementary Table 4 | 05 | `SupTable_4_cdEQTL_coefficients.csv` |
| Supplementary Table 5 | 13 (data from 11) | `SupTable_5_dryRvalidation.csv` |
| Supplementary Table 6 | 08 (data from 07) | `SupTable_6_cdEQTL_classification.csv` |
| Supplementary Tables 7-10 | 18 | `SupTable_{7,8,9,10}_*.csv` |
| Supplementary Tables 11-12 | 15 (data from 14) | `SupTable_{11,12}_*_cosineSimilarity.csv` |
| Supplementary Table 13 | 19 | `SupTable_13_PheWAS.csv` |
| Supplementary Table 14 | 21 (`21_chenComparison.R`) | `Results/published/chenComparison/chen_cdeQTL_merged.csv` |
| Supplementary Table 15 | 22 | `Results/published/chenPipeline/chenPipeline_54cdeQTL.csv` |
| Workbook | 21 | `Supplementary_Tables.xlsx` |

Figures 3, 4 and 5 were assembled or annotated by hand from these panels
(headers, colour boxes, category bars). Figure legends are in the manuscript.

