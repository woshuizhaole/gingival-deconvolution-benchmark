# Gingival deconvolution benchmark — analysis code

Code accompanying the manuscript *A tissue-matched single-cell reference and benchmark for deconvolving human gingival mucosa*.

The repository contains the analysis pipeline only: single-cell reference construction, the leave-out mixture benchmark, multi-solver comparison, cohort deconvolution, methylation anchoring, hierarchical recovery, spatial validation, and the summary tables. Figure rendering and manuscript-generation scripts are not part of this repository.

## Repository layout

```
code/pipeline/    reference construction and preprocessing
  01_load_qc_doublets.py        per-sample QC and doublet removal
  02_integrate_cluster.py       integration and clustering
  03_annotate.py                marker-based annotation
  04_refine_annotation.py       annotation refinement
  05_export_reference.py        export annotated atlas
  06_build_reference_artifacts.R  signature matrix and InstaPrism reference
  07_process_gse152042.py       held-out atlas for pseudo-bulk mixtures
  09_preprocess_bulk.R          bulk cohort preprocessing
  12_build_dwls_signature.R     DWLS signature matrix

code/analysis/    benchmarking and downstream analyses
  08_pseudobulk_benchmark.R     leave-out benchmark (CCC, Pearson, bias)
  09_multisolver_benchmark.R    InstaPrism, BayesPrism, MuSiC, DWLS, SVR, EPIC, NNLS
  10_deconvolve_and_composition.R  cohort composition and differential composition
  14_methylation_orthogonal.R   methylation-based immune and epithelial fractions
  16_synthesis.R                cross-cohort summary and direction re-estimation
  18_spatial_validation.R       per-section spatial composition
  E1_mrna_correction.R          read-to-cell scale correction
  E2_hierarchical_decollapse.R  hierarchical recovery of the immune compartment
  E3_anchor.R                   orthogonal anchoring of subtype fractions
  E4_multimethod.R              multi-method concordance
  E6_spatial_perspot.R          per-spot spatial deconvolution
  E7_flip_audit.R               direction audit across reference and cohort
  E9_detection_limit.R          single-type titration and detection limits

code/utils/       shared solver wrappers
```

## Data

All data are public. Accessions:

| Dataset | Role |
|---|---|
| GSE164241 | gingival single-cell reference (21 donors, 76,925 cells) |
| GSE152042 | independent gingival atlas for leave-out mixtures |
| GSE16134, GSE10334 | bulk microarray cohorts |
| GSE173078, GSE223924, GSE33774, GSE23586 | bulk RNA-seq and microarray cohorts |
| GSE173081 | genome-wide methylation, same donors as GSE173078 |
| GSE206621 | gingival Visium sections |

Download each series from GEO and place it under `data/raw/`. Processed objects are written under `data/processed/`, `data/interim/`, and `results/`.

## Environment

Two conda environments are provided:

```bash
conda env create -f env/environment_py.yml
conda env create -f env/environment_r.yml
```

The Python environment includes scanpy, scvi-tools, anndata, and scrublet; the R environment includes SimBu, InstaPrism, MuSiC, EPIC, EpiDISH, TOAST, and muscat. DWLS and BayesPrism are installed from their authors' repositories; `code/utils/deconv_utils.R` wraps the dampened weighted least-squares solver.

## Running

Each script sets a project root near the top (`PROJ` or `P`). Point it at your local copy of this repository before running. Scripts are run in this order:

```bash
python code/pipeline/01_load_qc_doublets.py
python code/pipeline/02_integrate_cluster.py
python code/pipeline/03_annotate.py
python code/pipeline/04_refine_annotation.py
python code/pipeline/05_export_reference.py
Rscript code/pipeline/06_build_reference_artifacts.R
python code/pipeline/07_process_gse152042.py
Rscript code/pipeline/09_preprocess_bulk.R
Rscript code/pipeline/12_build_dwls_signature.R

Rscript code/analysis/08_pseudobulk_benchmark.R
Rscript code/analysis/09_multisolver_benchmark.R
Rscript code/analysis/10_deconvolve_and_composition.R
Rscript code/analysis/E1_mrna_correction.R
Rscript code/analysis/E2_hierarchical_decollapse.R
Rscript code/analysis/E3_anchor.R
Rscript code/analysis/E4_multimethod.R
Rscript code/analysis/E6_spatial_perspot.R
Rscript code/analysis/E7_flip_audit.R
Rscript code/analysis/E9_detection_limit.R
Rscript code/analysis/14_methylation_orthogonal.R
Rscript code/analysis/18_spatial_validation.R
Rscript code/analysis/16_synthesis.R
```

Random seeds are fixed in the scripts that subsample or simulate (seed 0).

## Contact

Hailing Cui, Department of Stomatology, Characteristic Medical Center of the Chinese People's Armed Police Force, Tianjin, China. 13752561486@163.com
