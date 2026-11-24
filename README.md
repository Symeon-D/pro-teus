# Pro-teus

An R Shiny app for label-free proteomics: from a protein table to QC, filtering, normalization, imputation,
batch correction, differential expression, heatmaps and enrichment analysis. Deployed on Posit Connect Cloud.

## Workflow (tabs)
1. **Data** - expression table (Perseus / generic) + sample sheet, Group assignment, sample exclusion, QC plots, rMd outlier detection, Venn of detected proteins.
2. **Normalization** - NormalyzerDE-based evaluation and normalization.
3. **Imputation** - several methods with evaluation.
4. **Batch / SVA** - surrogate variable analysis (Group included by default).
5. **Differential expression** - limma (duplicateCorrelation, arrayWeights, F-tests), volcano plots, per-contrast summaries, Venn of significant proteins, protein plots.
6. **Heatmap** - normal / sideways / circular, rectangles / squares / circles, auto-fit to the window, Morpheus export, processed matrix.
7. **ORA / GSEA** - g:Profiler, clusterProfiler, GeneCodis, fgsea, rrvgo.

See [CHANGELOG.md](CHANGELOG.md) for what changed since commit `f522158`.

## Run locally
```r
Rscript install_packages.R     # once (some packages take a while)
shiny::runApp()
```

## Layout
```
app.R                 entry point
parts/                UI + server pieces, one pair per tab (the folder must NOT be named R/)
functions/            analysis and plotting functions
data/                 demo data
install_packages.R    package list
make_manifest.R       regenerates manifest.json for Posit Connect
```

## Deployment (Posit Connect Cloud)
Regenerate `manifest.json` with `make_manifest.R` after changing packages, commit, and redeploy from GitHub.

## Optional configuration
- `GITHUB_PAT` (token with `gist` scope): required for *Open in Morpheus* and *Open in Ora-cle*. Without it for Morpheus use *Download GCT* and open the file in
  <https://software.broadinstitute.org/morpheus/> (File > Open) and for Ora-cle download the enrichment file and upload the file in <https://symeon-d.github.io/ora-cle/> (Upload CSV / TSV).

## Notes
- Heavy jobs (ORA, imputation, NormalyzerDE) run in a background R process (`callr`), so *Stop process* works.
