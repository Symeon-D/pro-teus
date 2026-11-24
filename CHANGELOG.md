# Changelog

Everything below was added **after** commit `f522158`
("Inclusion of sva, multiple protein plots, protein search bars, fixes to volcano plots").

## Data tab
- **Perseus-style imports:** a `#!{Type}` annotation row is now read correctly (the "T" bug).
- **Sample exclusion** (Sample info tab): multi-select, *Remove from analysis* / *Restore* / *Restore all*.
  Removed samples are ignored by every tab; the filter, normalization, imputation and results reset; a warning appears if
  fewer than 2 groups remain; the status bar shows the number removed.
- **Outlier detection (rMd, pmartR):** new sub-tab with metric selection, p-value cutoff (numeric box, default 0.05),
  labelled rMd plot, flagged-samples table, per-sample boxplots of every metric (flagged sample marked with an X, one plot per row)
  and a *Remove detected outliers (n)* button. Legend title is "Group".
- **Venn of detected proteins:** region extraction (unique / common / discordant proteins) and a 3-sheet Excel download.
  The Venn functions no longer write PNG files to a temp folder unless a file name is given.
- Uploading new files no longer crashes the Sample-info tab (stale `group_vars`).

## Normalization
- Normalization is attempted first; if it cannot run the status bar reports *"Normalization skipped: too many missing values"*.

## Batch correction (SVA)
- The Group variable is included by default.

## Differential expression
- Volcano: Y axis defaults to P.Value; when fewer than 5 proteins are significant only those are labelled;
  fixed the endless refresh between label counts (label changes are debounced with an echo guard).
- "Differentially expressed proteins per contrast": a zero on the down side is drawn on the left.
- arrayWeights: *Show weights* button in the Volcano "Show data" style (Group / Sample / Weight rows, scrolls sideways).
- Protein plots: *Free y-axis* is off by default.

## Heatmap (largely rewritten)
- **Directions:** normal, sideways (samples as rows) and circular. **Shapes:** rectangles, square boxes, circles
  (circle drawing and linear 20-100 % size-by-value follow Morpheus.js).
- **Auto-fit:** every direction/shape is fitted into the visible window (no scrolling). A bigger box size or typed
  width/height makes the plot scrollable. Fonts (sample, group, gene), box size, figure size (cm) and Z-bar height show their
  auto-fitted value; each has a **reset** button that is only enabled when the value was changed.
- Size inputs apply at once on a single change and are throttled (0.5 s) while dragging.
- **Group boxes** are drawn with `decorate_annotation`, so they sit exactly on their samples; sample names can be rotated
  (squares); without sample names the colour boxes match the group boxes.
- **Legends:** separate Z-score bar and group legend (*Show group legend*), both centred on the heatmap body; bar height capped to the canvas.
- **Gene names:** position drop-down (normal: right/left, sideways: down/up); font auto-fits.
- **Circular:** one ring per sample (or per group when merged), no sample names, group names in the opening, small spacers
  between groups, opening size / start angle controls, auto thinning of rings when circlize runs out of space,
  Z-score bar right next to the circle. Options that do not apply are greyed out and unticked.
- Top N blank = all proteins, in every direction. Cluster column in the processed matrix only when rows are clustered.
- **Processed matrix:** scrollable DT table with CSV download, like the other tables.
- **Open in Morpheus:** exports the matrix as GCT 1.3 to a secret GitHub Gist (needs `GITHUB_PAT` with `gist` scope) and opens Morpheus;
  fallback *Download GCT*.

## Over-representation analysis
- Heatmap clusters are only offered as gene sets when a cluster column exists.

## Dependencies
- Added `pmartR` (`install_packages.R`).
