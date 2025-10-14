# BreastSubtypeR 1.1.9

## Highlights (from v1.1.3 onward)
- Paper published in *NAR Genomics and Bioinformatics* (2025), **Editor’s Choice** (DOI: 10.1093/nargab/lqaf131).
- Support for raw RNA-seq counts (requires gene lengths).
- **iBreastSubtypeR** refresh: cleaner UX, smarter AUTO guidance, consistent exports.

### Enhancements
- **SSPBC output now “full”:** `BS_sspbc()` and Shiny “sspbc” runs return a full metrics table (not calls-only).  
  - Exports map core label columns to the standard names (`Call_5class` / `Call_4class` when applicable).

- **Shiny: “Load example data…” button**  
  - One-click load of a small demo dataset from `inst/RshinyTest/` to explore the UI without uploads.  
  - Shows a notification on success; users can immediately run **Preprocess & map** and analyses.

### Documentation
- README/vignette: brief note on the example-data button and expected file locations.

### Compatibility Notes
- SSPBC “full” output keeps previous columns for calls; additional metrics may appear.  

## Upgrade Notes
- Raw RNA-seq counts are supported **from v1.1.3 onward** (requires gene lengths).
- If you previously parsed `BS` / `BS.Subtype`, switch to **`Call_5class` / `Call_4class`**.
- Package API unchanged.
