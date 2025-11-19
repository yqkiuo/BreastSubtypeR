# BreastSubtypeR 1.3.1

## Highlights (from v1.1.3 onward)

-   Paper published in *NAR Genomics and Bioinformatics* (2025), **Editor’s Choice** (DOI: 10.1093/nargab/lqaf131).
-   Support for raw RNA-seq counts (requires gene lengths).
-   **iBreastSubtypeR** refresh: cleaner UX, smarter AUTO guidance, consistent exports.

### Enhancements

-   **ssBC/ssBC.v2: singleton subgroup robustness:** Subgroups with `n=1` no longer error:
    -   Keeps matrix shape (`drop=FALSE`) and hardens dimnames/types.
    -   Primary path: original `sspPredict()`.\
        Fallback: nearest-centroid (Spearman) when needed.
    -   If there are **0 common PAM50 genes**, returns `NA` labels with shaped `distances`/`dist.RORSubtype` to avoid downstream errors.
    -   ROR computation guarded for incomplete inputs.
-   **SSPBC output now “full”:** `BS_sspbc()` and Shiny “sspbc” runs return a full metrics table (not calls-only).
    -   Exports map core label columns to the standard names (`Call_5class` / `Call_4class` when applicable).
-   **Shiny: “Load example data…” button**
    -   One-click load of a small demo dataset from `inst/RshinyTest/` to explore the UI without uploads.\
    -   Shows a notification on success; users can immediately run **Preprocess & map** and analyses.
-   **AUTO preflight UI (Shiny):** Now detects cohort kind (`TN`, `ER/HER2`, `ER`-only, `HER2`-only) and shows compact stats:
    -   ER/HER2 subgroups: **ER+/HER2−**, **ER−/HER2−**, **ER+/HER2+**, **ER−/HER2+**\
    -   TN cohorts: **TN** vs **nonTN**\
    -   Readiness uses the same minimums used by AUTO (sourced programmatically; no duplicated thresholds).
-   **Shorter notifications.**\
    -Routine toasts (e.g., “Step 1 complete. Proceed to Step 2.”) now auto-dismiss sooner to reduce UI clutter.

### Bug fixes

-   **TN cohorts + ssBC**: `BS_Multi()` now respects TN cohorts when methods are specified manually; `ssBC`/`ssBC.v2` switch to `s = "TN"` / `"TN.v2"` when a `TN` column indicates a TN cohort. Falls back to `s = "ER"` / `"ER.v2"` otherwise.
-   **AUTO internals**: fixed variable name typo (`samples_ERHER2.icd`).
-   **Mapping():** Robust ENTREZID coercion (from `as.character()` to `as.integer()` with suppressed warnings).

### Shiny

-   **Surface method warnings as toasts:**
    -   Runs are wrapped in a warning handler; package warnings (e.g., ssBC.v2 singleton fallbacks) appear as yellow notifications.\
    -   Warnings include subgroup, `n`, and example sample IDs for quick triage.
-   **Shiny preflight reset**:
    -   Fixed a stale cohort summary after switching data sources (manual uploads ↔ example). The preflight panel now revalidates once inputs change.

### Developer notes

-   Added lightweight internal logger `._msg()` and replaced scattered `message()` calls in AUTO to standardize package output without affecting CRAN/Bioc checks.

### Documentation

-   README/vignette: brief note on the example-data button and expected file locations.

### Compatibility Notes

-   SSPBC “full” output keeps previous columns for calls; additional metrics may appear.

## Upgrade Notes

-   Raw RNA-seq counts are supported **from v1.1.3 onward** (requires gene lengths).
-   If you previously parsed `BS` / `BS.Subtype`, switch to **`Call_5class` / `Call_4class`**.
-   Package API unchanged.
