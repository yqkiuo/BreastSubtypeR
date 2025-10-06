# BreastSubtypeR 1.1.5

## Highlights
- Paper published in *NAR Genomics and Bioinformatics* (2025), **Editor’s Choice** (DOI: 10.1093/nargab/lqaf131).
- Support for raw RNA-seq counts (requires gene lengths).
- **iBreastSubtypeR** refresh: cleaner UX, smarter AUTO guidance, consistent exports.
- Refined, data-driven thresholds for ER/HER2 skew detection in AUTO.
- Broader input-validation across the subtyping pipeline.

## Bug Fixes
- PAM50 variants no longer error with `calibration = "None"` or `"External"`.
- ssBC variants handle datasets with <50 PAM50 genes more robustly.
- Fixed `data.frame` issue (“`check.names` matched by multiple arguments”) via safe builders.
- Removed duplicate “Subtype == BS_*class” columns in downloads.
- Safer ROR merges on `PatientID` with clearer notifications.
- Eliminated `jsonlite` named-vector warning in plotting.

## Shiny (iBreastSubtypeR)
- New hero card + method chips (NC, SSP, ROR, AUTO) and “Why AUTO?” explainer modal.
- **AIMS:** 4-class toggle disabled (AIMS is 5-class only).
- **AUTO/AIMS/SSPBC:** Full-metrics export selector hidden.
- Clear per-method help + PAM50 calibration notes.
- Small UX polish: file dialogs + Mapping preserve scroll position.

**Exports**
- Two modes: **Calls only** and **Full metrics** (incl. ROR when available).
- Standard column names: **`Call_5class`** / **`Call_4class`** (internal `BS` / `BS.Subtype` mapped).
- **AUTO:** Calls = selected-k table (+ entropy); Full = selected-k merged with other-k (suffix `_4class` / `_5class`).
- Filenames encode method, class, and mode (e.g., `results-PAM50-5class-full.txt`).

**Behavior & Validation**
- **ROR** only for NC methods when `TSIZE`/`NODE` exist and are numeric (auto coercion + warnings).
- Stricter checks for **ssBC** subgroups (ER/HER2/TN presence and coding).
- **AUTO** preflight warns on missing or mis-coded ER/HER2.

## Core Features (unchanged)
- Unified NC (PAM50, cIHC, PCAPAM50, ssBC) and SSP (AIMS, SSPBC) subtyping under one API.
- `BS_Multi` with cohort-aware **AUTO** selection; `Mapping` for platform-agnostic preprocessing.
- Interactive Shiny app with Bioconductor-friendly exports.

## Upgrade Notes
- Raw RNA-seq counts are supported **from v1.1.3 onward** (requires gene lengths).
- If you previously parsed `BS` / `BS.Subtype`, switch to **`Call_5class` / `Call_4class`**.
- Package API unchanged.
