# BreastSubtypeR 1.1.2

## Highlights
- Our manuscript describing **BreastSubtypeR** has been accepted in *NAR Genomics and Bioinformatics* (2025) and selected as **Editor’s Choice**.  
  Citation: Yang Q, Hartman J, Sifakis EG. *BreastSubtypeR: A Unified R/Bioconductor Package for Intrinsic Molecular Subtyping in Breast Cancer Research.* NAR Genomics and Bioinformatics, 2025.  

## Bug Fixes and Maintenance
This patch release addresses minor issues identified after the initial 1.0.0 release while preserving all core functionality. Key fixes include:
- Improved robustness when handling different input data types, ensuring accurate subtyping for both normalized (e.g., TPM, FPKM) and unnormalized (e.g., raw count) gene expression matrices
- Enhanced handling of cohorts with extreme ER+ ratios in AUTO mode
- Updated documentation typos and parameter descriptions

## New Features and Bug Fixes
This release enhances stability, adds new input support, and improves documentation. Key updates include:
- **New:** Support for raw RNA-seq counts (requires user-supplied gene lengths)  
- Refined, data-driven thresholds for ER+ skew detection in AUTO mode  
- Fixed input validation edge cases in the subtyping pipeline  
- Minor documentation updates and typo corrections  

## Core Features (from v1.0.0, unchanged but highlighted for reference)
- **Comprehensive Intrinsic Subtyping:** Integrates multiple published algorithms (NC- and SSP-based), including PAM50 variants, AIMS, ssBC, sspbc, and others  
- **Unified Multi-Method Interface (`BS_Multi`):** Execute multiple classifiers through a consistent API and compare results side-by-side  
- **AUTO Mode (cohort-aware selection):** Evaluates cohort composition (e.g., receptor-status distribution, subtype purity, subgroup sizes) and automatically disables classifiers whose assumptions are violated, improving accuracy in skewed or small cohorts  
- **Standardised Input & Method-Specific Normalisation:** Supports FPKM and log₂-normalised microarray/nCounter data with method-appropriate transformations  
- **Optimised Probe/Gene Mapping:** Entrez ID–based mapping with conflict resolution to maximise marker coverage across platforms  
- **Interactive Shiny App (`iBreastSubtypeR`):** Local GUI for non-programmers that replicates core workflows while preserving data privacy  
- **Bioconductor Distribution & Reproducibility:** Unit tests, vignettes, and SummarizedExperiment compatibility for reproducible deployment  

## Upgrade Note
- Raw RNA-seq counts are supported **only from version 1.1.2 onwards**.  
- Users working with raw counts should ensure they are running **v1.1.2 or later**, as earlier versions do not implement this functionality.