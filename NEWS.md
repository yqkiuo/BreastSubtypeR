# BreastSubtypeR 1.1.1

## Bug Fixes and Maintenance
This patch release addresses minor issues identified after the initial 1.0.0 release while preserving all core functionality. Key fixes include:
- Resolved input validation edge cases in the subtyping pipeline
- Enhanced handling of cohorts with extreme ER+ ratios in AUTO mode
- Corrected data type handling for raw counts input
- Updated documentation typos and parameter descriptions

All core features remain fully functional:
- **Comprehensive Intrinsic Subtyping for Breast Cancer**: Integrates multiple published intrinsic subtyping methods, including NC-based approaches like the original PAM50 (Parker et al., J Clin Oncol, 2009) and SSP-based methods like AIMS (Paquet et al., J Natl Cancer Inst, 2015).
- **Multi-Method Subtyping Functionality**: Simultaneously predicts breast cancer intrinsic subtypes using a variety of validated methods for comparative analysis.
- **AUTO Mode Feature**: Evaluates the distribution of ER and HER2 status in the test cohort to automatically select subtyping methods that align with the cohort's characteristics, ensuring compatibility with method-specific assumptions for greater accuracy and reliability.
- **Optimized Gene Mapping**: Optimizes gene mapping using Entrez IDs to maximize the inclusion of genes across subtyping methods.
- **Streamlined Input and Output**: Provides standardized input/output formats to ensure smooth integration with other gene expression analysis tools.
- **User-Friendly Shiny App Interface**: Features a web-based GUI that runs entirely locally, ensuring data privacy with no online sharing, ideal for users who prefer a visual interface over R scripting.  

Users can safely upgrade from v1.0.0 for improved stability.