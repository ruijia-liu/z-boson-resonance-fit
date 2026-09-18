# Provenance and preparation

The original project is Ruijia Liu's 2025 PHYS10362 Z-boson coursework analysis. The owner has confirmed independent completion of the original project in the local project inventory. Original work included CSV handling, data selection, nonlinear resonance fitting, covariance errors, a lifetime calculation and chi-square plotting.

Two local submission variants were compared. The more descriptive variant with named exception handling was used as the behavioral reference. Both original scripts and old figures remain unchanged in the local archive; neither is distributed as a second runnable implementation here.

AI-assisted portfolio preparation on 17 September 2026 produced the refactored `analysis.py`, documentation and tests. This work:

- Preserved the original E²-numerator line shape and default input windows/one-pass residual cut.
- Replaced the ambiguously normalized partial-width amplitude with cross section at the mass parameter in nb.
- Corrected the extra c² in the lifetime conversion and propagated width uncertainty.
- Added positive fit bounds, explicit failures, row-level exclusion auditing and pre-cut fit results.
- Replaced fixed-amplitude contours and a grid minimum with profiled-amplitude contours relative to the fitted minimum.
- Added reproducible command-line paths, noninteractive plots, dependency records and physical/numerical checks.

These preparation changes should not be represented as features of the original submission. Student account identifiers and unrelated teaching materials are excluded. The supplied data are preserved byte-for-byte, with hashes in `results/fit_results.json`. Their original experimental/simulation source and redistribution license were not supplied; this repository does not assign a license to them or claim they are a particular experiment's measurements.
