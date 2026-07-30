# combatrefql 0.1.0

No `0.1.0` Git tag or GitHub release existed on 2026-07-30, so this follow-up
retains version 0.1.0 as a pre-release correction.

- Fixed within-batch reference designs to include an intercept and audit every fallback.
- Made automatic-reference tie resolution independent of factor-level order.
- Corrected QL fitting to use edgeR's v4 global common NB dispersion, never tagwise dispersion.
- Removed fabricated reference frequency and added honest selection fields.
- Derived convergence from edgeR failure/numerical checks; unavailable states are `NA`.
- Preserved integer-valued mapped counts as doubles above the 32-bit integer range.
- Added centralized, auditable lossless input repair and opt-in lossy policies.
- Added normal/verbose/quiet CLI reports, JSON diagnostics, gene status, logs,
  validation-only mode, flexible delimited files, timings, and stable exit codes.
