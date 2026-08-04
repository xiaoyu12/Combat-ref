# Validation report — 0.2.0

Unit tests cover the new API, S7 validation, typed failures, dispersion
separation, deterministic transport, zero preservation, and chunk equivalence.
The reproducible smoke simulation in `validation/simulation.R` checks dimensions,
integer-valued double output, and a fixed reference on synthetic NB counts.

This is an implementation smoke check, not broad statistical validation.
External real-data validation, predefined multi-scenario simulations, and
operating-characteristic assessment remain required. No production-readiness,
superiority, FDR, power, or biological-signal preservation claim is made.

On R 4.6.1 with seed 104729, 200 genes, 12 samples, and two batches, the raw
median B/A mean ratio was 1.7112 and the adjusted ratio was 1.6667; zero genes
were unsupported. This modest movement toward one is not evidence of general
correction quality.
