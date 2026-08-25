# INCL26_07a — spectral-function + Unified-CC QE + INCL FSI tune

Derived from `AR23_20i` (byte-identical copy, then reconfigured). Runnable as
tune id `INCL26_07a_00_000`.

## Physics content

 - **1p1h QE (CC):** `genie::HybridXSecAlgorithm/Unified-CC`, which delegates
   complex nuclei to the SF-Unified impulse-approximation model
   `genie::UnifiedQELPXSec/ZExp_lqcd` (LQCD Z-expansion axial form factor,
   `LwlynSmithFFCC/ZExp_lqcd`) and free nucleons to `LwlynSmithQELCCPXSec/Dipole`.
   Uses the standard `NewQELXSec` integrator — not the INCL-specific QE integrator.
 - **Ground state:** the 2024 Ankowski–Benhar–Sakuda spectral function for C12,
   `genie::SpectralFunc/pke12_2024` (data: `data/evgen/nucl/spectral_functions/pke12_2024.table`).
   Non-C12 targets fall back to the global `genie::LocalFGM/Default`.
 - **2p2h:** SuSAv2 MEC, retained unchanged from AR23_20i.
 - **FSI:** INCL intranuclear cascade, `genie::INCLCascadeIntranuke/Default`
   (`DeltaTransp-Enable` = true). Requires the INCL data environment from the
   `genie_inclxx` installation's `setup_env.sh`.
 - Pion-production parameters inherited from AR23_20i (originally G18_10a_02_11b).
 - **GENIE nuclear de-excitation photons OFF** (`NucDeExcitationSim`: `DoCarbon` = `DoArgon` = false,
   changed 2026-08-24; AR23_20i has both true — previous file kept as `NucDeExcitationSim.xml.bak-deex-on`).
   Nuclear de-excitation after the cascade is left entirely to INCL's ABLAXX.

## Files changed relative to AR23_20i

 - `ModelConfiguration.xml` — NuclearModel (C12 → SpectralFunc/pke12_2024),
   QEL-CC XSecModel (→ HybridXSecAlgorithm/Unified-CC), HadronTransp-Model
   (→ INCLCascadeIntranuke/Default) + DeltaTransp-Enable.
 - `EventGenerator.xml` — added (copied from SF26_22b). Overrides the QEL-CC
   thread to `genie::QELEventGeneratorINCL/Hybrid`, which prepares the QE vertex
   (struck nucleon positioned via the global `NucleusGenHybridStruck`) for hand-off
   to the INCL cascade. This is required by the INCL FSI choice — the master default
   `QELEventGeneratorSM/Default` does not feed the INCL cascade. `QELEventGeneratorINCL`
   and `NucleusGenHybridStruck` resolve from the global `$GENIE/config` (SF26_22b's
   local copies are byte-identical), so only `EventGenerator.xml` is needed locally.
   **2026-08-24:** the `QEL-CC-LAMBDA` thread (nubar CC quasi-elastic hyperon
   production, Lambda/Sigma0/Sigma-) now uses `genie::NucleusGenerator/Default`
   (= `NucleusGenHybridStruck`: INCL vertex + `INCLNucleus::reset` + `FermiMover/Default`,
   the same initial-state module as QEL-NC/RES/DIS) instead of the stock
   `VertexGenerator/Default` + `FermiMover/Default`. The stock pair never re-initialises
   the INCL nucleus, so every hyperon event ran the cascade on the previous event's stale
   nucleus and killed the job (`BaryonNumberConservation` exit(1), or a segfault on the
   dangling hit-nucleon pointer) — ~7.5e-4 per rockbox spill.
 - `TuneGeneratorList.xml` — dropped the two charm threads (`DIS-CC-CHARM`,
   `QEL-CC-CHARM`); `NGenerators` 18 -> 16, generators renumbered contiguously.
 - `NucDeExcitationSim.xml` — `DoCarbon`/`DoArgon` set to false (2026-08-24).
 - `CommonParam.xml`, `MECInteractionListGenerator.xml` — inherited from AR23_20i, unchanged.

## Naming note

FSI is INCL by explicit design (the `INCL26_` prefix). The `a` suffix here does
**not** follow the SF26_/LFG26_ series letter convention (where `a` = hA2018,
`b` = INCL FSI).
