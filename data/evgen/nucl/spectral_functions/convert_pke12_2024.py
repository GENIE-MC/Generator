#!/usr/bin/env python3
"""Convert pke12_2024.table.origin -> pke12_2024.table in GENIE's SpectralFunc format.

Source: Ankowski, Benhar & Sakuda (2024), "Determination of the proton spectral
function of 12C from (e,e'p) data" (suppl. material). P(|k|,E) in MeV^-4,
normalised to Z=6.

Why this conversion is needed
-----------------------------
GENIE's genie::SpectralFunc::LoadSFDataFile builds a *uniform* TH2D and reads a
3-line header (num_E num_p / E_min p_min / E_max p_max), like pke12_tot.data.
The 2024 table has a 2-line header and a NON-uniform energy grid:
  - fine   : 340  bins x 0.025 MeV over 13   < E < 21.5 MeV  (NIKHEF (e,e'p) fit)
  - coarse : 2785 bins x 0.1   MeV over 21.5 < E < 300  MeV  (Benhar model)
GENIE cannot read that directly. We resample energy onto a single UNIFORM grid.

Conversion (LOSSLESS, BIN = 0.025 MeV)
--------------------------------------
  - fine bins kept verbatim (already 0.025 MeV);
  - each coarse 0.1 MeV bin is SPLIT into 4 sub-bins of 0.025 MeV with the SAME
    density value (replicate, not average). Exact, since density is constant
    across the bin: density*0.1 == 4 * density*0.025, so int P dE is preserved.
  - momentum grid (40 bins x 20 MeV, centers 10..790) is already uniform -> kept.
Result: 11480 uniform 0.025-MeV energy bins over 13-300 MeV, centers 13.0125..299.9875.

P values are written unchanged: the table is already in GENIE's "N*P" convention
(int 4 pi k^2 P dk dE = Z = 6); GENIE divides by targetN=Z at read time.
"""
from pathlib import Path
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
IN_FILE  = SCRIPT_DIR / "pke12_2024.table.origin"
OUT_FILE = SCRIPT_DIR / "pke12_2024.table"

BIN = 0.025        # target uniform energy bin width [MeV]
SPLIT = 4          # 0.1 MeV coarse bin -> 4 x 0.025 MeV sub-bins
Z = 6              # C12 protons (file normalised to Z); used only for the self-test


def load_origin(path):
    """Parse the 2024 (non-uniform) table. Returns k, P[ik,iE], n1 (fine count)."""
    tok = np.fromstring(path.read_text(), sep=" ")
    n_k, dk = int(tok[0]), tok[1]
    n1, d1, n2, d2 = int(tok[2]), tok[3], int(tok[4]), tok[5]
    n_E = n1 + n2
    block = 1 + 2 * n_E
    body = tok[6:6 + n_k * block].reshape(n_k, block)
    k = body[:, 0]                  # MeV/c, centers 10..790
    E = body[0, 1::2]               # MeV (same grid every block)
    P = body[:, 2::2]               # MeV^-4 density (N*P convention)
    assert dk == 20.0 and d1 == 0.025 and d2 == 0.1, (dk, d1, d2)
    assert n_k == 40 and n1 == 340 and n2 == 2785, (n_k, n1, n2)
    return k, dk, E, P, n1


def resample_uniform(P, n1):
    """Fine part verbatim; each coarse bin split into SPLIT identical sub-bins."""
    P_fine = P[:, :n1]                              # [40, 340]
    P_coarse = P[:, n1:]                            # [40, 2785]
    P_new = np.concatenate([P_fine, np.repeat(P_coarse, SPLIT, axis=1)], axis=1)
    return P_new                                    # [40, 11480]


def write_genie(path, k, P_new, e_lo=13.0, e_hi=300.0, p_lo=0.0, p_hi=800.0):
    n_k, n_E = P_new.shape
    E_new = e_lo + (np.arange(n_E) + 0.5) * BIN     # bin centers 13.0125..299.9875
    lines = [f"{n_E} {n_k}", f"{int(e_lo)} {int(p_lo)}", f"{int(e_hi)} {int(p_hi)}"]
    for ik in range(n_k):
        lines.append(f"{k[ik]:.5E}")
        row = P_new[ik]
        # 4 (E, P) pairs per line, matching pke12_tot.data layout
        buf = []
        for ie in range(n_E):
            buf.append(f"{E_new[ie]:9.4f} {row[ie]:14.5E}")
            if len(buf) == 4:
                lines.append(" " + "  ".join(buf))
                buf = []
        if buf:
            lines.append(" " + "  ".join(buf))
    path.write_text("\n".join(lines) + "\n")
    return E_new


def selftest(out_path, k_src, P_src, n1_src):
    """Re-parse the OUTPUT with GENIE's uniform-grid reader logic and verify."""
    tok = np.fromstring(out_path.read_text(), sep=" ")
    num_E, num_p = int(tok[0]), int(tok[1])
    E_min, p_min, E_max, p_max = tok[2], tok[3], tok[4], tok[5]
    body = tok[6:6 + num_p * (1 + 2 * num_E)].reshape(num_p, 1 + 2 * num_E)
    k = body[:, 0]
    pairs = body[:, 1:].reshape(num_p, num_E, 2)
    E = pairs[0, :, 0]
    P = pairs[:, :, 1] / Z                          # GENIE divides by targetN=Z
    dk = (p_max - p_min) / num_p                     # 20 MeV
    dE = (E_max - E_min) / num_E                     # uniform 0.025 MeV

    w = 4.0 * np.pi * (k[:, None] ** 2) * P
    norm = float((w * dk * dE).sum())
    f_E = (w * dk).sum(axis=0)
    n_k = (w * dE).sum(axis=1)

    # ORIGINAL (non-uniform) table, per proton: norm + n(k) marginal
    dE_src = np.concatenate([np.full(n1_src, 0.025), np.full(P_src.shape[1] - n1_src, 0.1)])
    w_src = 4.0 * np.pi * (k_src[:, None] ** 2) * (P_src / Z)
    norm_src = float((w_src * dk * dE_src).sum())    # dk identical (20 MeV)
    nk_src = (w_src * dE_src).sum(axis=1)            # n(k) is invariant under the rebinning

    print(f"  output grid : {num_E} E-bins x {num_p} k-bins  "
          f"(E {E_min:.0f}-{E_max:.0f}, dE={dE:.4f} MeV ; k {p_min:.0f}-{p_max:.0f}, dk={dk:.0f})")
    print(f"  E centers   : {E[0]:.4f} .. {E[-1]:.4f}  (expect 13.0125 .. 299.9875)")
    print(f"  source norm   int 4pi k^2 (P/Z) dk dE = {norm_src:.6f}   (Ankowski et al. Z=6 precision)")
    print(f"  converted norm                        = {norm:.6f}")
    print(f"  norm preserved to {abs(norm - norm_src):.2e}   (lossless conversion)")
    print(f"  f(E) peak @ {E[f_E.argmax()]:.4f} MeV   <E> = {np.average(E, weights=f_E):.2f} MeV  (uniform grid)")
    print(f"  n(k) peak @ {k[n_k.argmax()]:.0f} MeV/c   n(k) max rel. diff vs source = "
          f"{np.max(np.abs(n_k - nk_src) / np.maximum(nk_src, 1e-30)):.2e}")
    assert abs(norm - norm_src) < 1e-9, (norm, norm_src)   # conversion is exact
    assert num_E == 11480 and num_p == 40
    print("  SELF-TEST PASSED  (uniform 0.025 MeV grid round-trips through GENIE's reader)")


def main():
    print(f"reading {IN_FILE.name}")
    k, dk, E, P, n1 = load_origin(IN_FILE)
    P_new = resample_uniform(P, n1)
    print(f"resampled to uniform {BIN} MeV: {P.shape[1]} -> {P_new.shape[1]} energy bins "
          f"({n1} fine verbatim + {P.shape[1]-n1}x{SPLIT} split)")
    E_new = write_genie(OUT_FILE, k, P_new)
    print(f"wrote {OUT_FILE.name}  ({OUT_FILE.stat().st_size/1e6:.1f} MB)")
    selftest(OUT_FILE, k, P, n1)


if __name__ == "__main__":
    main()
