#!/usr/bin/env python3
"""Convert a Benhar spectral-function table (``benhar-sf-*.data``) into the
``pke*_tot.data`` format read by GENIE's SpectralFunc nuclear model.

The two layouts hold the *same* P(|k|, E) data:

  benhar-sf-*.data : comment header (``#``) then one ``k  E  P`` triplet per
                     line (k, E in MeV; P the probability density).

  pke*_tot.data    : the format consumed by
                     ``src/Physics/NuclearState/SpectralFunc.cxx``::

                         num_E_bins  num_p_bins
                         E_min       p_min
                         E_max       p_max
                         <p_center_1>
                          E1 P1   E2 P2   E3 P3   E4 P4      (num_E_bins pairs,
                          ...                                  4 per line)
                         <p_center_2>
                          ...

                     The tabulated k and E values are bin *centers*; the header
                     gives the bin-edge ranges (centers +/- half a bin width).
                     Numbers use Fortran ``0.dddE+NN`` notation.

The conversion is lossless to the precision Benhar's tables carry (3 significant
figures), which round-trips exactly through the 3-decimal pke mantissa.

Usage:
    python convert_benhar_to_pke.py benhar-sf-56fe.data pke56_tot.data
    python convert_benhar_to_pke.py benhar-sf-56fe.data        # -> pke<A>_tot.data
"""
import argparse
import math
import re
import sys


def fort(value, decimals):
    """Format ``value`` as Fortran ``0.d..dE+NN`` (mantissa in [0.1, 1))."""
    if value == 0:
        return "0." + "0" * decimals + "E+00"
    exp = math.floor(math.log10(abs(value))) + 1
    mant = value / 10.0 ** exp
    s = f"{mant:.{decimals}f}"
    if abs(float(s)) >= 1.0:           # rounding pushed the mantissa up to 1.0
        exp += 1
        mant = value / 10.0 ** exp
        s = f"{mant:.{decimals}f}"
    return f"{s}E{exp:+03d}"


def fnum(x):
    """Render a header bound: integers without a trailing ``.0``."""
    return str(int(x)) if x == int(x) else repr(x)


def read_benhar(path):
    """Return ((sorted k centers), (sorted E centers), {(k, E): P})."""
    table = {}
    ks, es = [], []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        k, e, p = (float(x) for x in line.split())
        table[(k, e)] = p
        if k not in ks:
            ks.append(k)
        if e not in es:
            es.append(e)
    return sorted(ks), sorted(es), table


def convert(src, dst):
    ks, es, table = read_benhar(src)
    n_e, n_k = len(es), len(ks)
    if n_e < 2 or n_k < 2:
        sys.exit(f"{src}: need >=2 k and E values to infer bin widths")

    w_e = es[1] - es[0]
    w_k = ks[1] - ks[0]
    e_min, e_max = es[0] - w_e / 2, es[-1] + w_e / 2
    p_min, p_max = ks[0] - w_k / 2, ks[-1] + w_k / 2

    lines = [
        f"{n_e} {n_k}",
        f"{fnum(e_min)} {fnum(p_min)}",
        f"{fnum(e_max)} {fnum(p_max)}",
    ]
    for ip, k in enumerate(ks):
        mom = fort(k, 5)
        lines.append(mom if ip == 0 else " " + mom)   # mirror pke12: 1st block flush-left
        buf = ""
        for ie, e in enumerate(es):
            buf += f"{e:6.1f}   {fort(table[(k, e)], 3)}"
            if (ie + 1) % 4 == 0:
                lines.append(buf)
                buf = ""
        if buf:
            lines.append(buf)

    with open(dst, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {dst}  ({n_e} E bins x {n_k} k bins = {n_e * n_k} points)")


def default_output(src):
    """benhar-sf-56fe.data -> pke56_tot.data (use the mass number A)."""
    m = re.search(r"(\d+)", src.rsplit("/", 1)[-1])
    return f"pke{m.group(1)}_tot.data" if m else "pke_tot.data"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="benhar-sf-*.data file to convert")
    ap.add_argument("output", nargs="?", help="output pke*_tot.data (default: pke<A>_tot.data)")
    args = ap.parse_args()
    convert(args.input, args.output or default_output(args.input))


if __name__ == "__main__":
    main()
