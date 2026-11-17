"""Replace splines in a GENIE spline file with those of the same name from
other GENIE spline files (e.g. the smooth QE splines made by
construct_representative_spline.py).

    python replace_splines.py -b gxspl-NUbig.xml.gz -r smooth1.xml smooth2.xml ... -o new.xml.gz
    python replace_splines.py -b gxspl-NUbig.xml.gz -l smooth.list -o new.xml.gz

Splines are matched by their full name only.  For every match, the base
spline's knots are replaced by the revision spline's knots, and its nknots
attribute is set to the new number of knots (GENIE reads exactly nknots
knots).  Everything else in the base file is copied unchanged.

The base file is read line by line (it can be several GB, .xml or .xml.gz),
so memory use stays small; only the revision files are held in memory.  It
assumes GENIE's own layout: each <spline ...> start tag on its own line, and
</spline> on its own line.  The output is gzipped if its name ends in .gz.

A summary is printed: how many revision splines, how many were replaced,
how many weren't found, and how many changed their number of knots.  The
exit status is 1, and (unless --allow-missing) no output is written, if a
revision spline has no spline of that name in the base file, has no knots,
or is in two revision files, or if a revision file's tune differs from the
base file's (unless --ignore-tune).  With --allow-missing, splines missing
from the base file are only reported, and the status is 0 if nothing else
is wrong.

Only the Python standard library is used.
"""
import argparse
import gzip
import os
import re
import sys
import xml.etree.ElementTree as ET

NAME_RE = re.compile(r'name="([^"]*)"')
NKNOTS_RE = re.compile(r'nknots="[0-9]*"')
TUNE_RE = re.compile(r'<genie_tune\s+name="([^"]*)"')


def open_text(path, mode, gz=None):
    """open a plain or gzipped file as text (gzipped if gz, or by default if
    the name ends in .gz)"""
    if gz is None:
        gz = path.endswith(".gz")
    if gz:
        return gzip.open(path, mode + "t", encoding="ISO-8859-1")
    return open(path, mode, encoding="ISO-8859-1")


def read_revisions(files):
    """{spline name: (file, tune, [(E, xsec), ...])} from the revision files,
    plus a list of problems"""
    revs, problems = {}, []
    for f in files:
        with open_text(f, "r") as fin:
            root = ET.parse(fin).getroot()
        for tune in root.iter("genie_tune"):
            tname = tune.get("name")
            for sp in tune.iter("spline"):
                name = sp.get("name")
                knots = [((k.findtext("E") or "").strip(), (k.findtext("xsec") or "").strip())
                         for k in sp.findall("knot")]
                if name in revs:
                    problems.append(f"{name} is in {revs[name][0]} and {f}")
                if not knots or any(e == "" or x == "" for e, x in knots):
                    problems.append(f"{name} (in {f}) has no knots, or a knot without E or xsec")
                revs[name] = (f, tname, knots)
    return revs, problems


def main():
    p = argparse.ArgumentParser(
        description="Replace splines in a GENIE spline file with those of the same name from other "
                    "GENIE spline files (matched by name only). See the module docstring for details.")
    p.add_argument("-b", "--base", required=True, help="base file whose splines are replaced, e.g. the big file (.xml or .xml.gz)")
    p.add_argument("-r", "--revision", nargs="+", default=[], help="revision file(s) with the new splines")
    p.add_argument("-l", "--list", help="file listing revision files, one per line (added to -r)")
    p.add_argument("-o", "--outpath", default="./output.xml", help="output file (gzipped if it ends in .gz)")
    p.add_argument("--allow-missing", action="store_true",
                   help="write the output even if some revision splines aren't in the base file")
    p.add_argument("--ignore-tune", action="store_true",
                   help="don't require the revision files' tune to match the base file's")
    args = p.parse_args()

    files = list(args.revision)
    if args.list:
        with open(args.list) as fl:
            files += [l.strip() for l in fl if l.strip() and not l.lstrip().startswith("#")]
    if not files:
        p.error("no revision files: give -r and/or -l")
    if os.path.abspath(args.outpath) == os.path.abspath(args.base):
        p.error("the output file must not be the base file")

    revs, problems = read_revisions(files)
    print(f"{len(revs)} spline(s) in {len(files)} revision file(s)")

    found = {}           # name -> times found in the base file
    nk_changed = []      # (name, old nknots, new nknots)
    base_tune = None
    tmp = args.outpath + ".tmp"
    with open_text(args.base, "r") as fin, open_text(tmp, "w", gz=args.outpath.endswith(".gz")) as fout:
        skip, nbase, cur = False, 0, None
        for line in fin:
            if skip:
                if "</spline>" in line:
                    skip = False
                    if nbase != len(revs[cur][2]):
                        nk_changed.append((cur, nbase, len(revs[cur][2])))
                    fout.write(line)
                elif "<knot>" in line:
                    nbase += 1
                continue
            if base_tune is None and "<genie_tune" in line:
                m = TUNE_RE.search(line)
                base_tune = m.group(1) if m else ""
            if "<spline " in line:
                m = NAME_RE.search(line)
                if m and m.group(1) in revs:
                    cur = m.group(1)
                    knots = revs[cur][2]
                    fout.write(NKNOTS_RE.sub(f'nknots="{len(knots)}"', line))
                    for e, x in knots:
                        fout.write(f"\t<knot> <E> {e:>10} </E> <xsec> {x} </xsec> </knot>\n")
                    found[cur] = found.get(cur, 0) + 1
                    skip, nbase = True, 0
                    continue
            fout.write(line)
        if skip:
            problems.append(f"the base file ends inside spline {cur}")

    if not args.ignore_tune:
        for f in files:
            tunes = {t for (rf, t, _) in revs.values() if rf == f}
            for t in sorted(tunes - {base_tune}):
                problems.append(f"{f} is for tune {t}, the base file for {base_tune}")
    missing = [n for n in revs if n not in found]
    for n in missing[:5]:
        problems.append(f"{n} (in {revs[n][0]}) is not in the base file")
    if len(missing) > 5:
        problems.append(f"... and {len(missing) - 5} more not in the base file")
    for n, c in found.items():
        if c > 1:
            print(f"note: {n} appears {c} times in the base file (all replaced)")
    for n, a, b in nk_changed[:5]:
        print(f"note: {n}: {a} knots replaced by {b}")

    print(f"{sum(found.values())} replaced in the base file, {len(missing)} not found, "
          f"{len(nk_changed)} with a different number of knots")
    for pr in problems:
        print(f"problem: {pr}", file=sys.stderr)

    only_missing = all(pr.endswith("is not in the base file") or "more not in the base file" in pr
                       for pr in problems)
    if problems and not (args.allow_missing and only_missing):
        os.remove(tmp)
        print(f"no output written ({len(problems)} problem(s))", file=sys.stderr)
        return 1
    os.replace(tmp, args.outpath)
    print("Saved output XML to", args.outpath)
    return 0


if __name__ == "__main__":
    sys.exit(main())
