---
tags: [genie, splines]
---

<!--
Making a PDF of this note (this comment doesn't show in Obsidian, on GitHub or in the PDF).

Needs pandoc and a TeX distribution with xelatex (e.g. MacTeX or TeX Live;
on macOS: brew install pandoc, plus MacTeX):

    pandoc Running_GENIE_Splines_at_Fermilab.md -o Running_GENIE_Splines_at_Fermilab.pdf \
        --pdf-engine=xelatex -V mainfont="STIX Two Text" -V monofont="Menlo" \
        -V geometry:margin=2cm -V colorlinks=true -V linkcolor=blue

xelatex and a font with Greek letters and Unicode superscripts are needed for
the ν, μ, τ and "10⁻³⁸ cm²" in the text: STIX Two Text has them (free, from
the STIX project). pandoc's default font lacks them, and Times New Roman lacks
the superscripts. STIX Two Text has no "→", so the text avoids it.
Obsidian's [[wiki links]] and > [!warning] callouts come out as plain text,
and very long code lines may run off the page.
-->

# Running GENIE Splines at Fermilab

How to run the GENIE spline production scripts, and how they are put together.

**Contents**

- [[#Overview]]
- [[#The scripts]]
- [[#Quick start]]
    - [[#1. Initialize]]
    - [[#2. Submit the DAG]]
    - [[#3. Watch progress]]
    - [[#4. Re-run failed pieces]]
    - [[#5. Check the splines]]
- [[#The stages]]
    - [[#Stages 1 and 3 `gxsp_mkspl`|Stages 1 and 3: gxsp_mkspl]]
    - [[#Stage 2 `gxsp_add_freenuc`|Stage 2: gxsp_add_freenuc]]
    - [[#Stage 4 `gxsp_add_isotopes`|Stage 4: gxsp_add_isotopes]]
    - [[#QE averaging (`--average-qe`)|QE averaging (stages 3q, 3s, 4q)]]
    - [[#Tunes where QE averaging does nothing]]
    - [[#Stage 5 `gxsp_reduce_graph`|Stage 5: gxsp_reduce_graph]]
    - [[#Stage 6 `gxsp_mkups`|Stage 6: gxsp_mkups]]
    - [[#Shortening the splines `gxsp_cut_emax`|Shortening the splines: gxsp_cut_emax]]
- [[#Design]]
    - [[#Language and portability]]
    - [[#Configuration travels in `cfg.tar.gz`|Configuration travels in cfg.tar.gz]]
    - [[#How GENIE finds XML files]]
    - [[#Start-up code pasted into each script]]
    - [[#Access to the staging area (PNFS)]]
    - [[#Mapping `$PROCESS` to a pair|Mapping $PROCESS to a pair]]
    - [[#Free nucleons come from stage 1]]
    - [[#Naming]]
    - [[#Combining files with gspladd]]
    - [[#NaN, inf and negative values]]
    - [[#The reduction awk script]]
    - [[#Job resources (jobsub options)]]
    - [[#Exit codes]]
    - [[#Environment variables]]
- [[#Other Considerations]]
    - [[#Scratch space when running interactively]]
    - [[#Running a stage by hand doesn't use ifdh]]
- [[#Adapting to another site]]
    - [[#What a site has to provide]]
    - [[#Setting up GENIE]]
    - [[#The staging area and file access]]
    - [[#Without CVMFS]]
    - [[#The DAG and the batch system]]
    - [[#The UPS tarball]]
- [[#Open items]]

## Overview

GENIE doesn't compute total cross sections while it generates events. It reads them from precomputed **splines**: for each probe, target and process, the cross section at a set of energy knots (`--knots`, up to `--emax`). A spline production computes these for every probe (the six neutrino flavors or a chosen subset, or electrons) on every isotope in `isotopes.cfg`, for one tune, and packages them as a UPS `genie_xsec` product. Computing them all in one job would take far too long, so the work is split into stages that run as a grid DAG, with many independent jobs where possible.

The scripts are tailored to Fermilab's batch production environment:
- jobs are submitted with `jobsub` as a DAG
- the staging area is in dCache (`/pnfs/...`), reached with `ifdh` from the grid nodes
- software comes from CVMFS, with UPS on SL7 nodes (in the SL7 container) and spack elsewhere
- the result is packaged as a UPS `genie_xsec` product

Running elsewhere would need those parts replaced.

The idea behind the stages:

1. **Free nucleons first (stage 1).** Each probe on the free proton and the free neutron. These are small jobs, and their results are needed by everything after.
2. **Combine them (stage 2)** into one file, `gxspl-freenuc.xml`.
3. **Every probe x isotope (stage 3).** One job per pair (750 with the default neutrino list), all in parallel. Each reads `gxspl-freenuc.xml` as input, so splines already computed there are loaded rather than recomputed. This is the bulk of the CPU time.
4. **Combine everything (stage 4)** into the big file, `gxspl-NUbig.xml.gz` (`gxspl-ELECTRON.xml.gz` for electrons): the free-nucleon file plus all the stage 3 files.
5. **Reduce and convert to ROOT (stage 5).** Make the smaller `gxspl-NUsmall.xml` (only the isotopes marked `reduced`; neutrinos only) and `xsec_graphs.root`, the cross sections as ROOT TGraphs (`gspl2root`). Nothing is plotted: the graphs are there for users to draw or use from the ROOT file.
6. **Package (stage 6)** the spline files, the configuration and a UPS table into a `genie_xsec` tarball.

**Optional QE averaging (`--average-qe N1[,N2]`):** the QE cross section on a nucleus is computed by throwing nucleons at random, so each run has statistical noise. With this option, stage 3 also runs every pair more times with only the QE generator list and different seeds (stage 3q): N2 times for isotopes of carbon, oxygen and argon, N1 times for the rest. Stage 3s, in parallel jobs alongside stage 4, averages those runs per pair and fits a smooth representative spline; then stage 4q puts the smooth splines in place of the QE splines in the big file.

How it runs:

- **Initialize once, interactively:** `gxspinit` creates the staging area in PNFS, fixes every parameter in a configuration tarball that all the jobs use, and writes the DAG.
- **Submit the DAG:** each grid job fetches its inputs from the staging area and copies its products back, and the DAG makes each stage wait for the one before.
- **Watch and repair:** `gxsp_status` shows which products are missing. Jobs never delete or replace anything, and a job whose products already exist has nothing to do, so failed pieces can be re-run, or the whole DAG resubmitted, safely.

## The scripts

```
gxspinit             initialization (run interactively)
gxsp_mkspl           stages 1 and 3   (gmkspl, one probe x one target)
gxsp_add_freenuc     stage 2          (gspladd -> gxspl-freenuc.xml)
gxsp_add_isotopes    stage 4          (gspladd -> gxspl-NUbig.xml.gz / gxspl-ELECTRON.xml.gz)
gxsp_qe_smoothing    stage 3s         (only with --average-qe: average and smooth each pair's seeded QE files)
gxsp_qe_replace_big  stage 4q         (only with --average-qe: the smooth QE splines into the big file)
gxsp_reduce_graph    stage 5          (awk -> gxspl-NUsmall.xml, gspl2root -> xsec_graphs.root)
gxsp_mkups           stage 6          (UPS genie_xsec tarball)
gxsp_status          progress report (run interactively)
gxsp_check_lastknot  find splines whose last knot is zero after a non-zero one   (see Spline Checking)
gxsp_check_bounce    find splines that drop back to zero after being non-zero
gxsp_check_many0     find splines that start with several zero knots
gxsp_check_cubic     find where the cubic interpolation goes below zero or overshoots
gxsp_check_smooth    compare how smooth the splines of two files are (QE before / after 4q)
gxsp_cut_emax        shorten every spline to end at the first knot at or above an energy
                     (stage 6 puts it in the UPS tarball for users)
gxsplib.sh           shared functions
gxsp_bootstrap.sh    start-up code for grid jobs (pasted into each script by gxspinit)
```

A default isotope list comes with the scripts; `gxspinit` uses it unless `--fetch-isotopes` is given.

## Quick start

### 1. Initialize

Run this on an interactive node where PNFS is mounted:

```bash
gxspinit \
    --top /pnfs/genie/scratch/users/rhatcher/GXSP \
    --sversion v3_08_00 --tune AR23_20i_00_000 \
    --knots 250 --emax 1000 \
    --gversion v3_08_00 --gqual e26:prof
```

| Option | Meaning |
|------------|------------------------------------|
| `--top DIR` | parent of the staging area, e.g. `/pnfs/genie/scratch/users/rhatcher/GXSP` (required) |
| `--sversion V` | spline series version, e.g. `v3_08_00` (required) |
| `--tune T` | tune; must match `^[A-Za-z]+\d{2}_\d{2}[a-z]_[a-z0-9]{2}_[a-z0-9]{3}$`, e.g. `AR23_20i_00_000`, `G18_02a_02_11b`, `GHE19_00a_00_000`. The first two fields are the CMC (`AR23_20i`) (required) |
| `--knots N`, `--emax E` | number of spline knots, and maximum energy in GeV (required) |
| `--gversion V`, `--gqual Q` | GENIE code version and qualifiers to `setup` (required) |
| `--electron` | electron-probe production instead of the 6 neutrino flavors |
| `--nuflavors F1,F2,...` | only these neutrino flavors instead of all six; not with `--electron`. Names or PDG codes, mixed and in any order (`nue,nuebar`, `12,-12`, `numu,-12`); they're always used in the standard order `12 -12 14 -14 16 -16`. Duplicates, empty entries and non-neutrinos are errors. Stored as `NUFLAVORS` (PDG codes) in `define_cfg.sh`, and adds a suffix to the staging area and UPS names (see below) |
| `--genlist L` | event generator list (default `Default`) |
| `--fetch-isotopes FILE` | isotope list (default: the one that comes with the scripts); always stored as `cfg/isotopes.cfg` |
| `--fetch-tune-from DIR` | tune directory to ship in `cfg/` (found through `GXMLPATH` in jobs). It's named after the CMC (e.g. `AR23_20i/`), optionally with a full-tune directory below it (`AR23_20i/AR23_20i_ab_cde/`). Tunes of the form `XXX_YYY_00_000` only use the CMC configuration, so they don't need the subdirectory |
| `--average-qe N1[,N2]` | number of seeded QE stage-3 passes to average (default `0` = none): N2 for the isotopes whose Z is in `--qe-z` (default carbon, oxygen and argon), N1 for all others; N2 defaults to N1. A typical choice is `100,1000`. Stored as `NQE` and `NQEZ` in `define_cfg.sh`. See [[#QE averaging (`--average-qe`)]] |
| `--qe-z Z1,Z2,...` | the Z of the isotopes that get the N2 count of `--average-qe N1,N2` (default `6,8,18`: carbon, oxygen, argon). Positive integers up to 118, comma separated; any order, repeats dropped. Stored as `QEZ` (space separated, sorted) in `define_cfg.sh`. The README names them by element symbol |
| `--average-throws N` | only with `--average-qe` > 0: recorded as `NTHROWS` in `define_cfg.sh` for the seeded QE jobs (default `1`). It must be a positive integer. Given without `--average-qe`, it's ignored with a warning and `NTHROWS` is `0`. More throws make each seeded job slower |
| `--qe-chunk N` | seeds per seeded QE job (default `250`): each QE `-N` cluster in the DAG runs a range of up to `N` seeds. Only in the DAG (not in `define_cfg.sh`). See [[#QE averaging (`--average-qe`)]] for choosing it |
| `--qe-smooth-jobs N` | number of stage 3s jobs (default `25`, at most one per averaged pair); the pairs are shared out so each job fetches about the same number of seed files. Recorded as `QESMJOBS` in `define_cfg.sh` |
| `--qe-scripts DIR` | where stage 3s gets the QE spline smoothing scripts (default `/pnfs/genie/scratch/users/rhatcher/GXSPL/qe_spline_smoothing`); recorded as `QESCRIPTS` in `define_cfg.sh` |
| `--keep-negative` | keep negative cross section values in the splines. By default they're set to `0.0`, as NaN and inf always are, with the original file kept as `<file>-NEG` (see [[#NaN, inf and negative values]]). Recorded as `NEGZERO` (0 with this option, else 1) in `define_cfg.sh` |
| `--tgraph EMAX[,NPOINTS[,log]]` | stage 5's TGraphs (`xsec_graphs.root`) go up to `EMAX` GeV (default `125`), never more than `--emax`: a larger value is lowered to it, with a message. With `NPOINTS` (an integer, at least 2) each TGraph has that many points; without it, as many as the splines have knots. The points are evenly spaced in E, or in log E with `log` (`gspl2root -l`); `EMAX,,log` gives log spacing with the default number of points. Recorded as `TGEMAX`, `TGNPOINTS` (empty without it) and `TGLOG` (1 with `log`, else 0) in `define_cfg.sh` |
| `--jobsub-group G` | jobsub group in the DAG (default `$JOBSUB_GROUP`, else `genie`) |
| `--force` | redo `bin/` and `cfg/` of an already initialized staging area |
| `file.xml ...` | override XML files to ship in `cfg/` |

The staging area it creates is

```
STAGEDIR=$TOP/GXSPL-$SVERSION-$TUNESQ-k$KNOTS-e$EMAX[-electron | -<flavors>]
```

With `--nuflavors`, `-<flavors>` keeps a production with only some flavors apart from the full one (the UPS qualifier and tarball name get the same suffix, with `:` in the qualifier):

```
--nuflavors                         suffix                qualifier ends with
(not given, or all six)             (none)                :e1000
nue,nuebar                          -nue-nuebar           :e1000:nue:nuebar
nue,nuebar,numu,numubar (no tau)    -notau                :e1000:notau
numu,nue                            -nue-numu             :e1000:nue:numu
nutaubar                            -nutaubar             :e1000:nutaubar
```

The script prints the job counts and a suggested submit command when it finishes.

### 2. Submit the DAG

```bash
jobsub_submit_dag -G genie file://$STAGEDIR/cfg/genie_splines.dag
```

### 3. Watch progress

```bash
gxsp_status $STAGEDIR        # add -v to list every file, marked present / missing
                             #     -q for the stage lines only
                             #     -k (--keep-going) to check every stage
```

```
/pnfs/genie/scratch/users/rhatcher/GXSP/GXSPL-v3_08_00-AR2320i00000-k250-e1000
  AR23_20i_00_000  k250  e1000  Default  splines v3_08_00  GENIE v3_08_00 -q e26:prof
  probes: nue,nuebar,numu,numubar,nutau,nutaubar
  average-qe=0
stage 1   free nucleon splines           12 / 12      complete
stage 2   combine free nucleons           1 / 1       complete
stage 3   probe x isotope splines       748 / 750     MISSING 2:
    -p numu     -t C12     work-products/isotopes/gxspl-numu-C12.xml
    -p nutaubar -t Pb208   work-products/isotopes/gxspl-nutaubar-Pb208.xml
```

`gxsp_status` lists the staging area once (a single `find $STAGEDIR`) and checks every stage against that list; the lists of expected files are worked out from the configuration in one `awk` pass each. So apart from the `find`, a full `--average-qe 100,1000` production (about 124,000 files) takes a few seconds. Stages are checked in order, and checking stops at the first stage that has missing files. With `-k` (`--keep-going`) every stage is checked anyway, and the report ends with the stages that are incomplete, e.g. `2 stage(s) incomplete: 3 5`. This is useful when a later stage was run regardless of an earlier one's gaps. The exit status is 0 when every stage is complete and 1 otherwise, with or without `-k`. For stage 3, each missing file is listed with the `-p` / `-t` that remake it (see [[#4. Re-run failed pieces]]), in stage-3 order. At the very end, after the last stage line, the report notes any files whose NaN/inf or negative values were replaced (the `-NAN` and `-NEG` originals), whether or not every stage is complete (see [[#NaN, inf and negative values]]):

```
all stages complete
note: 6 file(s) had NaN/inf values replaced by 0 (originals kept as -NAN):
    work-products/isotopes/gxspl-nue-C12.xml-NAN
    ...
```

There's no note if there are none.

With `-v`, every file of each stage checked is listed, including the present ones in the incomplete stage. Each line starts with `present` or `missing`:

```
stage 2   combine free nucleons           1 / 1       complete
    present  work-products/gxspl-freenuc.xml
stage 3   probe x isotope splines       748 / 750     MISSING 2:
    present  -p nue      -t H2      work-products/isotopes/gxspl-nue-H2.xml
    ...
    missing  -p numu     -t C12     work-products/isotopes/gxspl-numu-C12.xml
    ...
```

Stage 3 and 3s lines give `-p` / `-t` on every line; stage 3q lines give `-p` / `-t` and `--seed=` (a range for consecutive seeds), sorted by pair, then seed (see [[#QE averaging (`--average-qe`)]]).

### 4. Re-run failed pieces

Stage 3 jobs can be re-run by hand with explicit arguments, the `-p` / `-t` that `gxsp_status` lists for the missing file:

```bash
jobsub_submit ...usual options... file://$STAGEDIR/bin/gxsp_mkspl \
    -p numu -t C12 --stagedir $STAGEDIR --input gxspl-freenuc.xml
```

For a seeded QE file add `--genlist CCQE` (`EMQE` for electrons) and the `--seed=N-M` from its line. `--seed` takes `N`, `N,M` or `N-M`, as `--seed X` or `--seed=X`. A missing stage 3s smooth file is remade with `gxsp_qe_smoothing --stagedir $STAGEDIR -p P -t T`, using the stage 3s jobsub options (see [[#QE averaging (`--average-qe`)]]).

Jobs never delete or replace anything in the staging area (see [[#Access to the staging area (PNFS)]]). So re-running is safe in a different way: a job whose products are already there has nothing to do and ends successfully, and one that is partly done only does what's missing.

- **A seed range cut short** (`--seed N,M`) can be submitted again as it is; only the seeds without an output run.
- **Resubmitting the whole DAG** is safe too: the finished jobs end right away, and only the missing pieces are made.
- **To redo something that exists**, remove it interactively first (PNFS is writable from an interactive node), then re-run its job.

After filling in the gaps, submit the remaining stages by hand, or resubmit the DAG.

### 5. Check the splines

`gxsp_check_lastknot`, `gxsp_check_bounce` and `gxsp_check_many0` look for zero cross sections in odd places in the big file (at the last knot, after being non-zero, at the start), and `gxsp_check_cubic` for where GENIE's cubic interpolation goes below zero or overshoots between non-zero knots. With `--average-qe`, `gxsp_check_smooth gxspl-NUbig-stdCCQE.xml.gz gxspl-NUbig.xml.gz` compares the smooth QE splines with the standard ones.

## The stages

Each stage after initialization runs on a grid node. Every grid job takes `--stagedir $STAGEDIR`, writes a `.log` file next to its output, copies its results back to `$STAGEDIR`, waits 25 s (`GXSP_SLEEP`) so PNFS shows the new files before the next DAG step starts, and then exits.

| Stage | DAG entry | Reads | Writes |
|--|---------------------------|---------------|-------------------------|
| 1 | 12 lines (2 for electrons, 2 per flavor with `--nuflavors`): `gxsp_mkspl -p P -t T --stagedir S`, with T `freep` or `freen` | — | `work-products/freenucs/gxspl-<probe>-<T>.xml` |
| 2 | `gxsp_add_freenuc` | all stage 1 files | `work-products/gxspl-freenuc.xml` |
| 3 | one line, `jobsub_submit -N <pairs>`: `gxsp_mkspl --input gxspl-freenuc.xml` | `gxspl-freenuc.xml` | `work-products/isotopes/gxspl-<probe>-<isotope>.xml` |
| 3q | only with `--average-qe`: one `-N <pairs>` line per `--qe-chunk` (250) seeds, `gxsp_mkspl --genlist CCQE --seed first,last`, plus `--pairs qez` or `--pairs other` for the extra seeds (`EMQE` for electrons) | — | `work-products/isotopes/gxspl-<probe>-<isotope>-CCQE/gxspl-<probe>-<isotope>-CCQE-seed<k>.xml` |
| 3s | only with `--average-qe`: `jobsub_submit -N <--qe-smooth-jobs>`: `gxsp_qe_smoothing`, alongside stage 4 | all stage 3q files | per averaged pair `-CCQE-average.*` and `-CCQE-smooth.*` |
| 4 | `gxsp_add_isotopes` | `gxspl-freenuc.xml`, all stage 3 files | `work-products/gxspl-NUbig.xml.gz` (`gxspl-ELECTRON.xml.gz`); with `--average-qe` instead `work-products/gxspl-NUbig-stdCCQE.xml.gz` (`gxspl-ELECTRON-stdEMQE.xml.gz`) |
| 4q | only with `--average-qe`: `gxsp_qe_replace_big` | all stage 3s smooth files, `gxspl-NUbig-stdCCQE.xml.gz` | `work-products/gxspl-NUbig.xml.gz` |
| 5 | `gxsp_reduce_graph` | `gxspl-NUbig.xml.gz` (uncompressed locally) | `work-products/gxspl-NUsmall.xml` (neutrinos only), `work-products/xsec_graphs.root` |
| 6 | `gxsp_mkups` | `gxspl-freenuc.xml`, `gxspl-NUbig.xml.gz`, stage 5 outputs, plus `cfg/` | `ups/genie_xsec-<X.YY.ZZ>-noarch-$TUNESQ-k$KNOTS-e$EMAX[-electron].tar.bz2`, `work-products/README.md` |

The DAG layout is:

```
<parallel>  stage 1 (12 lines)                                   </parallel>
<serial>    stage 2                                               </serial>
<parallel>  stage 3 Default cluster (-N pairs)
            + QE clusters (-N pairs, --seed ranges of <=250, 15h;
              seeds above min(N1,N2) with --pairs qez|other)       </parallel>
<serial>    stage 4, stage 5, stage 6                             </serial>
```

With `--average-qe`, the last block is instead:

```
<parallel>  stage 4
            + stage 3s (--qe-smooth-jobs lines, --job K, 12h)      </parallel>
<serial>    stage 4q, stage 5, stage 6                            </serial>
```

Stage 3 covers the production's isotopes except the two free nucleons, free-n (1000000010) and H1 (1000010010): their splines are exactly the stage 1 ones, already in `gxspl-freenuc.xml`, which stage 4 starts from (see [[#Free nucleons come from stage 1]]). With the default isotope list, a neutrino production has 127 isotopes, 125 of them in stage 3, so stage 3 is `-N 750`. An electron production has 10 isotopes (the lines marked `electron`), 9 in stage 3, so it is `-N 9`. If the list has nothing but the free nucleons, the DAG has no stage 3 at all. The jobsub resource options (memory, disk, expected lifetime, the SL7 singularity image) are set near the end of `gxspinit` (see [[#Job resources (jobsub options)]]).

### Stages 1 and 3: `gxsp_mkspl`

```
gxsp_mkspl [-p P -t T | --pairs G] --stagedir S [--input FILE] [--seed N | --seed N,M] [--genlist L]
```

- **`-p` / `--probe` and `-t` / `--isotope`:** each takes a name or a PDG code (`nue`/`12`, `Ar40`/`1000180400`, `freep`, `freen`). The script converts names to PDG codes for `gmkspl`. The short and long forms are interchangeable; the generated DAG uses `-p` and `-t`, before `--stagedir`.
- **`--input FILE`:** the script fetches `FILE` from `work-products/`, passes it to `gmkspl` as `--input-cross-sections`, and writes its output to `work-products/isotopes/`. Without `--input`, output goes to `work-products/freenucs/`, except for the seeded QE runs of stage 3q (below), which don't need an input file. In practice, `--input` marks a stage 3 job.
- **Neither `-p` nor `-t`:** the script picks the pair from the stage-3 list using `$PROCESS`, and either `--input` (stage 3) or `--genlist CCQE --seed ...` (stage 3q, in an `--average-qe` production) is required. This is how the stage 3 and 3q DAG lines run. Giving only one of `-p`/`-t` is an error.
- **`--pairs G`:** which list `$PROCESS` indexes: `all` (the default), `qez` (only the isotopes whose Z is in `QEZ`, from `--qe-z`; default carbon, oxygen and argon) or `other` (all but those). The `--average-qe N1,N2` clusters for the seeds above min(N1, N2) use it (see [[#QE averaging (`--average-qe`)]]). It can't be combined with `-p`/`-t`, and a `$PROCESS` beyond that list is an error (exit 11).
- **Outputs that already exist:** a run (or, with a seed range, a seed) whose `.xml` is already in the staging area is skipped, not redone. If all of them exist, the job ends at once without fetching its input. A re-run's log goes to `.log.1`, `.log.2`, … if the earlier log is still there.
- **`--seed N`:** an optional integer, passed unchanged to `gmkspl` as `--seed N`. Without it, no seed is passed.
- **`--seed N,M`:** runs `gmkspl` once for each seed from N to M in turn. `--seed N` is the same as `--seed N,N`, and `--seed N-M` the same as `--seed N,M`. Any of these can also be written with `=` instead of a space (`--seed=5`, `--seed=1,4`, `--seed=1-3`), the form `gxsp_status` prints.
  - **Copied back as it goes:** each seed's `.xml` and `.log` (and any `-NAN` or `-NEG` file) go to the staging area as soon as that seed finishes, and the local copy is then deleted. If the job is interrupted or times out partway through, the seeds already done are saved.
  - **Failures:** if one seed fails, the loop carries on with the rest. The job exits with the first failing seed's status and reports how many succeeded (e.g. `3 of 4 succeeded`).
  - **Done once per job:** the pair selection, the input-file fetch and (for seeded QE jobs) the `NewQELXSec.xml` copy. Each seed's log still records the `NewQELXSec.xml` source.
  - **Values:** seeds are integers of 0 or more (`0,2` gives `-seed0`, `-seed1`, `-seed2`). A negative seed, N greater than M, or anything other than `N` or `N,M`, is an error. The generated DAG uses seed ranges for the QE clusters (see [[#QE averaging (`--average-qe`)]]); otherwise add `--seed` to a DAG line or to a job you run by hand when you need it.
- **`--genlist L`:** uses `L` for `gmkspl --event-generator-list` in this job, in place of the `GENLIST` set at initialization. The job prints a message when it overrides the value, and its log header shows the list it used. It only affects this `gmkspl` run: stage 5's `gspl2root` and the UPS table and README still use the initialization `GENLIST`.
- **Output names for non-standard runs:** when `--genlist` is given and isn't `Default`, `-<genlist>` is added to the output name. When `--seed` is given, `-seed<N>` is added. The `.log` (and any `-NAN`) file names follow:

  ```
  options                        output file
  (neither)                      gxspl-numu-Ar40.xml
  --genlist CCQE                 gxspl-numu-Ar40-CCQE.xml
  --seed 12345                   gxspl-numu-Ar40-seed12345.xml
  --genlist CCQE --seed 12345    gxspl-numu-Ar40-CCQE-seed12345.xml
  --genlist Default --seed 7     gxspl-numu-Ar40-seed7.xml
  ```

  These files go in the same `work-products/` directory, next to the plain-named ones. Stages 2 and 4 and `gxsp_status` only look for the plain names, so these runs don't count toward or feed into the production.

The command it runs is:

```
gmkspl -e $EMAX -n $KNOTS --tune $TUNE -t <target> -p <probe> \
       --event-generator-list <--genlist or $GENLIST> [--input-cross-sections FILE] [--seed N] \
       -o gxspl-<probe>-<target>.xml
```

### Stage 2: `gxsp_add_freenuc`

```
gxsp_add_freenuc --stagedir S
```

1. If `work-products/gxspl-freenuc.xml` already exists, there's nothing to do.
2. Check that all stage 1 files (`work-products/freenucs/gxspl-<probe>-freep.xml` and `-freen.xml` for every probe) exist, re-listing and re-checking a few times if some aren't listed yet (see [[#Access to the staging area (PNFS)]]).
3. Copy them in and combine them with `gspladd` into `gxspl-freenuc.xml`, then upload it to `work-products/` with its log, `gxspl-freenuc.log`.

`gxspl-freenuc.xml` is the `--input` of every stage 3 job, the starting point of the big file in stage 4, and it goes in the UPS tarball in stage 6.

### Stage 4: `gxsp_add_isotopes`

1. Check that `gxspl-freenuc.xml` and all stage 3 files exist, using one listing of the staging area (re-listed and re-checked a few times if anything is missing; see [[#Access to the staging area (PNFS)]]).
2. For each stage-3 isotope, copy in its per-probe files, combine them into `gxspl-<isotope>.xml`, and delete the per-probe copies right away. For electrons there is only one probe, so the file is just renamed.
3. Combine `gxspl-freenuc.xml` (the free-n and H1 splines) and all the per-isotope files into the big file, gzip it, and upload it:
   - without `--average-qe`: as `work-products/gxspl-NUbig.xml.gz`, the big file itself
   - with `--average-qe` (`NQE` > 0): as `work-products/gxspl-NUbig-stdCCQE.xml.gz`, the big file with the standard QE splines (for electrons `gxspl-ELECTRON-stdEMQE.xml.gz`). Stage 4q reads it and writes `gxspl-NUbig.xml.gz` as a new file, so no job ever has to delete or replace the big file.

The big file is always kept gzipped in `work-products/`. Every version (`gxspl-NUbig.xml.gz`, `gxspl-NUbig-stdCCQE.xml.gz`, and the copy in the UPS tarball) is compressed with `gzip --best`, and it's only ever uncompressed in a worker node's scratch area. Stages 4q and 5 uncompress their local copy as needed. Stage 6 puts the `.gz` straight into the tarball for neutrinos; for electrons it uncompresses it, because `gxspl-ELECTRON.xml` is the default file named in the table.

### QE averaging (`--average-qe`)

With `gxspinit --average-qe N1[,N2] [--average-throws NTHROWS]`, and N1 or N2 > 0:

- **How many seeds:** isotopes whose Z is in `--qe-z` (default 6, 8, 18: carbon, oxygen and argon, any A) get N2 seeded runs per probe; all other stage-3 isotopes get N1. `--average-qe N` is the same as `N,N`. A typical choice is `--average-qe 100,1000`. A count of 0 leaves that group's QE splines as they come from stage 3: no seeded runs, no directories, no averaging. The Z comes from the PDG code `10LZZZAAAI`: Z = (pdg / 10000) % 1000, e.g. 18 for Ar40 (`1000180400`).

- **Which QE:** the seeded runs use the QE event generator list `<QE>`: `CCQE` for neutrinos, `EMQE` for electrons (`--electron`). `<QE>` also names the files and directories. The examples below are for neutrinos; for electrons read `EMQE` for `CCQE`, e.g. `gxspl-electron-C12-EMQE/gxspl-electron-C12-EMQE-seed1.xml`, `gxspl-ELECTRON-stdEMQE.xml.gz`.

- **Configuration:** `define_cfg.sh` records `NQE` (N1), `NQEZ` (N2), `QEZ` (`--qe-z`, default `"6 8 18"`: the Z values that get N2) and `NTHROWS` (default 1), so every grid job has them in its environment. The seeded QE jobs use `NTHROWS` to set `NumNucleonThrows` in their copy of `NewQELXSec.xml` (below). `gxspinit` and `gxsp_status` both show it.

- **Stage 3q (DAG):** the seeds are split into ranges of up to `--qe-chunk` seeds (default 250), and each range becomes one more cluster in the stage-3 `<parallel>` block, with `--expected-lifetime 15h`. Each runs `gxsp_mkspl --genlist CCQE --seed first,last`, with no `--input`: the QE splines on a nucleus don't need the free-nucleon ones. Each job loops over its range of seeds and copies each one back as it finishes. So that no job starts only to find it has nothing to do:
  - seeds 1..min(N1, N2) run for **every** pair, with the same `-N <probes x isotopes>` as the Default cluster;
  - the seeds above that, up to max(N1, N2), run only for the pairs that need them, with `--pairs qez` (the pairs whose isotope's Z is in `QEZ`, when N2 > N1) or `--pairs other` (the rest, when N1 > N2), and `-N` the number of those pairs. `--pairs` makes `gxsp_mkspl` map `$PROCESS` over that group only (see [[#Mapping `$PROCESS` to a pair]]).

  With the default neutrino list (125 stage-3 isotopes, 9 of them carbon, oxygen or argon: 750 pairs, 54 C/O/Ar, 696 others) and the default `--qe-chunk 250`:

  ```
  --average-qe            QE clusters
  1                       -N 750 --seed 1,1
  251                     -N 750 --seed 1,250
                          -N 750 --seed 251,251
  600                     -N 750 --seed 1,250
                          -N 750 --seed 251,500
                          -N 750 --seed 501,600
  20,150                  -N 750 --seed 1,20
                          -N 54  --seed 21,150     --pairs qez
  100,1000  (typical)     -N 750 --seed 1,100
                          -N 54  --seed 101,350    --pairs qez
                          -N 54  --seed 351,600    --pairs qez
                          -N 54  --seed 601,850    --pairs qez
                          -N 54  --seed 851,1000   --pairs qez
  100,1000 --qe-chunk 100 -N 750 --seed 1,100
                          -N 54  --seed 101,200    --pairs qez
                          ...    (nine qez clusters in all)
                          -N 54  --seed 901,1000   --pairs qez
  0,100                   -N 54  --seed 1,100      --pairs qez
  50,10                   -N 750 --seed 1,10
                          -N 696 --seed 11,50      --pairs other
  ```

  **Choosing `--qe-chunk`:** a job runs its seeds one after another, so its time is about the chunk size × the time per seed. A seeded CCQE `gmkspl` with 1 throw takes about 54 s (1000 of them in about 15 h), so 250 seeds take about 3.8 h of the 15 h `--expected-lifetime`. More `--average-throws` take longer. `gxspinit` prints its estimate, scaling the time in proportion to the throws (a rough guess), and warns if that is over 15 h. For example, `--average-throws 5` with the default chunk gives about 18.8 h and a warning, and `--qe-chunk 100` brings it down to about 7.5 h. A job that runs out of time keeps the seeds it finished; resubmitting the same line runs only the missing ones.

  The lifetime (`QET`) is set with the other jobsub options near the end of `gxspinit`. Their outputs follow the usual naming for non-standard runs, `gxspl-<probe>-<isotope>-CCQE-seed<k>.xml`, but go in a directory of their own per probe x isotope:
  ```
  work-products/isotopes/gxspl-numu-Ar40-CCQE/
      gxspl-numu-Ar40-CCQE-seed1.xml   gxspl-numu-Ar40-CCQE-seed1.log
      gxspl-numu-Ar40-CCQE-seed2.xml   gxspl-numu-Ar40-CCQE-seed2.log  ...
  ```
  `gxspinit` creates these directories for every pair whose count is > 0 (750 with the default list and both counts > 0: no free-n or H1). If one is missing, the job creates it (`ifdh mkdir` on the grid). The rule `gxsp_mkspl` uses is: a stage-3 isotope (not free-n or H1), `--genlist CCQE`, a `--seed`, and `--average-qe` in the configuration (N1 or N2 > 0), with or without `--input`. So a CCQE seeded run made by hand for that production also goes in the right directory. Any other run, including a seeded QE one in a production without `--average-qe`, still goes directly in `work-products/isotopes/`.
- **`NewQELXSec.xml` in the seeded jobs:** before running `gmkspl`, each seeded QE job copies `NewQELXSec.xml` into its working directory, which is first on `GXMLPATH`. It takes the first copy it finds, in this order (`<CMC>` is e.g. `AR23_20i` for `AR23_20i_00_000`):
  1. a `NewQELXSec.xml` given to `gxspinit` as an override file
  2. the shipped tune directory (`--fetch-tune-from`): `<CMC>/<TUNE>/NewQELXSec.xml`, then `<CMC>/NewQELXSec.xml`
  3. `$GENIE/config/<CMC>/<TUNE>/`, then `$GENIE/config/<CMC>/`, then `$GENIE/config/`

  The `<TUNE>/` levels are optional: a tune such as `AR23_20i_00_000` usually has only CMC-level configuration, and when there's no `<TUNE>/` directory the lookup goes on to the CMC level.

  Because the working directory is the first `$GXMLPATH` entry, GENIE uses this edited copy ahead of any tune or CMC copy (see [[#How GENIE finds XML files]]).

  It then sets the nucleon throws in the local copy to `NTHROWS`:
  ```bash
  awk -v n="${NTHROWS}" '{ sub(/NumNucleonThrows"> *[0-9]+/, "NumNucleonThrows\"> " n " ") } 1'
  ```
  The job log records which file was used and the resulting `NumNucleonThrows` line. The job fails with exit code 22 if there is no `NewQELXSec.xml` anywhere, or if the copy it found has no `NumNucleonThrows` value to set. Other jobs (Default stage 3, or a CCQE run without `--seed`) skip this step. `gxspinit` warns if the `--fetch-tune-from` directory isn't named after the tune's CMC, since the lookup (and GENIE) expect `<CMC>/<TUNE>/`.
- **`gxsp_status`:** after the Default stage-3 files are complete, it checks the seeded QE files as stage `3q`. It lists the missing ones sorted by probe x isotope (in stage-3 order), then by seed, with the `-p` / `-t` values for `gxsp_mkspl` and the seeds. Consecutive missing seeds of a pair go on one line, as `--seed=N-M` and a file name with `seed{N..M}`; after a gap in the seeds, a new line starts:
  ```
  stage 3q  QE seeded splines            64 / 72      MISSING 8:
      -p nue      -t C12     --seed=1-3       work-products/isotopes/gxspl-nue-C12-CCQE/gxspl-nue-C12-CCQE-seed{1..3}.xml
      -p nue      -t C12     --seed=5         work-products/isotopes/gxspl-nue-C12-CCQE/gxspl-nue-C12-CCQE-seed5.xml
      -p numu     -t N14     --seed=1         work-products/isotopes/gxspl-numu-N14-CCQE/gxspl-numu-N14-CCQE-seed1.xml
      -p numu     -t O16     --seed=2-4       work-products/isotopes/gxspl-numu-O16-CCQE/gxspl-numu-O16-CCQE-seed{2..4}.xml
  ```
  The `MISSING` count is of files, not lines. Each line can be redone by a make-up job with its `-p`, `-t` and `--seed=`, as printed (`gxsp_mkspl` takes `--seed` as `N`, `N,M` or `N-M`, with or without the `=`):
  ```bash
  gxsp_mkspl --stagedir $STAGEDIR --genlist CCQE -p nue -t C12 --seed=1-3
  ```
  With `-p` and `-t`, `gxsp_mkspl` doesn't use `$PROCESS`, and a seed range can cross the DAG's QE clusters: in the example, `numu-O16` seed 2 runs in the all-pairs cluster and seeds 3–4 in the `--pairs qez` one, and `--seed=2-4` redoes all three. If any Default stage-3 file is missing, it stops there and doesn't check the QE files. With `-v`, present files are grouped into ranges the same way.
- **Free-nucleon splines removed if `--input` is given:** `gmkspl` writes every spline of its `--input-cross-sections` file into its output, so a seeded QE file made with `--input gxspl-freenuc.xml` would also carry all of it, about 11 MB per file with AR23 at 250 knots, against a few tens of KB for the pair's own QE spline(s). The DAG doesn't pass `--input` to these jobs, but if one is given (by hand, or in an older DAG line), then after each seed's `gmkspl` the job removes the splines whose target is free-n or H1 (`tgt:1000000010;`, `tgt:1000010010;`) before copying the file back, and its log says how many went (e.g. `24 free-nucleon spline(s) from gxspl-freenuc.xml removed; 2 left`). Stage 3s doesn't use those splines (`average_splines.py` skips free nucleons). Default stage-3 files and CCQE runs without `--seed` keep them. Stage 3s works the same whether or not a seed file has them.
- **Stage 3s, `gxsp_qe_smoothing`:** averages each pair's seeded files and fits a smooth representative spline. It needs only the stage 3q files, so in the DAG it runs alongside stage 4, as `--qe-smooth-jobs` jobs (default 25, at most one per averaged pair):
  ```
  <parallel>  stage 4                                  (gxsp_add_isotopes)
              stage 3s --job 0                         (gxsp_qe_smoothing)
              ...
              stage 3s --job 24                        (gxsp_qe_smoothing)   </parallel>
  <serial>    stage 4q, stage 5, stage 6                                     </serial>
  ```
  **One DAG node per job, not a `-N` cluster:** when one job of a node's `-N` cluster fails, DAGMan removes the rest of that cluster, and removed jobs send back no output. A stage 3s job runs for hours and makes thousands of copies, so one failure would cost every job's unfinished work. As separate nodes, a failure stops only that job; the others run to the end. Stage 4q and the rest still wait for all of them, and don't run if any failed.
  **Sharing out the pairs:** every job works out the same assignment from the same lists (`gxsp_qe_smooth_assign`): the pairs go, largest seed count first, each to the job with the fewest seed files so far. So the jobs fetch about the same number of seed files. With the default neutrino list and `--average-qe 100,1000` (54 pairs × 1000 seeds plus 696 × 100, 123,600 seed files), each of the 25 jobs gets at most 5000. `gxspinit` prints the number of jobs and the largest share.

  **Resources:** 4 GB memory, 30 GB disk (`QEAVGD`: one pair's seed files at a time, up to about 11 GB for 1000 files if they were made with `--input` and hold the free-nucleon splines, plus the venv), `--expected-lifetime 12h` (`QEAVGT`). Like stage 4q, its DAG line has **no `--singularity-image`**: the job runs directly on the worker node, not in the SL7 container, where the smoothing scripts had trouble. It doesn't set up GENIE (it only needs `ifdh`).

  **Scripts and python:** it uses `average_splines.py` and `construct_representative_spline.py` from the directory given by `gxspinit --qe-scripts`. The default is `/pnfs/genie/scratch/users/rhatcher/GXSPL/qe_spline_smoothing`, and the directory is recorded as `QESCRIPTS` in `define_cfg.sh`. `gxspinit` warns if it can't see the two scripts there. The scripts run with `python3` (`GXSP_PYTHON` overrides). The `python3` on the node doesn't have the packages they need, so, unless `GXSP_PYTHON` is set, each job first builds a venv in its scratch area (`pyenv/`) and runs the scripts with it:
  ```bash
  python3 -m venv pyenv && source pyenv/bin/activate
  python3 -m pip install numpy
  python3 -m pip install pandas==2.2.3        # newer and meson complains
  python3 -m pip install fastparquet "pandas<2.3"
  python3 -m pip install tqdm colorama
  python3 -m pip install scikit-learn         # provides sklearn
  ```
  This needs the worker node to reach PyPI, and takes some time and disk at the start of every job. If it fails, the job stops with exit 24. `GXSP_VENV=0` turns the venv off: the job then warns and uses `python3` as is, which will only work if it has the packages.

  **What each job does:**
  1. Work out its pairs (from `--job K`, else `$PROCESS`, or the one pair given with `-p`/`-t`). If all their smooth files are already in the staging area, there's nothing to do. Otherwise check that all their seeded files exist, re-listing and re-checking a few times if some aren't listed yet (exit 17 if they're still missing), and copy the two scripts from `QESCRIPTS` (exit 23 if that fails).
  2. For each of its pairs, with `D` = `gxspl-<probe>-<isotope>-CCQE` (e.g. `gxspl-numu-Ar40-CCQE`):
     - copy its seeded files (N2 for the `--qe-z` isotopes; N1 for the others) from `work-products/isotopes/D/` into a local `./D/`
     - `python3 average_splines.py -b D -o D-average` makes `D-average.xml` and one parquet per spline it averages. Neutrino pairs have one spline, so that is `D-average.parquet`. Electron pairs have two, one per hit nucleon, so they get `D-average_nu11_tgt<isotope pdg>_N2212.parquet` and `..._N2112.parquet` instead, and no `D-average.parquet` (exit 20 if an expected parquet is missing; the error names the parquets that were made). The `.xml` and the parquets are copied back to `work-products/isotopes/D/`, and the local `.xml` and `./D/` are deleted, so only one pair's seed files are on the node at a time.
     - `python3 construct_representative_spline.py -i D-average.parquet -o D-smooth --tune $TUNE -p <probe pdg> --tgtpdg <isotope pdg>` makes `D-smooth.xml` and `.parquet`. NaN/inf (and, unless `--keep-negative`, negative) values in `D-smooth.xml` are set to 0.0, the original going back as `D-smooth.xml-NAN` or `-NEG`. Both files are copied back and the local copies deleted.
       **Electrons:** this runs twice, once per hit nucleon, each from that nucleon's parquet: `-i D-average_nu11_tgt<T>_N2212.parquet --nucleon 2212 -o D-smooth-p` and `-i D-average_nu11_tgt<T>_N2112.parquet --nucleon 2112 -o D-smooth-n`.
  3. Copy its log back as `work-products/stage3s-job<K>.log` (`stage3s-<probe>-<isotope>.log` with `-p`/`-t`).

  Free-n and H1 aren't averaged (their QE splines stay as they come from stage 1), and neither is a group with count 0. Nothing in the staging area is deleted or replaced. A pair that fails doesn't stop the job: a script failing or not making its expected files (20), a seed file that still can't be copied after the retries (13), or a result that can't be copied back (15). The job cleans up that pair's local files, says `D failed (status N); going on with the other pairs`, and does the rest. At the end it copies its log back and exits with the first failure's status, naming the failed pairs (`1 of 31 pair(s) failed: gxspl-numu-C12-CCQE`). Problems that would hit every pair (missing seeded files, the scripts or the venv) still stop the job at the start. A re-run picks up where the last attempt stopped: a pair whose smooth file(s) are there is skipped, and one whose average parquet(s) are there (those its missing smooth files need) skips the averaging.

  **`gxsp_status`** checks stage 3s after stage 3q: one product per averaged pair, `D-smooth.xml` (for electrons two, `D-smooth-p.xml` and `D-smooth-n.xml`). Missing ones are listed with the `-p` / `-t` that redo them, in stage-3 order:
  ```
  stage 3s  QE average + smooth           749 / 750     MISSING 1:
      -p numu     -t C12     work-products/isotopes/gxspl-numu-C12-CCQE/gxspl-numu-C12-CCQE-smooth.xml
  ```
  A make-up job for one pair:
  ```bash
  jobsub_submit ...stage 3s options... file://$STAGEDIR/bin/gxsp_qe_smoothing \
      --stagedir $STAGEDIR -p numu -t C12
  ```
  Or resubmit the stage 3s DAG line as it is: the jobs whose pairs are all done end at once.
- **Stage 4q, `gxsp_qe_replace_big`:** runs after stages 4 and 3s and before stage 5. It puts the smooth splines into the big file, in one pass, without python:
  1. If `work-products/gxspl-NUbig.xml.gz` already exists, there's nothing to do. Otherwise check that stage 4's `gxspl-NUbig-stdCCQE.xml.gz` and every stage 3s smooth file exist (exit 17 if some are still missing after the re-checks).
  2. Copy in all the smooth files.
  3. Read `gxspl-NUbig-stdCCQE.xml.gz` as a stream (`gzip -dc`), line by line, and write `gxspl-NUbig.xml` (`gxsp_qe_replace` in `gxsplib.sh`): a spline whose name matches one in a smooth file gets that file's knots, and its `nknots` is set to their number; every other line is copied as it is. Memory use is small whatever the size of the big file, and the pass takes about as long as reading and writing it once.
  4. Gzip the result and copy it in as the new file `work-products/gxspl-NUbig.xml.gz`.

  The log (`work-products/stage4q.log`) gives a summary: how many smooth splines, how many replaced, and any whose number of knots differed from the spline it replaced (the first 5 are named). If a smooth spline has no spline of the same name in the big file (or has no knots, or is in two files), the job stops with exit 20 and names it in the log; `gxspl-NUbig.xml.gz` doesn't exist until a run succeeds. Each run starts again from stage 4's untouched `gxspl-NUbig-stdCCQE.xml.gz`, so replacements are never applied twice.

  Resources: 4 GB memory, 10 GB disk (`BIGD`: the gzipped stage 4 file, the uncompressed new one, and its gzipped copy), 3h (`SHORTT`), and like stage 3s no SL7 container. Stages 5 and 6 then use the new big file. Without `--average-qe` (both counts 0), stages 3s and 4q aren't in the DAG, and the scripts exit 0 straight away if run anyway.

  For `gxsp_status`, stage 4 is complete when `gxspl-NUbig-stdCCQE.xml.gz` exists, and stage 4q when the new `gxspl-NUbig.xml.gz` is there.

  **Without these scripts:** GENIE's `replace_splines.py` (in the smoothing scripts' directory) does the same replacement for anyone running the smoothing by hand. Since 2026-11-17 it works the way stage 4q does: one pass over the big file, line by line (plain or `.gz`), any number of smooth files (`-r f1.xml f2.xml ...` or `-l list`), `nknots` set to the new knot count, and a summary. It writes no output, and exits 1, if a smooth spline isn't in the big file (`--allow-missing` to accept that), has no knots, is in two files, or is for another tune (`--ignore-tune`). Only the Python standard library is needed. On the G24_12a big file (1.3 GB uncompressed) it took about 14 s, and its output was the same as `gxsp_qe_replace`'s, line for line:

  ```
  python3 replace_splines.py -b gxspl-NUbig-stdCCQE.xml.gz -l smooth.list -o gxspl-NUbig.xml.gz
  ```

  The old version loaded the whole big file into ElementTree (far more than 4 GB of memory), took one smooth file per run, skipped silently any spline it couldn't match, and left `nknots` unchanged.

### Tunes where QE averaging does nothing

QE averaging only helps when the QE cross section comes from random nucleon throws. Some QE models don't throw nucleons at all: for them, every seed gives the same spline, `--average-throws` changes nothing, the average is just the standard spline, and the smoothing fit can only move it away from the model's exact value. Don't use `--average-qe` for those tunes: G24_12a, AR23/AR25 22i and 22m, and G21_11*.

- **G24_12a_00_000 (`genie::HybridXSecAlgorithm/Martini-QEL`): deterministic** (checked 2026-11-11). Nine seeded νμ Ar40 files (`--genlist CCQE`, 1 throw, seeds 1–9) were byte-for-byte identical to each other, and their knots identical to the standard spline (500 throws) in `gxspl-NUbig-stdCCQE.xml.gz`.
- **CRPA-QEL (AR23/AR25 22i and 22m) and SuSAv2-QEL (G21_11*): deterministic too**, from GENIE's source and config (no seeded files compared yet). `CRPA-QEL` gives nuclear targets to `SuSAv2QELPXSec/CRPA`, `SuSAv2-QEL` to `SuSAv2QELPXSec/SuSAv2`, and `Martini-QEL` to `MartiniQELPXSec/Default`. All three integrate with `genie::MECXSec/Default` (VEGAS over q0, q3 of a hadron-tensor table), not `NewQELXSec`, which is the only place `NumNucleonThrows` is used. `MECXSec` never calls GENIE's random number generator, so the seed doesn't reach it, and none of these tunes overrides the integrator. The G24_12a files above show this route gives the same spline for every seed.
- **Jitter that seeds can't average away:** VEGAS is itself a Monte Carlo integration (200,000 evaluations, 0.05% relative tolerance), but it starts the same way in every run. So these splines can still have some jitter, the same in every seed: averaging can't remove it, though smoothing could. Whether it is noticeable hasn't been looked at.

**To check the CRPA splines for jitter:**

- **Which tune:** any one of the four CRPA tunes (AR23/AR25 22i and 22m) will do. Their QE CC splines should be the same:
  - 22i and 22m differ only in the FSI model (`HAIntranuke2018` vs `HAIntranuke2025`), which doesn't enter the splines.
  - Where the tunes' nuclear models differ, it doesn't matter either: `SuSAv2QELPXSec` doesn't use the nuclear model, since the physics is in the hadron-tensor tables.

  Running a second CRPA tune and comparing would confirm this.
- **Which isotopes:** CRPA hadron-tensor tables exist only for C12, O16 and Ar40, for ν and ν̄ (`CRPAQELHadronTensorModel.xml`). Every other nucleus uses one of those tables, scaled by A and by the ratio of Fermi momenta. So C12, O16 and Ar40 are where the real CRPA splines are; the others inherit the shape, and any jitter, of their table.
- **The production:** an ordinary one, without `--average-qe`. For this check alone it can be cut down to `--genlist CCQE` and an isotope list (`--fetch-isotopes`) of C12, O16 and Ar40.
- **Measuring it:** `gxsp_check_smooth` compares two files, but given the same file twice with `-a` it reports each spline's jitter and flips per energy band (the "old" columns):

  ```
  gxsp_check_smooth -a -q -m 'Weak\[CC\],QES' gxspl-NUbig.xml.gz gxspl-NUbig.xml.gz
  ```

  `-m 'Weak\[CC\],QES'` is needed because the default `QELCC` doesn't match the CRPA spline names. Flips near 0.7 mean noise, near 0 a smooth curve.
- **Nieves baseline (2026-11-11):** `gxsp_check_smooth -m 'Weak\[CC\],QES' gxspl-NUbig-stdCCQE.xml.gz gxspl-NUbig.xml.gz` on an `--average-qe` production of a Nieves tune: 750 smoothed splines compared (2274 other matching splines identical in both files). Medians over the splines, standard run (500 nucleon throws) vs after averaging and smoothing:

  ```
  band      jitter std  jitter smooth  factor   flips std  flips smooth
  thr          0.326       0.118          2.2      0.67       0.00
  <0.1         0.209       0.0166        12        0.72       0.08
  0.1-1        0.035       0.00389        8.5      0.71       0.00
  1-10         0.00807     0.000329      22        0.69       0.00
  10-100       0.00608     5.6e-05      121        0.72       0.00
  >100         0.00555     1.5e-05      377        0.71       0.03
  ```

  So the standard Nieves splines are pure noise above 1 GeV at about 0.6–0.8% knot to knot (median). CRPA splines with flips near 0.7 at that level would be as noisy as unsmoothed Nieves ones; jitter near 1e-4 with flips near 0 would be as smooth as the smoothed ones.
- **Throw-based, averaging helps:** `NievesQELCCPXSec` (G18_10*, AR23/AR25 20i/20m and their sub-tunes, G24_20*, N24), `LwlynSmithQELCCPXSec` (G18_01*, G18_02*) and, for electrons, `UnifiedQELPXSec` (GEM26).

**To check a tune** before a long `--average-qe` production: make two seeded QE files for one pair, e.g. by hand with `gxsp_mkspl --stagedir $STAGEDIR --genlist CCQE -p numu -t Ar40 --seed=1-2` in a staging area of that tune, then compare them with `cmp`. If they are identical (and match the standard spline), averaging is pointless for that tune.

**Re-initializing without QE averaging** (a staging area already set up with `--average-qe`):

1. If the DAG is running, remove it first (`jobsub_rm` on its id).
2. Run `gxspinit` again with the same options, but without `--average-qe` (drop `--qe-z`, `--average-throws`, `--qe-chunk`, `--qe-smooth-jobs` and `--qe-scripts` too; `--average-throws` alone only gets a warning), and add `--force`. Give any override XML files again. `--force` rewrites only `bin/` and `cfg/` (with `NQE=0` in `define_cfg.sh` and a DAG without stages 3q, 3s and 4q); `work-products/` is left as it is.
3. If stage 4 has already run, its output is `work-products/gxspl-NUbig-stdCCQE.xml.gz`, but without `--average-qe` stage 4 writes `gxspl-NUbig.xml.gz`. The content is the same, so either copy it to the new name (`ifdh cp`, or `cp` where PNFS is mounted) or let stage 4 run again.
4. Submit the new DAG. Jobs whose products already exist end at once, so finished stage 1–3 work isn't redone.
5. The seeded directories `work-products/isotopes/gxspl-<probe>-<isotope>-CCQE/` are ignored with `NQE=0`; delete them if the space is needed.

### Stage 5: `gxsp_reduce_graph`

- **Neutrinos only:** runs `gawk -f reduce_gxspl.awk gxspl-NUbig.xml > gxspl-NUsmall.xml`.
- **Always:** runs `gspl2root` on the big file, for every probe and for the isotopes marked `root` in `isotopes.cfg` that are part of this production. If none are marked, it uses all of them. The energy limit is `-e $TGEMAX` (`gxspinit --tgraph`, default 125 GeV), or `$EMAX` if that is smaller. The lower limit is GENIE's `GVLD-Emin` for the tune (0.01 GeV for G18 and AR23). With `--tgraph EMAX,NPOINTS` it also passes `-n NPOINTS`; otherwise each TGraph has as many points as the longest spline has knots. The points are evenly spaced in E, or in log E with `--tgraph EMAX,[NPOINTS],log` (stage 5 then passes `-l`), and are the same for every TGraph in the file. Both the knot-count default and `-n` are new in GENIE 3.08.00 (September 2026); earlier versions always make 1000 points and have no `-n`, so `NPOINTS` needs 3.08.00 or later.

### Stage 6: `gxsp_mkups`

This stage builds the tree below, puts it in a `.tar.bz2`, and uploads it to `$STAGEDIR/ups`. It doesn't need GENIE set up, only `ifdh`.

```
genie_xsec/$SVERSION.version/NULL_$UPSQUALDASH
genie_xsec/$SVERSION/NULL/$UPSQUALDASH/ups/genie_xsec.table
genie_xsec/$SVERSION/NULL/$UPSQUALDASH/data/
    README.md  isotopes.cfg  xsec_graphs.root  gxspl-freenuc.xml.gz  gxsp_cut_emax
    neutrino:  gxspl-NUbig.xml.gz  gxspl-NUsmall.xml  reduce_gxspl.awk
               gxspl-FNALbig.xml.gz -> gxspl-NUbig.xml.gz   (old names)
               gxspl-FNALsmall.xml  -> gxspl-NUsmall.xml
    electron:  gxspl-ELECTRON.xml
    plus any tune directory / override XML files from cfg/
```

- **Tarball name:** `genie_xsec-3.08.00-noarch-AR2320i00000-k250-e1000.tar.bz2`. This is the only place the version is written with dots (`v3_08_00` becomes `3.08.00`), because the tarball has to be named that way to be processed properly. Everything else, including the directories inside the tarball and the version file, uses `$SVERSION` with underscores.
- **UPS qualifier:** `$TUNESQ:k$KNOTS:e$EMAX`, with `:electron` added for electron productions, or the `--nuflavors` suffix (`:nue:nuebar`, `:notau`, …) for a production with only some neutrino flavors. The directory names use dashes in place of colons (`UPSQUALDASH`).
- **Default spline file:** the table points `GENIEXSECFILE` at `gxspl-NUsmall.xml`, or at `gxspl-ELECTRON.xml` for electrons.
- **`GENIE_XSEC_EMAX`:** written with one decimal place (e.g. `1000.0`).
- **`gxsp_cut_emax`:** copied from `$STAGEDIR/bin` into `data/` (executable), for users who want a spline file that stops at a lower energy; see [[#Shortening the splines `gxsp_cut_emax`|Shortening the splines: gxsp_cut_emax]].
- **`README.md`:** a Markdown description of the product, in `data/`, with a copy put in `$STAGEDIR/work-products/README.md`. It has:
  - a summary table: product and qualifier, tune, generator list, probes (and `--nuflavors`), energy range and knots, isotope counts, QE averaging, and the GENIE version, date and user
  - the `setup` command and the variables the table sets
  - the files, how to make a smaller file with `reduce_gxspl.awk` (neutrinos), how to make a file that stops at a lower energy with `gxsp_cut_emax` (usage, and what changes below the cut)
  - a "TGraphs" section on `xsec_graphs.root`: its isotopes (the `root` ones, or all if none are marked, as in stage 5), the energies (from the tune's `GVLD-Emin` up to stage 5's `-e`; the number of points and whether they're spaced in E or log E, from `--tgraph`), the units (10⁻³⁸ cm²) and the directory and graph names
  - the probes, and the isotopes as a table with one row per element, `| Z | isotopes |`, e.g. `| 18 | Ar36 Ar39 Ar40 |` (those in the big file, and those in `gxspl-NUsmall.xml`)

  A `README.md` already in `work-products/` from an earlier attempt is left as it is.

### Shortening the splines: `gxsp_cut_emax`

`gxsp_cut_emax` makes a spline file that stops at a lower energy: every spline ends at its second knot at or above the cut `ECUT` (GeV), so the file still covers every energy up to `ECUT`. It is in `bin/` (and so in `$STAGEDIR/bin`) and in the UPS tarball's `data/`. Like `reduce_gxspl.awk`, it reads a file (plain or `.gz`) or a pipe, and writes to standard output:

```bash
gxsp_cut_emax 25 gxspl-NUbig.xml.gz > gxspl-NUbig-e25.xml
zcat gxspl-NUbig.xml.gz | gxsp_cut_emax 25 > gxspl-NUbig-e25.xml
gxsp_cut_emax -n 3 25 gxspl-NUbig.xml.gz | gzip --best > gxspl-NUbig-e25.xml.gz
```

```
ECUT         energy in GeV; each spline keeps its knots below it and the first N at or above it
-n N         knots to keep at or above ECUT (default 2)
-q           no summary line on stderr
file or -    input; standard input if no file is given, or for "-"
```

- **Each spline is cut on its own knots.** Splines with a threshold have their 5 knots below threshold on their own grid, so splines keep different numbers of knots. The `nknots="..."` on each `<spline>` line is set to the number kept: GENIE reads exactly that many.
- **Unchanged:** a spline that already ends below `ECUT`, and every line outside the splines. Every spline keeps at least 2 knots, even for a cut below its second knot.
- **Above the cut:** a cut file is only for energies up to `ECUT`. GENIE doesn't extrapolate: `Spline::Evaluate` gives 0 outside a spline's knot range, and it does so silently (its message is compiled in only with GENIE's low-level messages). `GMCJDriver` checks the flux's maximum energy only against `GEVGDriver::ValidEnergyRange()`, which comes from the event generators' configured validity, not from the splines. So with a flux that goes above `ECUT` there is no error or warning: neutrinos above the splines' last knots simply never interact. The UPS table's `GENIE_XSEC_EMAX` is still the full files' energy.
- **Just below the cut:** GENIE builds a cubic (`TSpline3`, not-a-knot ends) through the knots, so moving the last knot changes the cubic slightly in the last few intervals. That is why the default keeps 2 knots at or above `ECUT`: with only 1, the change at 25 GeV was nearly 3 times larger. On every 20th spline of an AR23 production, the largest change below the cut with the default was:

  ```
  cut        interval holding ECUT   next below   next     next     next
  10 GeV          1.7%                 0.50%       0.13%    0.033%   0.009%
  25 GeV          0.39%                0.11%       0.028%   0.007%   0.002%
  125 GeV         0.74%                0.25%       0.065%   0.017%   0.004%
  ```

  In the interval holding `ECUT`, only 27, 9 and 12 of the ~3,700 splines compared changed by more than 0.1% (at 10, 25 and 125 GeV). The largest changes are all in Nieves QES splines, whose knots are a little noisy. With `-n 1` at 25 GeV the largest change was 1.1%; with `-n 3`, 0.09%.
- **Exit status:** 0, or 2 for bad arguments or input that can't be read.

A cut of the AR23 big file takes about 45 s. The sizes, with the default `-n 2`, as `ls -h` gives them:

```
                gxspl-NUbig.xml               gxspl-NUsmall.xml
cut        uncompressed  share  gzipped    uncompressed  share  gzipped
full          1.3G       100%     -           518M       100%     -
10 GeV        646M        49%    65M          256M        49%    26M
25 GeV        780M        60%    85M          309M        60%    34M
125 GeV      1014M        78%   124M          402M        78%    49M
```

## Design

### Language and portability

Everything is bash. Grid nodes always have bash, the old tooling was bash, and Python versions on SL7 nodes are unreliable. The scripts also run under bash 3.2 (as on macOS): there is no `mapfile` and no associative arrays. The only Python is GENIE's own QE smoothing scripts in stage 3s.

### Configuration travels in `cfg.tar.gz`

`gxspinit` writes `$STAGEDIR/cfg/`, which holds:

- `define_cfg.sh`: all parameters (`TOPDIR`/`TOP`, `STAGEDIR`, `SVERSION`, `TUNE`, `TUNESQ`, `KNOTS`, `EMAX`, `GENLIST`, `ELECTRON`, `GVERSION`, `GQUAL`, `UPSQUAL`, `UPSQUALDASH`, ...). It also defines `gxsp_setup_genie` and sources `gxsplib.sh`.
- `gxsplib.sh`: the shared function library.
- `isotopes.cfg`
- `reduce_gxspl.awk` (neutrinos only)
- the optional tune directory and override XML files
- `cfg.tar.gz`: all of the above
- `genie_splines.dag`: kept outside the tarball

Every grid job fetches `cfg.tar.gz`, unpacks it in its scratch directory, and sources `define_cfg.sh`. That way all jobs use exactly the configuration chosen at initialization, and the shared library goes with it.

`gxsp_setup_genie` checks `/etc/system-release` for `Scientific Linux release 7.9`:

- **SL7:** sources `bootstrap_many_ups.sh` and `setup_genie_he_data.sh`, then runs `setup genie $GVERSION -q $GQUAL` (and `setup ifdhc` if `ifdh` isn't already available).
- **Anything else** (e.g. AlmaLinux 9): gets `ifdh` from spack:
  ```bash
  source /cvmfs/fermilab.opensciencegrid.org/packages/common/spack/current/NULL/share/spack/setup-env.sh
  spack load fife-utils@3.7.8%gcc@11.4.1
  ```
  Setting up GENIE itself there is still to be done.

After GENIE is set up, the job's working directory (where `cfg.tar.gz` was unrolled) is always put first on `GXMLPATH`. This happens whether or not a tune directory or override XML files were shipped, and it isn't added again if it's already there. The job prints the resulting `GXMLPATH`. Jobs that skip GENIE setup (stage 6, or `GXSP_NOSETUP=1`) don't change `GXMLPATH`.

### How GENIE finds XML files

GENIE always looks for an XML configuration file in this order:

1. directories in `$GXMLPATH`, in order of appearance
2. directories in `$GXMLPATHS` (deprecated)
3. directories in `$GALCONF` (deprecated)
4. the fully specified tune directory (`<CMC>/<TUNE>/`), except for `_00_000` tunes
5. the CMC directory (`<CMC>/`)
6. `GetXMLDefaultPath()`, i.e. `$GENIE/config`
7. `$GENIE_REWEIGHT/config`

What this means for the scripts:

- **Working-directory files win:** the job's working directory is the first `$GXMLPATH` entry, so any file there beats every tune, CMC and default copy. That covers the override XML files from `gxspinit` and the edited `NewQELXSec.xml` in the seeded QE jobs.
- **`_00_000` tunes:** these never use a tune-specific directory, so for them only the CMC level matters (as in the `NewQELXSec.xml` lookup).

### Start-up code pasted into each script

jobsub copies only the script file named after `file://` to the worker node. Before anything else, each job therefore has to make `ifdh` available and fetch `cfg.tar.gz`. That code lives once in `gxsp_bootstrap.sh`. Each stage script has a line

```bash
source "$(dirname "${BASH_SOURCE[0]}")/gxsp_bootstrap.sh"   # GXSP_INLINE
```

and when `gxspinit` copies the scripts to `$STAGEDIR/bin`, it replaces that line with the file's contents. The original copies still work when run directly, because they source the file next to them.

`gxsp_bootstrap` does the following:

1. Change to `$_CONDOR_SCRATCH_DIR`. When a script is run by hand without it, use a new temporary directory, so nothing is unpacked into your current directory.
2. Choose how to reach PNFS: `ifdh` when `GRID_USER` is set, otherwise plain Unix commands (`GXSP_IO` overrides this).
3. Set up `ifdh` if it's needed: on SL7 with UPS (`setup ifdhc`), elsewhere with spack (`spack load fife-utils@3.7.8%gcc@11.4.1`, as above).
4. Fetch `cfg.tar.gz`, unpack it, and source `define_cfg.sh`.
5. Call `gxsp_setup_genie`.

### Access to the staging area (PNFS)

Two PNFS rules apply everywhere, including during initialization when plain `ls`, `rm` and `cp` work:

- **No overwriting:** a file can't be copied over to update it; the existing copy has to be removed first.
- **No deleting from grid jobs:** a job's token may allow creating files in PNFS but not removing them, so grid jobs never delete or replace anything. Interactively you can still `rm` and copy files.
- **No `mv`:** use `cp` followed by `rm` instead.

`gxspinit` builds everything (the `bin/` scripts with the bootstrap inlined, and all the `cfg/` contents) in a local scratch directory. It then copies each item in with `stage_install`, which removes any existing file or directory of that name first. The `mv` calls in the grid scripts only rename files in the job's local scratch area, never in PNFS.

On the grid, all reads and writes go through a few helpers in `gxsplib.sh`, so the `ifdh` rules are handled in one place:

- **`gxsp_listing`:** runs `ifdh ls $STAGEDIR 99` (or `find` when PNFS is mounted) once per job. (Don't add `--debug=2`: with it, the spack `fife-utils` `ifdh` on AlmaLinux 9 aborts.) The output is cleaned up into one absolute path per line, with no doubled or trailing `/`, and cached in `stagedir.ls`. If the listing command fails (non-zero exit status), or the listing comes back empty, the job's stderr shows why: a warning with the command, the access mode and the exit status, then the command's whole stderr between two marker lines:
  ```
  ==================== BEGIN stderr of: ifdh ls /pnfs/.../GXSPL-... 99 ====================
  ...
  ==================== END stderr of: ifdh ls /pnfs/.../GXSPL-... 99 ======================
  ```
  That information would otherwise be lost with the scratch area. A listing whose command failed but still printed paths is used as it is.
- **`gxsp_have` / `gxsp_missing`:** look up exact paths in that cached listing. No repeated `ifdh ls` calls.
- **`gxsp_check_inputs` / `gxsp_require_inputs`:** the input check of stages 2, 4, 3s, 4q and 6. PNFS can be slow to show the files the previous DAG step has just written, and a listing can fail now and then. So if anything is missing, or no listing can be made, they take a fresh listing and check again, up to `GXSP_RETRIES` times (default 3), `GXSP_RETRY_WAIT` seconds apart (default 60). Only then do they fail: exit 17, listing what's still missing (e.g. `1 of 12 stage 1 files missing`), or exit 12 if no listing ever worked. They run in the job's main shell, so a failed listing can't slip past the check.
- **`gxsp_fetch REL [local]`:** checks the listing, then copies with full file names on both sides. A file not in the listing gets the same re-listings first, which covers the stage 3 and 5 jobs fetching what stages 2 and 4 just wrote. A file still missing after that is fatal. A failed copy is tried again, up to `GXSP_RETRIES` more times, `GXSP_RETRY_WAIT` seconds apart, with any partial local file removed first. That matters because `ifdh cp` (via `gfal-copy`) never overwrites a destination: a transfer that fails part way leaves a partial file, and `ifdh`'s own retry then fails with "File exists … overwrite is not set". Stage 3s makes thousands of copies per job, so a transient transfer error is likely sooner or later. After the last try the job exits 13 (`copy failed N times: REL`).
- **`gxsp_put LOCAL REL`:** copies a new file in. It never deletes or overwrites: if `REL` already exists it's left as it is, with a warning.
- **`gxsp_put_log LOCAL REL`:** for logs and `-NAN` / `-NEG` files. If `REL` exists from an earlier attempt, it uses the first free `REL.1`, `REL.2`, ….
- **`gxsp_skip_if_done REL ...`:** at the start of stages 2, 4, 3s, 4q, 5 and 6 (and per seed in `gxsp_mkspl`). If the job's products are all there already, it says so and ends successfully. Stage 5 also redoes only whichever of its two products is missing.

So the way to redo a product is to remove it interactively and re-run (or resubmit) its job.

The end-of-job sleep (`GXSP_SLEEP`, 25 s) also gives PNFS time to show a job's outputs before the next DAG step starts.

### Mapping `$PROCESS` to a pair

The order of stage-3 jobs is isotope-major and probe-minor:

```
PROCESS = isotope_index * NPROBE + probe_index        (both 0-based)
```

Isotopes are taken in the order they appear in `isotopes.cfg`: all active (uncommented) lines for neutrinos, only lines containing `electron` for electrons, leaving out the free nucleons (see below). Probes are ordered `12 -12 14 -14 16 -16`, or just `11`; with `--nuflavors`, only the chosen ones, in that same order (so `NPROBE` is the number of flavors). `gxsp_pair_from_index`, `gxsp_expected` and `gxsp_status` all use this same order, which is why the `PROCESS=` numbers in the status report are correct. With the default list, `PROCESS=0` is `gxspl-nue-H2.xml` and `749` is `gxspl-nutaubar-Pb208.xml`.

With `--pairs qez` or `--pairs other` (the `--average-qe N1,N2` clusters of the extra seeds), the same formula applies to a shorter isotope list: only the isotopes whose Z is in `QEZ` (6, 8 or 18 by default), or only the others, still in `isotopes.cfg` order. With the default list and `--qe-z`, `--pairs qez` has 9 isotopes (54 pairs), so `PROCESS=0` there is `nue` on the first carbon isotope.

### Free nucleons come from stage 1

The free neutron (`free-n`, 1000000010) and free proton (`H1`, 1000010010) are the stage 1 targets. Run again in stage 3 with `gxspl-freenuc.xml` as input, they only reproduce the stage 1 splines. In the AR23 v3_06_00 production, all 619 free-n and H1 splines in `gxspl-NUbig.xml.gz` are identical, knot for knot, to those in `gxspl-freenuc.xml.gz`. So:

- **Stage 3 and 3q:** no jobs on free-n or H1.
- **Stage 4:** starts the big file from `gxspl-freenuc.xml`, then adds the stage-3 isotopes.
- **Stage 3s:** doesn't average free-n or H1, and the seeded QE runs (stage 3q) don't read `gxspl-freenuc.xml`.
- **Everything else:** `isotopes.cfg` keeps its free-n and H1 lines, and they're used everywhere else. The reduction script still keeps them in `gxspl-NUsmall.xml` (both are marked `reduced`), the TGraph file still has them (`root`), and the README lists them.

Two consequences:

- **Free nucleons always in the big file:** the big file includes whatever `gxspl-freenuc.xml` has, free-n and H1 for every probe, even if `isotopes.cfg` leaves one of them out. For electrons this means the big file has the free-neutron splines, even though `free-n` isn't marked `electron` in `isotopes.cfg`. Those splines aren't in `README.md`'s isotope list, since that comes from `isotopes.cfg`.
- **Free nucleons only:** an isotope list with nothing but the free nucleons gives a DAG without stage 3, and a big file that is just `gxspl-freenuc.xml`.

> [!note] Editing `isotopes.cfg` after initialization
> Editing `cfg/isotopes.cfg` after initialization changes this mapping. Re-run `gxspinit --force` so that the DAG's `-N` and `cfg.tar.gz` stay consistent.

### Naming

File names use the probe and target names, never PDG codes: `gxspl-nue-freep.xml`, `gxspl-numubar-Ar40.xml`, `gxspl-Ar40.xml`. The free nucleons are `freep` (1000010010) and `freen` (1000000010) in the stage 1 file names. Other targets take their name from the second column of `isotopes.cfg`, in stage 3 and in the stage 4 per-isotope sums (`gxspl-Ar40.xml`).

### Combining files with gspladd

`gspladd` fails if its `-f` list is too long. `gxsp_spladd OUT files...` builds the comma-separated list, and if it's longer than `GXSP_MAXARGLEN` (default 20000 bytes), it:

1. combines the files in groups that fit,
2. combines the partial results the same way (recursing if needed), and
3. deletes the temporary `gxsp_part*.xml` files.

A single file is just copied. With the default limit, a normal production never needs the grouping: stage 4's list of 125 isotope files plus `gxspl-freenuc.xml` is about 1800 bytes. Set `GXSP_MAXARGLEN` lower if `gspladd` turns out to have a shorter limit.

### NaN, inf and negative values

A cross section should never be NaN, inf or negative. Every spline file a job makes is cleaned (`gxsp_scrub` in `gxsplib.sh`) before it's copied to the staging area:

- **NaN and inf:** always replaced by `0.0`.
- **Negative values:** replaced by `0.0` too, unless the production was set up with `gxspinit --keep-negative` (recorded as `NEGZERO` in `define_cfg.sh`: 1 by default, 0 with the option). `-0.0` is zero, not negative, and is left as it is. GENIE itself doesn't set a negative knot to zero: it treats it as non-zero and interpolates through it.
- **The original file** is uploaded next to the cleaned one: as `<file>-NAN` if it had any NaN or inf, otherwise as `<file>-NEG`. The job's log says how many values of each kind were replaced.

This happens in two places:
- **Stages 1, 3 and 3q:** `gxsp_mkspl`, right after `gmkspl`, for every output, including each seed of a seed range.
- **Stage 3s:** each smooth file made by `construct_representative_spline.py`, before it's copied back for stage 4q to put into the big file (`D-smooth.xml-NEG` holds the original).

Everything else, the big file included, is built from these cleaned files.

Only values inside `<xsec>` tags are looked at, so names like "Inelastic" aren't affected. `gxsp_status` lists any `-NAN` and `-NEG` files at the end of its report. Since jobs don't delete, such a file from an earlier attempt stays until you remove it, even if a later run of that job was clean.

### The reduction awk script

`reduce_gxspl.awk` lists every probe (`nukeep[...] = 1`) and every active isotope (`tgtkeep[...] = 1` if the line has `reduced`, otherwise `0`). The keep test is `tgtkeep[tgt] == 1 && nukeep[nu] == 1`, so entries set to `0` are dropped. With the default list it keeps 51 isotopes.

### Job resources (jobsub options)

The jobsub options for every job in the DAG are set in one block of shell variables near the end of `gxspinit`, just above where it writes `genie_splines.dag` (the lines starting `JSNOSING=`). `gxspinit` builds each line of the DAG from them; nothing else sets them.

| Variable | Setting | Used by |
|------------|-------------------------------|-----------------------|
| `JSNOSING` | `jobsub_submit`, group, usage model, CVMFS requirement (below) | every job |
| `JS` | `JSNOSING` plus the SL7 singularity image (below) | every stage except 3s and 4q |
| `BIGMEM` | `--memory=4095MB` | stages 4, 3s, 4q, 5, 6 |
| `SHORTT` | `--expected-lifetime 3h` | stages 2, 4, 4q, 5, 6 |
| `LONGT` | `--expected-lifetime 85200s` (about 23.7 h) | stage 1 |
| `SUPERLONGT` | `--expected-lifetime 170400s` (about 47.3 h) | stage 3 (Default) |
| `QET` | `--expected-lifetime 15h` (from `QETH=15`) | stage 3q (the QE seed clusters) |
| `QEAVGT` | `--expected-lifetime 12h` | stage 3s |
| `NORMALD` | `--disk 2GB` | stages 1, 2, 3, 3q |
| `BIGD` | `--disk 10GB` | stages 4, 4q, 5, 6 |
| `QEAVGD` | `--disk 30GB` | stage 3s: one pair's seed files at a time, plus the venv. This allows for 1000 seed files of about 11 MB, their size if they were made with `--input gxspl-freenuc.xml`; otherwise a seed file is a few tens of KB |

The full options in `JSNOSING` and `JS`:

```
JSNOSING  jobsub_submit -n --group ${JOBSUB_GROUP}
          --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE
          --append_condor_requirements='(TARGET.HAS_CVMFS_dune_opensciencegrid_org==true)'
JS        ${JSNOSING}
          --singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
```

The group comes from `gxspinit --jobsub-group` (default `genie`). Stages 1 to 3q don't set `--memory`, so they get jobsub's default.

To change a setting:

- **For new productions:** edit the variable in `gxspinit`. The next `gxspinit` run writes a DAG with the new value; `gxspinit --force` rewrites `cfg/genie_splines.dag` in an existing staging area. The QE lifetime is set in hours by `QETH` (`QET` is built from it), and `gxspinit`'s warning that the QE job time estimate is over the expected lifetime compares against the same `QETH`, so change only `QETH`. The QE job time estimate itself is explained in [[#QE averaging (`--average-qe`)]].
- **For one staging area only:** edit `$STAGEDIR/cfg/genie_splines.dag` before submitting it, or change the options on the line you use to resubmit a single stage. A later `gxspinit --force` replaces the edited DAG.

### Exit codes

The grid jobs' own codes are 10 to 24. That keeps them apart from what a GENIE program (`gmkspl`, `gspl2root`) can return and the job then passes on (1 to 5, 60 to 80, and 255 for -1), from the shell's 126 and 127 (command not executable / not found), and from 128 + N, a process killed by signal N (137 for `SIGKILL`, 143 for `SIGTERM`).

| Code | Meaning |
|----|------------------------------------|
| 10 | bad arguments (including an option given without its value), or cfg fetch/unpack failed |
| 11 | configuration problem (no isotopes selected, `$PROCESS` out of range for the `--pairs` list) |
| 12 | could not list the staging area (after the retries, where the stage has them) |
| 13 | required input missing, or copying it in failed (after the retries) |
| 14 | GENIE setup failed |
| 15 | copying output back failed (including when it already exists: jobs don't replace files) |
| 16 | `gmkspl` exited 0 but produced no file |
| 17 | inputs from an earlier stage still missing after the retries (stages 2, 4, 3s, 4q, 6) |
| 18 | `gspladd` failed (stages 2 and 4) |
| 19 | `gawk -f reduce_gxspl.awk` (stage 5) or `tar` (stage 6) failed |
| 20 | stage 3s: a smoothing script failed or didn't make its expected files; stage 4q: a smooth spline has no spline of the same name in the big file |
| 21 | compressing or uncompressing the big file failed (stages 4, 4q, 5, 6) |
| 22 | seeded QE job: no `NewQELXSec.xml` found, or no `NumNucleonThrows` in it |
| 23 | stage 3s: the smoothing scripts couldn't be copied from `QESCRIPTS` (`--qe-scripts`) |
| 24 | stage 3s: building the python venv (`python3 -m venv` or `pip install`) failed |
| other | the status of `gmkspl` (stages 1, 3, 3q) or `gspl2root` (stage 5), passed on |

The interactive tools (`gxspinit`, `gxsp_status`, the spline checks) aren't jobs and keep their own exit statuses: `gxspinit` exits 1 on an error, for example, and the spline checks exit with the number of splines they found.

Stage 3s goes on past a pair that fails and exits at the end with the first failed pair's status (13, 15 or 20).

On failure, a job prints its log to stdout so the log shows up in the job output.

### Environment variables

These change how the scripts behave; they are mostly useful for running stages by hand or off the grid.

| Variable | Effect |
|-------------|-----------------------------|
| `GXSP_IO=ifdh` or `local` | force how PNFS is reached (default: `ifdh` if `GRID_USER` is set) |
| `GXSP_NOSETUP=1` | skip `gxsp_setup_genie` |
| `GXSP_SLEEP=N` | end-of-job sleep (default 25 s) |
| `GXSP_MAXARGLEN=N` | gspladd `-f` length limit in bytes (default 20000) |
| `GXSP_LISTING=path` | where the cached staging-area listing goes (default `./stagedir.ls`) |
| `GXSP_RETRIES=N` | how many times to re-list and re-check missing inputs, or retry a failed copy from the staging area, before failing (default 3) |
| `GXSP_RETRY_WAIT=N` | seconds between those re-checks or retries (default 60) |
| `GXSP_PYTHON=cmd` | python used for the stage 3s smoothing scripts (default `python3`); when set, no venv is built |
| `GXSP_VENV=0` | don't build stage 3s's python venv (default: build it, on any node) |

## Other Considerations

### Scratch space when running interactively

On an interactive node `/tmp` is often small and can fill up, so set `TMPDIR` to somewhere with room before running the scripts by hand:

```bash
export TMPDIR=/path/with/space/tmp     # any directory with room; it must exist
mkdir -p $TMPDIR
```

Every scratch area the scripts make uses `$TMPDIR`, and falls back to `/tmp` only when `TMPDIR` isn't set:

- **`gxspinit`** builds `bin/` and `cfg/` in `$TMPDIR/gxspinit.XXXXXX`, prints the path (`scratch directory ...`), and removes it when it exits.
- **A stage script run by hand** (no `$_CONDOR_SCRATCH_DIR`) works in `$TMPDIR/<script>.XXXXXX` and prints the path (`scratch directory ... removed on success, kept on failure`). When the job succeeds (including when it has nothing to do), the directory is removed as it exits. When it fails, the directory is kept so you can look at what went wrong, and the last line says so: `exit status 17: scratch directory kept: ...`. Remove a kept one yourself when you're done: from stage 4 on it can hold big files (for the AR23 k250 e1000 production, `gxspl-NUbig.xml.gz` is 181 MB and uncompressed 1.4 GB).
- **`gxsp_status`** keeps its staging-area listing in a temporary file in `$TMPDIR`, removed when it exits.

If `TMPDIR` names a directory that doesn't exist (or can't be written), the script stops with `can't create a scratch directory in ...` and doesn't fall back to `/tmp` or the current directory.

On a grid node none of this applies: the jobs work in `$_CONDOR_SCRATCH_DIR`, which the batch system cleans up.

### Running a stage by hand doesn't use ifdh

How a script reaches the staging area depends on where it runs:

- **On a grid node** (`GRID_USER` set): every access goes through `ifdh`. The staging area is listed with `ifdh ls $STAGEDIR 99`, and files are copied with `ifdh cp`, with `ifdh mkdir` for a missing directory.
- **On an interactive node** (`GRID_USER` not set), where `/pnfs` is mounted: plain Unix commands. The listing comes from `find $STAGEDIR`, and files are copied with `cp` and `mkdir`. `ifdh` isn't used for any of it.

So a stage that works when run by hand hasn't necessarily been tested with `ifdh`. It also makes running a stage by hand a useful workaround when `ifdh` is the problem.

To make a stage run by hand use `ifdh`, as it would on the grid, set `GXSP_IO=ifdh`. `GXSP_IO=local` does the opposite. With `GXSP_IO=ifdh` you need `ifdh` set up first: on SL7 with `setup ifdhc`, elsewhere as below.

```bash
source /cvmfs/fermilab.opensciencegrid.org/packages/common/spack/current/NULL/share/spack/setup-env.sh
spack load fife-utils@3.7.8%gcc@11.4.1
GXSP_IO=ifdh $STAGEDIR/bin/gxsp_qe_smoothing --stagedir $STAGEDIR -p numu -t C12
```

The rules for the staging area are the same either way: a job never deletes or replaces anything there. Interactively you can still remove files yourself with `rm`.

## Adapting to another site

The scripts were written for Fermilab:

- GENIE comes from UPS on CVMFS;
- the staging area is in PNFS and reached from grid jobs with `ifdh`;
- the jobs are submitted as a jobsub DAG.

The production logic doesn't depend on any of that: the stages, file names, configuration, checks and `gxsp_status`. The Fermilab-specific parts are confined to a few places, listed below, so another site needs to change those and not much else. A good first step is to run a small production by hand on one machine with `GXSP_IO=local` (see [[#Running a stage by hand doesn't use ifdh]]), with an isotope list of a few lines (`--fetch-isotopes`), and then move to the batch system.

### What a site has to provide

- **bash** (3.2 or later), plus the usual tools `awk`, `gawk` (stage 5's reduction script), `gzip`, `bzip2`, `tar`, `sort`, `find`.
- **GENIE** with `gmkspl`, `gspladd` and `gspl2root` on the `PATH` after setup, and the environment GENIE needs (`GENIE`, ROOT, LHAPDF and so on).
- **A staging area** that interactive nodes and every batch job can reach, for reading and writing.
- **A batch system** that can run arrays of jobs (`-N` with `$PROCESS` = 0 … N−1) and make one group of jobs wait for another.
- **With `--average-qe` only:** `python3` (3.9 or later) for stage 3s, with internet access to PyPI from the batch nodes, or a ready python environment (below), and GENIE's QE smoothing scripts (`src/contrib/ccqe_spline_smoothing/` in the GENIE source).

### Setting up GENIE

GENIE is set up in one function, `gxsp_setup_genie`, which `gxspinit` writes into `define_cfg.sh` (look for `cat >> "${WORK}/define_cfg.sh"` in `gxspinit`). Every job calls it during start-up, unless `GXSP_NOSETUP=1`. At Fermilab, on SL7, it does this:

```bash
source /cvmfs/fermilab.opensciencegrid.org/products/genie/bootstrap_many_ups.sh
source /cvmfs/fermilab.opensciencegrid.org/products/genie/externals/pochoarus-genie_he_data/setup_genie_he_data.sh
setup genie "${GVERSION}" -q "${GQUAL}"
```

Elsewhere, replace those lines with the site's own setup, whatever form it takes: `spack load genie@…`, `module load …`, sourcing the environment script of a local GENIE build, or nothing if the jobs run in a container that already has GENIE. `--gversion` and `--gqual` are UPS terms, but they're only stored (`GVERSION`, `GQUAL`) and shown in logs, the README and the UPS table. So they can carry, for example, a spack spec or a build tag.

Keep the rest of the function as it is: it puts the working directory (where `cfg.tar.gz` is unpacked) first on `GXMLPATH`, which the shipped tune directory and the seeded QE jobs' `NewQELXSec.xml` rely on (see [[#How GENIE finds XML files]]). The function must return success only if `gmkspl` can be found. The `gxsp_is_sl7` test (from `/etc/system-release`) that picks the UPS path is also Fermilab's; a site with one kind of node can drop it.

### The staging area and file access

At Fermilab the staging area is in PNFS. Batch jobs can't see PNFS as a file system, so they reach it with `ifdh`, while interactive nodes use plain Unix commands (see [[#Access to the staging area (PNFS)]]). Every access goes through a few functions in `gxsplib.sh`:

```
gxsp_io_mode       which way: GXSP_IO=ifdh|local, else ifdh if GRID_USER is set (a Fermilab grid job), else local
gxsp_listing_try   one listing of the whole staging area: ifdh ls $STAGEDIR 99, or find $STAGEDIR
gxsp_cp            copy one file in or out: ifdh cp, or cp
gxsp_mkdir         make a directory: ifdh mkdir, or mkdir -p
```

There's one more copy in `gxsp_bootstrap.sh`: the fetch of `cfg/cfg.tar.gz` at start-up, before `gxsplib.sh` is available. The same file also sets up `ifdh` when it isn't there, from CVMFS (UPS `ifdhc` on SL7, spack `fife-utils` elsewhere).

- **A shared file system** (NFS, Lustre, GPFS, a mounted EOS, …) seen by both interactive and batch nodes is the easy case: set `GXSP_IO=local` in every job's environment, and the jobs use `cp`, `mkdir` and `find`. Nothing else needs to change. `GRID_USER` is only set on Fermilab's grid, so `local` is already the default elsewhere, but setting `GXSP_IO` explicitly is safer.
- **Remote storage only** (xrootd, gfal, S3, …): add a third mode to the four functions above and to the fetch in `gxsp_bootstrap.sh`. The listing must print every file and directory below the staging area, one full path per line (directories may end in `/`); the functions tidy up repeated and trailing slashes. Copies go both ways, between a local file and a path in the staging area.

Some of the design follows from PNFS and can be relaxed:

- **No deleting or replacing:** jobs never delete or overwrite anything in the staging area (Fermilab's grid tokens allowed creating files but not removing them). That's also what makes re-running safe, so it's worth keeping even where deleting is allowed.
- **Waiting for files to show up:** PNFS can be slow to list a file another job has just written. So each job sleeps `GXSP_SLEEP` (25 s) before it exits, and input checks re-list up to `GXSP_RETRIES` (3) times, `GXSP_RETRY_WAIT` (60 s) apart. On a file system that shows new files at once, `GXSP_SLEEP=0` and a short `GXSP_RETRY_WAIT` are fine (see [[#Environment variables]]).
- **Paths:** the example `--top` (`/pnfs/genie/scratch/users/rhatcher/GXSP`) and the default `--qe-scripts` directory are Fermilab paths. Give your own with `--top` and `--qe-scripts`.

### Without CVMFS

These are the scripts' only uses of `/cvmfs`:

```
/cvmfs/fermilab.opensciencegrid.org/products/genie/bootstrap_many_ups.sh
        UPS itself, for "setup genie" and "setup ifdhc" (gxsp_setup_genie, gxsp_bootstrap.sh)
/cvmfs/fermilab.opensciencegrid.org/products/genie/externals/pochoarus-genie_he_data/setup_genie_he_data.sh
        GENIE's extra data package genie_he_data, set up along with GENIE (gxsp_setup_genie)
/cvmfs/fermilab.opensciencegrid.org/packages/common/spack/current/NULL/share/spack/setup-env.sh
        spack, for "spack load fife-utils" (ifdh) on non-SL7 nodes
/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
        the SL7 container image, in the DAG's --singularity-image
```

The DAG also asks for nodes with the DUNE CVMFS repository mounted:

```
--append_condor_requirements='(TARGET.HAS_CVMFS_dune_opensciencegrid_org==true)'
```

Without CVMFS, all of these come out with the changes above:

- **GENIE and its data** come from the site's own installation (`gxsp_setup_genie`).
- **`ifdh`** isn't needed with a shared file system (`GXSP_IO=local`).
- **The SL7 container** is only there because Fermilab's GENIE builds need SL7. Run the GENIE stages wherever the site's GENIE runs.

Stages 3s and 4q don't use GENIE at all, only bash, awk, gzip and (stage 3s) python.

**The python environment for stage 3s:** each stage 3s job builds a venv with `pip` (numpy, pandas 2.2.3, fastparquet, tqdm, colorama, scikit-learn; see [[#QE averaging (`--average-qe`)]]). If the batch nodes can't reach PyPI, build that environment once in a place the jobs can see, and point `GXSP_PYTHON` at its `python3`. The jobs then use it as is and skip the venv.

### The DAG and the batch system

`gxspinit` writes `cfg/genie_splines.dag` for Fermilab's `jobsub_submit_dag`, a layer over HTCondor's DAGMan:

- `<serial>` and `<parallel>` blocks;
- one `jobsub_submit` command per line, with Fermilab options: `--group`, `--resource-provides=usage_model=…`, `--append_condor_requirements`, `--singularity-image`, and `-N` for job arrays.

It won't run anywhere else, but it is a complete, readable list of every job with its exact arguments, so it's the place to start from. Make one with the configuration you want, even on a machine without jobsub, and translate it.

What the scripts need from the batch system:

- **One script per job.** Each job runs one script from `$STAGEDIR/bin` with the arguments on its DAG line. `gxspinit` pastes the start-up code into each copy (see [[#Start-up code pasted into each script]]). So a job needs only that one file: everything else (`cfg.tar.gz`, inputs) it fetches from the staging area.
- **`$PROCESS` in job arrays.** The stage 3 and 3q lines are arrays (`-N n`), and each job finds its pair from `$PROCESS` = 0 … n−1 (see [[#Mapping `$PROCESS` to a pair]]). Stage 3s jobs take their share of pairs from `--job K` on their own DAG lines, so they need no array index (they also accept `$PROCESS` in place of `--job`). On a system that removes the rest of an array when one task fails, as DAGMan does with a `-N` cluster, separate jobs are the safer choice for long stages. On another system, set `PROCESS` from its own index in a small wrapper, e.g. `PROCESS=$SLURM_ARRAY_TASK_ID` for Slurm (array indices starting at 0), or `PROCESS=$(( SGE_TASK_ID - 1 ))` for Grid Engine.
- **A scratch directory.** The jobs work in `$_CONDOR_SCRATCH_DIR` (HTCondor's per-job directory). If that isn't set, they make a fresh directory in `$TMPDIR` (else `/tmp`), removed when the job succeeds and kept when it fails (see [[#Scratch space when running interactively]]). Point `TMPDIR` at the node's local disk, or set `_CONDOR_SCRATCH_DIR` in the wrapper. Stage 3s needs up to 30 GB there (see [[#Job resources (jobsub options)]]).
- **The order.** Each stage waits for the whole of the one before. The DAG layout is in [[#The stages]]; with `--average-qe`, stage 3s can run alongside stage 4, and stage 4q waits for both. Any dependency mechanism will do: Slurm's `--dependency=afterok:…`, HTCondor DAGMan directly, a workflow tool, or submitting each stage by hand once `gxsp_status` shows the one before complete. Since a job whose products exist ends at once, resubmitting a whole stage is always safe.
- **Resources.** Memory, disk and run time per stage are in [[#Job resources (jobsub options)]], as set for Fermilab's nodes. They're a starting point; the times depend on the CPUs. The Default stage 3 jobs have the longest limit (`--expected-lifetime` of about 47 h); a seeded QE job takes about 54 s per seed with one throw.

The DAG is written near the end of `gxspinit` (the block after `genie_splines.dag`). A site could change that block to write its own submission file instead. The job counts it uses come from `gxsp_njobs`, `gxsp_qe_clusters` and `gxsp_qe_smooth_njobs` in `gxsplib.sh`.

For a small production on a single machine, the mock DAG runner in the `test/` directory (`test/run_mock_dag.sh STAGEDIR`) runs a generated DAG line by line, expanding `-N n` into `PROCESS=0 … n−1`. It's written for the tests: it puts stand-in GENIE programs first on the `PATH` and skips the GENIE setup (`GXSP_NOSETUP=1`). A copy without those lines would run the real programs one job after another. That's slow, but enough for a few isotopes.

### The UPS tarball

Stage 6 (`gxsp_mkups`) packs the results as a Fermilab UPS product, `genie_xsec`: a version directory, a `.table` file, and the spline files and README under `data/`. A site that doesn't use UPS can skip stage 6 and take the files straight from `work-products/`:

- `gxspl-NUbig.xml.gz` (or `gxspl-ELECTRON.xml.gz`);
- `gxspl-NUsmall.xml`;
- `gxspl-freenuc.xml`;
- `xsec_graphs.root`;
- `README.md`.

Or it can use the tarball's `data/` directory, which doesn't depend on UPS.

## Open items

- [ ] Set up GENIE on non-SL7 nodes in `gxsp_setup_genie` (still a placeholder). `ifdh` there comes from spack's `fife-utils` (`spack load fife-utils@3.7.8%gcc@11.4.1`), as for stages 3s and 4q.
