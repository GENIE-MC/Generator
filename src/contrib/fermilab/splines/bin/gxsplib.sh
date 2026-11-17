# gxsplib.sh -- shared functions for the GENIE spline production scripts
#
# Intended to be sourced by bash.  Most functions expect the variables
# defined by define_cfg.sh (STAGEDIR, ELECTRON, ISOFILE, ...).
#
# Conventions
#   * "rel" paths are relative to $STAGEDIR (e.g. work-products/gxspl-freenuc.xml)
#   * job index ordering (stage 1 and 3) is target-major, probe-minor:
#        idx = itarget * NPROBE + iprobe       (both 0-based)
#
##############################################################################

# (sourced in an interactive login shell $0 is e.g. "-bash": drop the "-")
GXSP_ME=${GXSP_ME:-$(basename -- "${0#-}")}

GXSP_FREENUC_NAMES=(freep freen)
GXSP_FREENUC_PDGS=(1000010010 1000000010)

##############################################################################
# messages
##############################################################################
gxsp_msg()  { echo "[$(date '+%F %T')] ${GXSP_ME}: $*" ; }
gxsp_warn() { echo "[$(date '+%F %T')] ${GXSP_ME}: WARNING: $*" >&2 ; }
gxsp_die()  {
  local rc=$1 ; shift
  echo "[$(date '+%F %T')] ${GXSP_ME}: ERROR: $*" >&2
  exit "${rc}"
}

##############################################################################
# probes
##############################################################################
# gxsp_probe_pdg  <name|pdg>   -> pdg code   (return 1 if unknown)
gxsp_probe_pdg() {
  case "$1" in
    11|electron) echo 11  ;;  -11|positron) echo -11 ;;
    12|nue)      echo 12  ;;  -12|nuebar)   echo -12 ;;
    14|numu)     echo 14  ;;  -14|numubar)  echo -14 ;;
    16|nutau)    echo 16  ;;  -16|nutaubar) echo -16 ;;
    *) return 1 ;;
  esac
}
# gxsp_probe_name <name|pdg>   -> name       (return 1 if unknown)
gxsp_probe_name() {
  case "$1" in
    11|electron) echo electron ;;  -11|positron) echo positron ;;
    12|nue)      echo nue      ;;  -12|nuebar)   echo nuebar   ;;
    14|numu)     echo numu     ;;  -14|numubar)  echo numubar  ;;
    16|nutau)    echo nutau    ;;  -16|nutaubar) echo nutaubar ;;
    *) return 1 ;;
  esac
}
# all six neutrino flavors, in the order used everywhere (job index etc.)
GXSP_NUFLAVORS_ALL="12 -12 14 -14 16 -16"
# gxsp_probe_pdgs  -> space separated list of pdg codes for this production:
#   11 for electrons; otherwise NUFLAVORS (gxspinit --nuflavors), or all six
gxsp_probe_pdgs() {
  if [ "${ELECTRON:-0}" -eq 1 ]; then echo "11"
  elif [ -n "${NUFLAVORS:-}" ];  then echo "${NUFLAVORS}"
  else                               echo "${GXSP_NUFLAVORS_ALL}"
  fi
}
# gxsp_nuflavors_parse <comma list>  -> the pdg codes, in the standard order
#   (12 -12 14 -14 16 -16) whatever order they were given in.  Takes names
#   or pdg codes (nue / 12, nuebar / -12, ...).  Returns 1 (after saying
#   why) for an empty entry, something that isn't a neutrino, or a duplicate.
gxsp_nuflavors_parse() {
  local list=$1 f pdg seen=" " out="" p
  [[ "${list}" =~ ^[^,]+(,[^,]+)*$ ]] || { echo "empty entry in '${list}'" >&2 ; return 1 ; }
  for f in ${list//,/ } ; do
    pdg=$(gxsp_probe_pdg "${f}") || { echo "unknown neutrino flavor '${f}'" >&2 ; return 1 ; }
    case "${pdg}" in 11|-11) echo "'${f}' is not a neutrino flavor" >&2 ; return 1 ;; esac
    case "${seen}" in *" ${pdg} "*) echo "'${f}' given twice" >&2 ; return 1 ;; esac
    seen="${seen}${pdg} "
  done
  for p in ${GXSP_NUFLAVORS_ALL} ; do
    case "${seen}" in *" ${p} "*) out="${out:+${out} }${p}" ;; esac
  done
  echo "${out}"
}
# gxsp_nuflavors_tag <pdg list> <sep>  -> what is added to the staging area
#   and UPS qualifier names for a production with only some neutrino flavors:
#   nothing for all six, "notau" for nue nuebar numu numubar, otherwise the
#   flavor names joined by <sep>  (e.g. "nue-nuebar", "nue:nuebar")
gxsp_nuflavors_tag() {
  local pdgs=$1 sep=$2 p out=""
  if [ "${pdgs}" == "${GXSP_NUFLAVORS_ALL}" ]; then return 0 ; fi
  if [ "${pdgs}" == "12 -12 14 -14" ]; then echo "notau" ; return 0 ; fi
  for p in ${pdgs} ; do out="${out:+${out}${sep}}$(gxsp_probe_name "${p}")" ; done
  echo "${out}"
}

##############################################################################
# tunes:  AR23_20i_00_000 <-> AR2320i00000
##############################################################################
GXSP_TUNE_RE='^([A-Za-z]+)([0-9]{2})_([0-9]{2})([a-z])_([a-z0-9]{2})_([a-z0-9]{3})$'
GXSP_TUNESQ_RE='^([A-Za-z]+)([0-9]{2})([0-9]{2})([a-z])([a-z0-9]{2})([a-z0-9]{3})$'

gxsp_tune_valid()   { [[ "$1" =~ ${GXSP_TUNE_RE} ]] ; }
# comprehensive model configuration of a tune:  AR23_20i_00_000 -> AR23_20i
gxsp_tune_cmc() {
  [[ "$1" =~ ${GXSP_TUNE_RE} ]] || return 1
  echo "${BASH_REMATCH[1]}${BASH_REMATCH[2]}_${BASH_REMATCH[3]}${BASH_REMATCH[4]}"
}
gxsp_tune_squeeze() { echo "${1//_/}" ; }
gxsp_tune_unsqueeze() {
  [[ "$1" =~ ${GXSP_TUNESQ_RE} ]] || return 1
  local m=("${BASH_REMATCH[@]}")
  echo "${m[1]}${m[2]}_${m[3]}${m[4]}_${m[5]}_${m[6]}"
}

##############################################################################
# isotopes
#
#   isotopes.cfg lines:  100ZZZAAA0 name [abundance] [keywords ...]
#   text after '#' is ignored
##############################################################################
# gxsp_iso_select <file> [keyword]
#   print "pdg name" for every active line; if keyword is given only for
#   lines carrying that keyword (e.g. reduced, root, electron)
gxsp_iso_select() {
  awk -v key="${2:-}" '
    { sub(/#.*/, "") }
    NF < 2 || $1 !~ /^100[0-9]+$/ { next }
    key == "" { print $1, $2 ; next }
    { for (i = 3; i <= NF; i++) if ($i == key) { print $1, $2 ; break } }
  ' "$1"
}

# gxsp_load_isotopes
#   fill ISO_PDGS / ISO_NAMES arrays with all the targets of this production
#   (all active lines for neutrinos, "electron" lines otherwise); these are
#   what the big file covers (README, TGraphs, reduction)
gxsp_load_isotopes() {
  local key=""
  if [ "${ELECTRON:-0}" -eq 1 ]; then key="electron"; fi
  ISO_PDGS=() ; ISO_NAMES=()
  local pdg name
  while read -r pdg name ; do
    ISO_PDGS+=("${pdg}") ; ISO_NAMES+=("${name}")
  done < <(gxsp_iso_select "${ISOFILE}" "${key}")
  NISO=${#ISO_PDGS[@]}
  [ "${NISO}" -gt 0 ] || gxsp_die 11 "no isotopes selected from ${ISOFILE}"
}

# gxsp_graph_isotopes  (after gxsp_load_isotopes)
#   print "pdg name" for stage 5's TGraph isotopes: those marked "root" that
#   are in this production, or all of them (returning 1) if none are
gxsp_graph_isotopes() {
  local pdg name i n=0
  while read -r pdg name ; do
    for (( i = 0 ; i < NISO ; i++ )) ; do
      if [ "${ISO_PDGS[$i]}" == "${pdg}" ]; then echo "${pdg} ${name}" ; n=$(( n + 1 )) ; break ; fi
    done
  done < <(gxsp_iso_select "${ISOFILE}" root)
  [ ${n} -gt 0 ] && return 0
  for (( i = 0 ; i < NISO ; i++ )) ; do echo "${ISO_PDGS[$i]} ${ISO_NAMES[$i]}" ; done
  return 1
}

# gxsp_tgraph_emax  -> stage 5's TGraph upper energy: TGEMAX (gxspinit
#   --tgraph, default 125; a define_cfg.sh from before --tgraph has none),
#   but not beyond EMAX
gxsp_tgraph_emax() {
  local g=${TGEMAX:-125}
  if awk -v e="${EMAX}" -v g="${g}" 'BEGIN { exit !(e < g) }' ; then g=${EMAX} ; fi
  echo "${g}"
}

# gxsp_load_stage3 [group]
#   fill S3_PDGS / S3_NAMES with the stage-3 targets: the production's
#   isotopes without the free nucleons (free-n 1000000010, H1 1000010010),
#   whose splines come from stage 1 via gxspl-freenuc.xml.  NS3 may be 0.
#   group (for the --average-qe seed clusters, see gxsp_qe_n):
#     all    every stage-3 target (default)
#     qez    only those whose Z is in QEZ (gxspinit --qe-z)
#     other  only those whose Z is not
gxsp_load_stage3() {
  local group=${1:-all}
  gxsp_load_isotopes
  S3_PDGS=() ; S3_NAMES=()
  local i
  for (( i = 0 ; i < NISO ; i++ )) ; do
    case " ${GXSP_FREENUC_PDGS[*]} " in
      *" ${ISO_PDGS[$i]} "*) continue ;;
    esac
    case "${group}" in
      all)   ;;
      qez)   gxsp_is_qez "${ISO_PDGS[$i]}" || continue ;;
      other) gxsp_is_qez "${ISO_PDGS[$i]}" && continue ;;
      *)     gxsp_die 11 "gxsp_load_stage3: bad group '${group}'" ;;
    esac
    S3_PDGS+=("${ISO_PDGS[$i]}") ; S3_NAMES+=("${ISO_NAMES[$i]}")
  done
  NS3=${#S3_PDGS[@]}
}

# gxsp_pdg_z <pdg>  -> Z of a nucleus pdg code 10LZZZAAAI (e.g. 1000180400 -> 18)
gxsp_pdg_z() { echo $(( ($1 / 10000) % 1000 )) ; }
# gxsp_is_qez <pdg>  -> true if its Z is one of QEZ (the isotopes that get
#   the second --average-qe count; gxspinit --qe-z, default 6,8,18)
gxsp_is_qez() {
  case " ${QEZ} " in *" $(gxsp_pdg_z "$1") "*) return 0 ;; esac
  return 1
}
# gxsp_z_symbols "<Z> ..."  -> their element symbols, e.g. "6 8 18" -> "C, O, Ar"
gxsp_z_symbols() {
  echo "$1" | awk '
    BEGIN { n = split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn " \
      "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm " \
      "Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu " \
      "Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og", sym, " ") }
    { for (i = 1 ; i <= NF ; i++) out = out (i > 1 ? ", " : "") (($i in sym) ? sym[$i] : "Z" $i) }
    END { print out }'
}

# gxsp_target_pdg <name|pdg>   (freep/freen or anything in $ISOFILE)
gxsp_target_pdg() {
  local i
  for i in 0 1 ; do
    if [ "$1" == "${GXSP_FREENUC_NAMES[$i]}" ] || [ "$1" == "${GXSP_FREENUC_PDGS[$i]}" ]; then
      echo "${GXSP_FREENUC_PDGS[$i]}" ; return 0
    fi
  done
  gxsp_iso_select "${ISOFILE}" | awk -v t="$1" '$1 == t || $2 == t { print $1 ; f=1 ; exit } END { exit !f }'
}
# gxsp_target_name <name|pdg>
gxsp_target_name() {
  local i
  for i in 0 1 ; do
    if [ "$1" == "${GXSP_FREENUC_NAMES[$i]}" ] || [ "$1" == "${GXSP_FREENUC_PDGS[$i]}" ]; then
      echo "${GXSP_FREENUC_NAMES[$i]}" ; return 0
    fi
  done
  gxsp_iso_select "${ISOFILE}" | awk -v t="$1" '$1 == t || $2 == t { print $2 ; f=1 ; exit } END { exit !f }'
}

##############################################################################
# job index <-> (probe, target) mapping
##############################################################################
# gxsp_njobs <stage> [group]   (1 or 3; group as for gxsp_load_stage3)
gxsp_njobs() {
  local probes=($(gxsp_probe_pdgs))
  case "$1" in
    1) echo $(( ${#GXSP_FREENUC_PDGS[@]} * ${#probes[@]} )) ;;
    3) gxsp_load_stage3 "$2" ; echo $(( NS3 * ${#probes[@]} )) ;;
    *) return 1 ;;
  esac
}
# gxsp_pair_from_index <stage> <idx> [group]  -> sets PAIR_PROBE, PAIR_TARGET (pdg codes)
gxsp_pair_from_index() {
  local stage=$1 idx=$2 group=${3:-all}
  local probes=($(gxsp_probe_pdgs))
  local np=${#probes[@]}
  local targets
  case "${stage}" in
    1) targets=("${GXSP_FREENUC_PDGS[@]}") ;;
    3) gxsp_load_stage3 "${group}" ; targets=("${S3_PDGS[@]}") ;;
    *) gxsp_die 11 "gxsp_pair_from_index: bad stage ${stage}" ;;
  esac
  local n=$(( ${#targets[@]} * np ))
  if ! [[ "${idx}" =~ ^[0-9]+$ ]] || [ "${idx}" -ge "${n}" ]; then
    gxsp_die 11 "job index '${idx}' out of range [0,${n}) for stage ${stage} (${group} pairs)"
  fi
  PAIR_PROBE=${probes[$(( idx % np ))]}
  PAIR_TARGET=${targets[$(( idx / np ))]}
}
# gxsp_pair_file <probe> <target>  -> gxspl-<probename>-<targetname>.xml
gxsp_pair_file() {
  echo "gxspl-$(gxsp_probe_name "$1")-$(gxsp_target_name "$2").xml"
}

##############################################################################
# well known file names
##############################################################################
gxsp_freenuc_file() { echo "gxspl-freenuc.xml" ; }
gxsp_big_file() {
  if [ "${ELECTRON:-0}" -eq 1 ]; then echo "gxspl-ELECTRON.xml"
  else                               echo "gxspl-NUbig.xml"
  fi
}
# every version of the big file (gxspl-NUbig.xml, -stdCCQE) is kept in
# $STAGEDIR gzipped with --best, and uncompressed on the worker node for use
gxsp_big_gz() { echo "$(gxsp_big_file).gz" ; }
gxsp_small_file() { echo "gxspl-NUsmall.xml" ; }
gxsp_graph_file() { echo "xsec_graphs.root" ; }
# the tarball name (only) uses the dotted version: v3_08_00 -> 3.08.00
gxsp_sversion_dots() { local v=${SVERSION#v} ; echo "${v//_/.}" ; }
gxsp_ups_tarball() { echo "genie_xsec-$(gxsp_sversion_dots)-noarch-${UPSQUALDASH}.tar.bz2" ; }

# --average-qe products (gxsp_qe_on)
#   the QE event generator list of the seeded runs: CCQE for neutrinos,
#   EMQE for electrons.  It names the files and directories below (-<QE>).
gxsp_qe_genlist() {
  if [ "${ELECTRON:-0}" -eq 1 ]; then echo "EMQE" ; else echo "CCQE" ; fi
}
#   stage 3q: gxsp_mkspl --genlist <QE> --seed <k> output, k = 1..gxsp_qe_n
#             (same suffix rule as gxsp_mkspl: -<genlist>-seed<k>)
gxsp_qe_file() { local f=$(gxsp_pair_file "$1" "$2") ; echo "${f%.xml}-$(gxsp_qe_genlist)-seed$3.xml" ; }
#   ... kept in one directory per probe x isotope (created by gxspinit):
#   work-products/isotopes/gxspl-<probe>-<isotope>-<QE>/
gxsp_qe_dir() { local f=$(gxsp_pair_file "$1" "$2") ; echo "work-products/isotopes/${f%.xml}-$(gxsp_qe_genlist)" ; }
#   stage 3s, per probe x isotope (in that pair's -<QE> directory):
#     gxspl-<probe>-<isotope>-<QE>-average.{xml,parquet}   (average_splines.py)
#     gxspl-<probe>-<isotope>-<QE>-smooth.{xml,parquet}    (construct_representative_spline.py)
#   average_splines.py writes one parquet per spline it averages: with one
#   spline (neutrinos) it is <D>-average.parquet; with several (electrons:
#   one per hit nucleon) <D>-average_nu<probe>_tgt<target>_N<nucleon>.parquet
#   each, with pdg codes.  <D>-average.xml holds them all.
gxsp_qe_base()   { basename "$(gxsp_qe_dir "$1" "$2")" ; }
gxsp_qe_avg()    { echo "$(gxsp_qe_base "$1" "$2")-average" ; }
gxsp_qe_smooth() { echo "$(gxsp_qe_base "$1" "$2")-smooth" ; }
#   the smooth files of a pair: one for neutrinos; for electrons one per hit
#   nucleon (construct_representative_spline.py --nucleon 2212 / 2112), each
#   made from that nucleon's average parquet.
#   gxsp_qe_smooths <probe> <target>  (pdg codes)
#     -> lines "<smooth basename> <nucleon or -> <average parquet basename>"
gxsp_qe_smooths() {
  local b=$(gxsp_qe_smooth "$1" "$2") a=$(gxsp_qe_avg "$1" "$2")
  if [ "${ELECTRON:-0}" -eq 1 ]; then
    echo "${b}-p 2212 ${a}_nu$1_tgt$2_N2212"
    echo "${b}-n 2112 ${a}_nu$1_tgt$2_N2112"
  else
    echo "${b} - ${a}"
  fi
}
#   ... and the big file with the standard (un-averaged) QE splines:
#   stage 4's output when gxsp_qe_on, which stage 4q reads to make the big file
#   (putting the stage 3s smooth splines in place of its QE ones)
#   (gxspl-NUbig-stdCCQE.xml.gz, gxspl-ELECTRON-stdEMQE.xml.gz)
gxsp_big_std_file()  { local b=$(gxsp_big_file) ; echo "${b%.xml}-std$(gxsp_qe_genlist).xml.gz" ; }
#   how many seeded QE runs: --average-qe N1[,N2] gives NQE = N1 for most
#   isotopes and NQEZ = N2 (default N1) for those whose Z is in QEZ
#   (gxspinit --qe-z; default "6 8 18", carbon, oxygen, argon).  A count of 0 means that isotope's
#   QE splines are left as they are.
gxsp_qe_n() {
  if gxsp_is_qez "$1" ; then echo "${NQEZ:-0}" ; else echo "${NQE:-0}" ; fi
}
gxsp_qe_max() {
  local a=${NQE:-0} b=${NQEZ:-0}
  if [ "${a}" -gt "${b}" ]; then echo "${a}" ; else echo "${b}" ; fi
}
gxsp_qe_min() {
  local a=${NQE:-0} b=${NQEZ:-0}
  if [ "${a}" -lt "${b}" ]; then echo "${a}" ; else echo "${b}" ; fi
}
# gxsp_qe_on  -> true if this production averages QE (any count > 0)
gxsp_qe_on() { [ "$(gxsp_qe_max)" -gt 0 ] ; }
# gxsp_qe_seed_group <seed>  -> the pair group whose $PROCESS mapping the
#   cluster running that seed uses: seeds 1..min(N1,N2) are run for every
#   pair ("all"); the seeds above that only for the pairs that have more
#   ("qez" if N2 > N1, else "other"), so no job is left with nothing to do.
gxsp_qe_seed_group() {
  if [ "$1" -le "$(gxsp_qe_min)" ]; then echo all
  elif [ "${NQEZ:-0}" -gt "${NQE:-0}" ]; then echo qez
  else echo other
  fi
}
# gxsp_qe_clusters <chunk>  -> lines "<first> <last> <group>": the stage 3q
#   DAG clusters, seed ranges of up to <chunk> seeds; first those over all
#   pairs (1..min), then those over the group with the larger count
gxsp_qe_clusters() {
  local chunk=$1 lo=$(gxsp_qe_min) hi=$(gxsp_qe_max) s e
  for (( s = 1 ; s <= lo ; s += chunk )) ; do
    e=$(( s + chunk - 1 )) ; if [ ${e} -gt ${lo} ]; then e=${lo} ; fi
    echo "${s} ${e} all"
  done
  for (( s = lo + 1 ; s <= hi ; s += chunk )) ; do
    e=$(( s + chunk - 1 )) ; if [ ${e} -gt ${hi} ]; then e=${hi} ; fi
    echo "${s} ${e} $(gxsp_qe_seed_group ${s})"
  done
}
# gxsp_pair_table  -> lines "<probe pdg> <probe name> <target pdg> <target
#   name> <seeds>": every stage-3 pair, in stage-3 order (isotope, then
#   probe), with its --average-qe seed count (gxsp_qe_n; 0 without it).
#   The lists of expected files below are made from it in one awk pass:
#   calling the per-pair helpers (each a subshell, some reading isotopes.cfg)
#   for 750 pairs and 1000 seeds took over a minute in gxsp_status.
gxsp_pair_table() {
  local probes=($(gxsp_probe_pdgs)) names=() p i j z n
  for p in "${probes[@]}" ; do names+=("$(gxsp_probe_name "$p")") ; done
  gxsp_load_stage3
  for (( i = 0 ; i < NS3 ; i++ )) ; do
    z=$(( (${S3_PDGS[$i]} / 10000) % 1000 ))
    case " ${QEZ} " in *" ${z} "*) n=${NQEZ:-0} ;; *) n=${NQE:-0} ;; esac
    for (( j = 0 ; j < ${#probes[@]} ; j++ )) ; do
      echo "${probes[$j]} ${names[$j]} ${S3_PDGS[$i]} ${S3_NAMES[$i]} ${n}"
    done
  done
}

# gxsp_qe_entries  -> lines "<seed> <PROCESS> <rel path>": every stage 3q
#   product, seed-major, with the $PROCESS index (within its seed's group,
#   gxsp_qe_seed_group) of the job that makes it: seeds 1..min(N1,N2) over
#   every pair, the seeds above that over the pairs with the larger count
gxsp_qe_entries() {
  gxsp_pair_table | awk -v lo="$(gxsp_qe_min)" -v hi="$(gxsp_qe_max)" -v qe="$(gxsp_qe_genlist)" '
    { np++ ; f = "gxspl-" $2 "-" $4 "-" qe
      pre[np] = "work-products/isotopes/" f "/" f "-seed" ; n[np] = $5
      if (hi > lo && $5 == hi) gidx[np] = ng++ }
    END {
      for (s = 1 ; s <= hi ; s++)
        for (i = 1 ; i <= np ; i++) {
          if (s <= lo)        print s, i - 1, pre[i] s ".xml"
          else if (n[i] >= s) print s, gidx[i], pre[i] s ".xml"
        }
    }'
}

# gxsp_qe_pairs  -> lines "<probe pdg> <target pdg> <seeds>": the averaged
#   pairs (stage-3 order: isotope, then probe), skipping isotopes whose count is 0
gxsp_qe_pairs() {
  gxsp_pair_table | awk '$5 > 0 { print $1, $3, $5 }'
}
# gxsp_qe_smooth_njobs  -> how many stage 3s jobs: QESMJOBS (gxspinit
#   --qe-smooth-jobs, default 25), but no more than there are averaged pairs
gxsp_qe_smooth_njobs() {
  local np=$(gxsp_qe_pairs | wc -l | tr -d ' ') nj=${QESMJOBS:-25}
  if [ "${nj}" -gt "${np}" ]; then nj=${np} ; fi
  echo "${nj}"
}
# gxsp_qe_smooth_assign  -> lines "<job> <probe pdg> <target pdg> <seeds>",
#   in pair order: which stage 3s job ($PROCESS) does each pair.  Pairs go,
#   largest seed count first, to the job with the fewest seeds so far (ties:
#   the lowest job), so the jobs fetch about the same number of seed files.
#   Every job works this out the same way from the same lists.
gxsp_qe_smooth_assign() {
  local nj=$(gxsp_qe_smooth_njobs)
  gxsp_qe_pairs | awk '{ print NR, $0 }' | sort -k4,4nr -k1,1n |
  awk -v nj="${nj}" '
    { best = 0
      for (j = 1 ; j < nj ; j++) if (load[j] < load[best]) best = j
      load[best] += $4
      print $1, best, $2, $3, $4 }' |
  sort -k1,1n | awk '{ print $2, $3, $4, $5 }'
}

# gxsp_stage4_out  -> what stage 4 writes: the big file itself, or (when
#   stage 4q follows) the -std<QE> file.  Either way a new file, never a
#   replacement.
gxsp_stage4_out() {
  if gxsp_qe_on ; then gxsp_big_std_file ; else gxsp_big_gz ; fi
}

# gxsp_stages  -> the stages of this production, in order
gxsp_stages() {
  if gxsp_qe_on ; then echo "1 2 3 3q 3s 4 4q 5 6"
  else                 echo "1 2 3 4 5 6"
  fi
}

# gxsp_expected <stage>  -> rel paths of the products of that stage
#   (3q is seed-major, in the order of gxsp_qe_entries)
gxsp_expected() {
  local p t s
  case "$1" in
    3q) gxsp_qe_entries | awk '{ print $3 }' ;;
    3s) # a smooth file for every averaged pair (two for electrons: -p, -n;
        #   the same names as gxsp_qe_smooths)
        gxsp_pair_table | awk -v qe="$(gxsp_qe_genlist)" -v el="${ELECTRON:-0}" '
          $5 > 0 { f = "gxspl-" $2 "-" $4 "-" qe ; d = "work-products/isotopes/" f "/" f "-smooth"
                   if (el == 1) { print d "-p.xml" ; print d "-n.xml" } else print d ".xml" }' ;;
    4q) echo "work-products/$(gxsp_big_gz)" ;;
    1) for t in "${GXSP_FREENUC_PDGS[@]}" ; do
         for p in $(gxsp_probe_pdgs) ; do
           echo "work-products/freenucs/$(gxsp_pair_file "$p" "$t")"
         done
       done ;;
    2) echo "work-products/$(gxsp_freenuc_file)" ;;
    3) gxsp_pair_table | awk '{ print "work-products/isotopes/gxspl-" $2 "-" $4 ".xml" }' ;;
    4) echo "work-products/$(gxsp_stage4_out)" ;;
    5) if [ "${ELECTRON:-0}" -eq 0 ]; then echo "work-products/$(gxsp_small_file)" ; fi
       echo "work-products/$(gxsp_graph_file)" ;;
    6) echo "ups/$(gxsp_ups_tarball)" ;;
    *) return 1 ;;
  esac
}

##############################################################################
# staging area (PNFS) access
#
#   On grid nodes PNFS is only reachable through ifdh; interactively plain
#   Unix commands work.  GXSP_IO=ifdh|local overrides the guess.
##############################################################################
gxsp_io_mode() {
  if [ -n "${GXSP_IO:-}" ]; then echo "${GXSP_IO}"
  elif [ -n "${GRID_USER:-}" ]; then echo ifdh
  else echo local
  fi
}

GXSP_LISTING=${GXSP_LISTING:-${PWD}/stagedir.ls}

# gxsp_listing_try [--refresh]  -> create (once) a normalized recursive
#   listing of $STAGEDIR in $GXSP_LISTING; one absolute path per line, no
#   trailing /.  Returns 1 (after reporting why) if nothing could be listed.
gxsp_listing_try() {
  if [ "$1" == "--refresh" ]; then rm -f "${GXSP_LISTING}" ; fi
  if [ -s "${GXSP_LISTING}" ]; then return 0 ; fi
  local raw="${GXSP_LISTING}.raw" err="${GXSP_LISTING}.err" mode rc cmd
  mode=$(gxsp_io_mode)
  if [ "${mode}" == "ifdh" ]; then
    # (no --debug=2: with it, the spack fife-utils ifdh on AlmaLinux 9 aborts)
    cmd="ifdh ls ${STAGEDIR} 99"
    gxsp_msg "${cmd}"
    ifdh ls "${STAGEDIR}" 99 > "${raw}" 2> "${err}"
  else
    cmd="find ${STAGEDIR}"
    find "${STAGEDIR}" > "${raw}" 2> "${err}"
  fi
  rc=$?
  # the scratch area (and the .err) goes away, so a failed listing's stderr
  # is copied to the job's stderr
  if [ ${rc} -ne 0 ]; then
    gxsp_warn "listing failed: ${cmd} (io=${mode}, exit status ${rc}, GRID_USER=${GRID_USER:-unset})"
    gxsp_dump_err "${err}" "${cmd}"
  fi
  sed -e 's#//*#/#g' -e 's#/$##' "${raw}" | sort -u > "${GXSP_LISTING}"
  if [ -s "${GXSP_LISTING}" ]; then return 0 ; fi
  gxsp_warn "listing ${STAGEDIR} gave nothing (io=${mode}, exit status ${rc}, GRID_USER=${GRID_USER:-unset})"
  if [ ${rc} -eq 0 ]; then gxsp_dump_err "${err}" "${cmd}" ; fi
  rm -f "${GXSP_LISTING}"
  return 1
}
# gxsp_dump_err <errfile> <command>  -> copy a command's saved stderr to
#   stderr, between clearly marked delimiter lines
gxsp_dump_err() {
  echo "==================== BEGIN stderr of: $2 ====================" >&2
  if [ -s "$1" ]; then cat "$1" >&2 ; else echo "(stderr was empty)" >&2 ; fi
  echo "==================== END stderr of: $2 ======================" >&2
}
# gxsp_listing [--refresh]  -> as gxsp_listing_try, but exits (3) on failure
gxsp_listing() {
  gxsp_listing_try "$@" || gxsp_die 12 "could not list ${STAGEDIR}"
}

# gxsp_check_inputs <rel> ...
#   check that these files are in the staging area, taking a fresh listing
#   and checking again up to GXSP_RETRIES times (default 3), GXSP_RETRY_WAIT
#   seconds apart (default 60), before giving up -- PNFS can be slow to show
#   files the previous DAG step has just written, and a listing can fail.
#   Returns 0 if all are present; 3 if no listing could be made; otherwise 1,
#   printing the missing ones, with GXSP_NMISSING set to how many.
#   Runs in the current shell (never call it inside $( )).
gxsp_check_inputs() {
  local tries=$(( ${GXSP_RETRIES:-3} + 1 )) wait=${GXSP_RETRY_WAIT:-60}
  local i rc missing=""
  GXSP_NMISSING=0
  [ $# -gt 0 ] || return 0
  for (( i = 1 ; i <= tries ; i++ )) ; do
    if [ ${i} -gt 1 ]; then
      gxsp_msg "waiting ${wait} s, then listing ${STAGEDIR} again (check ${i} of ${tries})"
      sleep "${wait}"
    fi
    if ! gxsp_listing_try --refresh ; then rc=3 ; continue ; fi
    missing=$(gxsp_missing "$@")
    if [ -z "${missing}" ]; then return 0 ; fi
    rc=1
    GXSP_NMISSING=$(echo "${missing}" | wc -l | tr -d ' ')
    gxsp_msg "${GXSP_NMISSING} of $# inputs not listed (check ${i} of ${tries})"
  done
  if [ ${rc} -eq 1 ]; then
    echo "missing after ${tries} checks:"
    echo "${missing}"
  fi
  return ${rc}
}

# gxsp_require_inputs <what> <rel> ...   gxsp_check_inputs, exiting on failure
#   (8 = missing inputs, 3 = no listing) with "<n> of <m> <what> files missing"
gxsp_require_inputs() {
  local what=$1 rc ; shift
  gxsp_check_inputs "$@"
  rc=$?
  case ${rc} in
    0) return 0 ;;
    3) gxsp_die 12 "could not list ${STAGEDIR}" ;;
    *) gxsp_die 17 "${GXSP_NMISSING} of $# ${what} files missing" ;;
  esac
}

# gxsp_have <rel>   (uses the cached listing)
gxsp_have() {
  gxsp_listing
  grep -qxF "$(echo "${STAGEDIR}/$1" | sed -e 's#//*#/#g')" "${GXSP_LISTING}"
}

# gxsp_missing <rel> ...  -> print those not present; return 1 if any missing
gxsp_missing() {
  [ $# -gt 0 ] || return 0
  gxsp_listing
  # one pass over the listing rather than a grep per file
  printf '%s\n' "$@" | awk -v sd="${STAGEDIR}" '
    NR == FNR { have[$0] = 1 ; next }
    { p = sd "/" $0 ; gsub(/\/+/, "/", p) ; if (!(p in have)) { print ; bad = 1 } }
    END { exit bad }' "${GXSP_LISTING}" -
}

gxsp_cp() {
  if [ "$(gxsp_io_mode)" == "ifdh" ]; then ifdh cp "$1" "$2"
  else                                     cp "$1" "$2"
  fi
}
# gxsp_mkdir <rel>   create a directory in the staging area if it is missing
gxsp_mkdir() {
  gxsp_have "$1" && return 0
  gxsp_msg "creating ${1}/"
  if [ "$(gxsp_io_mode)" == "ifdh" ]; then ifdh mkdir "${STAGEDIR}/$1"
  else                                     mkdir -p "${STAGEDIR}/$1"
  fi
}

# NOTE: jobs never delete (or overwrite) anything in the staging area --
# grid tokens may allow creating files in PNFS but not removing them.  A
# product that already exists is kept; a job whose products all exist has
# nothing to do (gxsp_skip_if_done).  To redo something, remove it
# interactively first.

# gxsp_fetch <rel> [localname]   copy from the staging area; die on failure
#   (a file not in the cached listing gets gxsp_check_inputs' re-listings
#   first, in case PNFS hasn't shown it yet).  A failed copy is tried again
#   up to GXSP_RETRIES times, GXSP_RETRY_WAIT seconds apart: a transfer that
#   fails part way can leave a partial local file, which ifdh's own retry
#   then refuses to overwrite ("File exists"), so it's removed first.
gxsp_fetch() {
  local rel=$1 local=${2:-$(basename "$1")}
  local tries=$(( ${GXSP_RETRIES:-3} + 1 )) wait=${GXSP_RETRY_WAIT:-60} i
  if ! { gxsp_listing_try && gxsp_have "${rel}" ; } ; then
    gxsp_msg "${rel} not in the listing yet; listing again"
    gxsp_check_inputs "${rel}" || gxsp_die 13 "missing in staging area: ${rel}"
  fi
  for (( i = 1 ; i <= tries ; i++ )) ; do
    if [ ${i} -gt 1 ]; then
      gxsp_msg "copy of ${rel} failed; waiting ${wait} s, then trying again (copy ${i} of ${tries})"
      sleep "${wait}"
    fi
    rm -f "${local}"
    gxsp_cp "${STAGEDIR}/${rel}" "${local}" && return 0
  done
  rm -f "${local}"
  gxsp_die 13 "copy failed ${tries} times: ${rel}"
}

# gxsp_put <localfile> <rel>     copy a new file into the staging area.
#   Never deletes or overwrites: if <rel> already exists it is left as it is
#   and 2 is returned (1 if the copy fails).
gxsp_put() {
  local local=$1 rel=$2
  [ -f "${local}" ] || { gxsp_warn "no local ${local} to put" ; return 1 ; }
  if gxsp_have "${rel}" ; then
    gxsp_warn "${rel} already exists; not replaced (jobs don't delete in PNFS)"
    return 2
  fi
  gxsp_cp "${local}" "${STAGEDIR}/${rel}" || { gxsp_warn "put failed: ${rel}" ; return 1 ; }
  echo "$(echo "${STAGEDIR}/${rel}" | sed -e 's#//*#/#g')" >> "${GXSP_LISTING}"   # now listed
}

# gxsp_put_log <localfile> <rel>  copy a log (or other side file) in; if <rel>
#   exists already (an earlier attempt), use the first free <rel>.1, <rel>.2 ...
gxsp_put_log() {
  local local=$1 rel=$2 n=0 try=$2
  while gxsp_have "${try}" ; do n=$(( n + 1 )) ; try="${rel}.${n}" ; done
  gxsp_put "${local}" "${try}"
}

# gxsp_skip_if_done <rel> ...   if all these products are already in the
#   staging area, say so and end the job successfully (nothing to redo)
gxsp_skip_if_done() {
  [ $# -gt 0 ] || return 0
  gxsp_listing_try || return 0          # can't tell: carry on (puts never overwrite)
  [ -z "$(gxsp_missing "$@")" ] || return 0
  gxsp_msg "already in the staging area, nothing to do (remove interactively to redo):"
  printf '    %s\n' "$@"
  gxsp_finish 0
}

##############################################################################
# job helpers
##############################################################################
# gxsp_join <sep> <items...>
gxsp_join() { local IFS=$1 ; shift ; echo "$*" ; }

# gxsp_log_header <logfile> <description...>   start a fresh log
gxsp_log_header() {
  local log=$1 ; shift
  {
    echo "=== $* ==="
    echo "date      $(date '+%F %T %Z')"
    echo "host      $(uname -n)   $(uname -srm)"
    echo "os        $(cat /etc/system-release 2>/dev/null || echo unknown)"
    echo "pwd       ${PWD}"
    echo "stagedir  ${STAGEDIR}"
    echo "tune      ${TUNE}  genlist ${GENLIST}  knots ${KNOTS}  emax ${EMAX}"
    echo "genie     ${GVERSION} -q ${GQUAL}   splines ${SVERSION}  electron=${ELECTRON}"
    echo "process   ${PROCESS:-unset}   cluster ${CLUSTER:-unset}"
    echo
  } > "${log}"
}

# gxsp_run <logfile> <command...>   run and append output to the log,
#   returning the command's status
gxsp_run() {
  local log=$1 ; shift
  local t0=$(date '+%s')
  gxsp_msg "$*"
  echo "+ $*" >> "${log}"
  "$@" >> "${log}" 2>&1
  local rc=$?
  echo "+ status ${rc}, wall time $(( $(date '+%s') - t0 )) s" >> "${log}"
  return ${rc}
}

# gxsp_strip_freenuc <xml>   remove the free-nucleon splines (tgt free-n or
#   H1) from a gmkspl output, in place; sets STRIP_N (how many were removed).
#   gmkspl copies every spline of --input-cross-sections into its output, so
#   a seeded QE file would otherwise carry all of gxspl-freenuc.xml.
gxsp_strip_freenuc() {
  local xml=$1 out="$1.strip"
  awk '
    index($0, "<spline ") && (index($0, "tgt:1000000010;") || index($0, "tgt:1000010010;")) { skip = 1 ; n++ ; next }
    skip { if (index($0, "</spline>")) skip = 0 ; next }
    { print }
    END { print n + 0 > "/dev/stderr" }' "${xml}" > "${out}" 2> "${out}.n" || { rm -f "${out}" "${out}.n" ; return 1 ; }
  STRIP_N=$(cat "${out}.n")
  mv "${out}" "${xml}" && rm -f "${out}.n"
}

# gxsp_qe_replace <list>   (stdin: a big file; stdout: the same with splines
#   replaced)  For every spline in the files named in <list> (one per line:
#   the stage 3s smooth files), put its knots in place of those of the spline
#   with the same name in the big file, and set its nknots.  One pass, line
#   by line, so memory use is small whatever the size of the big file.  A
#   summary goes to stderr; exits 1 if a smooth spline matches nothing (or
#   has no knots, or two files have the same one), after writing the output.
gxsp_qe_replace() {
  awk -v list="$1" '
    BEGIN {
      while ((getline f < list) > 0) {
        nfile++ ; txt = ""
        while ((getline l < f) > 0) txt = txt " " l
        close(f)
        # tokens: each tag and each value on its own
        gsub(/</, " <", txt) ; gsub(/>/, "> ", txt)
        nt = split(txt, tok, /[ \t\r\n]+/)
        cur = "" ; intag = 0
        for (i = 1 ; i <= nt ; i++) {
          t = tok[i]
          if (t == "<spline") { intag = 1 ; continue }
          if (intag) {
            if (t ~ /^name="/) {
              cur = t ; sub(/^name="/, "", cur) ; sub(/".*$/, "", cur)
              if (cur in nk) { printf "gxsp_qe_replace: %s is in %s and %s\n", cur, src[cur], f > "/dev/stderr" ; bad = 1 }
              nk[cur] = 0 ; src[cur] = f ; nrep++
            }
            if (t ~ />$/) intag = 0
            continue
          }
          if (cur == "") continue
          if (t == "<E>")    e = tok[i + 1]
          if (t == "<xsec>") { k = ++nk[cur] ; KE[cur, k] = e ; KX[cur, k] = tok[i + 1] }
          if (t == "</spline>") cur = ""
        }
      }
      close(list)
      for (s in nk) if (nk[s] == 0) { printf "gxsp_qe_replace: %s (in %s) has no knots\n", s, src[s] > "/dev/stderr" ; bad = 1 }
    }
    skip && index($0, "</spline>") {
      skip = 0
      if (nbase != nk[rep]) { nmis++ ; if (nmis <= 5) printf "gxsp_qe_replace: %s: %d knots replaced by %d\n", rep, nbase, nk[rep] > "/dev/stderr" }
      print ; next
    }
    skip { if (index($0, "<knot>")) nbase++ ; next }
    index($0, "<spline ") && match($0, /name="[^"]*"/) {
      nm = substr($0, RSTART + 6, RLENGTH - 7)
      if (nm in nk) {
        line = $0 ; sub(/nknots="[0-9]*"/, "nknots=\"" nk[nm] "\"", line) ; print line
        for (k = 1 ; k <= nk[nm] ; k++) print "\t<knot> <E> " KE[nm, k] " </E> <xsec> " KX[nm, k] " </xsec> </knot>"
        found[nm]++ ; nbig++ ; skip = 1 ; rep = nm ; nbase = 0
        next
      }
    }
    { print }
    END {
      for (s in nk) if (!(s in found)) { nmiss++ ; if (nmiss <= 5) printf "gxsp_qe_replace: %s (in %s) is not in the big file\n", s, src[s] > "/dev/stderr" }
      for (s in found) if (found[s] > 1) printf "gxsp_qe_replace: %s appears %d times in the big file (all replaced)\n", s, found[s] > "/dev/stderr"
      printf "gxsp_qe_replace: %d spline(s) from %d file(s); %d replaced in the big file, %d not found; %d with a different number of knots\n", nrep, nfile, nbig, nmiss, nmis > "/dev/stderr"
      exit (bad || nmiss) ? 1 : 0
    }'
}

# gxsp_scrub <xml>   replace NaN/inf cross section values by 0.0, and (unless
#   NEGZERO=0, gxspinit --keep-negative) negative ones too ("-0.0" is zero
#   and is left alone).  If anything was replaced, the original is kept as
#   <xml>-NAN (it had NaN/inf) or else <xml>-NEG, and 1 is returned with
#   SCRUB_ORIG (that file), SCRUB_NNAN and SCRUB_NNEG (the counts) set.
gxsp_scrub() {
  local xml=$1 out="$1.scrub" counts kind
  SCRUB_ORIG="" SCRUB_NNAN=0 SCRUB_NNEG=0
  counts=$(awk -v zneg="${NEGZERO:-1}" -v out="${out}" '
    {
      if (match($0, /<xsec>[^<]*<\/xsec>/)) {
        v = substr($0, RSTART + 6, RLENGTH - 13) ; gsub(/[[:space:]]/, "", v)
        bad = 0
        if (v ~ /^[-+]?([nN][aA][nN]|[iI][nN][fF])/) { nnan++ ; bad = 1 }
        else if (zneg && v + 0 < 0)                 { nneg++ ; bad = 1 }
        if (bad) $0 = substr($0, 1, RSTART - 1) "<xsec> 0.0 </xsec>" substr($0, RSTART + RLENGTH)
      }
      print > out
    }
    END { print nnan + 0, nneg + 0 }' "${xml}") || { rm -f "${out}" ; return 0 ; }
  read -r SCRUB_NNAN SCRUB_NNEG <<< "${counts}"
  if [ "${SCRUB_NNAN}" -eq 0 ] && [ "${SCRUB_NNEG}" -eq 0 ]; then rm -f "${out}" ; return 0 ; fi
  if [ "${SCRUB_NNAN}" -gt 0 ]; then kind=NAN ; else kind=NEG ; fi
  SCRUB_ORIG="${xml}-${kind}"
  mv "${xml}" "${SCRUB_ORIG}" && mv "${out}" "${xml}"
  gxsp_warn "${xml}: ${SCRUB_NNAN} NaN/inf and ${SCRUB_NNEG} negative cross section value(s) replaced by 0.0 (original kept as $(basename "${SCRUB_ORIG}"))"
  return 1
}

# gxsp_spladd <out> <in...>   combine spline files with gspladd, splitting
#   into partial sums when the -f list would be too long for gspladd
GXSP_MAXARGLEN=${GXSP_MAXARGLEN:-20000}
GXSP_NPART=0
gxsp_spladd() {
  local out=$1 ; shift
  if [ $# -eq 0 ]; then gxsp_warn "gxsp_spladd: no inputs for ${out}" ; return 1 ; fi
  if [ $# -eq 1 ]; then cp "$1" "${out}" ; return $? ; fi
  local all=$(gxsp_join , "$@")
  if [ ${#all} -le ${GXSP_MAXARGLEN} ]; then
    gxsp_run "${GXSP_LOG:-/dev/null}" gspladd -f "${all}" -o "${out}"
    return $?
  fi
  local parts=() group=() glen=0 f p
  for f in "$@" ; do
    if [ ${#group[@]} -gt 0 ] && [ $(( glen + ${#f} + 1 )) -gt ${GXSP_MAXARGLEN} ]; then
      _gxsp_spladd_group "${group[@]}" || return 1
      parts+=("${GXSP_PART}") ; group=() ; glen=0
    fi
    group+=("${f}") ; glen=$(( glen + ${#f} + 1 ))
  done
  _gxsp_spladd_group "${group[@]}" || return 1
  parts+=("${GXSP_PART}")
  gxsp_spladd "${out}" "${parts[@]}"
  local rc=$?
  for p in "${parts[@]}" ; do
    case "${p}" in gxsp_part*.xml) rm -f "${p}" ;; esac
  done
  return ${rc}
}
# combine a group into a new partial-sum file whose name is left in GXSP_PART
_gxsp_spladd_group() {
  if [ $# -eq 1 ]; then GXSP_PART=$1 ; return 0 ; fi
  GXSP_NPART=$(( GXSP_NPART + 1 ))
  GXSP_PART="gxsp_part${GXSP_NPART}.xml"
  local part=${GXSP_PART}
  gxsp_spladd "${part}" "$@"
  local rc=$?
  GXSP_PART=${part}
  return ${rc}
}

# gxsp_finish <status>   common end of every grid job
# let PNFS "catch up"
gxsp_finish() {
  sleep "${GXSP_SLEEP:-25}"
  gxsp_msg "done, exit status $1"
  exit "$1"
}

# end-of-file gxsplib.sh
