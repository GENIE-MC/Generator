##############################################################################
# gxsp_bootstrap -- common start-up for the scripts run on grid nodes
#
#   jobsub ships only the script itself to the worker node, so gxspinit
#   inlines this file in place of the "source ... # GXSP_INLINE" line when
#   it copies the scripts into $STAGEDIR/bin.
#
#   gxsp_bootstrap <stagedir>
#     * cd to $_CONDOR_SCRATCH_DIR (if set; otherwise to a new directory in
#       $TMPDIR or /tmp, removed at exit if the job succeeds)
#     * make ifdh available (grid) and fetch + unroll cfg.tar.gz
#     * source define_cfg.sh (which also sources gxsplib.sh)
#     * gxsp_setup_genie   (skipped if GXSP_NOSETUP=1, for testing)
##############################################################################
gxsp_bootstrap() {
  local sd=$1
  if [ -z "${sd}" ]; then echo "${GXSP_ME}: --stagedir is required" >&2 ; exit 10 ; fi
  # work in the job scratch area; when run by hand use a fresh directory in
  # $TMPDIR (or /tmp if TMPDIR isn't set) so the unpacked cfg files don't
  # land in the caller's directory.  It is removed if the job succeeds and
  # kept (for a look at what went wrong) if it fails.
  if [ -n "${_CONDOR_SCRATCH_DIR:-}" ]; then
    cd "${_CONDOR_SCRATCH_DIR}" || exit 10
  else
    scratch=$(mktemp -d "${TMPDIR:-/tmp}/${GXSP_ME}.XXXXXX") && [ -d "${scratch}" ] ||
      { echo "${GXSP_ME}: can't create a scratch directory in ${TMPDIR:-/tmp}" >&2 ; exit 10 ; }
    GXSP_SCRATCH=${scratch} GXSP_CALLDIR=${PWD}
    cd "${GXSP_SCRATCH}" || exit 10
    echo "${GXSP_ME}: scratch directory ${GXSP_SCRATCH} (not on a grid node; in \$TMPDIR, else /tmp; removed on success, kept on failure)"
    trap 'gxsp_scratch_cleanup $?' EXIT
  fi

  local io=${GXSP_IO:-}
  if [ -z "${io}" ]; then
    if [ -n "${GRID_USER:-}" ]; then io=ifdh ; else io=local ; fi
  fi
  echo "${GXSP_ME}: bootstrap in ${PWD} (io=${io})"

  if [ "${io}" == "ifdh" ]; then
    if ! type ifdh > /dev/null 2>&1 ; then
      if grep -q 'Scientific Linux release 7.9' /etc/system-release 2>/dev/null ; then
        source /cvmfs/fermilab.opensciencegrid.org/products/genie/bootstrap_many_ups.sh
        setup ifdhc
        export IFDH_CP_MAXRETRIES=2
      else
        # non-SL7 (e.g. AlmaLinux 9): ifdh from spack's fife-utils
        source /cvmfs/fermilab.opensciencegrid.org/packages/common/spack/current/NULL/share/spack/setup-env.sh
        spack load fife-utils@3.7.8%gcc@11.4.1
      fi
    fi
    ifdh cp "${sd}/cfg/cfg.tar.gz" cfg.tar.gz
  else
    cp "${sd}/cfg/cfg.tar.gz" cfg.tar.gz
  fi
  if [ ! -f cfg.tar.gz ]; then
    echo "${GXSP_ME}: could not fetch ${sd}/cfg/cfg.tar.gz" >&2 ; exit 10
  fi
  tar xzf cfg.tar.gz || { echo "${GXSP_ME}: can't unroll cfg.tar.gz" >&2 ; exit 10 ; }
  source ./define_cfg.sh
  export GXSP_IO=${io}

  if [ "${STAGEDIR}" != "${sd}" ]; then
    gxsp_warn "--stagedir ${sd} differs from define_cfg.sh STAGEDIR ${STAGEDIR}"
  fi
  if [ "${GXSP_NOSETUP:-0}" != "1" ]; then
    gxsp_setup_genie || gxsp_die 14 "GENIE setup (${GVERSION} -q ${GQUAL}) failed"
  fi
}
# gxsp_scratch_cleanup <exit status>  (EXIT trap, only for a scratch
#   directory made by gxsp_bootstrap): remove it on success, keep it and say
#   where it is on failure
gxsp_scratch_cleanup() {
  if [ "$1" -eq 0 ]; then
    cd "${GXSP_CALLDIR:-/}" 2>/dev/null || cd /
    rm -rf "${GXSP_SCRATCH}"
  else
    echo "${GXSP_ME}: exit status $1: scratch directory kept: ${GXSP_SCRATCH}" >&2
  fi
}
# end-of-file gxsp_bootstrap.sh
