# Pythia6 Crash Fix for genie-dev

## Issue

genie-dev crashed with a segmentation violation during event generation, typically during hadronization or particle decay. The crash occurred randomly at different event numbers depending on the random seed.

**Stack trace:**
```
[libPythia6.dylib] pylogo_
[libPythia6.dylib] pylist_
[libPythia6.dylib] py1ent_ or py2ent_
[libGPhDcy] genie::Pythia6Decayer2023::Decay
```

## Investigation

1. **Initial hypothesis:** The crash appeared after commits on the `command-line-q2` branch that added EM-Q2-min configuration. This seemed to be the cause.

2. **Testing:** After reverting all code changes to match master, the crash persisted. This ruled out the source code changes as the cause.

3. **Library comparison:** Compared the Pythia6 libraries between production genie and genie-dev:

   ```bash
   md5 ~/opt/genie-dev/src/scripts/build/ext/v6_428/lib/libPythia6.dylib
   # 62c87c0054cec567701cf5133025960f

   md5 ~/opt/genie/lib/libPythia6.dylib
   # a59ec935d35865a2c13c2b125469ea3d (symlink to ROOTEGPythia6)
   ```

4. **Root cause identified:** Production genie has a symlink in `lib/` to a working Pythia6 library:
   ```
   ~/opt/genie/lib/libPythia6.dylib -> /Users/pbarham/opt/ROOTEGPythia6/lib/libPythia6.dylib
   ```

   genie-dev was missing this symlink, so `DYLD_LIBRARY_PATH` loaded the buggy library from `src/scripts/build/ext/v6_428/lib/` instead.

## Resolution

Create symlinks in genie-dev to use the correct Pythia6 libraries:

```bash
ln -sf /Users/pbarham/opt/ROOTEGPythia6/lib/libPythia6.dylib ~/opt/genie-dev/lib/libPythia6.dylib
ln -sf /Users/pbarham/opt/ROOTEGPythia6/lib/libEGPythia6.dylib ~/opt/genie-dev/lib/libEGPythia6.dylib
```

## Why This Works

The `genie_setup.sh` script sets `DYLD_LIBRARY_PATH` with `$GENIE/lib` before `$PYTHIA6`:

```bash
export DYLD_LIBRARY_PATH="$GENIE/lib:$PYTHIA6:$ROOTSYS/lib:$DYLD_LIBRARY_PATH"
```

With the symlinks in place, macOS loads the correct libraries from `$GENIE/lib` first. `libPythia6.dylib` provides the Fortran PYTHIA6 routines, while `libEGPythia6.dylib` provides ROOT's `TPythia6` wrapper library that rebuilt GENIE binaries may reference directly as `@rpath/libEGPythia6.dylib`.

## Prevention

When setting up a new GENIE development environment, ensure both Pythia6 library symlinks exist:

```bash
# Check if symlinks exist
ls -la $GENIE/lib/libPythia6.dylib $GENIE/lib/libEGPythia6.dylib

# If missing, create them (adjust path as needed for your system)
ln -sf /Users/pbarham/opt/ROOTEGPythia6/lib/libPythia6.dylib $GENIE/lib/libPythia6.dylib
ln -sf /Users/pbarham/opt/ROOTEGPythia6/lib/libEGPythia6.dylib $GENIE/lib/libEGPythia6.dylib
```

## Date

January 2026
