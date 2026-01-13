# Implementation Report: Configurable EM Q² Minimum

## Overview
Made the minimum Q² threshold for electromagnetic scattering events configurable via XML configuration files, replacing the previously hardcoded value of 0.02 GeV².

## Motivation
Electromagnetic scattering cross-sections diverge as Q² → 0. The minimum Q² cutoff is physics-dependent and experiments may need different thresholds. Previously, changing this value required modifying source code and recompiling. This implementation allows users to adjust the threshold via XML configuration.

## Implementation

### 1. XML Configuration (5 files)

Added new `Kinematics` parameter set to:
- `config/CommonParam.xml`
- `config/GEM21_11a/CommonParam.xml`
- `config/GEM21_11b/CommonParam.xml`
- `config/GEM21_11c/CommonParam.xml`
- `config/GEM21_11d/CommonParam.xml`

```xml
<param_set name="Kinematics">
  <param type="double" name="EM-Q2-min"> 0.02 </param>
</param_set>
```

Default value maintains backward compatibility (0.02 GeV²).

### 2. Framework Changes (4 files)

#### KineUtils.h/cxx
- Added `Q2min_cut` parameter with default values to electromagnetic kinematic limit functions:
  - `InelQ2Lim_W(El, ml, M, W, Q2min_cut = kMinQ2Limit)`
  - `Inelq2Lim_W(El, ml, M, W, q2min_cut = -1*kMinQ2Limit)`
  - `InelQ2Lim(El, ml, M, Q2min_cut = kMinQ2Limit)`
  - `Inelq2Lim(El, ml, M, q2min_cut = -1*kMinQ2Limit)`
- Updated implementations to use the parameter instead of hardcoded `kMinQ2Limit`
- Maintained `kMinQ2Limit = 0.02` as static default for backward compatibility

#### KPhaseSpace.h/cxx
- Added `GetQ2MinEM()` static method to load `EM-Q2-min` from configuration
- Uses lazy initialization with caching (loads once on first call)
- Updated 4 call sites to use `GetQ2MinEM()` instead of hardcoded value:
  - Line 573: Q2Lim_W() for general inelastic
  - Line 653: Q2Lim() for quasi-elastic
  - Line 683: Q2Lim() for MEC
  - Line 695: Q2Lim() for general inelastic

### 3. Physics Module Changes (8 files)

#### Issues Found and Fixed

**Critical Issue in DISKinematicsGenerator.cxx:**
- **Problem**: Generator samples in (x,y) space and converts to (W,Q²), but had **no check** to reject events with Q² < Q²min
- **Impact**: Events were leaking below the Q² threshold even with correct cross-sections
- **Fix**: Added Q²min check after computing Q² from (x,y):
  ```cpp
  double Q2min = interaction->ProcInfo().IsEM() ? KPhaseSpace::GetQ2MinEM() : controls::kMinQ2Limit;
  if(interaction->KinePtr()->Q2() < Q2min) continue;
  ```

**Hardcoded References in 7 Other Files:**
All used `utils::kinematics::electromagnetic::kMinQ2Limit` directly. Replaced with `KPhaseSpace::GetQ2MinEM()`:

1. **RESKinematicsGenerator.cxx** (line 290)
   - Used in max cross-section calculation for resonance production

2. **EmpiricalMECPXSec2015.cxx** (line 129)
   - Used in phase space validation for empirical MEC model

3. **MECUtils.cxx** (line 346)
   - Used in MEC tensor utility functions

4. **SuSAv2MECPXSec.cxx** (line 120)
   - Used in SuSAv2 MEC differential cross-section calculation

5. **MECGenerator.cxx** (line 868)
   - Used in MEC event generation kinematics

6. **SuSAv2QELPXSec.cxx** (line 68)
   - Used in SuSAv2 QEL differential cross-section calculation

7. **QELEventGeneratorSuSA.cxx** (line 103)
   - Used in SuSAv2 QEL event generation kinematics

## Usage

Users can now adjust the EM Q² minimum by editing the relevant `CommonParam.xml`:

```xml
<param type="double" name="EM-Q2-min"> 1.0 </param>  <!-- Changed from 0.02 to 1.0 -->
```

Changes take effect immediately after recompiling and regenerating cross-section splines.

## Important Notes

1. **Recompilation Required**: Source code changes require rebuilding GENIE
2. **Spline Regeneration Required**: Cross-section splines must be regenerated with the new Q²min value, as total cross-sections are computed by integrating over the allowed phase space
3. **Backward Compatibility**: Default value (0.02 GeV²) maintains existing behavior
4. **Scope**: Only affects electromagnetic scattering events (IsEM() == true); neutrino interactions continue using `controls::kMinQ2Limit = 1E-4 GeV²`

## Testing

Verified with:
```bash
# Set EM-Q2-min = 1.0 in config/GEM21_11a/CommonParam.xml
gmkspl -p 11 -t 1000180400 -e 8.0 -n 50 --event-generator-list EM -o spline.xml --tune GEM21_11a_00_000
gevgen -n 100 -p 11 -t 1000180400 -e 5.2 --cross-sections spline.xml --tune GEM21_11a_00_000 --event-generator-list EM
```

Confirmed all generated events have Q² ≥ 1.0 GeV².

## Files Modified

### Configuration (5 files)
- config/CommonParam.xml
- config/GEM21_11a/CommonParam.xml
- config/GEM21_11b/CommonParam.xml
- config/GEM21_11c/CommonParam.xml
- config/GEM21_11d/CommonParam.xml

### Framework (4 files)
- src/Framework/Utils/KineUtils.h
- src/Framework/Utils/KineUtils.cxx
- src/Framework/Interaction/KPhaseSpace.h
- src/Framework/Interaction/KPhaseSpace.cxx

### Physics Modules (8 files)
- src/Physics/DeepInelastic/EventGen/DISKinematicsGenerator.cxx
- src/Physics/Resonance/EventGen/RESKinematicsGenerator.cxx
- src/Physics/Multinucleon/XSection/EmpiricalMECPXSec2015.cxx
- src/Physics/Multinucleon/XSection/MECUtils.cxx
- src/Physics/Multinucleon/XSection/SuSAv2MECPXSec.cxx
- src/Physics/Multinucleon/EventGen/MECGenerator.cxx
- src/Physics/QuasiElastic/XSection/SuSAv2QELPXSec.cxx
- src/Physics/QuasiElastic/EventGen/QELEventGeneratorSuSA.cxx

**Total: 17 files modified**
