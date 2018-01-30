#include "INukeOset.h"
#include <cstddef>

// workaround to get access to last instance
INukeOset* osetUtils::currentInstance = NULL;

//! set up initial density and enegry values; set up pointer to current instance
INukeOset :: INukeOset () : fNuclearDensity (-1.0), fPionKineticEnergy (-1.0)
{
  osetUtils::currentInstance = this;
}

void INukeOset :: setCrossSections (const int &pionPDG, const double &protonFraction)
{  
  if (pionPDG == kPdgPi0)
  {
      fCexCrossSection = fCexCrossSections[2];
    fTotalCrossSection = fQelCrossSections[2] + fAbsorptionCrossSection;
  }
  else
  {
    // set channel for pi on proton and pi on neutron
    const int channelIndexOnProton  = (pionPDG == kPdgPiP); // 0 = pi-, 1 = pi+
    const int channelIndexOnNeutron = (pionPDG == kPdgPiM); // 0 = pi+, 1 = pi-

    // total xsec = (Z * xsec_proton + (A-Z) * xsec_neutron) / A
    fCexCrossSection = protonFraction * fCexCrossSections[channelIndexOnProton] +
                       (1.0 - protonFraction) * fCexCrossSections[channelIndexOnNeutron];

    fTotalCrossSection = protonFraction * fQelCrossSections[channelIndexOnProton] +
                        (1.0 - protonFraction) * fQelCrossSections[channelIndexOnNeutron] + fAbsorptionCrossSection;
  }
}  
