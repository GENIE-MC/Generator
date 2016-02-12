#include "INukeOset.h"
#include <cstddef>

// workaround to get access to last instance
INukeOset* osetUtils::currentInstance = NULL;

//! set up initial density and enegry values; set up pointer to current instance
INukeOset :: INukeOset () : nuclearDensity (-1.0), pionKineticEnergy (-1.0)
{
  osetUtils::currentInstance = this;
}

void INukeOset :: setCrossSections (const int &pionPDG, const double &protonFraction)
{  
  if (pionPDG == kPdgPi0)
  {
      cexCrossSection = cexCrossSections[2];
    totalCrossSection = qelCrossSections[2] +
                        absorptionCrossSection;
  }
  else
  {
    // set channel for pi on proton and pi on neutron
    const int channelIndexOnProton  = (pionPDG == kPdgPiP); // 0 = pi-, 1 = pi+
    const int channelIndexOnNeutron = (pionPDG == kPdgPiM); // 0 = pi+, 1 = pi-

    // total xsec = (Z * xsec_proton + (A-Z) * xsec_neutron) / A
    cexCrossSection = protonFraction * cexCrossSections[channelIndexOnProton] +
                      (1.0 - protonFraction) *
                      cexCrossSections[channelIndexOnNeutron];

    totalCrossSection = protonFraction *
                        qelCrossSections[channelIndexOnProton] +
                        (1.0 - protonFraction) *
                        qelCrossSections[channelIndexOnNeutron] +
                        + absorptionCrossSection;
  }
}  
