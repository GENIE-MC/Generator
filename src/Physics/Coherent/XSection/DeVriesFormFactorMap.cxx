//____________________________________________________________________________
/*
  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda
  University of Liverpool
  <mroda \at liverpool.ac.uk>

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <vector>

#include "Framework/Messenger/Messenger.h"

#include "Framework/Registry/RegistryItemTypeDef.h"
#include "Physics/Coherent/XSection/DeVriesFormFactorMap.h"

#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

DeVriesFormFactorMap::DeVriesFormFactorMap()
    : COHFormFactorI("genie::DeVriesFormFactorMap") {}
//____________________________________________________________________________
DeVriesFormFactorMap::DeVriesFormFactorMap(string config)
    : COHFormFactorI("genie::DeVriesFormFactorMap", config) {}
//____________________________________________________________________________
DeVriesFormFactorMap::DeVriesFormFactorMap(string name, string config)
    : COHFormFactorI(name, config) {}
//____________________________________________________________________________
DeVriesFormFactorMap::~DeVriesFormFactorMap() {}
//____________________________________________________________________________
double DeVriesFormFactorMap::ProtonFF(double Q, int pdg) const {

  const std::map<int, const genie::DeVriesFormFactor *>::const_iterator it =
      fNuclearFFs.find(pdg);

  if (it == fNuclearFFs.end())
    return 0.;

  return it->second->Calculator().FormFactor(Q);
}
//____________________________________________________________________________
double DeVriesFormFactorMap::NeutronFF(double Q, int pdg) const {

  // DeVries form factor are measured from EM interactions
  // so they are the EM charge distribution inside the nucleus
  // The neutron charge distribution is assumed to be the same
  // as the proton component, with a normalization as the total weak charge is
  // the number of neutron, not the number of protons.
  // In Fourier transform, the normalization is easy as
  // FF(Q=0) = total charge = Z for proton or (A-Z) for Neutron

  int z = pdg::IonPdgCodeToZ(pdg);
  double scale = (pdg::IonPdgCodeToA(pdg) - z) / (double)z;

  return scale * ProtonFF(Q, pdg);
}
//____________________________________________________________________________
genie::Range1D_t DeVriesFormFactorMap::QRange(int pdg) const {

  const std::map<int, const genie::DeVriesFormFactor *>::const_iterator it =
      fNuclearFFs.find(pdg);

  if (it == fNuclearFFs.end())
    return COHFormFactorI::QRange(pdg);

  return Range1D_t(it->second->Calculator().QMin(),
                   it->second->Calculator().QMax());
}
//____________________________________________________________________________
bool DeVriesFormFactorMap::HasNucleus(int pdg) const {

  return (fNuclearFFs.count(pdg) > 0);
}
//____________________________________________________________________________
void DeVriesFormFactorMap::LoadConfig(void) {

  fNuclearFFs.clear();

  bool good_configuration = true;

  // read the vector of algos from the xml file
  std::vector<RgKey> keys;
  GetParamVectKeys("COH-DV-FormFactor", keys);

  // Store pointers to subalgos in the local map
  for (unsigned int i = 0; i < keys.size(); ++i) {

    const Algorithm *algo = SubAlg(keys[i]);

    const DeVriesFormFactor *ff = dynamic_cast<const DeVriesFormFactor *>(algo);

    if (!ff) {
      good_configuration = false;
      LOG("DeVriesFormFactorMap", pERROR)
          << "SubAlgo with key " << keys[i] << " not retrieved";
    }

    if (fNuclearFFs.count(ff->NucleusPDG()) > 0) {

      good_configuration = false;
      LOG("DeVriesFormFactorMap", pERROR)
          << "Attempt to add a second DeVries form factor for PDG "
          << ff->NucleusPDG();
    }

    fNuclearFFs[ff->NucleusPDG()] = ff;

  } // loop over subalgo

  if (!good_configuration) {
    LOG("DeVriesFormFactorMap", pFATAL) << "Configuration not good, exiting";
    exit(78);
  }
}
//____________________________________________________________________________
