//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NeutralHeavyLepton/NHLDecayUtils.h"

using namespace genie;

//____________________________________________________________________________
string genie::utils::nhl::AsString(NHLDecayMode_t nhldm)
{
  switch(nhldm) {

  case (kNHLDcyNull):
    return "Invalid NHL decay mode!";
    break;

// add all other cases and specify a string for the decay channel
//
//
  case(kNHLDcyTEST):
  return "N -> e+e-";

  }
  return "Invalid HL decay mode!";
}
//____________________________________________________________________________
bool genie::utils::nhl::IsKinematicallyAllowed(NHLDecayMode_t nhldm, double M)
{
// Check the input mass of the NHL (M) against the sum of the masses of the
// particles to be produced in the input decay

  PDGCodeList decay_products = DecayProductList(nhldm);

  PDGLibrary * pdglib = PDGLibrary::Instance();

  double Msum = 0.;
  PDGCodeList::const_iterator it = decay_products.begin();
  for ( ; it != decay_products.end(); ++it)
  {
    int pdg_code = *it;
    TParticlePDG * p = pdglib->Find(pdg_code);
    if(p) {
       Msum += p->Mass();
    } else {
       LOG("NHL", pERROR)
        << "Decay list includes particle with unrecognised PDG code: "
        << pdg_code;
    }
  }

  return (M > Msum);
}
//____________________________________________________________________________
PDGCodeList genie::utils::nhl::DecayProductList(NHLDecayMode_t nhldm)
{
  bool allow_duplicate = true;
  PDGCodeList decay_products(allow_duplicate);

  switch(nhldm) {

  // Specify the final state particles in each NHL decay channel,
  // adding sections that look like
  //
  //  case (kNHLDcy...):
  //     decay_products.push_back(kPdgPositron);
  //     decay_products.push_back(kPdgElectron);
  //     decay_products.push_back(... some other particle PDG code);
  //     break;
  //
  case(kNHLDcyTEST):
    decay_products.push_back(kPdgPositron);
    decay_products.push_back(kPdgElectron);
    break;

  default :
    break;
  }

  return decay_products;
}
//____________________________________________________________________________
