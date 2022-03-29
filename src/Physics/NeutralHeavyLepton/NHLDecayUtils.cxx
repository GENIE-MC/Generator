//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
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
using namespace genie::NHL;

//____________________________________________________________________________
string genie::utils::nhl::AsString(NHLDecayMode_t nhldm)
{
  switch(nhldm) {

  case (kNHLDcyNull):
    return "Unspecified";
    break;

  case (kNHLDcyPiMu):
    return "N -> pi+- mu-+";
    break;

  case (kNHLDcyPiE):
    return "N -> pi+- e-+";
    break;

  case (kNHLDcyPi0Nu):
    return "N -> pi0 v";
    break;

  case (kNHLDcyNuNuNu):
    return "N -> v v v";
    break;

  case (kNHLDcyNuMuMu):
    return "N -> v mu+ mu-";
    break;

  case (kNHLDcyNuEE):
    return "N -> v e+ e-";
    break;

  case (kNHLDcyNuMuE):
    return "N -> v mu+- e-+";
    break;

  case (kNHLDcyPiPi0E):
    return "N -> pi+- pi0 e-+";
    break;

  case (kNHLDcyPiPi0Mu):
    return "N -> pi+- pi0 mu-+";
    break;

  case (kNHLDcyPi0Pi0Nu):
    return "N -> v pi0 pi0";
    break;

  case (kNHLDcyTEST):
    return "N -> e+ e-";
    break;

  default:
    return "Invalid NHL decay mode!";
  }
  //return "Invalid NHL decay mode!";
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

  case(kNHLDcyPiMu):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgMuon);
    break;

  case(kNHLDcyPiE):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgElectron);
    break;

  case(kNHLDcyPi0Nu):
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgNuMu); // could be nue or nutau too!
    break;

  case(kNHLDcyNuNuNu):
    decay_products.push_back(kPdgNuE);
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgNuTau); // again, any permutation of {e,mu,tau}^3 works
    break;

  case(kNHLDcyNuMuMu):
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgAntiMuon);
    break;

  case(kNHLDcyNuEE):
    decay_products.push_back(kPdgNuE);
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgPositron);
    break;

  case(kNHLDcyNuMuE):
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgPositron);
    break;

  case(kNHLDcyPiPi0E):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgElectron);
    break;

  case(kNHLDcyPiPi0Mu):
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgMuon);
    break;

  case(kNHLDcyPi0Pi0Nu):
    decay_products.push_back(kPdgNuMu);
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgPi0);
    break;

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
