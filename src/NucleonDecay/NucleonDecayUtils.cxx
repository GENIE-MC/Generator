//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 03, 2008 - CA
   First added in v2.7.1

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "GHEP/GHepParticle.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/NuclearUtils.h"
#include "Utils/PrintUtils.h"
#include "NucleonDecay/NucleonDecayUtils.h"

using namespace genie;
//using namespace genie::utils::nucleon_decay;

//____________________________________________________________________________
string genie::utils::nucleon_decay::AsString(NucleonDecayMode_t ndm)
{
  switch(ndm) {
    case (kNDp2epi0)     : return "p --> e+ pi0";      break;
    case (kNDp2mupi0)    : return "p --> mu+ pi0";     break;
    case (kNDp2eeta0)    : return "p --> e+ eta";      break;
    case (kNDp2mueta0)   : return "p --> mu+ eta";     break;
    case (kNDp2erho0)    : return "p --> e+ rho0";     break;
    case (kNDp2murho0)   : return "p --> mu+ pi0";     break;
    case (kNDp2eomega0)  : return "p --> e+ omega0";   break;
    case (kNDp2muomega0) : return "p --> mu+ omega0";  break;
    case (kNDn2epim)     : return "n --> e+ pi-";      break;
    case (kNDn2mupim)    : return "n --> mu+ pi-";     break;
    case (kNDp2nubarKp)  : return "p --> nubar K+";    break;
    default              : return "?";                 break;
  }
  return "??";
}
//____________________________________________________________________________
bool genie::utils::nucleon_decay::IsValidMode(NucleonDecayMode_t ndm)
{
  switch(ndm) {
    case (kNDp2epi0)     : 
    case (kNDp2mupi0)    : 
    case (kNDp2eeta0)    : 
    case (kNDp2mueta0)   : 
    case (kNDp2erho0)    : 
    case (kNDp2murho0)   : 
    case (kNDp2eomega0)  : 
    case (kNDp2muomega0) : 
    case (kNDn2epim)     : 
    case (kNDn2mupim)    : 
    case (kNDp2nubarKp)  : 
      return true;  
      break;
    default : 
      return false;  
      break;
  }
  return false;
}
//____________________________________________________________________________
int genie::utils::nucleon_decay::DecayedNucleonPdgCode(NucleonDecayMode_t ndm)
{
  switch(ndm) {
    case (kNDp2epi0)     : return kPdgProton;  break;
    case (kNDp2mupi0)    : return kPdgProton;  break;
    case (kNDp2eeta0)    : return kPdgProton;  break;
    case (kNDp2mueta0)   : return kPdgProton;  break;
    case (kNDp2erho0)    : return kPdgProton;  break;
    case (kNDp2murho0)   : return kPdgProton;  break;
    case (kNDp2eomega0)  : return kPdgProton;  break;
    case (kNDp2muomega0) : return kPdgProton;  break;
    case (kNDn2epim)     : return kPdgNeutron; break;
    case (kNDn2mupim)    : return kPdgNeutron; break;
    case (kNDp2nubarKp)  : return kPdgProton;  break;
    default              : return 0;           break;
  }
  return 0;
}
//____________________________________________________________________________
PDGCodeList genie::utils::nucleon_decay::DecayProductList(
  NucleonDecayMode_t ndm)
{
  bool allow_duplicate = true;
  PDGCodeList decay_products(allow_duplicate);

  switch(ndm) {
    case (kNDp2epi0) :
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgPi0); 
      break;
    case (kNDp2mupi0) : 
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgPi0); 
      break;
    case (kNDp2eeta0) : 
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgEta); 
      break;
    case (kNDp2mueta0) : 
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgEta); 
      break;
    case (kNDp2erho0) : 
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgRho0); 
      break;
    case (kNDp2murho0) : 
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgRho0); 
      break;
    case (kNDp2eomega0) : 
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgomega); 
      break;
    case (kNDp2muomega0) : 
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgomega); 
      break;
    case (kNDn2epim) : 
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgPiM); 
      break;
    case (kNDn2mupim) : 
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgPiM); 
      break;
    case (kNDp2nubarKp) : 
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgKP); 
      break;
    default : 
      break;
  }
  return decay_products;
}
//____________________________________________________________________________
GHepStatus_t genie::utils::nucleon_decay::DecayProductStatus(
  bool in_nucleus, int pdgc)
{
  if(in_nucleus) {
    if( pdgc == kPdgPi0   ||
        pdgc == kPdgPiM   ||
        pdgc == kPdgEta   ||
        pdgc == kPdgRho0  ||
        pdgc == kPdgomega ||
        pdgc == kPdgKP) 
    {
      return kIStHadronInTheNucleus;
    } 
  } 

  return kIStStableFinalState;
}
//____________________________________________________________________________
