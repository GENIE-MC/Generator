//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
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
#include "NeutronOsc/NeutronOscUtils.h"

using namespace genie;
//using namespace genie::utils::neutron_osc;

//____________________________________________________________________________
string genie::utils::neutron_osc::AsString(NeutronOscMode_t ndm)
{
  // this just maps the decay mode integers to string descriptors. replaced. -j
  switch(ndm) {
    case (kNORandom)          :    return "Random mode";
                                   break;
    case (kNOpto1pip1pi0)     :    return "p + nbar --> pi+ pi0";
                                   break;
    case (kNOpto1pip2pi0)     :    return "p + nbar --> pi+ 2pi0";
                                   break;
    case (kNOpto1pip3pi0)     :    return "p + nbar --> pi+ 3pi0";
                                   break;
    case (kNOpto2pip1pim1pi0) :    return "p + nbar --> 2pi+ pi- pi0";
                                   break;
    case (kNOpto2pip1pim2pi0) :    return "p + nbar --> 2pi+ pi- 2pi0";
                                   break;
    case (kNOpto2pip1pim2o)   :    return "p + nbar --> 2pi+ pi- 2omega";
                                   break;
    case (kNOpto3pip2pim1pi0) :    return "p + nbar --> 3pi+ 2pi- pi0";
                                   break;
    case (kNOnto1pip1pim)     :    return "n + nbar --> pi+ pi-";
                                   break;
    case (kNOnto2pi0)         :    return "n + nbar --> 2pi0";
                                   break;
    case (kNOnto1pip1pim1pi0) :    return "n + nbar --> pi+ pi- pi0";
                                   break;
    case (kNOnto1pip1pim2pi0) :    return "n + nbar --> pi+ pi- 2pi0";
                                   break;
    case (kNOnto1pip1pim3pi0) :    return "n + nbar --> pi+ pi- 3pi0";
                                   break;
    case (kNOnto2pip2pim)     :    return "n + nbar --> 2pi+ 2pi-";
                                   break;
    case (kNOnto2pip2pim1pi0) :    return "n + nbar --> 2pi+ 2pi- pi0";
                                   break;
    case (kNOnto1pip1pim1o)   :    return "n + nbar --> pi+ pi- omega";
                                   break;
    case (kNOnto2pip2pim2pi0) :    return "n + nbar --> 2pi+ 2pi- 2pi0";
                                   break;
    default                   :    return "?";
                                   break;
  }
  return "??";
}
//____________________________________________________________________________
bool genie::utils::neutron_osc::IsValidMode(NeutronOscMode_t ndm)
{
  // checks if a mode is valid. just straight replaced. -j
  switch(ndm) {
    case (kNORandom)          :
    case (kNOpto1pip1pi0)     :
    case (kNOpto1pip2pi0)     :
    case (kNOpto1pip3pi0)     :
    case (kNOpto2pip1pim1pi0) :
    case (kNOpto2pip1pim2pi0) :
    case (kNOpto2pip1pim2o)   :
    case (kNOpto3pip2pim1pi0) :
    case (kNOnto1pip1pim)     :
    case (kNOnto2pi0)         :
    case (kNOnto1pip1pim1pi0) :
    case (kNOnto1pip1pim2pi0) :
    case (kNOnto1pip1pim3pi0) :
    case (kNOnto2pip2pim)     :
    case (kNOnto2pip2pim1pi0) :
    case (kNOnto1pip1pim1o)   :
    case (kNOnto2pip2pim2pi0) :
      return true;  
      break;
    default : 
      return false;  
      break;
  }
  return false;
}
//____________________________________________________________________________
int genie::utils::neutron_osc::AnnihilatingNucleonPdgCode(NeutronOscMode_t ndm)
{
  // name isn't really accurate any more. instead of decayed nucleon, function
  // returns what particle the oscillated neutron annihilated with -j
  switch(ndm) {
    case (kNOpto1pip1pi0)     : return kPdgProton;  break;
    case (kNOpto1pip2pi0)     : return kPdgProton;  break;
    case (kNOpto1pip3pi0)     : return kPdgProton;  break;
    case (kNOpto2pip1pim1pi0) : return kPdgProton;  break;
    case (kNOpto2pip1pim2pi0) : return kPdgProton;  break;
    case (kNOpto2pip1pim2o)   : return kPdgProton;  break;
    case (kNOpto3pip2pim1pi0) : return kPdgProton;  break;
    case (kNOnto1pip1pim)     : return kPdgNeutron; break;
    case (kNOnto2pi0)         : return kPdgNeutron; break;
    case (kNOnto1pip1pim1pi0) : return kPdgNeutron; break;
    case (kNOnto1pip1pim2pi0) : return kPdgNeutron; break;
    case (kNOnto1pip1pim3pi0) : return kPdgNeutron; break;
    case (kNOnto2pip2pim)     : return kPdgNeutron; break;
    case (kNOnto2pip2pim1pi0) : return kPdgNeutron; break;
    case (kNOnto1pip1pim1o)   : return kPdgNeutron; break;
    case (kNOnto2pip2pim2pi0) : return kPdgNeutron; break;
    default                   : return 0;           break;
  }
  return 0;
}
//____________________________________________________________________________
PDGCodeList genie::utils::neutron_osc::DecayProductList(
  NeutronOscMode_t ndm)
{
  // ok so i think this is the first function where a straight replacement
  // isn't gonna cut it. all the nucleon decay modes are two-body, but that is
  // painfully untrue for nnbar. i just threw aaaaaaallll of the final state
  // particles into the vector, so let's just hope for the best -j

  // need to implement a lorentz boost into rest frame of two nucleons -j

  bool allow_duplicate = true;
  PDGCodeList decay_products(allow_duplicate);

  switch(ndm) {
    case (kNOpto1pip1pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOpto1pip2pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOpto1pip3pi0) :
      decay_products.push_back(kPdgPiP); 
      decay_products.push_back(kPdgPi0); 
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOpto2pip1pim1pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOpto2pip1pim2pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOpto2pip1pim2o) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgomega);
      decay_products.push_back(kPdgomega);
      break;
    case (kNOpto3pip2pim1pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiP); 
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOnto1pip1pim) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      break;
    case (kNOnto2pi0) :
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOnto1pip1pim1pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOnto1pip1pim2pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOnto1pip1pim3pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOnto2pip2pim) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPiM);
      break;
    case (kNOnto2pip2pim1pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      break;
    case (kNOnto1pip1pim1o) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgomega);
      break;
    case (kNOnto2pip2pim2pi0) :
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiP);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPiM);
      decay_products.push_back(kPdgPi0);
      decay_products.push_back(kPdgPi0);
      break;
    default :
      break;
  }
  return decay_products;
}
//____________________________________________________________________________
GHepStatus_t genie::utils::neutron_osc::DecayProductStatus(
  bool in_nucleus, int pdgc)
{
  // took out all the irrelevant particles -j
  if(in_nucleus) {
    if( pdgc == kPdgPi0   ||
        pdgc == kPdgPiM   ||
        pdgc == kPdgPiP   ||
        pdgc == kPdgomega)
    {
      return kIStHadronInTheNucleus;
    } 
  } 

  return kIStStableFinalState;
}
//____________________________________________________________________________
