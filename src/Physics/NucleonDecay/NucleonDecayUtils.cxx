//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
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

#include "Framework/Messenger/Messenger.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NucleonDecay/NucleonDecayUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
//using namespace genie::utils::nucleon_decay;

//____________________________________________________________________________
string genie::utils::nucleon_decay::AsString(NucleonDecayMode_t ndm,
					     int npdg)
{
  switch(ndm) {
  
  case (kNDNull):
    return "Invalid nucleon decay mode!";
    break;

  case (kNDN2eppi): 
    if (npdg == kPdgProton) {
      return "p --> e+ pi0";
    } else if (npdg == kPdgNeutron) {
      return "n --> e+ pi-";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDN2muppi): 
    if (npdg == kPdgProton) {
      return "p --> mu+ pi0";
    } else if (npdg == kPdgNeutron) {
      return "n --> mu+ pi-";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDN2nubarpi): 
    if (npdg == kPdgProton) {
      return "p --> nubar pi+";
    } else if (npdg == kPdgNeutron) {
      return "n --> nubar pi0";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDp2epeta):
    return "p --> e+ eta";
    break;
      
  case (kNDp2mupeta):
    return "p --> mu+ eta";
    break;

  case (kNDn2nubareta):
    return "n --> nubar eta";
    break;

  case (kNDN2eprho): 
    if (npdg == kPdgProton) {
      return "p --> e+ rho0";
    } else if (npdg == kPdgNeutron) {
      return "n --> e+ rho-";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDN2muprho): 
    if (npdg == kPdgProton) {
      return "p --> mu+ rho0";
    } else if (npdg == kPdgNeutron) {
      return "n --> mu+ rho-";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDN2nubarrho): 
    if (npdg == kPdgProton) {
      return "p --> nubar rho+";
    } else if (npdg == kPdgNeutron) {
      return "n --> nubar rho0";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDp2epomega):
    return "p --> e+ omega";
    break;

  case (kNDp2mupomega):
    return "p --> mu+ omega";
    break;

  case (kNDn2nubaromega):
    return "n --> nubar omega";
    break;

  case (kNDN2epK): 
    if (npdg == kPdgProton) {
      return "p --> e+ K0";
    } else if (npdg == kPdgNeutron) {
      return "n --> e+ K-";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDp2epK0s):
    return "p --> e+ K0s";
    break;

  case (kNDp2epK0l):
    return "p --> e+ K0l";
    break;

  case (kNDN2mupK): 
    if (npdg == kPdgProton) {
      return "p --> mu+ K0";
    } else if (npdg == kPdgNeutron) {
      return "n --> mu+ K-";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDp2mupK0s):
    return "p --> mu+ K0s";
    break;

  case (kNDp2mupK0l):
    return "p --> mu+ K0l";
    break;

  case (kNDN2nubarK):
    if (npdg == kPdgProton) {
      return "p --> nubar K+";
    } else if (npdg == kPdgNeutron) {
      return "n --> nubar K0";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDn2nubarK0s):
    return "n --> nubar K0s";
    break;

  case (kNDp2epKstar0):
    return "p --> e+ K*0";
    break;
      
  case (kNDN2nubarKstar):
    if (npdg == kPdgProton) {
      return "p --> nubar K*+";
    } else if (npdg == kPdgNeutron) {
      return "n --> nubar K*0";
    } else {
      return "Invalid nucleon decay mode!";
    }
    break;

  case (kNDp2eppippim):
    return "p --> e+ pi+ pi-";
    break;

  case (kNDp2eppi0pi0):
    return "p --> e+ pi0 pi0";
    break;

  case (kNDn2eppimpi0):
    return "n --> e+ pi- pi0";
    break;

  case (kNDp2muppippim):
    return "p --> mu+ pi+ pi-";
    break;

  case (kNDp2muppi0pi0):
    return "p --> mu+ pi0 pi0";
    break;

  case (kNDn2muppimpi0):
    return "n --> mu+ pi- pi0";
    break;

  case (kNDn2epK0pim):
    return "n --> e+ K0 pi-";
    break;

  case (kNDn2empip):
    return "n --> e- pi+";
    break;

  case (kNDn2mumpip):
    return "n --> mu- pi+";
    break;

  case (kNDn2emrhop):
    return "n --> e- rho+";
    break;

  case (kNDn2mumrhop):
    return "n --> mu- rho+";
    break;

  case (kNDn2emKp):
    return "n --> e- K+";
    break;

  case (kNDn2mumKp):
    return "n --> mu- K+";
    break;

  case (kNDp2empippip):
    return "p --> e- pi+ pi+";
    break;

  case (kNDn2empippi0):
    return "n --> e- pi+ pi0";
    break;

  case (kNDp2mumpippip):
    return "p --> mu- pi+ pi+";
    break;

  case (kNDn2mumpippi0):
    return "n --> mu- pi+ pi0";
    break;

  case (kNDp2empipKp):
    return "p --> e- pi+ K+";
    break;

  case (kNDp2mumpipKp):
    return "p --> mu- pi+ K+";
    break;

  case (kNDp2epgamma):
    return "p --> e+ gamma";
    break;

  case (kNDp2mupgamma):
    return "p --> mu+ gamma";
    break;

  case (kNDn2nubargamma):
    return "n --> nubar gamma";
    break;

  case (kNDp2epgammagamma):
    return "p --> e+ gamma gamma";
    break;

  case (kNDn2nubargammagamma):
    return "n --> nubar gamma gamma";
    break;

  case (kNDp2epepem):
    return "p --> e+ e+ e-";
    break;

  case (kNDp2epmupmum):
    return "p --> e+ mu+ mu-";
    break;

  case (kNDp2epnubarnu):
    return "p --> e+ nubar nu";
    break;

  case (kNDn2epemnubar):
    return "n --> e+ e- nubar";
    break;

  case (kNDn2mupemnubar):
    return "n --> mu+ e- nubar";
    break;

  case (kNDn2mupmumnubar):
    return "n --> mu+ mu- nubar";
    break;

  case (kNDp2mupepem):
    return "p --> mu+ e+ e-";
    break;

  case (kNDp2mupmupmum):
    return "p --> mu+ mu+ mu-";
    break;

  case (kNDp2mupnubarnu):
    return "p --> mu+ nubar nu";
    break;

  case (kNDp2emmupmup):
    return "p --> e- mu+ mu+";
    break;

  case (kNDn2threenus):
    return "n --> nubar nubar nu";
    break;

  case (kNDn2fivenus):
    return "n --> nubar nubar nubar nu nu";
    break;

  }
  return "Invalid nucleon decay mode!";
}
//____________________________________________________________________________
bool genie::utils::nucleon_decay::IsValidMode(NucleonDecayMode_t ndm,
					      int npdg)
{
  switch(ndm) {
   
  case (kNDN2eppi):
  case (kNDN2muppi):
  case (kNDN2nubarpi):
  case (kNDN2eprho):
  case (kNDN2muprho):
  case (kNDN2nubarrho):
  case (kNDN2epK):
  case (kNDN2mupK):
  case (kNDN2nubarK):
  case (kNDN2nubarKstar):
    if (npdg == kPdgProton || npdg == kPdgNeutron) {
      return true;  
    } else {
      return false;
    }
  break;

  case (kNDp2epeta):
  case (kNDp2mupeta):
  case (kNDp2epomega):
  case (kNDp2mupomega):
  case (kNDp2epK0s):
  case (kNDp2epK0l):
  case (kNDp2mupK0s):
  case (kNDp2mupK0l):
  case (kNDp2epKstar0):
  case (kNDp2eppippim):
  case (kNDp2eppi0pi0):
  case (kNDp2muppippim):
  case (kNDp2muppi0pi0):
  case (kNDp2empippip):
  case (kNDp2mumpippip):
  case (kNDp2empipKp):
  case (kNDp2mumpipKp):
  case (kNDp2epgamma):
  case (kNDp2mupgamma):
  case (kNDp2epgammagamma):
  case (kNDp2epepem):
  case (kNDp2epmupmum):
  case (kNDp2epnubarnu):
  case (kNDp2mupepem):
  case (kNDp2mupmupmum):
  case (kNDp2mupnubarnu):
  case (kNDp2emmupmup):

    if (npdg == kPdgProton || npdg == 0) {
      return true;  
    } else {
      return false;
    }
  break;

  case (kNDn2nubareta):
  case (kNDn2nubaromega):
  case (kNDn2nubarK0s):
  case (kNDn2eppimpi0):
  case (kNDn2muppimpi0):
  case (kNDn2epK0pim):
  case (kNDn2empip):
  case (kNDn2mumpip):
  case (kNDn2emrhop):
  case (kNDn2mumrhop):
  case (kNDn2emKp):
  case (kNDn2mumKp):
  case (kNDn2empippi0):
  case (kNDn2mumpippi0):
  case (kNDn2nubargamma):
  case (kNDn2nubargammagamma):
  case (kNDn2epemnubar):
  case (kNDn2mupemnubar):
  case (kNDn2mupmumnubar):
  case (kNDn2threenus):
  case (kNDn2fivenus):

    if (npdg == kPdgNeutron || npdg == 0) {
      return true;  
    } else {
      return false;
    }
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
    
  case (kNDp2epeta)           : return kPdgProton;  break;
  case (kNDp2mupeta)          : return kPdgProton;  break;
  case (kNDn2nubareta)        : return kPdgNeutron; break;
  case (kNDp2epomega)         : return kPdgProton;  break;
  case (kNDp2mupomega)        : return kPdgProton;  break;
  case (kNDn2nubaromega)      : return kPdgNeutron; break;
  case (kNDp2epK0s)           : return kPdgProton;  break;
  case (kNDp2epK0l)           : return kPdgProton;  break;
  case (kNDp2mupK0s)          : return kPdgProton;  break;
  case (kNDp2mupK0l)          : return kPdgProton;  break;
  case (kNDn2nubarK0s)        : return kPdgNeutron; break;
  case (kNDp2epKstar0)        : return kPdgProton;  break;
  case (kNDp2eppippim)        : return kPdgProton;  break;
  case (kNDp2eppi0pi0)        : return kPdgProton;  break;
  case (kNDn2eppimpi0)        : return kPdgNeutron; break;
  case (kNDp2muppippim)       : return kPdgProton;  break;
  case (kNDp2muppi0pi0)       : return kPdgProton;  break;
  case (kNDn2muppimpi0)       : return kPdgNeutron; break;
  case (kNDn2epK0pim)         : return kPdgNeutron; break;
  case (kNDn2empip)           : return kPdgNeutron; break;
  case (kNDn2mumpip)          : return kPdgNeutron; break;
  case (kNDn2emrhop)          : return kPdgNeutron; break;
  case (kNDn2mumrhop)         : return kPdgNeutron; break;
  case (kNDn2emKp)            : return kPdgNeutron; break;
  case (kNDn2mumKp)           : return kPdgNeutron; break;
  case (kNDp2empippip)        : return kPdgProton;  break;
  case (kNDn2empippi0)        : return kPdgNeutron; break;
  case (kNDp2mumpippip)       : return kPdgProton;  break;
  case (kNDn2mumpippi0)       : return kPdgNeutron; break;
  case (kNDp2empipKp)         : return kPdgProton;  break;
  case (kNDp2mumpipKp)        : return kPdgProton;  break;
  case (kNDp2epgamma)         : return kPdgProton;  break;
  case (kNDp2mupgamma)        : return kPdgProton;  break;
  case (kNDn2nubargamma)      : return kPdgNeutron; break;
  case (kNDp2epgammagamma)    : return kPdgProton;  break;
  case (kNDn2nubargammagamma) : return kPdgNeutron; break;
  case (kNDp2epepem)          : return kPdgProton;  break;
  case (kNDp2epmupmum)        : return kPdgProton;  break;
  case (kNDp2epnubarnu)       : return kPdgProton;  break;
  case (kNDn2epemnubar)       : return kPdgNeutron; break;
  case (kNDn2mupemnubar)      : return kPdgNeutron; break;
  case (kNDn2mupmumnubar)     : return kPdgNeutron; break;
  case (kNDp2mupepem)         : return kPdgProton;  break;
  case (kNDp2mupmupmum)       : return kPdgProton;  break;
  case (kNDp2mupnubarnu)      : return kPdgProton;  break;
  case (kNDp2emmupmup)        : return kPdgProton;  break;
  case (kNDn2threenus)        : return kPdgNeutron; break;
  case (kNDn2fivenus)         : return kPdgNeutron; break; 

  default               : return 0;           break;
  }
  return 0;
}

//____________________________________________________________________________
PDGCodeList genie::utils::nucleon_decay::DecayProductList(NucleonDecayMode_t ndm, 
							  int npdg)
{
  bool allow_duplicate = true;
  PDGCodeList decay_products(allow_duplicate);

  switch(ndm) {
  
  case (kNDN2eppi):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgPi0); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgPiM); 
    }
    break;  
    
  case (kNDN2muppi):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgPi0); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgPiM); 
    }
    break;  
    
  case (kNDN2nubarpi):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgPiP); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgPi0); 
    }
    break;  

  case (kNDp2epeta):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgEta); 
    break;    

  case (kNDp2mupeta):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgEta); 
    break;    

  case (kNDn2nubareta):
    decay_products.push_back(kPdgAntiNuE); 
    decay_products.push_back(kPdgEta); 
    break;    

  case (kNDN2eprho):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgRho0); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgRhoM); 
    }
    break;  

  case (kNDN2muprho):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgRho0); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgRhoM); 
    }
    break;  

  case (kNDN2nubarrho):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgRhoP); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgRho0); 
    }
    break;  

  case (kNDp2epomega):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgomega); 
    break;    

  case (kNDp2mupomega):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgomega); 
    break;    

  case (kNDn2nubaromega):
    decay_products.push_back(kPdgAntiNuE); 
    decay_products.push_back(kPdgomega); 
    break;    

  case (kNDN2epK):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgK0); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgPositron); 
      decay_products.push_back(kPdgKM); 
    }
    break;  

  case (kNDp2epK0s):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgK0S); 
    break;   

  case (kNDp2epK0l):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgK0L); 
    break;   

  case (kNDN2mupK):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgK0); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgAntiMuon); 
      decay_products.push_back(kPdgKM); 
    }
    break;  

  case (kNDp2mupK0s):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgK0S); 
    break;   

  case (kNDp2mupK0l):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgK0L); 
    break;   

  case (kNDN2nubarK):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgKP); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgK0); 
    }
    break;  

  case (kNDn2nubarK0s):
    decay_products.push_back(kPdgAntiNuE); 
    decay_products.push_back(kPdgK0S); 
    break;   

  case (kNDp2epKstar0):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgKStar0); 
    break;   

  case (kNDN2nubarKstar):
    if (npdg == kPdgProton) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgKStarP); 
    } else if (npdg == kPdgNeutron) {
      decay_products.push_back(kPdgAntiNuE); 
      decay_products.push_back(kPdgKStar0); 
    }
    break;  

  case (kNDp2eppippim):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPiM);
    break;   

  case (kNDp2eppi0pi0):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgPi0);
    break;   

  case (kNDn2eppimpi0):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgPiM);
    decay_products.push_back(kPdgPi0);
    break;   

  case (kNDp2muppippim):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPiM);
    break;   

  case (kNDp2muppi0pi0):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgPi0);
    decay_products.push_back(kPdgPi0);
    break;   

  case (kNDn2muppimpi0):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgPiM);
    decay_products.push_back(kPdgPi0);
    break;   

  case (kNDn2epK0pim):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgK0);
    decay_products.push_back(kPdgPiM);
    break;   

  case (kNDn2empip):
    decay_products.push_back(kPdgElectron); 
    decay_products.push_back(kPdgPiP);
    break;   

  case (kNDn2mumpip):
    decay_products.push_back(kPdgMuon); 
    decay_products.push_back(kPdgPiP);
    break;   

  case (kNDn2emrhop):
    decay_products.push_back(kPdgElectron); 
    decay_products.push_back(kPdgRhoP);
    break;   

  case (kNDn2mumrhop):
    decay_products.push_back(kPdgMuon); 
    decay_products.push_back(kPdgRhoP);
    break;   

  case (kNDn2emKp):
    decay_products.push_back(kPdgElectron); 
    decay_products.push_back(kPdgKP);
    break;   

  case (kNDn2mumKp):
    decay_products.push_back(kPdgMuon); 
    decay_products.push_back(kPdgKP);
    break;   

  case (kNDp2empippip):
    decay_products.push_back(kPdgElectron); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPiP);
    break;   

  case (kNDn2empippi0):
    decay_products.push_back(kPdgElectron); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPi0);
    break;   

  case (kNDp2mumpippip):
    decay_products.push_back(kPdgMuon); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPiP);
    break;   

  case (kNDn2mumpippi0):
    decay_products.push_back(kPdgMuon); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgPi0);
    break;   

  case (kNDp2empipKp):
    decay_products.push_back(kPdgElectron); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgKP);
    break;   

  case (kNDp2mumpipKp):
    decay_products.push_back(kPdgMuon); 
    decay_products.push_back(kPdgPiP);
    decay_products.push_back(kPdgKP);
    break;   

  case (kNDp2epgamma):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgGamma);
    break;  

  case (kNDp2mupgamma):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgGamma);
    break;   

  case (kNDn2nubargamma):
    decay_products.push_back(kPdgAntiNuE); 
    decay_products.push_back(kPdgGamma);
    break;  

  case (kNDp2epgammagamma):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgGamma);
    decay_products.push_back(kPdgGamma);
    break;  

  case (kNDn2nubargammagamma):
    decay_products.push_back(kPdgAntiNuE); 
    decay_products.push_back(kPdgGamma);
    decay_products.push_back(kPdgGamma);
    break;  

  case (kNDp2epepem):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgPositron);
    decay_products.push_back(kPdgElectron);
    break; 

  case (kNDp2epmupmum):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgAntiMuon);
    decay_products.push_back(kPdgMuon);
    break; 

  case (kNDp2epnubarnu):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgAntiNuE);
    decay_products.push_back(kPdgNuE);
    break;

  case (kNDn2epemnubar):
    decay_products.push_back(kPdgPositron); 
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgAntiNuE);
    break;

  case (kNDn2mupemnubar):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgElectron);
    decay_products.push_back(kPdgAntiNuE);
    break;

  case (kNDn2mupmumnubar):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgMuon);
    decay_products.push_back(kPdgAntiNuE);
    break;

  case (kNDp2mupepem):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgPositron);
    decay_products.push_back(kPdgElectron);
    break;

  case (kNDp2mupmupmum):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgAntiMuon);
    decay_products.push_back(kPdgMuon);
    break;

  case (kNDp2mupnubarnu):
    decay_products.push_back(kPdgAntiMuon); 
    decay_products.push_back(kPdgAntiNuE);
    decay_products.push_back(kPdgNuE);
    break;

  case (kNDp2emmupmup):
    decay_products.push_back(kPdgElectron); 
    decay_products.push_back(kPdgAntiMuon);
    decay_products.push_back(kPdgAntiMuon);
    break;

  case (kNDn2threenus):
    decay_products.push_back(kPdgAntiNuE); 
    decay_products.push_back(kPdgAntiNuE);
    decay_products.push_back(kPdgNuE);
    break;

  case (kNDn2fivenus):
    decay_products.push_back(kPdgAntiNuE); 
    decay_products.push_back(kPdgAntiNuE);
    decay_products.push_back(kPdgAntiNuE);
    decay_products.push_back(kPdgNuE);
    decay_products.push_back(kPdgNuE);
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
    if( pdg::IsHadron(pdgc) )
    {
      return kIStHadronInTheNucleus;
    } 
  } 

  return kIStStableFinalState;
}
//____________________________________________________________________________
