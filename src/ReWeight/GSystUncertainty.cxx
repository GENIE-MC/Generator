//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "ReWeight/GSystUncertainty.h"
#include "ReWeight/GSystConst.h"

using namespace genie;
using namespace genie::rew;

GSystUncertainty * GSystUncertainty::fInstance = 0;
//____________________________________________________________________________
GSystUncertainty::GSystUncertainty()
{
//  fInstance = 0;
}
//____________________________________________________________________________
GSystUncertainty::~GSystUncertainty()
{
  fInstance = 0;
}
//____________________________________________________________________________
GSystUncertainty * GSystUncertainty::Instance()
{
  if(fInstance == 0) {
    LOG("ReW", pINFO) << "GSystUncertainty late initialization";
    static GSystUncertainty::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new GSystUncertainty;
  }
  return fInstance;
}
//____________________________________________________________________________
double GSystUncertainty::OneSigmaErr(GSyst_t syst) const
{
   switch(syst) {
     case ( kSystNuXSec_MaQEL ) : 
       return kNuXSec_MaQEL_1SigmaErr;
       break;
     case ( kSystNuXSec_MvQEL ) : 
       return kNuXSec_MvQEL_1SigmaErr;
       break;
     case ( kSystNuXSec_MaRES ) : 
       return kNuXSec_MaRES_1SigmaErr;
       break;
     case ( kSystNuXSec_MvRES ) : 
       return kNuXSec_MvRES_1SigmaErr;
       break;
     case ( kSystNuXSec_MaCOHPi ) : 
       return kNuXSec_MaCOHPi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvpCC1pi ) : 
       return kNuXSec_RvpCC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvpCC2pi ) : 
       return kNuXSec_RvpCC2pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvpNC1pi ) : 
       return kNuXSec_RvpNC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvpNC2pi ) : 
       return kNuXSec_RvpNC2pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvnCC1pi ) : 
       return kNuXSec_RvnCC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvnCC2pi ) : 
       return kNuXSec_RvnCC2pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvnNC1pi ) : 
       return kNuXSec_RvnNC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvnNC2pi ) : 
       return kNuXSec_RvnNC2pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarpCC1pi ) : 
       return kNuXSec_RvbarpCC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarpCC2pi ) : 
       return kNuXSec_RvbarpCC2pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarpNC1pi ) : 
       return kNuXSec_RvbarpNC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarpNC2pi ) : 
       return kNuXSec_RvbarpNC2pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarnCC1pi ) : 
       return kNuXSec_RvbarnCC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarnCC2pi ) : 
       return kNuXSec_RvbarnCC2pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarnNC1pi ) : 
       return kNuXSec_RvbarnNC1pi_1SigmaErr;
       break;
     case ( kSystNuXSec_RvbarnNC2pi ) : 
       return kNuXSec_RvbarnNC2pi_1SigmaErr;
       break;
     case ( kSystINuke_MFPTwk_pi ) : 
       return kINuke_MFPTwk_pi_1SigmaErr;
       break;
     case ( kSystINuke_MFPTwk_N  ) : 
       return kINuke_MFPTwk_N_1SigmaErr;
       break;
     case ( kSystINuke_CExTwk_pi ) : 
       return kINuke_CExTwk_pi_1SigmaErr;
       break;
     case ( kSystINuke_ElTwk_pi ) : 
       return kINuke_ElTwk_pi_1SigmaErr;
       break;
     case ( kSystINuke_InelTwk_pi ) : 
       return kINuke_InelTwk_pi_1SigmaErr;
       break;
     case ( kSystINuke_AbsTwk_pi ) : 
       return kINuke_AbsTwk_pi_1SigmaErr;
       break;
     case ( kSystINuke_PiProdTwk_pi ) : 
       return kINuke_PiProdTwk_pi_1SigmaErr;
       break;
     case ( kSystINuke_CExTwk_N ) : 
       return kINuke_CExTwk_N_1SigmaErr;
       break;
     case ( kSystINuke_ElTwk_N ) : 
       return kINuke_ElTwk_N_1SigmaErr;
       break;
     case ( kSystINuke_InelTwk_N ) : 
       return kINuke_InelTwk_N_1SigmaErr;
       break;
     case ( kSystINuke_AbsTwk_N ) : 
       return kINuke_AbsTwk_N_1SigmaErr;
       break;
     case ( kSystINuke_PiProdTwk_N ) : 
       return kINuke_PiProdTwk_N_1SigmaErr;
       break;
     default: 
        return 0.;
        break;
  }
  return 0.;
}
//____________________________________________________________________________

