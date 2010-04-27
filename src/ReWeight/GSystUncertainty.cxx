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
 @ Apr 27, 2010 - CA
   Included new parameters in preparation for the Summer 2010 T2K analyses.
   Added option to override the default 1\sigma errors.
*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "ReWeight/GSystUncertainty.h"

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
double GSystUncertainty::OneSigmaErr(GSyst_t s) const
{
  map<GSyst_t,double>::const_iterator it = fOneSigmeErrMap.find(s);
  if(it != fOneSigmeErrMap.end()) return it->second;
  return 0;
}
//____________________________________________________________________________
void GSystUncertainty::OverrideDefaultUncertainty(GSyst_t s, double onesigerr)
{
  fOneSigmeErrMap[s] = onesigerr;
}
//____________________________________________________________________________
void GSystUncertainty::SetDefaults(void)
{
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_NormCCQE,       0.15));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MaCCQEshape,    0.10));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MaCCQE,         0.15));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MvCCQE,         0.05));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_NormCCRES,      0.20));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MaCCRESshape,   0.10));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MvCCRESshape,   0.05));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MaCCRES,        0.20));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MvCCRES,        0.10));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_MaCOHPi,        0.40));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvpCC1pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvpCC2pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvpNC1pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvpNC2pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvnCC1pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvnCC2pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvnNC1pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvnNC2pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarpCC1pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarpCC2pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarpNC1pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarpNC2pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarnCC1pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarnCC2pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarnNC1pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_RvbarnNC2pi,    0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystNuXSec_NormCCSafeDIS,  0.05));

  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystHadrnz_FormZone,       0.50));

  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_MFPTwk_pi,       0.20));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_MFPTwk_N,        0.20));

  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_CExTwk_pi,       0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_ElTwk_pi,        0.10));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_InelTwk_pi,      0.40));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_AbsTwk_pi,       0.30));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_PiProdTwk_pi,    0.20));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_CExTwk_N,        0.50));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_ElTwk_N,         0.30));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_InelTwk_N,       0.40));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_AbsTwk_N,        0.20));
  fOneSigmeErrMap.insert(map<GSyst_t,double>::value_type(kSystINuke_PiProdTwk_N,     0.20));
}
//____________________________________________________________________________
