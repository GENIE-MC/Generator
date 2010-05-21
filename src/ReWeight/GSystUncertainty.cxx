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
    fInstance->SetDefaults();
  }
  return fInstance;
}
//____________________________________________________________________________
double GSystUncertainty::OneSigmaErr(GSyst_t s) const
{
  map<GSyst_t,double>::const_iterator it = fOneSigErrMap.find(s);
  if(it != fOneSigErrMap.end()) return it->second;
  return 0;
}
//____________________________________________________________________________
void GSystUncertainty::OverrideDefaultUncertainty(GSyst_t s, double onesigerr)
{
  fOneSigErrMap[s] = onesigerr;
}
//____________________________________________________________________________
void GSystUncertainty::SetDefaults(void)
{
  map<GSyst_t, double> & m = fOneSigErrMap;

  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_NormCCQE,       0.15));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_MaCCQEshape,    0.10));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_MaCCQE,         0.15));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_NormCCRES,      0.20));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_MaCCRESshape,   0.10));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_MvCCRESshape,   0.05));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_MaCCRES,        0.20));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_MvCCRES,        0.10));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_MaCOHpi,        0.40));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_R0COHpi,        0.20));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvpCC1pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvpCC2pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvpNC1pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvpNC2pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvnCC1pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvnCC2pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvnNC1pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvnNC2pi,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarpCC1pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarpCC2pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarpNC1pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarpNC2pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarnCC1pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarnCC2pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarnNC1pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_RvbarnNC2pi,    0.50));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_NormCCSafeDIS,  0.05));
  m.insert(map<GSyst_t,double>::value_type(kXSecTwkDial_DISNuclMod,     1.00));
  m.insert(map<GSyst_t,double>::value_type(kSystNucl_CCQEPauliSupViaKF, 0.05));
  m.insert(map<GSyst_t,double>::value_type(kHadrAGKYTwkDial_xF1pi,      0.20));
  m.insert(map<GSyst_t,double>::value_type(kHadrAGKYTwkDial_pT1pi,      0.03));
  m.insert(map<GSyst_t,double>::value_type(kHadroTwkDial_FormZone,      0.50));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_MFP_pi,        0.20));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_MFP_N,         0.20));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrCEx_pi,      0.50));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrElas_pi,     0.10));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrInel_pi,     0.40));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrAbs_pi,      0.30));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrPiProd_pi,   0.20));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrCEx_N,       0.50));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrElas_N,      0.30));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrInel_N,      0.40));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrAbs_N,       0.20));
  m.insert(map<GSyst_t,double>::value_type(kINukeTwkDial_FrPiProd_N,    0.20));
}
//____________________________________________________________________________
