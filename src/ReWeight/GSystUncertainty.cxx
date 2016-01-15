//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Apr 27, 2010 - CA
   Included new parameters in preparation for the Summer 2010 T2K analyses.
   Added option to override the default 1\sigma errors.
 @ Nov 25, 2010 - CA
   Allow for asymmetric 1 sigma fractional errors.
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
double GSystUncertainty::OneSigmaErr(GSyst_t s, int sign) const
{
  if(sign > 0) {
    map<GSyst_t,double>::const_iterator it = fOneSigPlusErrMap.find(s);
    if(it != fOneSigPlusErrMap.end()) return it->second;
    return 0;
  } 
  else 
  if(sign < 0) {
    map<GSyst_t,double>::const_iterator it = fOneSigMnusErrMap.find(s);
    if(it != fOneSigMnusErrMap.end()) return it->second;
    return 0;
  } 
  else {
    // Handle default argument (sign=0)
    // Case added for compatibility purposes since most existing weight 
    // calcutators call GSystUncertainty::OneSigmaErr(GSyst_t) and the error 
    // on most GSyst_t params is symmetric.
    double err = 0.5 * (
        this->OneSigmaErr(s, +1) + this->OneSigmaErr(s, -1));
    return err;
  }
}
//____________________________________________________________________________
void GSystUncertainty::SetUncertainty(
   GSyst_t s, double plus_err, double minus_err)
{
  //fOneSigPlusErrMap.insert( map<GSyst_t,double>::value_type(s, plus_err ) );
  //fOneSigMnusErrMap.insert( map<GSyst_t,double>::value_type(s, minus_err) );
  fOneSigPlusErrMap[s] = plus_err;
  fOneSigMnusErrMap[s] =  minus_err;
}
//____________________________________________________________________________
void GSystUncertainty::SetDefaults(void)
{
  this->SetUncertainty( kXSecTwkDial_MaNCEL,         0.25, 0.25);
  this->SetUncertainty( kXSecTwkDial_EtaNCEL,        0.30, 0.30);
  this->SetUncertainty( kXSecTwkDial_NormCCQE,       0.20, 0.15);
  this->SetUncertainty( kXSecTwkDial_MaCCQEshape,    0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_MaCCQE,         0.25, 0.15);
  this->SetUncertainty( kXSecTwkDial_NormCCRES,      0.20, 0.20);
  this->SetUncertainty( kXSecTwkDial_MaCCRESshape,   0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_MvCCRESshape,   0.05, 0.05);
  this->SetUncertainty( kXSecTwkDial_MaCCRES,        0.20, 0.20);
  this->SetUncertainty( kXSecTwkDial_MvCCRES,        0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_NormNCRES,      0.20, 0.20);
  this->SetUncertainty( kXSecTwkDial_MaNCRESshape,   0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_MvNCRESshape,   0.05, 0.05);
  this->SetUncertainty( kXSecTwkDial_MaNCRES,        0.20, 0.20);
  this->SetUncertainty( kXSecTwkDial_MvNCRES,        0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_MaCOHpi,        0.40, 0.40);
  this->SetUncertainty( kXSecTwkDial_R0COHpi,        0.10, 0.10);
  this->SetUncertainty( kXSecTwkDial_RvpCC1pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvpCC2pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvpNC1pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvpNC2pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvnCC1pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvnCC2pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvnNC1pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvnNC2pi,       0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarpCC1pi,    0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarpCC2pi,    0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarpNC1pi,    0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarpNC2pi,    0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarnCC1pi,    0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarnCC2pi,    0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarnNC1pi,    0.50, 0.50);
  this->SetUncertainty( kXSecTwkDial_RvbarnNC2pi,    0.50, 0.50);

  // From Debdatta's thesis: 
  //   Aht  = 0.538 +/- 0.134	
  //   Bht  = 0.305 +/- 0.076
  //   CV1u = 0.291 +/- 0.087
  //   CV2u = 0.189 +/- 0.076

  this->SetUncertainty( kXSecTwkDial_AhtBY,          0.25, 0.25);
  this->SetUncertainty( kXSecTwkDial_BhtBY,          0.25, 0.25);
  this->SetUncertainty( kXSecTwkDial_CV1uBY,         0.30, 0.30);
  this->SetUncertainty( kXSecTwkDial_CV2uBY,         0.40, 0.40);

  this->SetUncertainty( kXSecTwkDial_AhtBYshape,     0.25, 0.25);
  this->SetUncertainty( kXSecTwkDial_BhtBYshape,     0.25, 0.25);
  this->SetUncertainty( kXSecTwkDial_CV1uBYshape,    0.30, 0.30);
  this->SetUncertainty( kXSecTwkDial_CV2uBYshape,    0.40, 0.40);

  this->SetUncertainty( kXSecTwkDial_DISNuclMod,     1.00, 1.00);
  this->SetUncertainty( kSystNucl_CCQEPauliSupViaKF, 0.30, 0.30);
  this->SetUncertainty( kHadrAGKYTwkDial_xF1pi,      0.20, 0.20);
  this->SetUncertainty( kHadrAGKYTwkDial_pT1pi,      0.03, 0.03);
  this->SetUncertainty( kHadrNuclTwkDial_FormZone,   0.50, 0.50);

  // From INTRANUKE pi+A and N+A mode comparisons with hadron scattering data:
  //
  this->SetUncertainty( kINukeTwkDial_MFP_pi,        0.20, 0.20);
  this->SetUncertainty( kINukeTwkDial_MFP_N,         0.20, 0.20);
  this->SetUncertainty( kINukeTwkDial_FrCEx_pi,      0.50, 0.50);
  this->SetUncertainty( kINukeTwkDial_FrElas_pi,     0.10, 0.10);
  this->SetUncertainty( kINukeTwkDial_FrInel_pi,     0.40, 0.40);
  this->SetUncertainty( kINukeTwkDial_FrAbs_pi,      0.30, 0.30);
  this->SetUncertainty( kINukeTwkDial_FrPiProd_pi,   0.20, 0.20);
  this->SetUncertainty( kINukeTwkDial_FrCEx_N,       0.50, 0.50);
  this->SetUncertainty( kINukeTwkDial_FrElas_N,      0.30, 0.30);
  this->SetUncertainty( kINukeTwkDial_FrInel_N,      0.40, 0.40);
  this->SetUncertainty( kINukeTwkDial_FrAbs_N,       0.20, 0.20);
  this->SetUncertainty( kINukeTwkDial_FrPiProd_N,    0.20, 0.20);

  this->SetUncertainty( kRDcyTwkDial_BR1gamma,       0.50, 0.50);
  this->SetUncertainty( kRDcyTwkDial_BR1eta,         0.50, 0.50);
}
//____________________________________________________________________________
