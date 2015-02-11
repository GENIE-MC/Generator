//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCCQEvec

\brief    Reweighting vector form factors in GENIE CCQE neutrino cross
          section calculations.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  May 24, 2010

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_VEC_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_VEC_H_

#include <map>
#include <string>

#include "ReWeight/GReWeightI.h"

using std::map;
using std::string;

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;

namespace rew   {

 class GReWeightNuXSecCCQEvec : public GReWeightI 
 {
 public:
   GReWeightNuXSecCCQEvec();
  ~GReWeightNuXSecCCQEvec();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

   // various config options
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }

 private:

   void Init (void);

   XSecAlgorithmI * fXSecModel_bba;  ///< CCQE model with BBA05  f/f (default)
   XSecAlgorithmI * fXSecModel_dpl;  ///< CCQE model with dipole f/f ("maximally" tweaked)

   double fFFTwkDial;    ///< tweaking dial (0: bba/default, +1: dipole)

   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?

   TFile *    fTestFile;
   TNtupleD * fTestNtp;
 };

} // rew   namespace
} // genie namespace

#endif

