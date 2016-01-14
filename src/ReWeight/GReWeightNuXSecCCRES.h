//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCCRES

\brief    Reweight GENIE CC resonance neutrino-production

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCRES_H_
#define _G_REWEIGHT_NU_XSEC_CCRES_H_

#include <map>
#include <string>

using std::map;
using std::string;

#include "ReWeight/GReWeightI.h"

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew   {

 class GReWeightNuXSecCCRES : public GReWeightI 
 {
 public:
   static const int kModeMaMv             = 0;
   static const int kModeNormAndMaMvShape = 1;

   GReWeightNuXSecCCRES();
  ~GReWeightNuXSecCCRES();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

   // various config options
   void SetMode     (int mode) { fMode       = mode; }
   void SetMaPath   (string p) { fMaPath     = p;    }
   void SetMvPath   (string p) { fMvPath     = p;    }
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }

 private:

   void   Init                (void);
   double CalcWeightNorm      (const EventRecord & event);
   double CalcWeightMaMvShape (const EventRecord & event);
   double CalcWeightMaMv      (const EventRecord & event);

   XSecAlgorithmI * fXSecModelDef;    ///< default model
   XSecAlgorithmI * fXSecModel;       ///< tweaked model
   Registry *       fXSecModelConfig; ///< config in tweaked model

   int    fMode;         ///< 0: Ma/Mv, 1: Norm and MaShape/MvShape
   string fMaPath;       ///< M_{A} path in configuration
   string fMvPath;       ///< M_{V} path in configuration
   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?
   double fNormTwkDial;  ///<
   double fNormDef;      ///<
   double fNormCurr;     ///<
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<
   double fMvTwkDial;    ///<
   double fMvDef;        ///<
   double fMvCurr;       ///<

   TFile *    fTestFile;
   TNtupleD * fTestNtp;
 };

} // rew   namespace
} // genie namespace

#endif

