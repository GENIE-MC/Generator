//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecNCRES

\brief    Reweight GENIE NC resonance neutrino-production cross section.
          Basically a clone of the corresponding CC code.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_NCRES_H_
#define _G_REWEIGHT_NU_XSEC_NCRES_H_

#include <map>
#include <string>

using std::map;
using std::string;

#include "ReWeight/GReWeightI.h"

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew   {

 class GReWeightNuXSecNCRES : public GReWeightI 
 {
 public:
   static const int kModeMaMv             = 0;
   static const int kModeNormAndMaMvShape = 1;

   GReWeightNuXSecNCRES();
  ~GReWeightNuXSecNCRES();

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
   bool   fRewNue;       ///< reweight nu_e NC?
   bool   fRewNuebar;    ///< reweight nu_e_bar NC?
   bool   fRewNumu;      ///< reweight nu_mu NC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar NC?
   double fNormTwkDial;  ///<
   double fNormDef;      ///<
   double fNormCurr;     ///<
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<
   double fMvTwkDial;    ///<
   double fMvDef;        ///<
   double fMvCurr;       ///<

 };

} // rew   namespace
} // genie namespace

#endif

