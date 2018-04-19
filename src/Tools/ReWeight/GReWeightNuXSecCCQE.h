//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCCQE

\brief    Reweighting CCQE GENIE neutrino cross sections

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_CCQE_H_
#define _G_REWEIGHT_NU_XSEC_CCQE_H_

//#define _G_REWEIGHT_CCQE_DEBUG_

#include <string>

#include "Tools/ReWeight/GReWeightModel.h"

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew   {

 class GReWeightNuXSecCCQE : public GReWeightModel 
 {
 public:
   static const int kModeMa               = 0;
   static const int kModeNormAndMaShape   = 1;
   static const int kModeZExp             = 2;

   static const int fZExpMaxSyst          = 4; ///< maximum number of systematics

   GReWeightNuXSecCCQE();
   GReWeightNuXSecCCQE(std::string model, std::string type);
  ~GReWeightNuXSecCCQE();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // various config options
   void SetMode     (int mode) { fMode       = mode; }
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }
   void SetMaPath   (string p) { fMaPath     = p;    }
   // z-expansion specific options
   void SetZExpPath    (string p){ fZExpPath    = p;   }

 private:

   void   Init              (void);
   double CalcWeightNorm    (const EventRecord & event);
   double CalcWeightMaShape (const EventRecord & event);
   double CalcWeightMa      (const EventRecord & event);
   double CalcWeightZExp    (const EventRecord & event);

   XSecAlgorithmI * fXSecModelDef;    ///< default model
   XSecAlgorithmI * fXSecModel;       ///< tweaked model
   Registry *       fXSecModelConfig; ///< config in tweaked model
   string fFFModel; ///< String name of form factor model
   bool fModelIsDipole; ///< Using dipole form factors?
   bool fModelIsZExp;   ///< Using Zexp form factors?
   std::string fManualModelName; ///< If using a tweaked model that isn't the same as default, name
   std::string fManualModelType; ///< If using a tweaked model that isn't the same as default, type

   int    fMode;         ///< 0: Ma, 1: Norm and MaShape, 2: Z-Expansion
   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?
   string fMaPath;       ///< M_{A} path in config Registry
   double fNormTwkDial;  ///<
   double fNormDef;      ///<
   double fNormCurr;     ///<
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<

   // unused // int     fZExpCurrIdx; ///< current coefficient index
   int     fZExpMaxCoef; ///< max number of coefficients to use
   string  fZExpPath;    ///< algorithm path to get coefficients
   double  fZExpTwkDial[fZExpMaxSyst]; ///< 
   double  fZExpDef    [fZExpMaxSyst]; ///<
   double  fZExpCurr   [fZExpMaxSyst]; ///< array of current parameter values

#ifdef _G_REWEIGHT_CCQE_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif
 };

} // rew   namespace
} // genie namespace

#endif

