//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecNCEL

\brief    Reweighting NCEL GENIE neutrino cross sections

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Nov 25, 2010

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_NU_XSEC_NCEL_H_
#define _G_REWEIGHT_NU_XSEC_NCEL_H_

//#define _G_REWEIGHT_NCEL_DEBUG_

#include <string>

#include "Tools/ReWeight/GReWeightModel.h"

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew   {

 class GReWeightNuXSecNCEL : public GReWeightModel 
 {
 public:
   GReWeightNuXSecNCEL();
   GReWeightNuXSecNCEL(std::string model, std::string type);
  ~GReWeightNuXSecNCEL();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // various config options
   void RewNue      (bool tf ) { fRewNue     = tf;   }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf;   }
   void RewNumu     (bool tf ) { fRewNumu    = tf;   }
   void RewNumubar  (bool tf ) { fRewNumubar = tf;   }
   void SetMaPath   (string p) { fMaPath     = p;    }
   void SetEtaPath  (string p) { fEtaPath    = p;    }

 private:

   void Init(void);

   XSecAlgorithmI * fXSecModelDef;    ///< default model
   XSecAlgorithmI * fXSecModel;       ///< tweaked model
   Registry *       fXSecModelConfig; ///< config in tweaked model

   bool   fRewNue;       ///< reweight nu_e CC?
   bool   fRewNuebar;    ///< reweight nu_e_bar CC?
   bool   fRewNumu;      ///< reweight nu_mu CC?
   bool   fRewNumubar;   ///< reweight nu_mu_bar CC?
   string fMaPath;       ///< M_{A} path in config Registry
   string fEtaPath;      ///< eta path in config Registry
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<
   double fEtaTwkDial;    ///<
   double fEtaDef;        ///<
   double fEtaCurr;       ///<
   
   std::string fManualModelName; ///< If using a tweaked model that isn't the same as default, name
   std::string fManualModelType; ///< If using a tweaked model that isn't the same as default, type

#ifdef _G_REWEIGHT_NCEL_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif
 };

} // rew   namespace
} // genie namespace

#endif

