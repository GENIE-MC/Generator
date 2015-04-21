//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecCOH

\brief    Reweighting GENIE coherent neutrino-nucleus cross sections

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

#ifndef _G_REWEIGHT_NU_XSEC_COH_H_
#define _G_REWEIGHT_NU_XSEC_COH_H_

#include <map>
#include <string>

#include "ReWeight/GReWeightI.h"

using std::map;
using std::string;

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew   {

 class GReWeightNuXSecCOH : public GReWeightI 
 {
 public:
   GReWeightNuXSecCOH();
  ~GReWeightNuXSecCOH();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

   // various config options
   void RewNue      (bool tf ) { fRewNue     = tf; }
   void RewNuebar   (bool tf ) { fRewNuebar  = tf; }
   void RewNumu     (bool tf ) { fRewNumu    = tf; }
   void RewNumubar  (bool tf ) { fRewNumubar = tf; }
   void RewCC       (bool tf ) { fRewCC      = tf; }
   void RewNC       (bool tf ) { fRewNC      = tf; }
   void SetMaPath   (string p) { fMaPath     = p;  }
   void SetR0Path   (string p) { fR0Path     = p;  }

 private:

   void Init (void);

   XSecAlgorithmI * fXSecModel;       ///<
   Registry *       fXSecModelConfig; ///<

   bool   fRewNue;       ///< reweight nu_e?
   bool   fRewNuebar;    ///< reweight nu_e_bar?
   bool   fRewNumu;      ///< reweight nu_mu?
   bool   fRewNumubar;   ///< reweight nu_mu_bar?
   bool   fRewCC;        ///< reweight CC?
   bool   fRewNC;        ///< reweight NC?
   string fMaPath;       ///< M_{A} path in config Registry
   string fR0Path;       ///< R_{0} path in config Registry
   double fMaTwkDial;    ///<
   double fMaDef;        ///<
   double fMaCurr;       ///<
   double fR0TwkDial;    ///<
   double fR0Def;        ///<
   double fR0Curr;       ///<

 };

} // rew   namespace
} // genie namespace

#endif

