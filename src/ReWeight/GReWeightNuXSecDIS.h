//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightNuXSecDIS

\brief    Reweighting GENIE DIS neutrino-nucleus cross sections

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

#ifndef _G_REWEIGHT_NU_XSEC_DIS_H_
#define _G_REWEIGHT_NU_XSEC_DIS_H_

#include "ReWeight/GReWeightI.h"

class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;
class Registry;

namespace rew   {

 class GReWeightNuXSecDIS : public GReWeightI 
 {
 public:
   static const int kModeABCV12u      = 0;
   static const int kModeABCV12uShape = 1;

   GReWeightNuXSecDIS();
  ~GReWeightNuXSecDIS();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

   // various config options
   void SetMode      (int    m )  { fMode       = m;  }
   void RewNue       (bool   tf)  { fRewNue     = tf; }
   void RewNuebar    (bool   tf)  { fRewNuebar  = tf; }
   void RewNumu      (bool   tf)  { fRewNumu    = tf; }
   void RewNumubar   (bool   tf)  { fRewNumubar = tf; }
   void RewCC        (bool   tf)  { fRewCC      = tf; }
   void RewNC        (bool   tf)  { fRewNC      = tf; }
   void SetAhtBYPath (string p )  { fAhtBYPath  = p;  }
   void SetBhtBYPath (string p )  { fBhtBYPath  = p;  }
   void SetCV1uBYPath(string p )  { fCV1uBYPath = p;  }
   void SetCV2uBYPath(string p )  { fCV2uBYPath = p;  }
   void SetWminCut   (double W )  { fWmin       = W;  }
   void SetQ2minCut  (double Q2)  { fQ2min      = Q2; }

 private:

   void   Init                   (void);
   double CalcWeightABCV12u      (const genie::EventRecord & event); ///< rew. Aht,Bht,CV1u,CV2u
   double CalcWeightABCV12uShape (const genie::EventRecord & event); ///< rew. AhtShape,BhtShape,CV1uShape,CV2uShape

   XSecAlgorithmI * fXSecModelDef;    ///< default model
   XSecAlgorithmI * fXSecModel;       ///< tweaked model
   Registry *       fXSecModelConfig; ///< config in tweaked model

   bool   fRewNue;               ///< reweight nu_e?
   bool   fRewNuebar;            ///< reweight nu_e_bar?
   bool   fRewNumu;              ///< reweight nu_mu?
   bool   fRewNumubar;           ///< reweight nu_mu_bar?
   bool   fRewCC;                ///< reweight CC?
   bool   fRewNC;                ///< reweight NC?

   int    fMode;            ///< 0: Aht,Bht,CV1u,CV2u, 1:AhtShape,BhtShape,CV1uShape,CV2uShape
   double fWmin;            ///< W_{min}  cut. Reweight only events with W  > W_{min}.
   double fQ2min;           ///< Q2_{min} cut. Reweight only events with Q2 > Q2_{min}.
   double fAhtBYTwkDial;    ///< tweak dial for BY parameter: Aht
   double fBhtBYTwkDial;    ///< tweak dial for BY parameter: Bht
   double fCV1uBYTwkDial;   ///< tweak dial for BY parameter: CV1u
   double fCV2uBYTwkDial;   ///< tweak dial for BY parameter: CV2u
   double fAhtBYDef;        ///<
   double fBhtBYDef;        ///<
   double fCV1uBYDef;       ///<
   double fCV2uBYDef;       ///<
   double fAhtBYCur;        ///<
   double fBhtBYCur;        ///<
   double fCV1uBYCur;       ///<
   double fCV2uBYCur;       ///<
   string fAhtBYPath;       ///<
   string fBhtBYPath;       ///<
   string fCV1uBYPath;      ///<
   string fCV2uBYPath;      ///<

   TFile *    fTestFile;
   TNtupleD * fTestNtp;
 };

} // rew   namespace
} // genie namespace

#endif

