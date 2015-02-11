//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightAGKY

\brief    Reweighting the GENIE AGKY (free-nucleon) hadronization model

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_AGKY_H_
#define _G_REWEIGHT_AGKY_H_

#include "ReWeight/GReWeightI.h"

using namespace genie::rew;
using namespace genie;

class TF1;
class TFile;
class TNtupleD;

namespace genie {
namespace rew   {

 class GReWeightAGKY : public GReWeightI 
 {
 public:
   GReWeightAGKY();
  ~GReWeightAGKY();

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

 private:

   void   Init       (void);
   double RewxFpT1pi (const EventRecord & event);

   bool   fRewNue;              ///< reweight nu_e?
   bool   fRewNuebar;           ///< reweight nu_e_bar?
   bool   fRewNumu;             ///< reweight nu_mu?
   bool   fRewNumubar;          ///< reweight nu_mu_bar?
   bool   fRewCC;               ///< reweight CC?
   bool   fRewNC;               ///< reweight NC?
   double fXFmin;               ///<
   double fXFmax;               ///<
   double fPT2min;              ///<
   double fPT2max;              ///<
   TF1 *  fBaryonXFpdf;         ///<
   TF1 *  fBaryonPT2pdf;        ///<
   TF1 *  fBaryonXFpdfTwk;      ///<
   TF1 *  fBaryonPT2pdfTwk;     ///<
   double fDefPeakBaryonXF;     ///<
   double fDefAvgPT2;           ///<
   double fPeakBaryonXFTwkDial; ///<
   double fAvgPT2TwkDial;       ///<
   double fI0XFpdf;             ///<
   double fI0PT2pdf;            ///<

   TFile *    fTestFile;
   TNtupleD * fTestNtp;
 };

} // rew
} // genie

#endif

