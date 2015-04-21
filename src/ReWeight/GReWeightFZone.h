//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightFZone

\brief    Reweighting the formation zone model

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Sep 20, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_FZONE_H_
#define _G_REWEIGHT_FZONE_H_

#include "ReWeight/GReWeightI.h"

using namespace genie::rew;
using namespace genie;

namespace genie {
namespace rew   {

 class GReWeightFZone : public GReWeightI 
 {
 public:
   GReWeightFZone();
  ~GReWeightFZone();

   // implement the GReWeightI interface
   bool   IsHandled      (GSyst_t syst);
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);
   double CalcChisq      (void);

   // other config options
   // set to match values used at event generation
   void SetR0 (double R0 )        { fR0         = R0;  }
   void SetNR (double NR )        { fNR         = NR;  }
   void SetCT0Pion(double ct0)    { fct0pion    = ct0; }
   void SetCT0Nucleon(double ct0) { fct0nucleon = ct0; }
   void SetK  (double k  )        { fK          = k;   }

 private:

   void Init(void);

   double fFZoneTwkDial; ///< formation zone tweaking dial

   double fNR;         ///<
   double fR0;         ///<
   double fct0pion;    ///<
   double fct0nucleon; ///<
   double fK;          ///<

 };

} // rew
} // genie

#endif

