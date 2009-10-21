//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightAGKY

\brief    Reweighting the GENIE AGKY (free-nucleon) hadronization model

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep 10, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_AGKY_H_
#define _G_REWEIGHT_AGKY_H_

#include "ReWeight/GReWeightI.h"

using namespace genie::rew;
using namespace genie;

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
 };

} // rew
} // genie

#endif

