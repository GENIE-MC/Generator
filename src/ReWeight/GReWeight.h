//____________________________________________________________________________
/*!

\class    genie::rew::GReWeight

\brief    Interface to the GENIE event reweighting engines

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_H_
#define _G_REWEIGHT_H_

#include "ReWeight/GSystSet.h"
#include "ReWeight/GReWeightI.h"

namespace genie {

class EventRecord;

namespace rew   {

 class GReWeight
 {
 public:
   GReWeight();
  ~GReWeight();

   GSystSet & Systematics (void);                             ///<
   void       Reconfigure (void);                             ///<
   double     CalcWeight  (const genie::EventRecord & event); ///<
   double     CalcChisq   (void);                             ///<
   void       Print       (void);                             ///<       

  private:

   void Init    (void);
   void CleanUp (void);

   GSystSet     fSystSet;         ///< set of enabled nuisance parameters
   GReWeightI * fReWeightNuXSec;  ///< handles all enabled GENIE cross section params
   GReWeightI * fReWeightAGKY;    ///< handles all enabled GENIE hadronization params
   GReWeightI * fReWeightINuke;   ///< handles all enabled GENIE intranuclear rescattering params
 };

} // rew   namespace
} // genie namespace

#endif

