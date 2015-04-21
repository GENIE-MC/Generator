//____________________________________________________________________________
/*!

\class    genie::rew::GReWeight

\brief    Interface to the GENIE event reweighting engines

\author   Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_H_
#define _G_REWEIGHT_H_

#include <string>
#include <map>

#include "ReWeight/GSystSet.h"
#include "ReWeight/GReWeightI.h"

using std::string;
using std::map;

namespace genie {

class EventRecord;

namespace rew   {

 class GReWeight
 {
 public:
   GReWeight();
  ~GReWeight();

   void        AdoptWghtCalc (string name, GReWeightI* wcalc);   ///< add concrete weight calculator, transfers ownership
   GReWeightI* WghtCalc      (string name);                      ///< access a weight calculator by name
   GSystSet &  Systematics   (void);                             ///< set of enabled systematic params & values
   void        Reconfigure   (void);                             ///< reconfigure weight calculators with new params
   double      CalcWeight    (const genie::EventRecord & event); ///< calculate weight for input event
   double      CalcChisq     (void);                             ///< calculate penalty chisq for current values of tweaking dials
   void        Print         (void);                             ///< print

  private:

   void CleanUp (void);

   GSystSet                  fSystSet;   ///< set of enabled nuisance parameters
   map<string, GReWeightI *> fWghtCalc;  ///< concrete weight calculators
 };

} // rew   namespace
} // genie namespace

#endif

