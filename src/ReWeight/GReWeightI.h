//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightI

\brief    GENIE event reweighting engine ABC 

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_ABC_H_
#define _G_REWEIGHT_ABC_H_

#include "ReWeight/GSyst.h"
#include "ReWeight/GSystSet.h"

namespace genie {

class EventRecord;

namespace rew   {

 class GReWeightI 
 {
 public:
  virtual ~GReWeightI() 
  { 

  } 

  //
  // define the GReWeightI interface
  //

  //! does the current weight calculator handle the input nuisance param?
  virtual bool IsHandled (GSyst_t syst) = 0; 

  //! update the value for the specified nuisance param
  virtual void SetSystematic (GSyst_t syst, double val) = 0; 

  //!  set all nuisance parameters to default values
  virtual void Reset (void) = 0;            

  //! propagate updated nuisance parameter values to actual MC, etc
  virtual void Reconfigure (void) = 0;            
  
  //! calculate a weight for the input event using the current nuisance param values
  virtual double CalcWeight (const genie::EventRecord & event) = 0;  

  //! calculate penalty factors
  virtual double CalcChisq (void) = 0;        

 protected:

   GReWeightI() 
   { 
   
   }
 };

} // rew   namespace
} // genie namespace

#endif

