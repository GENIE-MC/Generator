//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightModel

\brief    GENIE event reweighting engine common parameters

\author   Steve Dennis <s.r.dennis \at liverpool.ac.uk>
          University of Liverpool

\created  April 2018

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_MODEL_BASE_H_
#define _G_REWEIGHT_MODEL_BASE_H_

#include "Tools/ReWeight/GReWeightI.h"

namespace genie {

class EventRecord;

namespace rew   {

 class GReWeightModel : public GReWeightI
 {
 public:
  GReWeightModel(std::string name);
  ~GReWeightModel();

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
  
  //! Should we calculate the old weight ourselves, or use the one from the input tree? Default on.
  virtual void UseOldWeightFromFile(bool);
  
  //! If using the weight from the file, how many times should we check by calculating it ourself?
  virtual void SetNWeightChecks(int);

 protected:
   bool fUseOldWeightFromFile;
   int  fNWeightChecksToDo;
   int  fNWeightChecksDone;
   bool fFailedWeightCheck;
   
   std::string fName;
 };

} // rew   namespace
} // genie namespace

#endif

