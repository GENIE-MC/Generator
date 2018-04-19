//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________


#include "Framework/Messenger/Messenger.h"
#include "Tools/ReWeight/GReWeightModel.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightModel::GReWeightModel(std::string name) :
GReWeightI(),
fUseOldWeightFromFile(true),
fNWeightChecksToDo(20),
fNWeightChecksDone(0),
fFailedWeightCheck(false),
fName(name)
{
  
}
//_______________________________________________________________________________________
GReWeightModel::~GReWeightModel()
{
  if (fFailedWeightCheck && fUseOldWeightFromFile) {
    LOG("ReW",pWARN) << fName<< ": You used the weights from the files but the"
      <<" check against the calculated weights failed. Your weights are probably wrong!";
  }
}
//_______________________________________________________________________________________
void GReWeightModel::SetNWeightChecks(int n)
{
  fNWeightChecksToDo = n;
}
//_______________________________________________________________________________________
void GReWeightModel::UseOldWeightFromFile(bool should_we)
{
  fUseOldWeightFromFile = should_we;
}
//_______________________________________________________________________________________
