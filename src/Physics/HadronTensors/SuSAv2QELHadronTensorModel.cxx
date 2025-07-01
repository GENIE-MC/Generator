//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

// GENIE includes
#include "Physics/HadronTensors/SuSAv2QELHadronTensorModel.h"
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"
#include "Physics/HadronTensors/TabulatedLabFrameHadronTensor.h"

//____________________________________________________________________________
genie::SuSAv2QELHadronTensorModel::SuSAv2QELHadronTensorModel()
  : genie::TabulatedHadronTensorModelI("genie::SuSAv2QELHadronTensorModel")
{

}

//____________________________________________________________________________
genie::SuSAv2QELHadronTensorModel::SuSAv2QELHadronTensorModel(std::string config)
  : genie::TabulatedHadronTensorModelI("genie::SuSAv2QELHadronTensorModel", config)
{

}

//____________________________________________________________________________
genie::SuSAv2QELHadronTensorModel::~SuSAv2QELHadronTensorModel()
{

}

//____________________________________________________________________________
genie::HadronTensorI* genie::SuSAv2QELHadronTensorModel::ParseTensorFile(
  const std::string& full_file_name) const
{
  return new TabulatedLabFrameHadronTensor( full_file_name );
}
