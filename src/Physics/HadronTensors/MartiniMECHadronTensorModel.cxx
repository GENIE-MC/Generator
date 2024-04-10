//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

// GENIE includes
#include "Physics/HadronTensors/MartiniMECHadronTensorModel.h"
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"
#include "Physics/HadronTensors/TabulatedLabFrameHadronTensor.h"

//____________________________________________________________________________
genie::MartiniMECHadronTensorModel::MartiniMECHadronTensorModel()
  : genie::TabulatedHadronTensorModelI("genie::MartiniMECHadronTensorModel")
{

}

//____________________________________________________________________________
genie::MartiniMECHadronTensorModel::MartiniMECHadronTensorModel(std::string config)
  : genie::TabulatedHadronTensorModelI("genie::MartiniMECHadronTensorModel", config)
{

}

//____________________________________________________________________________
genie::MartiniMECHadronTensorModel::~MartiniMECHadronTensorModel()
{

}

//____________________________________________________________________________
genie::HadronTensorI* genie::MartiniMECHadronTensorModel::ParseTensorFile(
  const std::string& full_file_name) const
{
  return new TabulatedLabFrameHadronTensor( full_file_name );
}
