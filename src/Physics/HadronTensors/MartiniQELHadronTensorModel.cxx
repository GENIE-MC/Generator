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
#include "Physics/HadronTensors/MartiniQELHadronTensorModel.h"
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"
#include "Physics/HadronTensors/TabulatedLabFrameHadronTensor.h"

//____________________________________________________________________________
genie::MartiniQELHadronTensorModel::MartiniQELHadronTensorModel()
  : genie::TabulatedHadronTensorModelI("genie::MartiniQELHadronTensorModel")
{

}

//____________________________________________________________________________
genie::MartiniQELHadronTensorModel::MartiniQELHadronTensorModel(std::string config)
  : genie::TabulatedHadronTensorModelI("genie::MartiniQELHadronTensorModel", config)
{

}

//____________________________________________________________________________
genie::MartiniQELHadronTensorModel::~MartiniQELHadronTensorModel()
{

}

//____________________________________________________________________________
genie::HadronTensorI* genie::MartiniQELHadronTensorModel::ParseTensorFile(
  const std::string& full_file_name) const
{
  return new TabulatedLabFrameHadronTensor( full_file_name );
}
