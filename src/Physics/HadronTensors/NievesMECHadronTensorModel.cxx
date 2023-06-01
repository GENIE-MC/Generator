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
#include "Physics/HadronTensors/NievesMECHadronTensorModel.h"
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"
#include "Physics/HadronTensors/TabulatedLabFrameHadronTensor.h"

//____________________________________________________________________________
genie::NievesMECHadronTensorModel::NievesMECHadronTensorModel()
  : genie::TabulatedHadronTensorModelI("genie::NievesMECHadronTensorModel")
{

}

//____________________________________________________________________________
genie::NievesMECHadronTensorModel::NievesMECHadronTensorModel(std::string config)
  : genie::TabulatedHadronTensorModelI("genie::NievesMECHadronTensorModel", config)
{

}

//____________________________________________________________________________
genie::NievesMECHadronTensorModel::~NievesMECHadronTensorModel()
{

}

//____________________________________________________________________________
genie::HadronTensorI* genie::NievesMECHadronTensorModel::ParseTensorFile(
  const std::string& full_file_name) const
{
  return new TabulatedLabFrameHadronTensor( full_file_name );
}
