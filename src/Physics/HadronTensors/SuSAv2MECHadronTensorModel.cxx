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
#include "Physics/HadronTensors/SuSAv2MECHadronTensorModel.h"
#include "Physics/HadronTensors/TabulatedHadronTensorModelI.h"
#include "Physics/HadronTensors/TabulatedLabFrameHadronTensor.h"

//____________________________________________________________________________
genie::SuSAv2MECHadronTensorModel::SuSAv2MECHadronTensorModel()
  : genie::TabulatedHadronTensorModelI("genie::SuSAv2MECHadronTensorModel")
{

}

//____________________________________________________________________________
genie::SuSAv2MECHadronTensorModel::SuSAv2MECHadronTensorModel(std::string config)
  : genie::TabulatedHadronTensorModelI("genie::SuSAv2MECHadronTensorModel", config)
{

}

//____________________________________________________________________________
genie::SuSAv2MECHadronTensorModel::~SuSAv2MECHadronTensorModel()
{

}

//____________________________________________________________________________
genie::HadronTensorI* genie::SuSAv2MECHadronTensorModel::ParseTensorFile(
  const std::string& full_file_name) const
{
  return new TabulatedLabFrameHadronTensor( full_file_name );
}
