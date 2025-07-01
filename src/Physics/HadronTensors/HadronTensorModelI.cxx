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

// standard library includes
#include <string>

// GENIE includes
#include "Physics/HadronTensors/HadronTensorModelI.h"


//____________________________________________________________________________
genie::HadronTensorModelI::HadronTensorModelI()
  : genie::Algorithm()
{

}

//____________________________________________________________________________
genie::HadronTensorModelI::HadronTensorModelI(std::string name)
  : genie::Algorithm( name )
{

}

//____________________________________________________________________________
genie::HadronTensorModelI::HadronTensorModelI(std::string name,
  std::string config) : genie::Algorithm(name, config)
{

}

//____________________________________________________________________________
genie::HadronTensorModelI::~HadronTensorModelI()
{

}
