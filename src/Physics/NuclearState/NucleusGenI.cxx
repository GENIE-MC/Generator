//____________________________________________________________________________
/*!

  \class    genie::NucleusGenerator

  \brief    It visits the event record & computes a Fermi motion momentum for
  initial state nucleons bound in nuclei.
  Is a concrete implementation of the EventRecordVisitorI interface.

  \author   Liang Liu <liangliu \at fnal.gov>
  Fermi National Accelerator Laboratory

  \created  October 17, 2024

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#include "Physics/NuclearState/NucleusGenI.h"

namespace genie {

  NucleusGenI::NucleusGenI(string name ) :
    EventRecordVisitorI( name )
  {

  }
  //___________________________________________________________________________
  NucleusGenI::NucleusGenI(string name, string config) :
    EventRecordVisitorI( name, config)
  {

  }
  //___________________________________________________________________________
  NucleusGenI::~NucleusGenI()
  {

  }
  //____________________________________________________________________________
  void NucleusGenI::LoadConfig(void)
  {

  }
  //____________________________________________________________________________
}

