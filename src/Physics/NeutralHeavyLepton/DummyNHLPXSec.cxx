//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Physics/NeutralHeavyLepton/DummyNHLPXSec.h"

using namespace genie;

//____________________________________________________________________________
DummyNHLPXSec::DummyNHLPXSec() :
XSecAlgorithmI("genie::DummyNHLPXSec")
{

}
//____________________________________________________________________________
DummyNHLPXSec::DummyNHLPXSec(string config) :
XSecAlgorithmI("genie::DummyNHLPXSec", config)
{

}
//____________________________________________________________________________
DummyNHLPXSec::~DummyNHLPXSec()
{

}
//____________________________________________________________________________
double DummyNHLPXSec::XSec(const Interaction * , KinePhaseSpace_t ) const
{
  return 0;
}
//____________________________________________________________________________
double DummyNHLPXSec::Integral(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
bool DummyNHLPXSec::ValidProcess(const Interaction * ) const
{
  return true;
}
//____________________________________________________________________________
