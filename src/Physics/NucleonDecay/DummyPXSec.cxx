//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/NucleonDecay/DummyPXSec.h"

using namespace genie;

//____________________________________________________________________________
DummyPXSec::DummyPXSec() :
XSecAlgorithmI("genie::DummyPXSec")
{

}
//____________________________________________________________________________
DummyPXSec::DummyPXSec(string config) :
XSecAlgorithmI("genie::DummyPXSec", config)
{

}
//____________________________________________________________________________
DummyPXSec::~DummyPXSec()
{

}
//____________________________________________________________________________
double DummyPXSec::XSec(const Interaction * , KinePhaseSpace_t ) const
{
  return 0;
}
//____________________________________________________________________________
double DummyPXSec::Integral(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
bool DummyPXSec::ValidProcess(const Interaction * ) const
{
  return true;
}
//____________________________________________________________________________
