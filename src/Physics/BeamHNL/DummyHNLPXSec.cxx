//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/BeamHNL/DummyHNLPXSec.h"

using namespace genie;

//____________________________________________________________________________
DummyHNLPXSec::DummyHNLPXSec() :
XSecAlgorithmI("genie::DummyHNLPXSec")
{

}
//____________________________________________________________________________
DummyHNLPXSec::DummyHNLPXSec(string config) :
XSecAlgorithmI("genie::DummyHNLPXSec", config)
{

}
//____________________________________________________________________________
DummyHNLPXSec::~DummyHNLPXSec()
{

}
//____________________________________________________________________________
double DummyHNLPXSec::XSec(const Interaction * , KinePhaseSpace_t ) const
{
  return 0;
}
//____________________________________________________________________________
double DummyHNLPXSec::Integral(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
bool DummyHNLPXSec::ValidProcess(const Interaction * ) const
{
  return true;
}
//____________________________________________________________________________
