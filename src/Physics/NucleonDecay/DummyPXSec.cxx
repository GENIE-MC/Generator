//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 10, 2011 - CA
   First added in v2.7.1.

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
