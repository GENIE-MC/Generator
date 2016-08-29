//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
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

#include "NeutronOsc/NOscDummyPXSec.h"

using namespace genie;

//____________________________________________________________________________
NOscDummyPXSec::NOscDummyPXSec() :
XSecAlgorithmI("genie::NoscDummyPXSec")
{

}
//____________________________________________________________________________
NOscDummyPXSec::NOscDummyPXSec(string config) :
XSecAlgorithmI("genie::NOscDummyPXSec", config)
{

}
//____________________________________________________________________________
NOscDummyPXSec::~NOscDummyPXSec()
{

}
//____________________________________________________________________________
double NOscDummyPXSec::XSec(const Interaction * , KinePhaseSpace_t ) const
{
  return 0;
}
//____________________________________________________________________________
double NOscDummyPXSec::Integral(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
bool NOscDummyPXSec::ValidProcess(const Interaction * ) const
{
  return true;
}
//____________________________________________________________________________
