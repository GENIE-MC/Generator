//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Nuclear/NuclearModelI.h"

using namespace genie;

//____________________________________________________________________________
NuclearModelI::NuclearModelI() :
Algorithm()
{

}
//____________________________________________________________________________
NuclearModelI::NuclearModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
NuclearModelI::NuclearModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
NuclearModelI::~NuclearModelI()
{

}
//____________________________________________________________________________
double NuclearModelI::RemovalEnergy(void) const
{
  return fCurrRemovalEnergy;
}
//____________________________________________________________________________
double NuclearModelI::Momentum(void) const
{
  return fCurrMomentum.Mag();
}
//____________________________________________________________________________
TVector3 NuclearModelI::Momentum3(void) const
{
  return fCurrMomentum;
}
//____________________________________________________________________________

