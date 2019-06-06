//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 02, 2009 - CA
   Extended the decayer interface with the `UnInhibitDecay(int,TDecayChannel*) 
   const' and `InhibitDecay(int,TDecayChannel*) const' pure virtual methods.

*/
//____________________________________________________________________________

#include "Physics/Decay/DecayModelI.h"

using namespace genie;

//____________________________________________________________________________
DecayModelI::DecayModelI() :
Algorithm()
{

}
//____________________________________________________________________________
DecayModelI::DecayModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
DecayModelI::DecayModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
DecayModelI::~DecayModelI() 
{ 

}
//____________________________________________________________________________



