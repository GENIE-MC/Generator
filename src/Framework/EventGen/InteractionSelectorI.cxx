//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - December 05, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionSelectorI.h"
#include "Framework/Interaction/Interaction.h"

using namespace genie;

//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI() :
Algorithm()
{

}
//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
InteractionSelectorI::~InteractionSelectorI()
{

}
//___________________________________________________________________________
