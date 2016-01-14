//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TClonesArray.h>
#include <TH1D.h>

#include "Fragmentation/HadronizationModelI.h"
#include "Interaction/Interaction.h"
#include "PDG/PDGCodeList.h"

using namespace genie;

//____________________________________________________________________________
HadronizationModelI::HadronizationModelI() :
Algorithm()
{

}
//____________________________________________________________________________
HadronizationModelI::HadronizationModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
HadronizationModelI::HadronizationModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
HadronizationModelI::~HadronizationModelI() 
{ 

}
//____________________________________________________________________________



