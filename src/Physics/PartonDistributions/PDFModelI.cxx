//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Physics/PartonDistributions/PDFModelI.h"

using namespace genie;

//____________________________________________________________________________
PDFModelI::PDFModelI() :
Algorithm()
{

}
//____________________________________________________________________________
PDFModelI::PDFModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
PDFModelI::PDFModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
PDFModelI::~PDFModelI()
{

}
//____________________________________________________________________________
