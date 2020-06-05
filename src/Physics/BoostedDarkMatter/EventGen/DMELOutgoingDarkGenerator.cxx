//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Author: Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison

         Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 
*/
//____________________________________________________________________________

#include "Framework/GHEP/GHepRecord.h"
#include "Physics/BoostedDarkMatter/EventGen/DMELOutgoingDarkGenerator.h"

using namespace genie;

//___________________________________________________________________________
DMELOutgoingDarkGenerator::DMELOutgoingDarkGenerator() :
OutgoingDarkGenerator("genie::DMELOutgoingDarkGenerator")
{

}
//___________________________________________________________________________
DMELOutgoingDarkGenerator::DMELOutgoingDarkGenerator(string config):
OutgoingDarkGenerator("genie::DMELOutgoingDarkGenerator", config)
{

}
//___________________________________________________________________________
DMELOutgoingDarkGenerator::~DMELOutgoingDarkGenerator()
{

}
//___________________________________________________________________________
void DMELOutgoingDarkGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  // We need to do a modified version of the standard procedure in this case
  // Some of the pieces are identical, but the kinematics are very different
  OutgoingDarkGenerator::ProcessEventRecord(evrec);  
}
//___________________________________________________________________________
