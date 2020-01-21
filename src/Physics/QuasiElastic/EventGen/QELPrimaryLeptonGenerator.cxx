//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 
*/
//____________________________________________________________________________

#include "Framework/GHEP/GHepRecord.h"
#include "Physics/QuasiElastic/EventGen/QELPrimaryLeptonGenerator.h"

using namespace genie;

//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::QELPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::QELPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
QELPrimaryLeptonGenerator::~QELPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void QELPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in QEL events

  // no modification is required to the std implementation
  PrimaryLeptonGenerator::ProcessEventRecord(evrec);
}
//___________________________________________________________________________
