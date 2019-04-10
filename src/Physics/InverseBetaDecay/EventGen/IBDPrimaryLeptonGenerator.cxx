//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Corey Reed <cjreed \at nikhef.nl> - October 29, 2009
         using code from the QELKinematicGenerator written by
         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 03, 2004

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/GHEP/GHepRecord.h"
#include "Physics/InverseBetaDecay/EventGen/IBDPrimaryLeptonGenerator.h"

using namespace genie;

//___________________________________________________________________________
IBDPrimaryLeptonGenerator::IBDPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::IBDPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
IBDPrimaryLeptonGenerator::IBDPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::IBDPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
IBDPrimaryLeptonGenerator::~IBDPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void IBDPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in IBD events

  // no modification is required to the std implementation
  PrimaryLeptonGenerator::ProcessEventRecord(evrec);
}
//___________________________________________________________________________
