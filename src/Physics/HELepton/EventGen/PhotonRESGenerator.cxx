//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/PhotonRESGenerator.h"

//___________________________________________________________________________
PhotonRESGenerator::PhotonRESGenerator() :
EventRecordVisitorI("genie::PhotonRESGenerator")
{

}
//___________________________________________________________________________
PhotonRESGenerator::PhotonRESGenerator(string config) :
EventRecordVisitorI("genie::PhotonRESGenerator", config)
{

}
//___________________________________________________________________________
PhotonRESGenerator::~PhotonRESGenerator()
{

}
//___________________________________________________________________________
void PhotonRESGenerator::ProcessEventRecord(GHepRecord * event) const
{

  fWDecayer->ProcessEventRecord(event);

}
//____________________________________________________________________________
void PhotonRESGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESGenerator::LoadConfig(void)
{

  fWDecayer = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("WDecayer"));
  assert(fWDecayer);

}
//____________________________________________________________________________
