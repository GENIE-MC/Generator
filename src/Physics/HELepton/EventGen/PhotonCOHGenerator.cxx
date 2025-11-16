//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/PhotonCOHGenerator.h"

//___________________________________________________________________________
PhotonCOHGenerator::PhotonCOHGenerator() :
EventRecordVisitorI("genie::PhotonCOHGenerator")
{
}
//___________________________________________________________________________
PhotonCOHGenerator::PhotonCOHGenerator(string config) :
EventRecordVisitorI("genie::PhotonCOHGenerator", config)
{
}
//___________________________________________________________________________
PhotonCOHGenerator::~PhotonCOHGenerator()
{

}
//___________________________________________________________________________
void PhotonCOHGenerator::ProcessEventRecord(GHepRecord * event) const
{

  fWDecayer->ProcessEventRecord(event);

}
//____________________________________________________________________________
void PhotonCOHGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHGenerator::LoadConfig(void)
{

  fWDecayer = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("WDecayer"));
  assert(fWDecayer);

}
//____________________________________________________________________________
