//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 John Plows <komninos-john.plows \at physics.ox.ac.uk>
 University of Oxford

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Physics/BeamHNL/HNLGeomRecordVisitorI.h"

using namespace genie;
using namespace genie::hnl;

//____________________________________________________________________________
GeomRecordVisitorI::GeomRecordVisitorI() :
  EventRecordVisitorI()
{

}
//____________________________________________________________________________
GeomRecordVisitorI::GeomRecordVisitorI(string name) :
  EventRecordVisitorI(name)
{

}
//____________________________________________________________________________
GeomRecordVisitorI::GeomRecordVisitorI(string name, string config) :
  EventRecordVisitorI(name, config)
{

}
//____________________________________________________________________________
GeomRecordVisitorI::~GeomRecordVisitorI()
{

}
//____________________________________________________________________________
