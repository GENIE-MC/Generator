//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI() :
Algorithm()
{

}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
XSecAlgorithmI::~XSecAlgorithmI()
{

}
//___________________________________________________________________________
bool XSecAlgorithmI::ValidKinematics(const Interaction* interaction) const
{
// can offer common implementation for all concrete x-section models because
// the input interaction is aware of its kinematic limits

  if ( interaction->TestBit(kISkipKinematicChk) ) return true;

  const KPhaseSpace& kps = interaction->PhaseSpace();

  if ( ! kps.IsAboveThreshold() ) {
     LOG("XSecBase", pINFO)  << "*** Below energy threshold";
     return false;
  }
  if ( ! kps.IsAllowed() ) {
     LOG("XSecBase", pINFO)  << "*** Not in allowed kinematical space";
     return false;
  }
  return true;
}
//___________________________________________________________________________

