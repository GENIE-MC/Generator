//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Pauli Kehayias (Tufts Univ.)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 13, 2009 - PK
   This file was added in v2.5.1

*/
//____________________________________________________________________________

#include <iostream>

#include "validation/Merenyi/MerenyiNuEvent.h"

using namespace genie;
using namespace genie::vld_merenyi;

//____________________________________________________________________________
MerenyiNuEvent::MerenyiNuEvent()
{
  this->zeroAll();
}
//____________________________________________________________________________
MerenyiNuEvent::MerenyiNuEvent(int numHadrons)
{
  this->zeroAll();
  numHad = numHadrons;
}
//____________________________________________________________________________
MerenyiNuEvent::MerenyiNuEvent(
   int numHadrons, int interactionType, double nuEnergy,
   double lepXMom, double lepYMom, double lepZMom, double lepEnergy)
{
  this->zeroAll();
  
  numHad  = numHadrons;
  intType = interactionType;
  nuE     = nuEnergy;
  
  lepton.SetMomentum(lepXMom, lepYMom, lepZMom, lepEnergy);
}
//____________________________________________________________________________
void MerenyiNuEvent::zeroAll()
{
  hadList.clear();
  
  lepton.SetMomentum(0,0,0,0);

  nuE     = 0;
  numGam  = 0;
  numPi0  = 0;
  numPiP  = 0;
  numPiM  = 0;
  numN    = 0;
  numP    = 0;                                                         
  numStr0 = 0;
  numStrP = 0;
  numStrM = 0;
}
//____________________________________________________________________________
void MerenyiNuEvent::addParticle(
            double particleID, double px, double py, double pz, double E)
{
  TParticle hadron;
  hadron.SetPdgCode((int) particleID);
  hadron.SetMomentum(px, py, pz, E);
  hadList.push_back(hadron);
}
//____________________________________________________________________________
double MerenyiNuEvent::getLepAngCos()
{
  double lepPz  = lepton.Pz();
  double pTot   = lepton.P();
  double cosAng = lepPz / pTot;
  return cosAng;
}
//____________________________________________________________________________

