//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 09, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "EVGModules/PrimaryVtxGenerator.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::units;
using namespace genie::utils;

//___________________________________________________________________________
PrimaryVtxGenerator::PrimaryVtxGenerator() :
EventRecordVisitorI("genie::PrimaryVtxGenerator")
{

}
//___________________________________________________________________________
PrimaryVtxGenerator::PrimaryVtxGenerator(string config) :
EventRecordVisitorI("genie::PrimaryVtxGenerator", config)
{

}
//___________________________________________________________________________
PrimaryVtxGenerator::~PrimaryVtxGenerator()
{

}
//___________________________________________________________________________
void PrimaryVtxGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("VtxGenerator", pDEBUG) << "Generating an interaction vertex";

  Interaction * interaction = evrec->Summary();

  const InitialState & init_state = interaction->InitState();
  const Target & target = init_state.Tgt();

  if (!target.IsNucleus()) {
    LOG("VtxGenerator", pINFO) << "No nuclear target found - Vtx = (0,0,0)";
    return;
  }

  // compute the target radius
  double Rn = nuclear::Radius(target); // = R*A^(1/3), GeV
  Rn  /= units::m; // GeV -> m

  LOG("VtxGenerator", pINFO) << "Rnucl = " << Rn << " m";

  // generate a random position inside a spherical nucleus with radius R
  RandomGen * rnd = RandomGen::Instance();

  double R      =     Rn*(rnd->RndGen().Rndm());  // [0,Rn]
  double phi    =  2*kPi*(rnd->RndGen().Rndm());  // [0,2pi]
  double cos8   = -1 + 2*(rnd->RndGen().Rndm());  // [-1,1]
  double sin8   = TMath::Sqrt(1-cos8*cos8);
  double cosphi = TMath::Cos(phi);
  double sinphi = TMath::Sin(phi);
  double x      = Rn * sin8 * cosphi;
  double y      = Rn * sin8 * sinphi;
  double z      = Rn * cos8;

  LOG("VtxGenerator", pINFO)
            << "R = " << R/m << " m, phi = " << phi << ", cos8 = " << cos8;
  LOG("VtxGenerator", pINFO)
        << "x = " << x/m << " m, y = " << y/m << " m, z = " << z/m << " m";

  GHepParticle * probe = evrec->Probe();
  assert(probe);
  probe->SetPosition(x, y, z, 0.);

  GHepParticle * nucleon = evrec->HitNucleon();
  assert(nucleon);
  nucleon->SetPosition(x, y, z, 0.);
}
//___________________________________________________________________________
