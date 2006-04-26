//____________________________________________________________________________
/*!

\program testMuELoss

\brief   test program for the MuELoss utility package

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created March 10, 2006

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Units.h"
#include "MuELoss/MuELossI.h"
#include "MuELoss/MuELMaterial.h"
#include "MuELoss/MuELProcess.h"
#include "Messenger/Messenger.h"

using std::string;
using namespace genie;
using namespace genie::mueloss;

int main(int argc, char ** argv)
{
  const int N = 14;
  double E[N] = {1,5,10,15,20,30,50,100,200,500, 1000,2000, 5000, 9000}; //GeV

  MuELMaterial_t mt = eMuIron;

  AlgFactory * algf = AlgFactory::Instance();

  const MuELossI * betheBloch = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                         "genie::mueloss::BetheBlochModel","Default"));

  const MuELossI * petrukhinShestakov = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                 "genie::mueloss::PetrukhinShestakovModel","Default"));

  const MuELossI * kokoulinPetroukhin = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                  "genie::mueloss::KokoulinPetrukhinModel","Default"));

  const MuELossI * bezroukovBugaev = 
         dynamic_cast<const MuELossI *> (algf->GetAlgorithm(
                     "genie::mueloss::BezrukovBugaevModel","Default"));

  assert ( betheBloch         );
  assert ( petrukhinShestakov );
  assert ( kokoulinPetroukhin );
  assert ( bezroukovBugaev    );

  double myunits_conversion = units::GeV/(units::g/units::cm2); 
  string myunits_name       = " GeV/(gr/cm^2)";

  LOG("Main", pINFO) 
     << "---------- Computing/Printing muon energy losses in "
                       << MuELMaterial::AsString(mt) << " ----------";

  for(int i=0; i<N; i++)  {
      
     //-------- due to: ionization 
     LOG("Main", pINFO) 
       << "Process: " << MuELProcess::AsString(betheBloch->Process())
       << ", Model: " << betheBloch->Id().Key() 
       << " : \n -dE/dx(E=" << E[i] << ") = " 
       << betheBloch->dE_dx(E[i],mt) / myunits_conversion << myunits_name;

     //-------- due to: bremsstrahlung
     LOG("Main", pINFO) 
       << "Process: " << MuELProcess::AsString(petrukhinShestakov->Process())
       << ", Model: " << petrukhinShestakov->Id().Key() 
       << " : \n -dE/dx(E=" << E[i] << ") = " 
       << petrukhinShestakov->dE_dx(E[i],mt) / myunits_conversion << myunits_name;

     //-------- due to: e-e+ pair production
     LOG("Main", pINFO) 
       << "Process: " << MuELProcess::AsString(kokoulinPetroukhin->Process())
       << ", Model: " << kokoulinPetroukhin->Id().Key() 
       << " : \n -dE/dx(E=" << E[i] << ") = " 
       << kokoulinPetroukhin->dE_dx(E[i],mt) / myunits_conversion << myunits_name;

     //-------- due to: photonuclear interactions
     LOG("Main", pINFO) 
       << "Process: " << MuELProcess::AsString(bezroukovBugaev->Process())
       << ", Model: " << bezroukovBugaev->Id().Key() 
       << " : \n -dE/dx(E=" << E[i] << ") = " 
       << bezroukovBugaev->dE_dx(E[i],mt) / myunits_conversion << myunits_name
       << "\n\n";
  }
  return 0;
}

