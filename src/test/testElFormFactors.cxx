//____________________________________________________________________________
/*!

\program testElFormFactors

\brief   program used for testing / debugging the elastic form factors

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <string>

#include <TFile.h>
#include <TNtupleD.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Base/ELFormFactors.h"
#include "Base/ELFormFactorsModelI.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using std::string;

//__________________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  TNtupleD * elffnt = new TNtupleD("elffnt","","Q2:mod:Gep:Gmp:Gen:Gmn");

  AlgFactory * algf = AlgFactory::Instance();

  const ELFormFactorsModelI * dipole =
      dynamic_cast<const ELFormFactorsModelI *> (
        algf->GetAlgorithm("genie::DipoleELFormFactorsModel", "Default"));
  const ELFormFactorsModelI * bba2003 =
      dynamic_cast<const ELFormFactorsModelI *> (
        algf->GetAlgorithm("genie::BBA03ELFormFactorsModel", "Default"));
  const ELFormFactorsModelI * bba2005 =
      dynamic_cast<const ELFormFactorsModelI *> (
        algf->GetAlgorithm("genie::BBA05ELFormFactorsModel", "Default"));

  ELFormFactors elff;

  Interaction * interaction = 
                   Interaction::QELCC(1056026000,kPdgProton,kPdgNuMu,10);

  for(int iq=0; iq<100; iq++) {

   double Q2 = iq*0.01 + 0.01;
   interaction->KinePtr()->SetQ2(Q2);

   elff.SetModel(dipole);
   elff.Calculate(interaction);
   elffnt->Fill(Q2,0,elff.Gep(),elff.Gmp(),elff.Gen(),elff.Gmn());

   elff.SetModel(bba2003);
   elff.Calculate(interaction);
   elffnt->Fill(Q2,1,elff.Gep(),elff.Gmp(),elff.Gen(),elff.Gmn());

   elff.SetModel(bba2005);
   elff.Calculate(interaction);
   elffnt->Fill(Q2,2,elff.Gep(),elff.Gmp(),elff.Gen(),elff.Gmn());
  }

  TFile f("./elff.root","recreate");
  elffnt->Write();
  f.Close();

  return 0;
}
//__________________________________________________________________________
