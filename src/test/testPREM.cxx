//____________________________________________________________________________
/*!

\program gtestPREM

\brief   Test tehe PREM model

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created Aug 25, 2009

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>

#include "Conventions/Constants.h"
#include "Utils/PREM.h"

using namespace genie;

//____________________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  TNtuple * earth_density = new TNtuple("earth_density","","r:rho");

  const double dr   = 1; // km
  const double rmax = constants::kREarth + 10; // km

  double r = 0;
  while(r < rmax) {
     double rho = utils::prem::Density(r*units::km) / units::g_cm3; 
     earth_density->Fill(r, rho);
     r += dr;
  }

  TFile f("./prem.root","recreate");
  earth_density->Write();
  f.Close();

  return 0;
}
//____________________________________________________________________________
