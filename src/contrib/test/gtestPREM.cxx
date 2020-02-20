//____________________________________________________________________________
/*!

\program gtestPREM

\brief   Test tehe PREM model

\author  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

\created Aug 25, 2009

\cpright Copyright (c) 2003-2020, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtuple.h>

#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Utils/PREM.h"

using namespace genie;

//____________________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  TNtuple * earth_density = new TNtuple("earth_density","","r:rho");

  const double dr   = 1. * units::km; 
  const double rmax = constants::kREarth;

  double r = 0;
  while(r < rmax) {
     double rho = utils::prem::Density(r); 
     earth_density->Fill(r/units::km, rho/units::g_cm3);
     r += dr;
  }

  TFile f("./prem.root","recreate");
  earth_density->Write();
  f.Close();

  return 0;
}
//____________________________________________________________________________
