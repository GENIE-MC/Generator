//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 25, 2009 - CA
   Was first added in the development version 2.5.1

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Utils/PREM.h"

//___________________________________________________________________________
double genie::utils::prem::Density(double r)
{
// Return the Earth density according to the PREM model
// Inputs: 
//   r, Distance from the centre of the Earth (in km)
// Outputs: 
//   rho, Earth density (in gr/cm3)
//

  r = TMath::Max(0., r);

  double rho = 0.;
  double x   = r/constants::kREarth;

  if (r <= 1221.5 ) 
  { 
    rho = 13.0885 - 8.8381*x*x; 
  }
  else if (r >  1221.5 && r <= 3480.0 ) 
  { 
    rho = 12.5815 - 1.2638*x - 3.6426*x*x - 5.5281*x*x*x; 
  }
  else if (r >  3480.0 && r <= 5701.0 ) 
  { 
    rho = 7.9565  - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x; 
  } 
  else if (r >  5701.0 && r <= 5771.0 ) 
  { 
    rho = 5.3197  - 1.4836*x; 
  } 
  else if (r >  5771.0 && r <= 5971.0 ) 
  { 
    rho = 11.2494 - 8.0298*x; 
  } 
  else if (r >  5971.0 && r <= 6151.0 ) 
  { 
    rho = 7.1089  - 3.8045*x; 
  } 
  else if (r >  6151.0 && r <= 6346.6 ) 
  { 
    rho = 2.691   + 0.6924*x; 
  } 
  else if (r >  6346.6 && r <= 6356.0 ) 
  { 
    rho = 2.90; 
  } 
  else if (r >  6356.0 && r <= 6368.0 ) 
  { 
    rho = 2.60; 
  } 
  else if (r >  6368.0 && r <= constants::kREarth) 
  { 
    rho = 1.02; 
  } 

  return rho;
}
//___________________________________________________________________________
