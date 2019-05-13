//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

\author   Aaron Meyer <asmeyer2012 \at uchicago.edu>
          based off DipoleELFormFactorsModel by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <sstream>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/QuasiElastic/XSection/ZExpAxialFormFactorModel.h"
#include "Framework/Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
ZExpAxialFormFactorModel::ZExpAxialFormFactorModel() :
AxialFormFactorModelI("genie::ZExpAxialFormFactorModel")
{

}
//____________________________________________________________________________
ZExpAxialFormFactorModel::ZExpAxialFormFactorModel(string config) :
AxialFormFactorModelI("genie::ZExpAxialFormFactorModel", config)
{

}
//____________________________________________________________________________
ZExpAxialFormFactorModel::~ZExpAxialFormFactorModel()
{
  delete[] fZ_An;
}
//____________________________________________________________________________
double ZExpAxialFormFactorModel::FA(const Interaction * interaction) const
{
  // calculate and return FA
  double q2 = interaction->KinePtr()->q2();
  double zparam = this->CalculateZ(q2);
  if (zparam != zparam) // checks for nan
  {
    LOG("ZExpAxialFormFactorModel",pWARN) << "Undefined expansion parameter";
    return 0.;
  }
  double fa = 0.;
  for (int ki=0;ki<=fKmax+(fQ4limit ? 4 : 0);ki++)
  {
    fa = fa + TMath::Power(zparam,ki) * fZ_An[ki];
  }

  return fa;
}
//____________________________________________________________________________
double ZExpAxialFormFactorModel::CalculateZ(double q2) const
{

  // calculate z expansion parameter
  double znum  = TMath::Sqrt(fTcut - q2) - TMath::Sqrt(fTcut - fT0);
  double zden  = TMath::Sqrt(fTcut - q2) + TMath::Sqrt(fTcut - fT0);

  return znum/zden;
}
//____________________________________________________________________________
void ZExpAxialFormFactorModel::FixCoeffs(void)
{
  //if      (fKmax < 1 ) { fKmax = 1;  }
  //else if (fKmax > 10) { fKmax = 10; }

  if (fQ4limit) this->FixQ4Limit();
  else          this->FixA0();
}
//____________________________________________________________________________
void ZExpAxialFormFactorModel::FixA0(void)
{
  // Function to fix form factor such that FA(q2=0) = gA
  // For T0 = 0, this will set A0 = gA
  double zparam = this->CalculateZ(0.);
  double fa = 0.;
  for (int ki=1;ki<=fKmax;ki++)                    
  {
    fa = fa + TMath::Power(zparam,ki) * fZ_An[ki];
  }
  fZ_An[0] = fFA0 - fa; 

}
//____________________________________________________________________________
void ZExpAxialFormFactorModel::FixQ4Limit(void)
{
  // fixes parameters such that the model limits to 1/Q^4 at large t
  // -- requires at least 5 parameters to do so --
  // 4 parameters for Q^4 behavior, 1 for normalization to FA(q2=0)=gA

  // will use A_0 and A_Kmax through A_Kmax-3 to do the fitting
  // calculate some useful numbers (redundancy for readability)
  double kp4 = (double)fKmax+4;
  double kp3 = (double)fKmax+3;
  double kp2 = (double)fKmax+2;
  double kp1 = (double)fKmax+1;
  double kp0 = (double)fKmax  ;
  //double km5 = (double)fKmax-5;
  double z0   = this->CalculateZ(0.);
  double zkp4 = TMath::Power(z0,(int)kp4);
  double zkp3 = TMath::Power(z0,(int)kp3);
  double zkp2 = TMath::Power(z0,(int)kp2);
  double zkp1 = TMath::Power(z0,(int)kp1);

  // denominator
  double denom = \
  6. -    kp4*kp3*kp2*zkp1 + 3.*kp4*kp3*kp1*zkp2 \
     - 3.*kp4*kp2*kp1*zkp3 +    kp3*kp2*kp1*zkp4;

  // extra parameters (effectively constants)
  // number refers to the number of derivatives
  double b0  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    b0 = b0 + fZ_An[ki];
  }

  double b0z = -fFA0;
  for (int ki = 1;ki <= fKmax;ki++) 
  {
    b0z = b0z + fZ_An[ki]*TMath::Power(z0,ki);
  }

  double b1  = 0.;
  for (int ki = 1;ki <= fKmax;ki++) 
  {
    b1 = b1 + ki*fZ_An[ki];
  }
  
  double b2  = 0.;
  for (int ki = 1;ki <= fKmax;ki++) 
  {
    b2 = b2 + ki*(ki-1)*fZ_An[ki];
  }

  double b3  = 0.;
  for (int ki = 1;ki <= fKmax;ki++) 
  {
    b3 = b3 + ki*(ki-1)*(ki-2)*fZ_An[ki];
  }

  // Assign new parameters
  fZ_An[(int)kp4] = (1./denom) *                           \
  ( (b0-b0z)*kp3*kp2*kp1                                   \
  + b3*( -1. + .5*kp3*kp2*zkp1 - kp3*kp1*zkp2              \
         +.5*kp2*kp1*zkp3                               )  \
  + b2*(  3.*kp1 - kp3*kp2*kp1*zkp1                        \
         +kp3*kp1*(2*fKmax+1)*zkp2 - kp2*kp1*kp0*zkp3   )  \
  + b1*( -3.*kp2*kp1 + .5*kp3*kp2*kp2*kp1*zkp1             \
         -kp3*kp2*kp1*kp0*zkp2 + .5*kp2*kp1*kp1*kp0*zkp3)  );

  fZ_An[(int)kp3] = (1./denom) *                           \
  ( -3.*(b0-b0z)*kp4*kp2*kp1                               \
  + b3*(  3. - kp4*kp2*zkp1 + (3./2.)*kp4*kp1*zkp2         \
         -.5*kp2*kp1*zkp4                               )  \
  + b2*( -3.*(3*fKmax+4) + kp4*kp2*(2*fKmax+3)*zkp1        \
         -3.*kp4*kp1*kp1*zkp2 + kp2*kp1*kp0*zkp4        )  \
  + b1*(  3.*kp1*(3*fKmax+8) - kp4*kp3*kp2*kp1*zkp1        \
  +(3./2.)*kp4*kp3*kp1*kp0*zkp2 - .5*kp2*kp1*kp1*kp0*zkp4) );

  fZ_An[(int)kp2] = (1./denom) *                           \
  ( 3.*(b0-b0z)*kp4*kp3*kp1                                \
  + b3*( -3. + .5*kp4*kp3*zkp1 - (3./2.)*kp4*kp1*zkp3      \
         +kp3*kp1*zkp4                                  )  \
  + b2*(  3.*(3*fKmax+5) - kp4*kp3*kp2*zkp1                \
         +3.*kp4*kp1*kp1*zkp3 - kp3*kp1*(2*fKmax+1)*zkp4)  \
  + b1*( -3.*kp3*(3*fKmax+4) + .5*kp4*kp3*kp3*kp2*zkp1     \
    -(3./2.)*kp4*kp3*kp1*kp0*zkp3 + kp3*kp2*kp1*kp0*zkp4)  );

  fZ_An[(int)kp1] = (1./denom) *                           \
  ( -(b0-b0z)*kp4*kp3*kp2                                  \
  + b3*(  1. - .5*kp4*kp3*zkp2 + kp4*kp2*zkp3              \
         -.5*kp3*kp2*zkp4                               )  \
  + b2*( -3.*kp2 + kp4*kp3*kp2*zkp2                        \
         -kp4*kp2*(2*fKmax+3)*zkp3 + kp3*kp2*kp1*zkp4)     \
  + b1*(  3.*kp3*kp2 - .5*kp4*kp3*kp3*kp2*zkp2             \
         +kp4*kp3*kp2*kp1*zkp3 - .5*kp3*kp2*kp2*kp1*zkp4)  );

  fZ_An[0] = (1./denom) *                                  \
  ( -6.*b0z                                                \
  + b0*(  kp4*kp3*kp2*zkp1 - 3.*kp4*kp3*kp1*zkp2           \
         +3.*kp4*kp2*kp1*zkp3 - kp3*kp2*kp1*zkp4        )  \
  + b3*( -zkp1 + 3.*zkp2 - 3.*zkp3 + zkp4               )  \
  + b2*(  3.*kp2*zkp1 - 3.*(3*fKmax+5)*zkp2                \
         +3.*(3*fKmax+4)*zkp3 - 3.*kp1*zkp4             )  \
  + b1*( -3.*kp3*kp2*zkp1 + 3.*kp3*(3*fKmax+4)*zkp2        \
         -3.*kp1*(3*fKmax+8)*zkp3 + 3.*kp2*kp1*zkp4     )  );
}
//____________________________________________________________________________
void ZExpAxialFormFactorModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ZExpAxialFormFactorModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void ZExpAxialFormFactorModel::LoadConfig(void)
{
// get config options from the configuration registry or set defaults 
// from the global parameter list

  GetParam( "QEL-Q4limit", fQ4limit ) ;
  GetParam( "QEL-Kmax", fKmax ) ;

  GetParam( "QEL-T0", fT0 ) ;
  GetParam( "QEL-T0", fT0 ) ;
  GetParam( "QEL-Tcut", fTcut ) ;

  GetParam( "QEL-FA0", fFA0 ) ;
  assert(fKmax > 0);

  // z expansion coefficients
  if (fQ4limit) fZ_An = new double [fKmax+5];
  else          fZ_An = new double [fKmax+1];

  // load the user-defined coefficient values
  // -- A0 and An for n<fKmax are calculated from other means
  for (int ip=1;ip<fKmax+1;ip++) {
    ostringstream alg_key;
    alg_key << "QEL-Z_A" << ip;
    GetParam( alg_key.str(), fZ_An[ip] ) ;
  }

  this->FixCoeffs();
}
//____________________________________________________________________________

