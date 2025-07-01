//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

\author   Kaushik Borah <kaushik.borah99 \at uky.edu>
          University of Kentucky, Lexington, KY 40506, USA
          based off ZExpAxialFormFactorModel by
          Aaron Meyer <asmeyer2012 \at uchicago.edu>
          University of Chicago, Chicago, Illinois, 60637, USA

 For the class documentation see the corresponding header file.



*/
//____________________________________________________________________________

#include <TMath.h>
#include <sstream>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/QuasiElastic/XSection/ZExpELFormFactorModel.h"
#include "Framework/Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
ZExpELFormFactorModel::ZExpELFormFactorModel() :
ELFormFactorsModelI("genie::ZExpELFormFactorModel")
{

}
//____________________________________________________________________________
ZExpELFormFactorModel::ZExpELFormFactorModel(string config) :
ELFormFactorsModelI("genie::ZExpELFormFactorModel", config)
{

}
//____________________________________________________________________________
ZExpELFormFactorModel::~ZExpELFormFactorModel()
{
}
//____________________________________________________________________________
double ZExpELFormFactorModel::Gep(const Interaction * interaction) const
{
  // calculate and return Gep
  double q2 = interaction->KinePtr()->q2();
  double zparam = this->CalculateZ(q2);
  if (zparam != zparam) // checks for nan
  {
    LOG("ZExpELFormFactorModel",pWARN) << "Undefined expansion parameter";
    return 0.;
  }
  double gep = 0.;
  for (int ki=0;ki<=fKmax+(fQ4limit ? 4 : 0);ki++)
  {
    gep = gep + TMath::Power(zparam,ki) * fZ_APn[ki];
  }

  return gep;
}
//____________________________________________________________________________
double ZExpELFormFactorModel::Gmp(const Interaction * interaction) const
{
    // calculate and return Gmp
  double q2 = interaction->KinePtr()->q2();
  double zparam = this->CalculateZ(q2);
  if (zparam != zparam) // checks for nan
  {
    LOG("ZExpELFormFactorModel",pWARN) << "Undefined expansion parameter";
    return 0.;
  }
  double gmp = 0.;
  for (int ki=0;ki<=fKmax+(fQ4limit ? 4 : 0);ki++)
  {
    gmp = gmp + TMath::Power(zparam,ki) * fZ_BPn[ki];
  }

  return gmp;
}

//____________________________________________________________________________
double ZExpELFormFactorModel::Gen(const Interaction * interaction) const
{
  // calculate and return Gen
  double q2 = interaction->KinePtr()->q2();
  double zparam = this->CalculateZ(q2);
  if (zparam != zparam) // checks for nan
  {
    LOG("ZExpELFormFactorModel",pWARN) << "Undefined expansion parameter";
    return 0.;
  }
  double gen = 0.;
  for (int ki=0;ki<=fKmax+(fQ4limit ? 4 : 0);ki++)
  {
    gen = gen + TMath::Power(zparam,ki) * fZ_ANn[ki];
  }

  return gen;
}
//____________________________________________________________________________
double ZExpELFormFactorModel::Gmn(const Interaction * interaction) const
{
  // calculate and return Gmn
  double q2 = interaction->KinePtr()->q2();
  double zparam = this->CalculateZ(q2);
  if (zparam != zparam) // checks for nan
  {
    LOG("ZExpELFormFactorModel",pWARN) << "Undefined expansion parameter";
    return 0.;
  }
  double gmn = 0.;
  for (int ki=0;ki<=fKmax+(fQ4limit ? 4 : 0);ki++)
  {
    gmn = gmn + TMath::Power(zparam,ki) * fZ_BNn[ki];
  }

  return gmn;
}
//____________________________________________________________________________
double ZExpELFormFactorModel::CalculateZ(double q2) const
{

  // calculate z expansion parameter
  double znum  = TMath::Sqrt(fTcut - q2) - TMath::Sqrt(fTcut - fT0);
  double zden  = TMath::Sqrt(fTcut - q2) + TMath::Sqrt(fTcut - fT0);

  return znum/zden;
}
//____________________________________________________________________________
void ZExpELFormFactorModel::FixCoeffs(void)
{
  //if      (fKmax < 1 ) { fKmax = 1;  }
  //else if (fKmax > 10) { fKmax = 10; }

  if (fQ4limit) this->FixQ4Limit();
  else          this->FixEL0();
}
//____________________________________________________________________________
void ZExpELFormFactorModel::FixEL0(void)
{
  // Function to fix form factor such that Gep(q2=0) = 1
  // For T0 = 0, this will set Gep0 = 1
  double zparam = this->CalculateZ(0.);
  double gep = 0.;
  for (int ki=1;ki<=fKmax;ki++)
  {
    gep = gep + TMath::Power(zparam,ki) * fZ_APn[ki];
  }
  fZ_APn[0] = fGep0 - gep;

  // Function to fix form factor such that Gmp(q2=0) = 2.7928
  // For T0 = 0, this will set Gmp0 = 2.7928

  double gmp = 0.;
  for (int ki=1;ki<=fKmax;ki++)
  {
    gmp = gmp + TMath::Power(zparam,ki) * fZ_BPn[ki];
  }
  fZ_BPn[0] = fGmp0 - gmp;

  // Function to fix form factor such that Gen(q2=0) = 0
  // For T0 = 0, this will set Gen0 = 0

  double gen = 0.;
  for (int ki=1;ki<=fKmax;ki++)
  {
    gen = gen + TMath::Power(zparam,ki) * fZ_ANn[ki];
  }
  fZ_ANn[0] = fGen0 - gen;

  // Function to fix form factor such that Gmn(q2=0) = -1.9130
  // For T0 = 0, this will set Gmn0 = -1.9130

  double gmn = 0.;
  for (int ki=1;ki<=fKmax;ki++)
  {
    gmn = gmn + TMath::Power(zparam,ki) * fZ_BNn[ki];
  }
  fZ_BNn[0] = fGmn0 - gmn;

}
//____________________________________________________________________________
void ZExpELFormFactorModel::FixQ4Limit(void)
{
  // fixes parameters such that the model limits to 1/Q^4 at large t
  // -- requires at least 5 parameters to do so --
  // 4 parameters for Q^4 behavior

  // will use AP_0 and AP_Kmax through AP_Kmax-3 to do the fitting
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
    b0 = b0 + fZ_APn[ki];
  }

  double b0z = -fGep0;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    b0z = b0z + fZ_APn[ki]*TMath::Power(z0,ki);
  }

  double b1  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    b1 = b1 + ki*fZ_APn[ki];
  }

  double b2  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    b2 = b2 + ki*(ki-1)*fZ_APn[ki];
  }

  double b3  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    b3 = b3 + ki*(ki-1)*(ki-2)*fZ_APn[ki];
  }

  // Assign new parameters
  fZ_APn[(int)kp4] = (1./denom) *                           \
  ( (b0-b0z)*kp3*kp2*kp1                                   \
  + b3*( -1. + .5*kp3*kp2*zkp1 - kp3*kp1*zkp2              \
         +.5*kp2*kp1*zkp3                               )  \
  + b2*(  3.*kp1 - kp3*kp2*kp1*zkp1                        \
         +kp3*kp1*(2*fKmax+1)*zkp2 - kp2*kp1*kp0*zkp3   )  \
  + b1*( -3.*kp2*kp1 + .5*kp3*kp2*kp2*kp1*zkp1             \
         -kp3*kp2*kp1*kp0*zkp2 + .5*kp2*kp1*kp1*kp0*zkp3)  );

  fZ_APn[(int)kp3] = (1./denom) *                           \
  ( -3.*(b0-b0z)*kp4*kp2*kp1                               \
  + b3*(  3. - kp4*kp2*zkp1 + (3./2.)*kp4*kp1*zkp2         \
         -.5*kp2*kp1*zkp4                               )  \
  + b2*( -3.*(3*fKmax+4) + kp4*kp2*(2*fKmax+3)*zkp1        \
         -3.*kp4*kp1*kp1*zkp2 + kp2*kp1*kp0*zkp4        )  \
  + b1*(  3.*kp1*(3*fKmax+8) - kp4*kp3*kp2*kp1*zkp1        \
  +(3./2.)*kp4*kp3*kp1*kp0*zkp2 - .5*kp2*kp1*kp1*kp0*zkp4) );

  fZ_APn[(int)kp2] = (1./denom) *                           \
  ( 3.*(b0-b0z)*kp4*kp3*kp1                                \
  + b3*( -3. + .5*kp4*kp3*zkp1 - (3./2.)*kp4*kp1*zkp3      \
         +kp3*kp1*zkp4                                  )  \
  + b2*(  3.*(3*fKmax+5) - kp4*kp3*kp2*zkp1                \
         +3.*kp4*kp1*kp1*zkp3 - kp3*kp1*(2*fKmax+1)*zkp4)  \
  + b1*( -3.*kp3*(3*fKmax+4) + .5*kp4*kp3*kp3*kp2*zkp1     \
    -(3./2.)*kp4*kp3*kp1*kp0*zkp3 + kp3*kp2*kp1*kp0*zkp4)  );

  fZ_APn[(int)kp1] = (1./denom) *                           \
  ( -(b0-b0z)*kp4*kp3*kp2                                  \
  + b3*(  1. - .5*kp4*kp3*zkp2 + kp4*kp2*zkp3              \
         -.5*kp3*kp2*zkp4                               )  \
  + b2*( -3.*kp2 + kp4*kp3*kp2*zkp2                        \
         -kp4*kp2*(2*fKmax+3)*zkp3 + kp3*kp2*kp1*zkp4)     \
  + b1*(  3.*kp3*kp2 - .5*kp4*kp3*kp3*kp2*zkp2             \
         +kp4*kp3*kp2*kp1*zkp3 - .5*kp3*kp2*kp2*kp1*zkp4)  );

  fZ_APn[0] = (1./denom) *                                  \
  ( -6.*b0z                                                \
  + b0*(  kp4*kp3*kp2*zkp1 - 3.*kp4*kp3*kp1*zkp2           \
         +3.*kp4*kp2*kp1*zkp3 - kp3*kp2*kp1*zkp4        )  \
  + b3*( -zkp1 + 3.*zkp2 - 3.*zkp3 + zkp4               )  \
  + b2*(  3.*kp2*zkp1 - 3.*(3*fKmax+5)*zkp2                \
         +3.*(3*fKmax+4)*zkp3 - 3.*kp1*zkp4             )  \
  + b1*( -3.*kp3*kp2*zkp1 + 3.*kp3*(3*fKmax+4)*zkp2        \
         -3.*kp1*(3*fKmax+8)*zkp3 + 3.*kp2*kp1*zkp4     )  );






  // fixes parameters such that the model limits to 1/Q^4 at large t
  // -- requires at least 5 parameters to do so --
  // 4 parameters for Q^4 behavior

  // will use BP_0 and BP_Kmax through BP_Kmax-3 to do the fitting
  // calculate some useful numbers (redundancy for readability)
  double kkp4 = (double)fKmax+4;
  double kkp3 = (double)fKmax+3;
  double kkp2 = (double)fKmax+2;
  double kkp1 = (double)fKmax+1;
  double kkp0 = (double)fKmax  ;
  //double km5 = (double)fKmax-5;
  double kz0   = this->CalculateZ(0.);
  double zkkp4 = TMath::Power(kz0,(int)kkp4);
  double zkkp3 = TMath::Power(kz0,(int)kkp3);
  double zkkp2 = TMath::Power(kz0,(int)kkp2);
  double zkkp1 = TMath::Power(kz0,(int)kkp1);

  // denominator
  double kdenom = \
  6. -    kkp4*kkp3*kkp2*zkkp1 + 3.*kkp4*kkp3*kkp1*zkkp2 \
     - 3.*kkp4*kkp2*kkp1*zkkp3 +    kkp3*kkp2*kkp1*zkkp4;

  // extra parameters (effectively constants)
  // number refers to the number of derivatives
  double kb0  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    kb0 = kb0 + fZ_BPn[ki];
  }

  double kb0z = -fGmp0;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    kb0z = kb0z + fZ_BPn[ki]*TMath::Power(kz0,ki);
  }

  double kb1  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    kb1 = kb1 + ki*fZ_BPn[ki];
  }

  double kb2  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    kb2 = kb2 + ki*(ki-1)*fZ_BPn[ki];
  }

  double kb3  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    kb3 = kb3 + ki*(ki-1)*(ki-2)*fZ_BPn[ki];
  }

  // Assign new parameters
  fZ_BPn[(int)kkp4] = (1./kdenom) *                           \
  ( (kb0-kb0z)*kkp3*kkp2*kkp1                                   \
  + kb3*( -1. + .5*kkp3*kkp2*zkkp1 - kkp3*kkp1*zkkp2              \
         +.5*kkp2*kkp1*zkkp3                               )  \
  + kb2*(  3.*kkp1 - kkp3*kkp2*kkp1*zkkp1                        \
         +kkp3*kkp1*(2*fKmax+1)*zkkp2 - kkp2*kkp1*kkp0*zkkp3   )  \
  + kb1*( -3.*kkp2*kkp1 + .5*kkp3*kkp2*kkp2*kkp1*zkkp1             \
         -kkp3*kkp2*kkp1*kkp0*zkkp2 + .5*kkp2*kkp1*kkp1*kkp0*zkkp3)  );

  fZ_BPn[(int)kkp3] = (1./kdenom) *                           \
  ( -3.*(kb0-kb0z)*kkp4*kkp2*kkp1                               \
  + kb3*(  3. - kkp4*kkp2*zkkp1 + (3./2.)*kkp4*kkp1*zkkp2         \
         -.5*kkp2*kkp1*zkkp4                               )  \
  + kb2*( -3.*(3*fKmax+4) + kkp4*kkp2*(2*fKmax+3)*zkkp1        \
         -3.*kkp4*kkp1*kkp1*zkkp2 + kkp2*kkp1*kkp0*zkkp4        )  \
  + kb1*(  3.*kkp1*(3*fKmax+8) - kkp4*kkp3*kkp2*kkp1*zkkp1        \
  +(3./2.)*kkp4*kkp3*kkp1*kkp0*zkkp2 - .5*kkp2*kkp1*kkp1*kkp0*zkkp4) );

  fZ_BPn[(int)kkp2] = (1./kdenom) *                           \
  ( 3.*(kb0-kb0z)*kkp4*kkp3*kkp1                                \
  + kb3*( -3. + .5*kkp4*kkp3*zkkp1 - (3./2.)*kkp4*kkp1*zkkp3      \
         +kkp3*kkp1*zkkp4                                  )  \
  + kb2*(  3.*(3*fKmax+5) - kkp4*kkp3*kkp2*zkkp1                \
         +3.*kkp4*kkp1*kkp1*zkkp3 - kkp3*kkp1*(2*fKmax+1)*zkkp4)  \
  + kb1*( -3.*kkp3*(3*fKmax+4) + .5*kkp4*kkp3*kkp3*kkp2*zkkp1     \
    -(3./2.)*kkp4*kkp3*kkp1*kkp0*zkkp3 + kkp3*kkp2*kkp1*kkp0*zkkp4)  );

  fZ_BPn[(int)kkp1] = (1./kdenom) *                           \
  ( -(kb0-kb0z)*kkp4*kkp3*kkp2                                  \
  + kb3*(  1. - .5*kkp4*kkp3*zkkp2 + kkp4*kkp2*zkkp3              \
         -.5*kkp3*kkp2*zkkp4                               )  \
  + kb2*( -3.*kkp2 + kkp4*kkp3*kkp2*zkkp2                        \
         -kkp4*kkp2*(2*fKmax+3)*zkkp3 + kkp3*kkp2*kkp1*zkkp4)     \
  + kb1*(  3.*kkp3*kkp2 - .5*kkp4*kkp3*kkp3*kkp2*zkkp2             \
         +kkp4*kkp3*kkp2*kkp1*zkkp3 - .5*kkp3*kkp2*kkp2*kkp1*zkkp4)  );

  fZ_BPn[0] = (1./kdenom) *                                  \
  ( -6.*kb0z                                                \
  + kb0*(  kkp4*kkp3*kkp2*zkkp1 - 3.*kkp4*kkp3*kkp1*zkkp2           \
         +3.*kkp4*kkp2*kkp1*zkkp3 - kkp3*kkp2*kkp1*zkkp4        )  \
  + kb3*( -zkkp1 + 3.*zkkp2 - 3.*zkkp3 + zkkp4               )  \
  + kb2*(  3.*kkp2*zkkp1 - 3.*(3*fKmax+5)*zkkp2                \
         +3.*(3*fKmax+4)*zkkp3 - 3.*kkp1*zkkp4             )  \
  + kb1*( -3.*kkp3*kkp2*zkkp1 + 3.*kkp3*(3*fKmax+4)*zkkp2        \
         -3.*kkp1*(3*fKmax+8)*zkkp3 + 3.*kkp2*kkp1*zkkp4     )  );






  // fixes parameters such that the model limits to 1/Q^4 at large t
  // -- requires at least 5 parameters to do so --
  // 4 parameters for Q^4 behavior

  // will use AN_0 and AN_Kmax through AN_Kmax-3 to do the fitting
  // calculate some useful numbers (redundancy for readability)
  double lp4 = (double)fKmax+4;
  double lp3 = (double)fKmax+3;
  double lp2 = (double)fKmax+2;
  double lp1 = (double)fKmax+1;
  double lp0 = (double)fKmax  ;
  //double km5 = (double)fKmax-5;
  double lz0   = this->CalculateZ(0.);
  double zlp4 = TMath::Power(lz0,(int)lp4);
  double zlp3 = TMath::Power(lz0,(int)lp3);
  double zlp2 = TMath::Power(lz0,(int)lp2);
  double zlp1 = TMath::Power(lz0,(int)lp1);

  // ldenominator
  double ldenom = \
  6. -    lp4*lp3*lp2*zlp1 + 3.*lp4*lp3*lp1*zlp2 \
     - 3.*lp4*lp2*lp1*zlp3 +    lp3*lp2*lp1*zlp4;

  // extra parameters (effectively constants)
  // number refers to the number of derivatives
  double lb0  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    lb0 = lb0 + fZ_ANn[ki];
  }

  double lb0z = -fGen0;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    lb0z = lb0z + fZ_ANn[ki]*TMath::Power(lz0,ki);
  }

  double lb1  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    lb1 = lb1 + ki*fZ_ANn[ki];
  }

  double lb2  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    lb2 = lb2 + ki*(ki-1)*fZ_ANn[ki];
  }

  double lb3  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    lb3 = lb3 + ki*(ki-1)*(ki-2)*fZ_ANn[ki];
  }

  // Assign new parameters
  fZ_ANn[(int)lp4] = (1./ldenom) *                           \
  ( (lb0-lb0z)*lp3*lp2*lp1                                   \
  + lb3*( -1. + .5*lp3*lp2*zlp1 - lp3*lp1*zlp2              \
         +.5*lp2*lp1*zlp3                               )  \
  + lb2*(  3.*lp1 - lp3*lp2*lp1*zlp1                        \
         +lp3*lp1*(2*fKmax+1)*zlp2 - lp2*lp1*lp0*zlp3   )  \
  + lb1*( -3.*lp2*lp1 + .5*lp3*lp2*lp2*lp1*zlp1             \
         -lp3*lp2*lp1*lp0*zlp2 + .5*lp2*lp1*lp1*lp0*zlp3)  );

  fZ_ANn[(int)lp3] = (1./ldenom) *                           \
  ( -3.*(lb0-lb0z)*lp4*lp2*lp1                               \
  + lb3*(  3. - lp4*lp2*zlp1 + (3./2.)*lp4*lp1*zlp2         \
         -.5*lp2*lp1*zlp4                               )  \
  + lb2*( -3.*(3*fKmax+4) + lp4*lp2*(2*fKmax+3)*zlp1        \
         -3.*lp4*lp1*lp1*zlp2 + lp2*lp1*lp0*zlp4        )  \
  + lb1*(  3.*lp1*(3*fKmax+8) - lp4*lp3*lp2*lp1*zlp1        \
  +(3./2.)*lp4*lp3*lp1*lp0*zlp2 - .5*lp2*lp1*lp1*lp0*zlp4) );

  fZ_ANn[(int)lp2] = (1./ldenom) *                           \
  ( 3.*(lb0-lb0z)*lp4*lp3*lp1                                \
  + lb3*( -3. + .5*lp4*lp3*zlp1 - (3./2.)*lp4*lp1*zlp3      \
         +lp3*lp1*zlp4                                  )  \
  + lb2*(  3.*(3*fKmax+5) - lp4*lp3*lp2*zlp1                \
         +3.*lp4*lp1*lp1*zlp3 - lp3*lp1*(2*fKmax+1)*zlp4)  \
  + lb1*( -3.*lp3*(3*fKmax+4) + .5*lp4*lp3*lp3*lp2*zlp1     \
    -(3./2.)*lp4*lp3*lp1*lp0*zlp3 + lp3*lp2*lp1*lp0*zlp4)  );

  fZ_ANn[(int)lp1] = (1./ldenom) *                           \
  ( -(lb0-lb0z)*lp4*lp3*lp2                                  \
  + lb3*(  1. - .5*lp4*lp3*zlp2 + lp4*lp2*zlp3              \
         -.5*lp3*lp2*zlp4                               )  \
  + lb2*( -3.*lp2 + lp4*lp3*lp2*zlp2                        \
         -lp4*lp2*(2*fKmax+3)*zlp3 + lp3*lp2*lp1*zlp4)     \
  + lb1*(  3.*lp3*lp2 - .5*lp4*lp3*lp3*lp2*zlp2             \
         +lp4*lp3*lp2*lp1*zlp3 - .5*lp3*lp2*lp2*lp1*zlp4)  );

  fZ_ANn[0] = (1./ldenom) *                                  \
  ( -6.*lb0z                                                \
  + lb0*(  lp4*lp3*lp2*zlp1 - 3.*lp4*lp3*lp1*zlp2           \
         +3.*lp4*lp2*lp1*zlp3 - lp3*lp2*lp1*zlp4        )  \
  + lb3*( -zlp1 + 3.*zlp2 - 3.*zlp3 + zlp4               )  \
  + lb2*(  3.*lp2*zlp1 - 3.*(3*fKmax+5)*zlp2                \
         +3.*(3*fKmax+4)*zlp3 - 3.*lp1*zlp4             )  \
  + lb1*( -3.*lp3*lp2*zlp1 + 3.*lp3*(3*fKmax+4)*zlp2        \
         -3.*lp1*(3*fKmax+8)*zlp3 + 3.*lp2*lp1*zlp4     )  );







  // fixes parameters such that the model limits to 1/Q^4 at large t
  // -- requires at least 5 parameters to do so --
  // 4 parameters for Q^4 behavior

  // will use BN_0 and BN_Kmax through BN_Kmax-3 to do the fitting
  // calculate some useful numbers (redundancy for readability)
  double llp4 = (double)fKmax+4;
  double llp3 = (double)fKmax+3;
  double llp2 = (double)fKmax+2;
  double llp1 = (double)fKmax+1;
  double llp0 = (double)fKmax  ;
  //double km5 = (double)fKmax-5;
  double llz0   = this->CalculateZ(0.);
  double zllp4 = TMath::Power(llz0,(int)llp4);
  double zllp3 = TMath::Power(llz0,(int)llp3);
  double zllp2 = TMath::Power(llz0,(int)llp2);
  double zllp1 = TMath::Power(llz0,(int)llp1);

  // denominator
  double lldenom = \
  6. -    llp4*llp3*llp2*zllp1 + 3.*llp4*llp3*llp1*zllp2 \
     - 3.*llp4*llp2*llp1*zllp3 +    llp3*llp2*llp1*zllp4;

  // extra parameters (effectively constants)
  // number refers to the number of derivatives
  double llb0  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    llb0 = llb0 + fZ_BNn[ki];
  }

  double llb0z = -fGmn0;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    llb0z = llb0z + fZ_BNn[ki]*TMath::Power(llz0,ki);
  }

  double llb1  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    llb1 = llb1 + ki*fZ_BNn[ki];
  }

  double llb2  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    llb2 = llb2 + ki*(ki-1)*fZ_BNn[ki];
  }

  double llb3  = 0.;
  for (int ki = 1;ki <= fKmax;ki++)
  {
    llb3 = llb3 + ki*(ki-1)*(ki-2)*fZ_BNn[ki];
  }

  // Assign new parameters
  fZ_BNn[(int)llp4] = (1./lldenom) *                           \
  ( (llb0-llb0z)*llp3*llp2*llp1                                   \
  + llb3*( -1. + .5*llp3*llp2*zllp1 - llp3*llp1*zllp2              \
         +.5*llp2*llp1*zllp3                               )  \
  + llb2*(  3.*llp1 - llp3*llp2*llp1*zllp1                        \
         +llp3*llp1*(2*fKmax+1)*zllp2 - llp2*llp1*llp0*zllp3   )  \
  + llb1*( -3.*llp2*llp1 + .5*llp3*llp2*llp2*llp1*zllp1             \
         -llp3*llp2*llp1*llp0*zllp2 + .5*llp2*llp1*llp1*llp0*zllp3)  );

  fZ_BNn[(int)llp3] = (1./lldenom) *                           \
  ( -3.*(llb0-llb0z)*llp4*llp2*llp1                               \
  + llb3*(  3. - llp4*llp2*zllp1 + (3./2.)*llp4*llp1*zllp2         \
         -.5*llp2*llp1*zllp4                               )  \
  + llb2*( -3.*(3*fKmax+4) + llp4*llp2*(2*fKmax+3)*zllp1        \
         -3.*llp4*llp1*llp1*zllp2 + llp2*llp1*llp0*zllp4        )  \
  + llb1*(  3.*llp1*(3*fKmax+8) - llp4*llp3*llp2*llp1*zllp1        \
  +(3./2.)*llp4*llp3*llp1*llp0*zllp2 - .5*llp2*llp1*llp1*llp0*zllp4) );

  fZ_BNn[(int)llp2] = (1./lldenom) *                           \
  ( 3.*(llb0-llb0z)*llp4*llp3*llp1                                \
  + llb3*( -3. + .5*llp4*llp3*zllp1 - (3./2.)*llp4*llp1*zllp3      \
         +llp3*llp1*zllp4                                  )  \
  + llb2*(  3.*(3*fKmax+5) - llp4*llp3*llp2*zllp1                \
         +3.*llp4*llp1*llp1*zllp3 - llp3*llp1*(2*fKmax+1)*zllp4)  \
  + llb1*( -3.*llp3*(3*fKmax+4) + .5*llp4*llp3*llp3*llp2*zllp1     \
    -(3./2.)*llp4*llp3*llp1*llp0*zllp3 + llp3*llp2*llp1*llp0*zllp4)  );

  fZ_BNn[(int)llp1] = (1./lldenom) *                           \
  ( -(llb0-llb0z)*llp4*llp3*llp2                                  \
  + llb3*(  1. - .5*llp4*llp3*zllp2 + llp4*llp2*zllp3              \
         -.5*llp3*llp2*zllp4                               )  \
  + llb2*( -3.*llp2 + llp4*llp3*llp2*zllp2                        \
         -llp4*llp2*(2*fKmax+3)*zllp3 + llp3*llp2*llp1*zllp4)     \
  + llb1*(  3.*llp3*llp2 - .5*llp4*llp3*llp3*llp2*zllp2             \
         +llp4*llp3*llp2*llp1*zllp3 - .5*llp3*llp2*llp2*llp1*zllp4)  );

  fZ_BNn[0] = (1./lldenom) *                                  \
  ( -6.*llb0z                                                \
  + llb0*(  llp4*llp3*llp2*zllp1 - 3.*llp4*llp3*llp1*zllp2           \
         +3.*llp4*llp2*llp1*zllp3 - llp3*llp2*llp1*zllp4        )  \
  + llb3*( -zllp1 + 3.*zllp2 - 3.*zllp3 + zllp4               )  \
  + llb2*(  3.*llp2*zllp1 - 3.*(3*fKmax+5)*zllp2                \
         +3.*(3*fKmax+4)*zllp3 - 3.*llp1*zllp4             )  \
  + llb1*( -3.*llp3*llp2*zllp1 + 3.*llp3*(3*fKmax+4)*zllp2        \
         -3.*llp1*(3*fKmax+8)*zllp3 + 3.*llp2*llp1*zllp4     )  );




}
//____________________________________________________________________________
void ZExpELFormFactorModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ZExpELFormFactorModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void ZExpELFormFactorModel::LoadConfig(void )
{
// get config options from the configuration registry or set defaults
// from the global parameter list

  GetParam( "QEL-Q4limit", fQ4limit ) ;
  GetParam( "QEL-Kmax", fKmax ) ;

  GetParam( "QEL-T0", fT0 ) ;
  GetParam( "QEL-Tcut", fTcut ) ;

  GetParam( "QEL-Gep0", fGep0 ) ;
  GetParam( "QEL-Gmp0", fGmp0 ) ;
  GetParam( "QEL-Gen0", fGen0 ) ;
  GetParam( "QEL-Gmn0", fGmn0 ) ;
  assert(fKmax > 0);

  // z expansion coefficients
  std::vector<double> tmp_fZ_APn, tmp_fZ_BPn, tmp_fZ_ANn, tmp_fZ_BNn;
  if(fQ4limit){
    fZ_APn.resize(fKmax+5);
    fZ_BPn.resize(fKmax+5);
    fZ_ANn.resize(fKmax+5);
    fZ_BNn.resize(fKmax+5);
  }
  else{
    fZ_APn.resize(fKmax+1);
    fZ_BPn.resize(fKmax+1);
    fZ_ANn.resize(fKmax+1);
    fZ_BNn.resize(fKmax+1);
  }
  if(this->GetParamVect("QEL-Z_AP", tmp_fZ_APn) != fKmax){
    LOG("ZExpELFormFactorModel",pERROR) << "Wrong size of AP coefficients " << tmp_fZ_APn.size();
    exit(1);
  }
  if(this->GetParamVect("QEL-Z_BP", tmp_fZ_BPn) != fKmax){
    LOG("ZExpELFormFactorModel",pERROR) << "Wrong size of BP coefficients " << tmp_fZ_BPn.size();
    exit(1);
  }
  if(this->GetParamVect("QEL-Z_AN", tmp_fZ_ANn) != fKmax){
    LOG("ZExpELFormFactorModel",pERROR) << "Wrong size of AN coefficients " << tmp_fZ_ANn.size();
    exit(1);
  }
  if(this->GetParamVect("QEL-Z_BN", tmp_fZ_BNn) != fKmax){
    LOG("ZExpELFormFactorModel",pERROR) << "Wrong size of BN coefficients " << tmp_fZ_BNn.size();
    exit(1);
  }

  // load the user-defined coefficient values
  // -- AP0 and APn for n<fKmax are calculated from other means
  for (int ip=1;ip<fKmax+1;ip++) {
    fZ_APn[ip] = tmp_fZ_APn[ip-1];
    fZ_BPn[ip] = tmp_fZ_BPn[ip-1];
    fZ_ANn[ip] = tmp_fZ_ANn[ip-1];
    fZ_BNn[ip] = tmp_fZ_BNn[ip-1];
  }

  this->FixCoeffs();
  Interaction * interaction = new Interaction();
  for (int i=0;i<10;i++)  {
    double Q2 = i*0.1;
    interaction->KinePtr()->SetQ2( Q2);
    LOG("ZExpELFormFactorModel",pNOTICE)
      << "Q2=" <<Q2
      <<" : Gep(Q2)= " <<this->Gep( interaction);
  }
  for (int i=0;i<10;i++)  {
    double Q2 = i*0.1;
    interaction->KinePtr()->SetQ2( Q2);
    LOG("ZExpELFormFactorModel",pNOTICE)
      << "Q2=" <<Q2
      <<" : Gmp(Q2)= " <<this->Gmp( interaction);
  }
  for (int i=0;i<10;i++)  {
    double Q2 = i*0.1;
    interaction->KinePtr()->SetQ2( Q2);
    LOG("ZExpELFormFactorModel",pNOTICE)
      << "Q2=" <<Q2
      <<" : Gen(Q2)= " <<this->Gen( interaction);
  }
  for (int i=0;i<10;i++)  {
    double Q2 = i*0.1;
    interaction->KinePtr()->SetQ2( Q2);
    LOG("ZExpELFormFactorModel",pNOTICE)
      << "Q2=" <<Q2
      <<" : Gmn(Q2)= " <<this->Gmn( interaction);
  }
  delete interaction;
}
//____________________________________________________________________________

