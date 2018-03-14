//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          adapted from  fortran code provided by
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>,
          Joint Institute for Nuclear Research,  Institute for Theoretical and Experimental Physics
          Vladimir Lyubushkin,
          Joint Institute for Nuclear Research
          Vadim Naumov <vnaumov@theor.jinr.ru>,
          Joint Institute for Nuclear Research
          based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/QuasiElastic/XSection/BBA07ELFormFactorsModel.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Conventions/Constants.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BBA07ELFormFactorsModel::BBA07ELFormFactorsModel() :
ELFormFactorsModelI("genie::BBA07ELFormFactorsModel")
{

}
//____________________________________________________________________________
BBA07ELFormFactorsModel::BBA07ELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::BBA07ELFormFactorsModel", config)
{

}
//____________________________________________________________________________
BBA07ELFormFactorsModel::~BBA07ELFormFactorsModel()
{

}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gep(const Interaction * interaction) const
{
  const Kinematics & kine   = interaction->Kine();
  double q2 = kine.q2();
  double M2 = kProtonMass*kProtonMass;
  double tp = -q2/(4*M2);
  double xp  = 2.0/(1.0+TMath::Sqrt(1.0+1.0/tp));
  double GEp = (1.0+fGep.a1*tp)/(1.0+tp*(fGep.b1+tp*(fGep.b2+fGep.b3*tp)));
  double gep = AN(xp,fGep.p1,fGep.p2,fGep.p3,fGep.p4,fGep.p5,fGep.p6,fGep.p7)*GEp;
  return gep;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gmp(const Interaction * interaction) const
{
  const Kinematics & kine   = interaction->Kine();
  double q2 = kine.q2();
  double M2 = kProtonMass*kProtonMass;
  double tp = -q2/(4*M2);
  double xp  = 2.0/(1.0+TMath::Sqrt(1.0+1.0/tp));
  double GMp = (1.0+fGmp.a1*tp)/(1.0+tp*(fGmp.b1+tp*(fGmp.b2+fGmp.b3*tp)));
  double gmp = AN(xp,fGmp.p1,fGmp.p2,fGmp.p3,fGmp.p4,fGmp.p5,fGmp.p6,fGmp.p7)*GMp;
  gmp *= fMuP;
  return gmp;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gen(const Interaction * interaction) const
{
  const Kinematics & kine   = interaction->Kine();
  double q2 = kine.q2();
  double M2 = kNeutronMass*kNeutronMass;
  double tn = -q2/(4*M2);
  double xn  = 2.0/(1.0+TMath::Sqrt(1.0+1.0/tn));
  double gep   = this->Gep(interaction);
  double gen = AN(xn,fGen.p1,fGen.p2,fGen.p3,fGen.p4,fGen.p5,fGen.p6,fGen.p7)*gep*1.7*tn/(1+3.3*tn);
  return gen;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::Gmn(const Interaction * interaction) const
{
  const Kinematics & kine   = interaction->Kine();
  double q2 = kine.q2();
  double M2 = kNeutronMass*kNeutronMass;
  double tn = -q2/(4*M2);
  double xn  = 2.0/(1.0+TMath::Sqrt(1.0+1.0/tn));
  double gmp   = this->Gmp(interaction);
  double gmn = AN(xn,fGmn.p1,fGmn.p2,fGmn.p3,fGmn.p4,fGmn.p5,fGmn.p6,fGmn.p7)*gmp;
  gmn *= fMuN/fMuP;
  return gmn;
}
//____________________________________________________________________________
void BBA07ELFormFactorsModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA07ELFormFactorsModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BBA07ELFormFactorsModel::LoadConfig(void)
{
  
  //-- load the BBA2007 fit coefficients
  GetParam( "BBA07-Gep-a1", fGep.a1) ;
  GetParam( "BBA07-Gep-b1", fGep.b1) ;
  GetParam( "BBA07-Gep-b2", fGep.b2) ;
  GetParam( "BBA07-Gep-b3", fGep.b3) ;
  GetParam( "BBA07-Gmp-a1", fGmp.a1) ;
  GetParam( "BBA07-Gmp-b1", fGmp.b1) ;
  GetParam( "BBA07-Gmp-b2", fGmp.b2) ;
  GetParam( "BBA07-Gmp-b3", fGmp.b3) ;
  GetParam( "BBA07-Gep-p1", fGep.p1) ;
  GetParam( "BBA07-Gep-p2", fGep.p2) ;
  GetParam( "BBA07-Gep-p3", fGep.p3) ;
  GetParam( "BBA07-Gep-p4", fGep.p4) ;
  GetParam( "BBA07-Gep-p5", fGep.p5) ;
  GetParam( "BBA07-Gep-p6", fGep.p6) ;
  GetParam( "BBA07-Gep-p7", fGep.p7) ;
  GetParam( "BBA07-Gen-p1", fGen.p1) ;
  GetParam( "BBA07-Gen-p2", fGen.p2) ;
  GetParam( "BBA07-Gen-p3", fGen.p3) ;
  GetParam( "BBA07-Gen-p4", fGen.p4) ;
  GetParam( "BBA07-Gen-p5", fGen.p5) ;
  GetParam( "BBA07-Gen-p6", fGen.p6) ;
  GetParam( "BBA07-Gen-p7", fGen.p7) ;
  GetParam( "BBA07-Gmp-p1", fGmp.p1) ;
  GetParam( "BBA07-Gmp-p2", fGmp.p2) ;
  GetParam( "BBA07-Gmp-p3", fGmp.p3) ;
  GetParam( "BBA07-Gmp-p4", fGmp.p4) ;
  GetParam( "BBA07-Gmp-p5", fGmp.p5) ;
  GetParam( "BBA07-Gmp-p6", fGmp.p6) ;
  GetParam( "BBA07-Gmp-p7", fGmp.p7) ;
  GetParam( "BBA07-Gmn-p1", fGmn.p1) ;
  GetParam( "BBA07-Gmn-p2", fGmn.p2) ;
  GetParam( "BBA07-Gmn-p3", fGmn.p3) ;
  GetParam( "BBA07-Gmn-p4", fGmn.p4) ;
  GetParam( "BBA07-Gmn-p5", fGmn.p5) ;
  GetParam( "BBA07-Gmn-p6", fGmn.p6) ;
  GetParam( "BBA07-Gmn-p7", fGmn.p7) ;

   //-- anomalous magnetic moments
  GetParam( "AnomMagnMoment-P", fMuP ) ;
  GetParam( "AnomMagnMoment-N", fMuN ) ;
}
//____________________________________________________________________________
double BBA07ELFormFactorsModel::AN (double x,double c1, double c2, double c3,double c4,double c5, double c6, double c7) const
{
	 const double d1  = (0.0-1.0/6)*(0.0-2.0/6)*(0.0-3.0/6)*(0.0-4.0/6)*(0.0-5.0/6)*(0.0-1.0);
     const double d2  = (1.0/6-0.0)*(1.0/6-2.0/6)*(1.0/6-3.0/6)*(1.0/6-4.0/6)*(1.0/6-5.0/6)*(1.0/6-1.0);
     const double d3  = (2.0/6-0.0)*(2.0/6-1.0/6)*(2.0/6-3.0/6)*(2.0/6-4.0/6)*(2.0/6-5.0/6)*(2.0/6-1.0);
     const double d4  = (3.0/6-0.0)*(3.0/6-1.0/6)*(3.0/6-2.0/6)*(3.0/6-4.0/6)*(3.0/6-5.0/6)*(3.0/6-1.0);
     const double d5  = (4.0/6-0.0)*(4.0/6-1.0/6)*(4.0/6-2.0/6)*(4.0/6-3.0/6)*(4.0/6-5.0/6)*(4.0/6-1.0);
     const double d6  = (5.0/6-0.0)*(5.0/6-1.0/6)*(5.0/6-2.0/6)*(5.0/6-3.0/6)*(5.0/6-4.0/6)*(5.0/6-1.0);
     const double d7  = (1.0-0.0)*(1.0-1.0/6)*(1.0-2.0/6)*(1.0-3.0/6)*(1.0-4.0/6)*(1.0-5.0/6);
     
     return c1*        (x-1.0/6)*(x-2.0/6)*(x-3.0/6)*(x-4.0/6)*(x-5.0/6)*(x-1.0)/d1+
            c2*(x-0.0)*          (x-2.0/6)*(x-3.0/6)*(x-4.0/6)*(x-5.0/6)*(x-1.0)/d2+
            c3*(x-0.0)*(x-1.0/6)*          (x-3.0/6)*(x-4.0/6)*(x-5.0/6)*(x-1.0)/d3+
            c4*(x-0.0)*(x-1.0/6)*(x-2.0/6)*          (x-4.0/6)*(x-5.0/6)*(x-1.0)/d4+
            c5*(x-0.0)*(x-1.0/6)*(x-2.0/6)*(x-3.0/6)*          (x-5.0/6)*(x-1.0)/d5+
            c6*(x-0.0)*(x-1.0/6)*(x-2.0/6)*(x-3.0/6)*(x-4.0/6)*          (x-1.0)/d6+
            c7*(x-0.0)*(x-1.0/6)*(x-2.0/6)*(x-3.0/6)*(x-4.0/6)*(x-5.0/6)          /d7;
}
//____________________________________________________________________________
