//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Hugh Gallagher <hugh.gallagher@tufts.edu>

 From code provided by:
 Igor Kakorin <idkakorin@gmail.com>
 Joint Institute for Nuclear Research, Dubna
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/QuasiElastic/XSection/MArunAxialFormFactorModel.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
MArunAxialFormFactorModel::MArunAxialFormFactorModel() :
AxialFormFactorModelI("genie::MArunAxialFormFactorModel")
{

}
//____________________________________________________________________________
MArunAxialFormFactorModel::MArunAxialFormFactorModel(string config) :
AxialFormFactorModelI("genie::MArunAxialFormFactorModel", config)
{

}
//____________________________________________________________________________
MArunAxialFormFactorModel::~MArunAxialFormFactorModel()
{

}
//____________________________________________________________________________
double MArunAxialFormFactorModel::FA(const Interaction * interaction) const
{
  const InitialState & init_state = interaction->InitState();
  const Target & target = init_state.Tgt();
  // get scattering parameters
  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  double dn;
  if (target.A()>2)
  {
  double E = init_state.ProbeE(kRfLab);
  dn = TMath::Power(1.-q2/TMath::Power(fMa*(1+fE0/E), 2), 2);
  }
  else
  dn = TMath::Power(1.-q2/fMa2, 2);
  double fa = fFA0/dn;

  LOG("MArunAxialFormFactorModel", pDEBUG) << "FA(q^2 = " << q2 << ") = " << fa;
  return fa;
}
//____________________________________________________________________________
void MArunAxialFormFactorModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void MArunAxialFormFactorModel::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void MArunAxialFormFactorModel::LoadConfig(void)
{
  // axial mass
	GetParam( "QEL-Ma", fMa ) ;
	fMa2 = TMath::Power(fMa,2);

  // E0 for calculating running axial mass: Ma*(1+E0/Enu)
	GetParam( "QEL-E0", fE0 ) ;

	// FA(q2 = 0)
	GetParam( "QEL-FA0", fFA0 ) ;

}
//____________________________________________________________________________
