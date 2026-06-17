//____________________________________________________________________________
/*
RSHelicityAmplModelDMI: the interface for the Dark Matter helicity Amplitudes

Zachary W. Orr, Colorado State University
*/
//____________________________________________________________________________
#include "Physics/BoostedDarkMatter/XSection/RSHelicityAmplModelDMI.h" //Specific to DM

using namespace genie;

//____________________________________________________________________________
RSHelicityAmplModelDMI::RSHelicityAmplModelDMI() :
Algorithm()
{

}
//____________________________________________________________________________
RSHelicityAmplModelDMI::RSHelicityAmplModelDMI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
RSHelicityAmplModelDMI::RSHelicityAmplModelDMI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelDMI::~RSHelicityAmplModelDMI()
{

}
//____________________________________________________________________________

//___________added for quark charges
//____________________________________________________________________________
void RSHelicityAmplModelDMI::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSHelicityAmplModelDMI::LoadConfig(void)
{

  // quark couplings to mediator
  double QuL, QuR, QdL, QdR;
  this->GetParam( "UpLeftCharge", QuL ) ;
  this->GetParam( "UpRightCharge", QuR ) ;
  this->GetParam( "DownLeftCharge", QdL ) ;
  this->GetParam( "DownRightCharge", QdR ) ;
  fQuV = 0.5*(QuL + QuR);
  fQuA = 0.5*(- QuL + QuR);
  fQdV = 0.5*(QdL + QdR);
  fQdA = 0.5*(- QdL + QdR);
}
//_____________
