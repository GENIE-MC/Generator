//____________________________________________________________________________
/*!

\class    genie::ScatteringParams

\brief    A Registry subclass to hold scattering parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 08, 2004

*/ 
//____________________________________________________________________________

#include "Interaction/ScatteringParams.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
ScatteringParams::ScatteringParams() :
Registry()
{
  SetName("ScatteringParams");
}
//____________________________________________________________________________
ScatteringParams::ScatteringParams(const Registry & params) :
Registry(params)
{
  SetName("ScatteringParams");
}
//____________________________________________________________________________
ScatteringParams::~ScatteringParams()
{

}
//____________________________________________________________________________
double ScatteringParams::x(void) const
{
  double x = 0;
  
  if( Exists("x") ) Get("x",x);
  else LOG("Interaction", pWARN) << "scattering parameter x does not exist";

  return x;
}
//____________________________________________________________________________
double ScatteringParams::y(void) const
{
  double y = 0;

  if( Exists("y") ) Get("y",y);
  else LOG("Interaction", pWARN) << "scattering parameter y does not exist";

  return y;
}
//____________________________________________________________________________
double ScatteringParams::Q2(void) const
{
  if( Exists("Q2") ) { 

     double Q2 = 0;  
     Get("Q2",Q2);  
     return  Q2; 

  } else if( Exists("q2") ) { 

     double q2 = 0;  
     Get("q2",q2);  
     return -q2; 

  } else
       LOG("Interaction", pWARN) << "scattering parameter Q^2 does not exist";

  return 0;
}
//____________________________________________________________________________
double ScatteringParams::q2(void) const
{
  if( Exists("q2") ) { 

     double q2 = 0;  
     Get("q2",q2);  
     return  q2; 

  } else if( Exists("Q2") ) { 

     double Q2 = 0;  
     Get("Q2",Q2);  
     return -Q2; 

  } else
       LOG("Interaction", pWARN) << "scattering parameter q^2 does not exist";

  return 0;
}
//____________________________________________________________________________
double ScatteringParams::W(void) const
{
  double W = 0;

  if( Exists("W") ) Get("W",W);
  else LOG("Interaction", pWARN) << "scattering parameter W does not exist";

  return W;
}
//____________________________________________________________________________
double ScatteringParams::lnQ2(void) const
{
  if( Exists("Q2") ) {

    double Q2 = 0;   
    Get("Q2",Q2);
    return TMath::Log(Q2);

  } else 
       LOG("Interaction", pWARN) << "scattering parameter Q^2 does not exist";

  return 0;
}
//____________________________________________________________________________
double ScatteringParams::lnW(void) const
{
  if( Exists("W") ) {

    double W = 0;
    Get("W",W);
    return TMath::Log(W);

  } else 
       LOG("Interaction", pWARN) << "scattering parameter Q^2 does not exist";

  return 0;
}
//____________________________________________________________________________




