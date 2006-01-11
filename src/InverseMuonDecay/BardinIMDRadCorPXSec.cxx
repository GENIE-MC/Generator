//____________________________________________________________________________
/*!

\class    genie::BardinIMDRadCorPXSec

\brief    Computes the Inverse Muon Decay (IMD) diff. cross section, dxsec/dy,
          where y is the interaction inelasticity, using the Bardin -
          Dokuchaeva model which includes all 1-loop radiative corrections. \n

          This is a 'trully' inclusive IMD cross section, i.e. the brem. cross
          section (dxsec_brem/dy)|w>w0 [see Bardin paper, cited below] is not
          subtracted from the IMD cross section and therefore it is not suitable
          for experimental situations where a photon energy trigger threshold
          is applied.

          BardinIMDRadCorPXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      D.Yu.Bardin and V.A.Dokuchaeva, Nucl.Phys.B287:839 (1987)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "InverseMuonDecay/BardinIMDRadCorPXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BardinIMDRadCorPXSec::BardinIMDRadCorPXSec() :
XSecAlgorithmI("genie::BardinIMDRadCorPXSec")
{

}
//____________________________________________________________________________
BardinIMDRadCorPXSec::BardinIMDRadCorPXSec(string config) :
XSecAlgorithmI("genie::BardinIMDRadCorPXSec", config)
{

}
//____________________________________________________________________________
BardinIMDRadCorPXSec::~BardinIMDRadCorPXSec()
{

}
//____________________________________________________________________________
double BardinIMDRadCorPXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> GetInitialState();

  double E    = init_state.GetProbeE(kRfLab);
  double sig0 = kGF_2 * kElectronMass * E / kPi;
  double re   = 0.5 * kElectronMass / E;
  double r    = (kMuonMass_2 / kElectronMass_2) * re;
  double y    = interaction->GetKinematics().y();

  y = 1-y;  //Note: y = (Ev-El)/Ev but in Bardin's paper y=El/Ev.

  double ymin = r + re;
  double ymax = 1 + re + r*re / (1+re);

  ymax = TMath::Min(ymax,0.9999999); // avoid ymax=1, due to a log(1-y)

  LOG("InverseMuDecay", pDEBUG)
                << "sig0 = " << sig0 << ", r = " << r << ", re = " << re;
  LOG("InverseMuDecay", pDEBUG)
                        << "allowed y: [" << ymin << ", " << ymax << "]";

  if(y<ymin || y>ymax) return 0;

  double dsig_dy = 2 * sig0 * ( 1 - r + (kAem/kPi) * Fa(re,r,y) );

  LOG("InverseMuDecay", pINFO)
     << "dxsec[1-loop]/dy (Ev = " << E << ", y = " << y << ") = " << dsig_dy;

  return dsig_dy;
}
//____________________________________________________________________________
bool BardinIMDRadCorPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return true;
}
//____________________________________________________________________________
bool BardinIMDRadCorPXSec::ValidKinematics(
                                        const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();
  double E = init_state.GetProbeE(kRfLab);
  double s = kElectronMass_2 + 2*kElectronMass*E;

  //-- check if it is kinematically allowed
  if(s < kMuonMass_2) {
     LOG("InverseMuDecay", pINFO)
        << "Ev = " << E << " (s = " << s << ") is below threshold (s-min = "
        << kMuonMass_2 << ") for IMD";
     return false;
  }
  return true;
}
//____________________________________________________________________________
double BardinIMDRadCorPXSec::Fa(double re, double r, double y) const
{
  double y2  = y * y;
  double rre = r * re;
  double r_y = r/y;
  double y_r = y/r;

  double fa = 0;

  fa = (1-r) *       ( TMath::Log(y2/rre) * TMath::Log(1-r_y) +
                       TMath::Log(y_r) * TMath::Log(1-y) -
                       this->Li2( r ) +
                       this->Li2( y ) +
                       this->Li2( (r-y) / (1-y) ) +
                       1.5 * (1-r) * TMath::Log(1-r)
                     )
                     +

       0.5*(1+3*r) * ( this->Li2( (1-r_y) / (1-r) ) -
                       this->Li2( (y-r)   / (1-r) ) -
                       TMath::Log(y_r) * TMath::Log( (y-r) / (1-r) )
                     )
                     +

       this->P(1,r,y)                   -
       this->P(2,r,y) * TMath::Log(r)   -
       this->P(3,r,y) * TMath::Log(re)  +
       this->P(4,r,y) * TMath::Log(y)   +
       this->P(5,r,y) * TMath::Log(1-y) +
       this->P(6,r,y) * (1 - r_y) * TMath::Log(1-r_y);

  return fa;
}
//____________________________________________________________________________
double BardinIMDRadCorPXSec::P(int i, double r, double y) const
{
  int kmin = -3;
  int kmax =  2;

  double p = 0;

  for(int k = kmin; k <= kmax; k++) {

     double c  = this->C(i,k,r);
     double yk = TMath::Power(y,k);

     p += (c*yk);
  }
  return p;
}
//____________________________________________________________________________
double BardinIMDRadCorPXSec::Li2(double z) const
{
  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * alg_base = algf->GetAlgorithm("genie::Simpson1D");

  const IntegratorI * integrator =
                              dynamic_cast<const IntegratorI *> (alg_base);

  const int    nsteps  = 101;
  const double epsilon = 1e-2;
  const double min     = epsilon;
  const double max     = 1.0 - epsilon;
  const double step    = (max-min)/(nsteps-1);

  UnifGrid grid;

  grid.AddDimension(nsteps, min, max);

  FunctionMap fmap(grid);

  for(int i = 0; i < nsteps; i++) {

    double t  = min + i * step;
    double f  = TMath::Log(1-z*t)/t;

    fmap.AddPoint(f, i);
  }

  double li2 = integrator->Integrate(fmap);

  return li2;
}
//____________________________________________________________________________
double BardinIMDRadCorPXSec::C(int i, int k, double r) const
{
  if        ( i == 1 ) {

      if      (k == -3) return -0.19444444*pow(r,3.);
      else if (k == -2) return (0.083333333+0.29166667*r)*pow(r,2.);
      else if (k == -1) return -0.58333333*r - 0.5*pow(r,2.) - pow(r,3.)/6.;
      else if (k ==  0) return -1.30555560 + 3.125*r + 0.375*pow(r,2.);
      else if (k ==  1) return -0.91666667 - 0.25*r;
      else if (k ==  2) return 0.041666667;
      else              return 0.;

  } else if ( i == 2 ) {

      if      (k == -3) return 0.;
      else if (k == -2) return 0.5*pow(r,2.);
      else if (k == -1) return 0.5*r - 2*pow(r,2.);
      else if (k ==  0) return 0.25 - 0.75*r + 1.5*pow(r,2);
      else if (k ==  1) return 0.5;
      else if (k ==  2) return 0.;
      else              return 0.;

  } else if ( i == 3 ) {

      if      (k == -3) return 0.16666667*pow(r,3.);
      else if (k == -2) return 0.25*pow(r,2.)*(1-r);
      else if (k == -1) return r-0.5*pow(r,2.);
      else if (k ==  0) return 0.66666667;
      else if (k ==  1) return 0.;
      else if (k ==  2) return 0.;
      else              return 0.;

  } else if ( i == 4 ) {

      if      (k == -3) return 0.;
      else if (k == -2) return pow(r,2.);
      else if (k == -1) return r*(1-4.*r);
      else if (k ==  0) return 1.5*pow(r,2.);
      else if (k ==  1) return 1.;
      else if (k ==  2) return 0.;
      else              return 0.;

  } else if ( i == 5 ) {

      if      (k == -3) return 0.16666667*pow(r,3.);
      else if (k == -2) return -0.25*pow(r,2.)*(1+r);
      else if (k == -1) return 0.5*r*(1+3*r);
      else if (k ==  0) return -1.9166667+2.25*r-1.5*pow(r,2);
      else if (k ==  1) return -0.5;
      else if (k ==  2) return 0.;
      else              return 0.;

  } else if ( i == 6 ) {

      if      (k == -3) return 0.;
      else if (k == -2) return 0.16666667*pow(r,2.);
      else if (k == -1) return -0.25*r*(r+0.33333333);
      else if (k ==  0) return 1.25*(r+0.33333333);
      else if (k ==  1) return 0.5;
      else if (k ==  2) return 0.;
      else              return 0.;

  } else return 0.;
}
//____________________________________________________________________________
