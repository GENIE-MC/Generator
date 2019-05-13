//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <sstream>

#include <TSystem.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/SpectralFunc1d.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
SpectralFunc1d::SpectralFunc1d() :
NuclearModelI("genie::SpectralFunc1d")
{

}
//____________________________________________________________________________
SpectralFunc1d::SpectralFunc1d(string config) :
NuclearModelI("genie::SpectralFunc1d", config)
{

}
//____________________________________________________________________________
SpectralFunc1d::~SpectralFunc1d()
{
  this->CleanUp();
}
//____________________________________________________________________________
bool SpectralFunc1d::GenerateNucleon(const Target & target) const
{
  RandomGen * rnd = RandomGen::Instance();
  int Z = target.Z();

  map<int, double>::const_iterator  dbl_it;
  map<int, Spline*>::const_iterator spl_it;

  // Select fermi momentum from the integrated (over removal energies) s/f.
  //
  spl_it = fSFk.find(Z);
  dbl_it = fMaxProb.find(Z);
  if(spl_it == fSFk.end() || dbl_it == fMaxProb.end()) {
    fCurrRemovalEnergy = 0.;
    fCurrMomentum.SetXYZ(0.,0.,0.);
    return false;
  }

  double prob_max = dbl_it->second;
  LOG("SpectralFunc1", pDEBUG) << "Max probability = " << prob_max;

  double p = 0;
  unsigned int niter = 0;
  while(1) {
    if(niter > kRjMaxIterations) {
       LOG("SpectralFunc1", pWARN)
        << "Couldn't generate a hit nucleon after " << niter << " iterations";
       return false;
    }
    niter++;

    if(fUseRFGMomentumCutoff) p = fPCutOff * rnd->RndGen().Rndm();
    else p = rnd->RndGen().Rndm();

    double prob  = spl_it->second->Evaluate(p);
    double probg = prob_max * rnd->RndGen().Rndm();
    LOG("SpectralFunc1", pDEBUG) << "Trying p = " << p << " / prob = " << prob;

    bool accept = (probg<prob);
    if(accept) break;
  }

  LOG("SpectralFunc1", pINFO) << "|p,nucleon| = " << p;

  double costheta = -1. + 2. * rnd->RndGen().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->RndGen().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  double px = p*sintheta*cosfi;
  double py = p*sintheta*sinfi;
  double pz = p*costheta;

  fCurrMomentum.SetXYZ(px,py,pz);

  // Set removal energy
  // Do it either in the same way as in the FG model or by using the average
  // removal energy for the seleced pF as calculated from the s/f itself
  //
  if(fUseRFGRemovalE) {
    dbl_it = fNucRmvE.find(Z);
    if(dbl_it != fNucRmvE.end()) fCurrRemovalEnergy = dbl_it->second;
    else fCurrRemovalEnergy = nuclear::BindEnergyPerNucleon(target);
  } else {
    spl_it = fSFw.find(Z);
    if(spl_it==fSFw.end()) {
       fCurrRemovalEnergy = 0.;
       fCurrMomentum.SetXYZ(0.,0.,0.);
       return false;
    } else fCurrRemovalEnergy = spl_it->second->Evaluate(p);
  }

  return true;
}
//____________________________________________________________________________
double SpectralFunc1d::Prob(
                  double p, double /*w*/, const Target & target) const
{
  if(fUseRFGMomentumCutoff) { if(p > fPCutOff) return 0; }

  int Z = target.Z();
  map<int, Spline*>::const_iterator spl_it = fSFk.find(Z);

  if(spl_it == fSFk.end()) return 0;
  else return TMath::Max(0.,spl_it->second->Evaluate(p));
}
//____________________________________________________________________________
void SpectralFunc1d::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc1d::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void SpectralFunc1d::LoadConfig(void)
{
  LOG("SpectralFunc1", pDEBUG) << "Loading coonfiguration for SpectralFunc1d";

  this->CleanUp();

  // Load spectral function data.
  // Hopefully analytical expressions will be available soon.
  // Currently I have spectral functions for C12 and Fe56 only.
  //
  string data_dir =
      string(gSystem->Getenv("GENIE")) + string("/data/evgen/nucl/spectral_functions/");

  string c12_sf1dk_file  = data_dir + "benhar-sf1dk-12c.data";
  string fe56_sf1dk_file = data_dir + "benhar-sf1dk-56fe.data";
  string c12_sf1dw_file  = data_dir + "benhar-sf1dw-12c.data";
  string fe56_sf1dw_file = data_dir + "benhar-sf1dw-56fe.data";

  Spline * spl = 0;

  spl = new Spline(c12_sf1dk_file);
  fSFk.insert(map<int, Spline*>::value_type(6,spl));
  spl = new Spline(fe56_sf1dk_file);
  fSFk.insert(map<int, Spline*>::value_type(26,spl));

  spl = new Spline(c12_sf1dw_file);
  fSFw.insert(map<int, Spline*>::value_type(6,spl));
  spl = new Spline(fe56_sf1dw_file);
  fSFw.insert(map<int, Spline*>::value_type(26,spl));

  // scan for max.
  map<int, Spline*>::const_iterator spliter;
  int    n    = 100;
  double pmin =  0.000;
  double dp   =  0.010;
  for(spliter = fSFk.begin(); spliter != fSFk.end(); ++spliter) {
    double prob_max = 0;
    int Z = spliter->first;
    spl = spliter->second;
    for(int i=0; i<n; i++) {
       double p = pmin + i*dp;
       prob_max = TMath::Max(prob_max, spl->Evaluate(p));
    }
    fMaxProb.insert(map<int,double>::value_type(Z,prob_max));
  }

  // Check whether to use the same removal energies as in the FG model or
  // to use the average removal energy for the selected fermi momentum
  // (computed from the spectral function itself)
  GetParam( "UseRFGRemovalE", fUseRFGRemovalE ) ;

  // Check whether to use the same momentum cutoff as in the FG model
  GetParam("UseRFGMomentumCutoff", fUseRFGMomentumCutoff ) ;

  //Get the momentum cutoff
  GetParam( "RFG-MomentumCutOff", fPCutOff ) ;

  // Removal energies as used in the FG model
  // Load removal energy for specific nuclei from either the algorithm's
  // configuration file or the UserPhysicsOptions file.
  // If none is used use Wapstra's semi-empirical formula.
  for(int Z=1; Z<140; Z++) {
    for(int A=Z; A<3*Z; A++) {
      ostringstream key;
      int pdgc = pdg::IonPdgCode(A,Z);
      key   << "RFG-NucRemovalE@Pdg="     << pdgc;
      RgKey rgkey   = key.str();
      double eb ;
      if ( GetParam( rgkey, eb, false ) ) {
    	eb = TMath::Max(eb, 0.);
        LOG("BodekRitchie", pNOTICE)
          << "Nucleus: " << pdgc << " -> using Eb =  " << eb << " GeV";
        fNucRmvE.insert(map<int,double>::value_type(Z,eb));
      }
    }
  }
}
//____________________________________________________________________________
void SpectralFunc1d::CleanUp(void)
{
  map<int, Spline*>::iterator spliter;

  for(spliter = fSFk.begin(); spliter != fSFk.end(); ++spliter) {
    Spline * spl = spliter->second;
    if(spl) delete spl;
  }
  for(spliter = fSFw.begin(); spliter != fSFw.end(); ++spliter) {
    Spline * spl = spliter->second;
    if(spl) delete spl;
  }
  fSFk.clear();
  fSFw.clear();
  fNucRmvE.clear();
  fMaxProb.clear();
}
//____________________________________________________________________________

