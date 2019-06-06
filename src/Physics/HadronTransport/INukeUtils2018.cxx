//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Jim Dobson <j.dobson07 \at imperial.ac.uk>
         Imperial College London

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Aaron Meyer <asm58 \at pitt.edu>
         Pittsburgh University

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 04, 2009 - JD
   Was first added in v2.5.1. Adapted from the T2K GENIE reweighting tool.
 @ Mar 05, 2009 - CA
   Modified ReconstructHadronFateHA() to work with hadron+A event files in
   addition to neutrino event files.
 @ Sep 10, 2009 - CA
   Added MeanFreePath(), Dist2Exit(), Dist2ExitMFP()
 @ Sep 30, 2009 - CA
   Added StepParticle() from Intranuke.cxx
 @ Oct 02, 2009 - CA
   Added test MeanFreePath_Delta().
 @ Jul 15, 2010 - AM
   Added common utility functions used by both hA and hN mode. Updated
   MeanFreePath to separate proton and neutron cross sections. Added general
   utility functions.
 @ Jan 9, 2015 - SD, NG, TG
   Added 2014 version of INTRANUKE codes for v2.9.0.  Uses INukeHadroData2014,
   but no changes to mean free path.
*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/HadronTransport/INukeException.h"
#include "Physics/HadronTransport/INukeUtils2018.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Registry/Registry.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/HadronTransport/INukeOset.h"
#include "Physics/HadronTransport/INukeOsetTable.h"
#include "Physics/HadronTransport/INukeOsetFormula.h"
#include "TComplex.h"

using std::ostringstream;
using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
double genie::utils::intranuke2018::MeanFreePath(
   int pdgc, const TLorentzVector & x4, const TLorentzVector & p4,
   double A, double Z, double nRpi, double nRnuc, const bool useOset, const bool altOset, const bool xsecNNCorr, string INukeMode)
{
// Calculate the mean free path (in fm) for a pions and nucleons in a nucleus
//
// Inputs
//  pdgc : Hadron PDG code
//  x4   : Hadron 4-position in the nucleus coordinate system (units: fm)
//  p4   : Hadron 4-momentum (units: GeV)
//  A    : Nucleus atomic mass number
//  nRpi : Controls the pion ring size in terms of de-Broglie wavelengths
//  nRnuc: Controls the nuclepn ring size in terms of de-Broglie wavelengths
//
  bool is_pion    = pdgc == kPdgPiP || pdgc == kPdgPi0 || pdgc == kPdgPiM;
  bool is_nucleon = pdgc == kPdgProton || pdgc == kPdgNeutron;
  bool is_kaon    = pdgc == kPdgKP;
  bool is_gamma   = pdgc == kPdgGamma;

  if(!is_pion && !is_nucleon && !is_kaon && !is_gamma) return 0.;

  // before getting the nuclear density at the current position
  // check whether the nucleus has to become larger by const times the
  // de Broglie wavelength -- that is somewhat empirical, but this
  // is what is needed to get piA total cross sections right.
  // The ring size is different for light nuclei (using gaus density) /
  // heavy nuclei (using woods-saxon density).
  // The ring size is different for pions / nucleons.
  //
  double momentum = p4.Vect().Mag(); // hadron momentum in GeV
  double ring = (momentum>0) ? 1.240/momentum : 0; // de-Broglie wavelength

  if(A<=20) { ring /= 2.; }

  /*
  if      (is_pion               ) { ring *= nRpi;  }
  else if (is_nucleon            ) { ring *= nRnuc; }
  else if (is_gamma || is_kaon || useOset) { ring = 0.;     }
  */
  if(INukeMode=="hN2018")
    {
      if      (is_pion               ) { ring *= nRpi;  }
      else if (is_nucleon            ) { ring *= nRnuc; }
      else if (is_gamma || is_kaon || useOset) { ring = 0.;}
    }
  else
    {
      if      (is_pion    || is_kaon ) { ring *= nRpi;  }
      else if (is_nucleon            ) { ring *= nRnuc; }
      else if (is_gamma              ) { ring = 0.;     }
    }

  // get the nuclear density at the current position
  double rnow = x4.Vect().Mag();
  double rho  = A * utils::nuclear::Density(rnow,(int) A,ring);

  // the hadron+nucleon cross section will be evaluated within the range
  // of the input spline and assumed to be const outside that range
  //
  double ke = (p4.Energy() - p4.M()) / units::MeV;  // kinetic energy in MeV
  ke = TMath::Max(INukeHadroData2018::fMinKinEnergy,   ke);
  ke = TMath::Min(INukeHadroData2018::fMaxKinEnergyHN, ke);

  // get total xsection for the incident hadron at its current
  // kinetic energy
  double sigtot = 0;
  double ppcnt = (double) Z/ (double) A; // % of protons remaining
  INukeHadroData2018 * fHadroData2018 = INukeHadroData2018::Instance();

  if (is_pion and (INukeMode == "hN2018") and useOset and ke < 350.0)
    sigtot = sigmaTotalOset (ke, rho, pdgc, ppcnt, altOset);
  else if (pdgc == kPdgPiP)
    { sigtot = fHadroData2018 -> XSecPipp_Tot() -> Evaluate(ke)*ppcnt;
      sigtot+= fHadroData2018 -> XSecPipn_Tot() -> Evaluate(ke)*(1-ppcnt);}
  else if (pdgc == kPdgPi0)
    { sigtot = fHadroData2018 -> XSecPi0p_Tot() -> Evaluate(ke)*ppcnt;
      sigtot+= fHadroData2018 -> XSecPi0n_Tot() -> Evaluate(ke)*(1-ppcnt);}
  else if (pdgc == kPdgPiM)
    { sigtot = fHadroData2018 -> XSecPipn_Tot() -> Evaluate(ke)*ppcnt;
      sigtot+= fHadroData2018 -> XSecPipp_Tot() -> Evaluate(ke)*(1-ppcnt);}
  else if (pdgc == kPdgProton)
    {
      sigtot = fHadroData2018 -> XSecPp_Tot()   -> Evaluate(ke)*ppcnt;
      //sigtot+= fHadroData2018 -> XSecPn_Tot()   -> Evaluate(ke)*(1-ppcnt);

      PDGLibrary * pLib = PDGLibrary::Instance();
      double hc = 197.327;
      double R0 = 1.25 * TMath::Power(A,1./3.) + 2.0 * 0.65; // should all be in units of fm
      double Mp = pLib->Find(2212)->Mass();
      double M  = pLib->Find(pdgc)->Mass();
      //double E  = (p4.Energy() - Mp) * 1000.; // Convert GeV to MeV.
      double E = ke;
      if (Z*hc/137./x4.Vect().Mag() > E)  // Coulomb correction (Cohen, Concepts of Nuclear Physics, pg. 259-260)
        {
          double z  = 1.0; // charge for single proton
          double Bc = z*Z*hc/137./R0;
          double x  = E/Bc;
          double f  = TMath::ACos(TMath::Power(x,0.5)) - TMath::Power(x*(1-x),0.5);
          double B  = 0.63*z*Z*TMath::Power((M/Mp)/E,0.5);
          double Pc = TMath::Exp(-B*f);
          sigtot *= Pc;
        }
      sigtot+= fHadroData2018 -> XSecPn_Tot()   -> Evaluate(ke)*(1-ppcnt);

      double E0 = TMath::Power(A,0.2)*12.;
      if (INukeMode=="hN2018"){if(ke<E0){sigtot=0.0;}}  //empirical - needed to cut off large number of low energy nucleons
      //      LOG("INukeUtils",pDEBUG) "sigtot for proton= " << sigtot << "; KE= " << ke << "; E0= " << E0;
    }
  else if (pdgc == kPdgNeutron)
    {
      sigtot = fHadroData2018 -> XSecPn_Tot()   -> Evaluate(ke)*ppcnt;
      sigtot+= fHadroData2018 -> XSecNn_Tot()   -> Evaluate(ke)*(1-ppcnt);
      double E0 = TMath::Power(A,0.2)*12.;
      if (INukeMode=="hN2018"){if(ke<E0){sigtot=0.0;}}  //empirical - needed to cut off large number of low energy nucleons
      //      LOG("INukeUtils",pDEBUG) "sigtot for neutron= " << sigtot << "; KE= " << ke;
    }
  else if (pdgc == kPdgKP)
    { sigtot = fHadroData2018 -> XSecKpN_Tot()  -> Evaluate(ke);
      // this factor is used to empirically get agreement with tot xs data, justified historically.
      sigtot*=1.1;}
  else if (pdgc == kPdgGamma)
    { sigtot = fHadroData2018 -> XSecGamp_fs()  -> Evaluate(ke)*ppcnt;
      sigtot+= fHadroData2018 -> XSecGamn_fs()  -> Evaluate(ke)*(1-ppcnt);}
  else {
     return 0;
  }

  // the xsection splines in INukeHadroData return the hadron x-section in
  // mb -> convert to fm^2
  sigtot *= (units::mb / units::fm2);

  // avoid defective error handling
  //if(sigtot<1E-6){sigtot=1E-6;}

  if (xsecNNCorr and is_nucleon)
    sigtot *= INukeNucleonCorr::getInstance()->
      getAvgCorrection (rho, A, p4.E() - PDGLibrary::Instance()->Find(pdgc)->Mass());   //uses lookup tables

  // avoid defective error handling
  if(sigtot<1E-6){sigtot=1E-6;}

  // compute the mean free path
  double lamda = 1. / (rho * sigtot);

  // exits if lamda is InF (if cross section is 0)
  if( ! TMath::Finite(lamda) ) {
     return -1;
  }

/*
  LOG("INukeUtils", pDEBUG)
     << "sig_total = " << sigtot << " fm^2, rho = " << rho
     << " fm^-3  => mfp = " << lamda << " fm.";
*/
  return lamda;
}
//____________________________________________________________________________
double genie::utils::intranuke2018::MeanFreePath_Delta(
   int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A)
{
//
// **test**
//

// Calculate the mean free path (in fm) for Delta's in a nucleus
//
// Inputs
//  pdgc : Hadron PDG code
//  x4   : Hadron 4-position in the nucleus coordinate system (units: fm)
//  p4   : Hadron 4-momentum (units: GeV)
//  A    : Nucleus atomic mass number
//
  bool is_deltapp = (pdgc==kPdgP33m1232_DeltaPP);
  if(!is_deltapp) return 0.;

  // get the nuclear density at the current position
  double rnow = x4.Vect().Mag();
  double rho  = A * utils::nuclear::Density(rnow,(int) A);

  // the Delta+N->N+N cross section will be evaluated within the range
  // of the input spline and assumed to be const outside that range
  double ke = (p4.Energy() - p4.M()) / units::MeV;  // kinetic energy in MeV
  ke = TMath::Max(1500., ke);
  ke = TMath::Min(   0., ke);

  // get the Delta+N->N+N cross section
  double sig = 0;
  if      (ke< 500) sig=20;
  else if (ke<1000) sig=40;
  else              sig=30;

  // value is in mb -> convert to fm^2
  sig *= (units::mb / units::fm2);

  // compute the mean free path
  double lamda = 1. / (rho * sig);

  // exits if lamda is InF (if cross section is 0)
  if( ! TMath::Finite(lamda) ) {
     return -1;
  }

  return lamda;
}
//____________________________________________________________________________
double genie::utils::intranuke2018::ProbSurvival(
  int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A, double Z,
  double mfp_scale_factor, double nRpi, double nRnuc, double NR, double R0)
{
// Calculate the survival probability for a hadron inside a nucleus
//
// Inputs
//  pdgc : Hadron PDG code
//  x4 : Hadron 4-position in the nucleus coordinate system (units: fm)
//  p4 : Hadron 4-momentum (units: GeV)
//  A : Nucleus atomic mass number
//  mfp_scale_factor: Tweaks the mean free path (mfp -> mfp*scale). Def: 1.0
//  nRpi: Controls the pion ring size in terms of de-Broglie wavelengths
//  nRnuc: Controls the nuclepn ring size in terms of de-Broglie wavelengths
//  NR: How far away to track the hadron, in terms of the corresponding
//      nuclear radius. Def: 3
//  R0: R0 in R=R0*A^1/3 (units:fm). Def. 1.4

   double prob = 1.0;

   double step = 0.05; // fermi
   double R    = NR * R0 * TMath::Power(A, 1./3.);

   TVector3 dr3 = p4.Vect().Unit();  // unit vector along its direction
   TLorentzVector dr4(dr3,0);

   LOG("INukeUtils", pDEBUG)
     << "Calculating survival probability for hadron with PDG code = " << pdgc
     << " and momentum = " << p4.P() << " GeV";
   LOG("INukeUtils", pDEBUG)
     << "mfp scale = " << mfp_scale_factor
     << ", nRpi = " << nRpi << ", nRnuc = " << nRnuc << ", NR = " << NR
     << ", R0 = " << R0 << " fm";

   TLorentzVector x4_curr(x4); // current position

   while(1) {
     double rnow  = x4_curr.Vect().Mag();
     if (rnow > (R+step)) break;

     x4_curr += (step*dr4);
     rnow = x4_curr.Vect().Mag();
     double mfp =
       genie::utils::intranuke2018::MeanFreePath(pdgc,x4_curr,p4,A,Z,nRpi,nRnuc);
     double mfp_twk = mfp * mfp_scale_factor;

     double dprob = (mfp_twk>0) ? TMath::Exp(-step/mfp_twk) : 0.;
     prob*=dprob;
/*
     LOG("INukeUtils", pDEBUG)
       << "+ step size = " << step << " fm, |r| = " << rnow << " fm, "
       << "mfp = " << mfp_twk << "fm (nominal mfp = " << mfp << " fm): "
       << "dPsurv = " << dprob << ", current Psurv = " << prob;
*/
   }

   LOG("INukeUtils", pDEBUG) << "Psurv = " << prob;

   return prob;
}
//____________________________________________________________________________
double genie::utils::intranuke2018::Dist2Exit(
   const TLorentzVector & x4, const TLorentzVector & p4,
   double A, double NR, double R0)
{
// Calculate distance within a nucleus (units: fm) before we stop tracking
// the hadron.
// See previous functions for a description of inputs.
//
   double R    = NR * R0 * TMath::Power(A, 1./3.);
   double step = 0.05; // fermi

   TVector3 dr3 = p4.Vect().Unit();  // unit vector along its direction
   TLorentzVector dr4(dr3,0);

   TLorentzVector x4_curr(x4); // current position

   double d=0;
   while(1) {
        double rnow  = x4_curr.Vect().Mag();
        x4_curr += (step*dr4);
        d+=step;
        rnow = x4_curr.Vect().Mag();
        if (rnow > R) break;
   }
   return d;
}
//____________________________________________________________________________
double genie::utils::intranuke2018::Dist2ExitMFP(
   int pdgc, const TLorentzVector & x4, const TLorentzVector & p4,
   double A, double Z, double NR, double R0)
{
// Calculate distance within a nucleus (expressed in terms of 'mean free
// paths') before we stop tracking the hadron.
// See previous functions for a description of inputs.
//

// distance before exiting in mean free path lengths
//
   double R    = NR * R0 * TMath::Power(A, 1./3.);
   double step = 0.05; // fermi

   TVector3 dr3 = p4.Vect().Unit();  // unit vector along its direction
   TLorentzVector dr4(dr3,0);

   TLorentzVector x4_curr(x4); // current position

   double d=0;
   double d_mfp=0;
   while(1) {
        double rnow  = x4_curr.Vect().Mag();
        x4_curr += (step*dr4);
        d+=step;
        rnow = x4_curr.Vect().Mag();

        double lambda = genie::utils::intranuke2018::MeanFreePath(pdgc,x4_curr,p4,A,Z);
        d_mfp += (step/lambda);

        if (rnow > R) break;
   }
   return d_mfp;
}
//____________________________________________________________________________
void genie::utils::intranuke2018::StepParticle(
           GHepParticle * p, double step, double nuclear_radius)
{
// Steps a particle starting from its current position (in fm) and moving
// along the direction of its current momentum by the input step (in fm).
// The particle is stepped in a straight line.
// If a nuclear radius is set then the following check is performed:
// If the step is too large and takes the the particle far away from the
// nucleus then its position is scaled back so that the escaped particles are
// always within a ~1fm from the "outer nucleus surface"

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils", pDEBUG)
      << "Stepping particle [" << p->Name() << "] by dr = " << step << " fm";
#endif

  // Step particle
  TVector3 dr = p->P4()->Vect().Unit();      // unit vector along its direction
  dr.SetMag(step);                           // spatial step size
  double dt = 0;                             // temporal step:
  TLorentzVector dx4(dr,dt);                 // 4-vector step
  TLorentzVector x4new = *(p->X4()) + dx4;   // new position

  if(nuclear_radius > 0.) {
    // Check position against nuclear boundary. If the particle was stepped
    // too far away outside the nuclear boundary bring it back to within
    // 1fm from that boundary
    double epsilon = 1; // fm
    double r       = x4new.Vect().Mag(); // fm
    double rmax    = nuclear_radius+epsilon;
    if(r > rmax) {
       LOG("INukeUtils", pINFO)
         << "Particle was stepped too far away (r = " << r << " fm)";
       LOG("INukeUtils", pINFO)
         << "Placing it " << epsilon
         << " fm outside the nucleus (r' = " << rmax << " fm)";
       double scale = rmax/r;
       x4new *= scale;
    }//r>rmax
  }//nucl radius set

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils", pDEBUG)
      << "\n Init direction = " << print::Vec3AsString(&dr)
      << "\n Init position (in fm,nsec) = " << print::X4AsString(p->X4())
      << "\n Fin  position (in fm,nsec) = " << print::X4AsString(&x4new);
#endif

  p->SetPosition(x4new);
}


//___________________________________________________________________________
// Method to handle compound nucleus considerations, preequilibrium
//    and equilibrium
//    Alex Bell -- 6/17/2008
void genie::utils::intranuke2018::PreEquilibrium(
  GHepRecord * ev, GHepParticle * p,
  int &RemnA, int &RemnZ, TLorentzVector &RemnP4,
  bool /* DoFermi */, double /* FermiFac */,
  const NuclearModelI* /* Nuclmodel */, double NucRmvE, EINukeMode mode)
{

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils", pDEBUG)
    << "PreEquilibrium() is invoked for a : " << p->Name()
    << " whose kinetic energy is : " << p->KinE();
#endif

  // Random number generator
  RandomGen * rnd = RandomGen::Instance();
  //unused PDGLibrary * pLib = PDGLibrary::Instance();

  bool allow_dup = true;
  PDGCodeList list(allow_dup); // list of final state particles

  double ppcnt = (double) RemnZ / (double) RemnA; // % of protons left

  // figure out the final state conditions

  if(p->Pdg()==kPdgProton) list.push_back(kPdgProton);
  else if(p->Pdg()==kPdgNeutron) list.push_back(kPdgNeutron);

  for(int i=0;i<3;i++)
    {
      if(rnd->RndFsi().Rndm()<ppcnt)
        {
          list.push_back(kPdgProton);
          RemnZ--;
        }
      else list.push_back(kPdgNeutron);

      RemnA--;

      ppcnt = (double) RemnZ / (double) RemnA;
    }

  // Add the fermi energy of the three nucleons to the phase space
  /*
  if(DoFermi)
    {
      Target target(ev->TargetNucleus()->Pdg());
      TVector3 pBuf = p->P4()->Vect();
      double mBuf = p->Mass();
      double eBuf = TMath::Sqrt(pBuf.Mag2() + mBuf*mBuf);
      TLorentzVector tSum(pBuf,eBuf);
      double mSum = 0.0;
      vector<int>::const_iterator pdg_iter;
      for(pdg_iter=++(list.begin());pdg_iter!=list.end();++pdg_iter)
        {
          target.SetHitNucPdg(*pdg_iter);
          Nuclmodel->GenerateNucleon(target);
          mBuf = pLib->Find(*pdg_iter)->Mass();
          mSum += mBuf;
          pBuf = FermiFac * Nuclmodel->Momentum3();
          eBuf = TMath::Sqrt(pBuf.Mag2() + mBuf*mBuf);
          tSum += TLorentzVector(pBuf,eBuf);
          RemnP4 -= TLorentzVector(pBuf,eBuf-mBuf);
        }
      TLorentzVector dP4 = tSum + TLorentzVector(TVector3(0,0,0),-mSum);
      p->SetMomentum(dP4);
      }
  */
  // do the phase space decay & save all f/s particles to the event record
  bool success = genie::utils::intranuke2018::PhaseSpaceDecay(ev,p,list,RemnP4,NucRmvE,mode);
  if(success)  LOG("INukeUtils2018",pINFO) << "Successful phase space decay for pre-equilibrium nucleus FSI event";
  else
    {
      exceptions::INukeException exception;
      exception.SetReason("Phase space generation of pre-equilibrium nucleus final state failed, details above");
      throw exception;
    }

  int p_loc = 0;
  while(p_loc<ev->GetEntries())
    {
      GHepParticle * p_ref = ev->Particle(p_loc);
      if(!p->ComparePdgCodes(p_ref)) p_loc++;
      else
        {
          if(!p->CompareStatusCodes(p_ref)) p_loc++;
          else
            {
              if(!p->CompareMomentum(p_ref)) p_loc++;
              else break;
            }
        }
     }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils", pDEBUG)
    << "Particle at: " << p_loc;
#endif

 // find the appropriate daughter
  vector<int> * descendants = ev->GetStableDescendants(p_loc);

  int loc       = p_loc + 1;
  int f_loc     = p_loc + 1;
  double energy = ev->Particle(loc)->E();

/*  // (1) least energetic
  double min_en = energy;

  for(unsigned int j=0;j<descendants->size();j++)
    {
      loc = (*descendants)[j];
      energy = ev->Particle(loc)->E();
      if(energy<min_en)
        {
          f_loc = loc;
          min_en = energy;
        }
    }
*/
  // (2) most energetic
  double max_en = energy;

  for(unsigned int j=0;j<descendants->size();j++)
    {
      loc = (*descendants)[j];
      energy = ev->Particle(loc)->E();
      if(energy>max_en)
        {
          f_loc = loc;
          max_en = energy;
        }
    }

  // (3) 1st particle
  // ...just use the defaulted f_loc

  delete descendants;

  // change particle status for decaying particle - take out as test
  //ev->Particle(f_loc)->SetStatus(kIStIntermediateState);
  // decay a clone particle
  GHepParticle * t = new GHepParticle(*(ev->Particle(f_loc)));
  t->SetFirstMother(f_loc);
  //next statement was in Alex Bell's original code - PreEq, then Equilibrium using particle with highest energy.  Note it gets IST=kIStIntermediateState.
  //genie::utils::intranuke2018::Equilibrium(ev,t,RemnA,RemnZ,RemnP4,DoFermi,FermiFac,Nuclmodel,NucRmvE,mode);

  delete t;
}
//___________________________________________________________________________
// Method to handle Equilibrium reaction
// Alex Bell -- 6/17/2008
void genie::utils::intranuke2018::Equilibrium(
  GHepRecord * ev, GHepParticle * p,
  int &RemnA, int &RemnZ, TLorentzVector &RemnP4,
  bool /* DoFermi */, double /* FermiFac */,
  const NuclearModelI* /* Nuclmodel */, double NucRmvE, EINukeMode mode)
{

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils", pDEBUG)
    << "Equilibrium() is invoked for a : " << p->Name()
    << " whose kinetic energy is : " << p->KinE();
#endif

  // Random number generator
  RandomGen * rnd = RandomGen::Instance();
  //usused PDGLibrary * pLib = PDGLibrary::Instance();

  bool allow_dup = true;
  PDGCodeList list(allow_dup); // list of final state particles

  // % of protons left
  double ppcnt = (double) RemnZ / (double) RemnA;

  // figure out the final state conditions

  if(p->Pdg()==kPdgProton) list.push_back(kPdgProton);
  else if(p->Pdg()==kPdgNeutron) list.push_back(kPdgNeutron);

  //add additonal particles to stack
  for(int i=0;i<4;i++)
    {
      if(rnd->RndFsi().Rndm()<ppcnt)
        {
          list.push_back(kPdgProton);
          RemnZ--;
        }
      else list.push_back(kPdgNeutron);

      RemnA--;

      ppcnt = (double) RemnZ / (double) RemnA;
    }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("INukeUtils2018", pDEBUG)
    << "Remnant nucleus (A,Z) = (" << RemnA << ", " << RemnZ << ")";
#endif

  // Add the fermi energy of the three nucleons to the phase space
  /*  if(DoFermi)
    {
      Target target(ev->TargetNucleus()->Pdg());
      TVector3 pBuf = p->P4()->Vect();
      double mBuf = p->Mass();
      double eBuf = TMath::Sqrt(pBuf.Mag2() + mBuf*mBuf);
      TLorentzVector tSum(pBuf,eBuf);
      double mSum = 0.0;
      vector<int>::const_iterator pdg_iter;
      for(pdg_iter=++(list.begin());pdg_iter!=list.end();++pdg_iter)
        {
          target.SetHitNucPdg(*pdg_iter);
          Nuclmodel->GenerateNucleon(target);
          mBuf = pLib->Find(*pdg_iter)->Mass();
          mSum += mBuf;
          pBuf = FermiFac * Nuclmodel->Momentum3();
          eBuf = TMath::Sqrt(pBuf.Mag2() + mBuf*mBuf);
          tSum += TLorentzVector(pBuf,eBuf);
          RemnP4 -= TLorentzVector(pBuf,eBuf-mBuf);
        }
      TLorentzVector dP4 = tSum + TLorentzVector(TVector3(0,0,0),-mSum);
      p->SetMomentum(dP4);
    }
  */
  // do the phase space decay & save all f/s particles to the record
  bool success = genie::utils::intranuke2018::PhaseSpaceDecay(ev,p,list,RemnP4,NucRmvE,mode);
  if (success) LOG("INukeUtils",pINFO) << "successful equilibrium interaction";
  else
    {
      exceptions::INukeException exception;
      exception.SetReason("Phase space generation of compound nucleus final state failed, details above");
      throw exception;
    }

}


//___________________________________________________________________________
bool genie::utils::intranuke2018::TwoBodyCollision(
  GHepRecord* ev, int pcode, int tcode, int scode, int s2code, double C3CM,
  GHepParticle* p, GHepParticle* t, int &RemnA, int &RemnZ, TLorentzVector &RemnP4, EINukeMode mode)
{
  // Aaron Meyer (10/29/09)
  // Adapted from kinematics in other function calls
  //
  // C3CM is the cosine of the scattering angle, calculated before calling
  // p and t are the output particles, must be allocated before calling
  // pcode,tcode,scode,s2code are initial and final particle PDG codes in scattering
  // return value used for error checking

  // Kinematic variables

  double M1, /* M2, */ M3, M4; // rest energies, in GeV
  double E3L, P3L, E4L, P4L;
  TVector3 tP1L, tPtot, tbeta, tbetadir, tTrans, tVect;
  TVector3 tP1zCM, tP2zCM, tP3L, tP4L;

  // Library instance for reference
  PDGLibrary * pLib = PDGLibrary::Instance();

  // random number generator
  //RandomGen * rnd = RandomGen::Instance();

  // handle fermi target
  Target target(ev->TargetNucleus()->Pdg());

  // get mass for particles
  M1 = pLib->Find(pcode)->Mass();
  // usused // M2 = pLib->Find(tcode)->Mass();
  M3 = pLib->Find(scode)->Mass();
  M4 = pLib->Find(s2code)->Mass();

  // get lab energy and momenta and assign to 4 vectors
  TLorentzVector t4P1L = *p->P4();
  TLorentzVector t4P2L = *t->P4();

  // binding energy
  double bindE = 0.025; // empirical choice, might need to be improved
  //double bindE = 0.0;

  LOG("TwoBodyCollision",pNOTICE) << "M1 = " << t4P1L.M() << " , M2 = " << t4P2L.M();
  LOG("TwoBodyCollision",pNOTICE) << "E1 = " << t4P1L.E() << " , E2 = " << t4P2L.E();

  if ( (p->Energy()-p->Mass()) < bindE ){bindE = 0.;} // if the probe's energy is less than the binding energy, set the binding energy to 0.

  // don't use BE unless kinetic energy >> BE.
  if((pcode==2112||pcode==2212)&&(t4P1L.E()-M1)<.1) bindE = 0.0;
  if((pcode==211||pcode==-211||pcode==111)&&(t4P1L.E()-M1)<.08) bindE = 0.0;
  if((pcode==321)&&(t4P1L.E()-M1)<.1) bindE = 0.0;

  // carry out scattering
  TLorentzVector t4P3L, t4P4L;
  if (!TwoBodyKinematics(M3,M4,t4P1L,t4P2L,t4P3L,t4P4L,C3CM,RemnP4,bindE))
    {
      P3L = t4P3L.Vect().Mag();
      P4L = t4P4L.Vect().Mag();
      E3L = t4P3L.E();
      E4L = t4P4L.E();

      LOG("TwoBodyCollision",pNOTICE)
        << "TwoBodyKinematics fails: C3CM, P3 = " << C3CM << "  "
        << P3L << "   " << E3L << "\n" << "             P4 = "
        << P4L << "   " << E4L ;
      return false;     //covers all possiblities for now
    }

  // error checking
  P3L = t4P3L.Vect().Mag();
  P4L = t4P4L.Vect().Mag();
  E3L = t4P3L.E();
  E4L = t4P4L.E();
  LOG("INukeUtils",pINFO)
    << "C3CM, P3 = " << C3CM << "  "
    << P3L << "   " << E3L << "\n" << "             P4 = "
    << P4L << "   " << E4L ;

  // handle very low momentum particles
  if(!(TMath::Finite(P3L)) || P3L<.001)
    {
      LOG("INukeUtils",pINFO)
        << "Particle 3 momentum small or non-finite: " << P3L
        << "\n" << "--> Assigning .001 as new momentum";
      P3L = .001;
      E3L = TMath::Sqrt(P3L*P3L + M3*M3);
    }
  if(!(TMath::Finite(P4L)) || P4L<.001)
    {
      LOG("INukeUtils",pINFO)
        << "Particle 4 momentum small or non-finite: " << P4L
        << "\n" << "--> Assigning .001 as new momentum";
      P4L = .001;
      E4L = TMath::Sqrt(P4L*P4L + M4*M4);
    }

  // if this is going to be on in the future, remember to not apply PB for Oset
  /*// pauli blocking  turn off for now to better match data
  //  if(P3L<fFermiMomentum && pdg::IsNeutronOrProton(scode) ||
  //     P4L<fFermiMomentum && pdg::IsNeutronOrProton(s2code) )
  if(P3L<.25 && pdg::IsNeutronOrProton(scode) ||
     P4L<.25 && pdg::IsNeutronOrProton(s2code) )
    {
      LOG("INukeUtils",pNOTICE)<< "TwoBodyCollision failed: Pauli blocking";
      p->SetStatus(kIStHadronInTheNucleus);
      RemnP4 -= TLorentzVector(0,0,0,bindE);
      return false;
      }*/

  // update remnant nucleus
  RemnP4 -= t4P2L;
  LOG("INukeUtils",pINFO)
    << "t4P2L= " << t4P2L.E() << "  " << t4P2L.Z()
    << "  RemnP4= " << RemnP4.E() << "   " << RemnP4.Z()  ;
  if (tcode==kPdgProton) {RemnZ--;RemnA--;}
  else if(tcode==kPdgNeutron) RemnA--;

  // create t particle w/ appropriate momenta, code, and status
  // Set target's mom to be the mom of the hadron that was cloned
  t->SetFirstMother(p->FirstMother());
  t->SetLastMother(p->LastMother());

  // adjust p to reflect scattering
  p->SetPdgCode(scode);
  p->SetMomentum(t4P3L);

  t->SetPdgCode(s2code);
  t->SetMomentum(t4P4L);

  if (mode==kIMdHN)
  {
    p->SetStatus(kIStHadronInTheNucleus);
    t->SetStatus(kIStHadronInTheNucleus);
  }
  else
  {
    p->SetStatus(kIStStableFinalState);
    t->SetStatus(kIStStableFinalState);
  }
  LOG("INukeUtils",pINFO) << "Successful 2 body collision";
  return true;

}
//___________________________________________________________________________
bool genie::utils::intranuke2018::TwoBodyKinematics(
  double M3, double M4, TLorentzVector t4P1L, TLorentzVector t4P2L,
  TLorentzVector &t4P3L, TLorentzVector &t4P4L, double C3CM, TLorentzVector &RemnP4, double bindE)
{
  // Aaron Meyer (05/17/10)
  // Adapted from kinematics in other function calls
  //
  // Outgoing particle masses M3,M4
  // Scatters particles according to normal two body collisions
  //
  // bindE is the binding energy (GeV) of a particle removed from the nucleus (default 0)
  // For nonzero binding energy, remove the binding energy from the total energy,
  //   then put both of the particles back on mass shell by shifting momentum/energy
  //   of target
  // Momentum only shifted in the direction parallel to the probe's motion
  //
  // Rotates final transverse component of momentum by a random angle from 0->2pi
  // Return value for error checking
  // Gives outgoing 4-momenta of particles 3 and 4 (t4P3L, t4P4L respectively)
  //
  // All 4-momenta should be on mass shell

  double E1L, E2L, P1L, P2L, E3L, P3L;
  double beta, gm; // speed and gamma for CM frame in Lab
  double S3CM; // sin of scattering angle
  double PHI3;
  double E1CM, E2CM, E3CM, P3CM;//, E4CM, P4CM;
  double P3zL, P3tL;//, P4zL, P4tL;
  double Et;
  double theta1, theta2, theta5, P1zL, P2zL, P1tL, P2tL;
  TVector3 tbeta, tbetadir, tTrans, tVect;
  TVector3 tP1zCM, tP2zCM, vP3L;
  TLorentzVector t4P1buf, t4P2buf, t4Ptot;

  // Library instance for reference
  //PDGLibrary * pLib = PDGLibrary::Instance();

  // random number generator
  RandomGen * rnd = RandomGen::Instance();

  // error checking
  if (C3CM < -1. || C3CM > 1.) return false;

  // calculate sine from scattering angle
  S3CM = TMath::Sqrt(1.0 - C3CM*C3CM);

  // fill buffers
  t4P1buf = t4P1L;
  t4P2buf = t4P2L;

  // get lab energy and momenta
  E1L = t4P1buf.E();
  P1L = t4P1buf.P();
  E2L = t4P2buf.E();
  P2L = t4P2buf.P();
  t4Ptot = t4P1buf + t4P2buf;

  LOG("INukeUtils",pINFO) <<"M1   "<<t4P1L.M()<<  ", M2    "<<t4P2L.M();
  LOG("INukeUtils",pINFO) <<"E1L  "<<E1L<< ", E1CM  "<<E1CM;
  LOG("INukeUtils",pINFO) <<"bindE = " << bindE;

  // binding energy
  if (bindE!=0)
    {

      E1L -= bindE;

      if (E1L+E2L < M3+M4)
        {
          LOG("INukeUtils",pNOTICE) <<"TwoBodyKinematics Failed: Forbidden by binding energy";
          LOG("INukeUtils",pNOTICE) <<"E1L, E2L, M3, M4 : "<< E1L <<", "<< E2L <<", "<< M3 <<", "<< M4;
          t4P3L.SetPxPyPzE(0,0,0,0);
          t4P4L.SetPxPyPzE(0,0,0,0);
          return false;
          }
    }

  // calculate beta and gamma
  tbeta = t4Ptot.Vect() * (1.0 / (E1L + E2L));
  tbetadir = tbeta.Unit();
  beta = tbeta.Mag();
  gm = 1.0 / TMath::Sqrt(1.0 - beta*beta);

  // get angle and component info
  theta1 = t4P1buf.Angle(tbeta);
  theta2 = t4P2buf.Angle(tbeta);
  P1zL = P1L*TMath::Cos(theta1);
  P2zL = P2L*TMath::Cos(theta2);
  P1tL = TMath::Sqrt(P1L*P1L - P1zL*P1zL);
  P2tL = -TMath::Sqrt(P2L*P2L - P2zL*P2zL);
  tVect.SetXYZ(1,0,0);
  if(TMath::Abs((tVect - tbetadir).Mag())<.01) tVect.SetXYZ(0,1,0);
  theta5 = tVect.Angle(tbetadir);
  tTrans = (tVect - TMath::Cos(theta5)*tbetadir).Unit();

  // boost to CM frame to get scattered particle momenta
  E1CM = gm*E1L - gm*beta*P1zL;
  tP1zCM = gm*P1zL*tbetadir - gm*tbeta*E1L;
  E2CM = gm*E2L - gm*beta*P2zL;
  tP2zCM = gm*P2zL*tbetadir - gm*tbeta*E2L;
  Et = E1CM + E2CM;
//-------

  LOG("INukeUtils",pINFO) <<"M1   "<<t4P1L.M()<<  ", M2    "<<t4P2L.M();
  LOG("INukeUtils",pINFO) <<"E1L  "<<E1L<< ", E1CM  "<<E1CM;
  LOG("INukeUtils",pINFO) <<"P1zL "<<P1zL<<", P1zCM "<<tP1zCM.Mag()<<", P1tL "<<P1tL;
  LOG("INukeUtils",pINFO) <<"E2L  "<<E2L<< ", E2CM  "<<E2CM;
  LOG("INukeUtils",pINFO) <<"P2zL "<<P2zL<<", P2zCM "<<tP2zCM.Mag()<<", P2tL "<<P2tL;
  LOG("INukeUtils",pINFO) <<"C3CM "<<C3CM;

//-------
  E3CM = (Et*Et + M3*M3 - M4*M4) / (2.0 * Et);

  // check to see if decay is viable
  if(E3CM*E3CM - M3*M3<0 || E3CM<0 || Et<0)
  {
    if (Et<0) LOG("INukeUtils",pNOTICE) <<"TwoBodyKinematics Failed: Total energy is negative";
    if (E3CM<M3) LOG("INukeUtils",pNOTICE) <<"TwoBodyKinematics Failed: Scattered Particle 3 CM energy is too small";
    if (E3CM*E3CM - M3*M3<0) LOG("INukeUtils",pNOTICE) <<"TwoBodyKinematics Failed: Scattered Particle 3 CM momentum is nonreal";
    t4P3L.SetPxPyPzE(0,0,0,0);
    t4P4L.SetPxPyPzE(0,0,0,0);
    return false;
  }
  P3CM = TMath::Sqrt(E3CM*E3CM - M3*M3);

  // boost back to lab
  P3zL = gm*beta*E3CM + gm*P3CM*C3CM;
  P3tL = P3CM*S3CM;

  P3L = TMath::Sqrt(P3zL*P3zL + P3tL*P3tL);
  E3L = TMath::Sqrt(P3L*P3L + M3*M3);

  //-------

  double E4CM = Et-E3CM;
  double P4zL = gm*beta*E4CM - gm*P3CM*C3CM;
  double P4tL = -1.*P3tL;
  double P4L = TMath::Sqrt(P4zL*P4zL + P4tL*P4tL);
  double E4L = TMath::Sqrt(P4L*P4L + M4*M4);

  LOG("INukeUtils",pINFO) <<"M3   "<< M3 <<  ", M4    "<< M4;
  LOG("INukeUtils",pINFO) <<"E3L   "<<E3L<< ", E3CM "<<E3CM;
  LOG("INukeUtils",pINFO) <<"P3zL  "<<P3zL<<", P3tL "<<P3tL;
  LOG("INukeUtils",pINFO) <<"C3L   "<<P3zL/P3L;
  LOG("INukeUtils",pINFO) <<"Check:";
  LOG("INukeUtils",pINFO) <<"E4L   "<<E4L<< ", E4CM "<<E4CM;
  LOG("INukeUtils",pINFO) <<"P4zL  "<<P4zL<<", P4tL "<<P4tL;
  LOG("INukeUtils",pINFO) <<"P4L   "<<P4L;
  LOG("INukeUtils",pINFO) <<"C4L   "<<P4zL/P4L;

  double echeck = E1L + E2L - (E3L + E4L);
  double pzcheck = P1zL+ P2zL - (P3zL + P4zL);
  double ptcheck = P1tL+ P2tL - (P3tL + P4tL);

  LOG("INukeUtils",pINFO) <<"Check 4-momentum conservation -  Energy  "<<echeck<<", z momentum "<<pzcheck << ",    transverse momentum  " << ptcheck ;

  // -------

  // handle very low momentum particles
  if(!(TMath::Finite(P3L)) || P3L<.001)
    {
      LOG("INukeUtils",pINFO)
        << "Particle 3 momentum small or non-finite: " << P3L
        << "\n" << "--> Assigning .001 as new momentum";
      P3tL = 0;
      P3zL = .001;
      P3L = .001;
      E3L = TMath::Sqrt(P3L*P3L + M3*M3);
    }

  // get random phi angle, distributed uniformally in 360 deg
  PHI3 = 2 * kPi * rnd->RndFsi().Rndm();

  vP3L = P3zL*tbetadir + P3tL*tTrans;
  vP3L.Rotate(PHI3,tbetadir);

  t4P3L.SetVect(vP3L);
  t4P3L.SetE(E3L);

  t4P4L = t4P1buf + t4P2buf - t4P3L;
  t4P4L-= TLorentzVector(0,0,0,bindE);
  /*LOG("INukeUtils",pINFO) <<"GENIE:";
  LOG("INukeUtils",pINFO) <<"E4L   "<<t4P4L.E();
  LOG("INukeUtils",pINFO) <<"P4zL  "<<t4P4L.Vect()*tbetadir<<", P4tL "<<-1.*TMath::Sqrt(t4P4L.Vect().Mag2()-TMath::Power(t4P4L.Vect()*tbetadir,2.));
  LOG("INukeUtils",pINFO) <<"P4L   "<<t4P4L.Vect().Mag();
  LOG("INukeUtils",pINFO) <<"C4L   "<<t4P4L.Vect()*tbetadir/t4P4L.Vect().Mag();  */

  if(t4P4L.Mag2()<0 || t4P4L.E()<0)
  {
    LOG("INukeUtils",pNOTICE)<<"TwoBodyKinematics Failed: Target mass or energy is negative";
    t4P3L.SetPxPyPzE(0,0,0,0);
    t4P4L.SetPxPyPzE(0,0,0,0);
    return false;
  }

  if (bindE!=0) RemnP4 += TLorentzVector(0,0,0,bindE);
  return true;
}
//___________________________________________________________________________
bool genie::utils::intranuke2018::ThreeBodyKinematics(
  GHepRecord* ev, GHepParticle* p, int tcode, GHepParticle* s1, GHepParticle* s2, GHepParticle* s3,
  bool DoFermi, double FermiFac, double FermiMomentum, const NuclearModelI* Nuclmodel)
{

  // Aaron Meyer (7/15/10)
  //
  // Kinematics used in utils::intranuke2018::PionProduction
  // Handles the kinematics of three body scattering
  //
  // s1,s2,and s3 are pointers to particles with the PDG code that needs to be scattered
  // the last four variables are for Fermi momentum and pauli blocking, default will use none of them

  // kinematic variables
  double M1, M2, M3, M4, M5; // rest energies, in GeV
  double P1L, P2L, P3L, P4L, P5L;
  double E1L, E2L, E3L, E4L, E5L;
  double E1CM, E2CM, P3tL;
  double PizL, PitL, PiL, EiL;
  double EiCM, P4CM2, E4CM2, E5CM2, P3CM, E3CM;
  double beta, gm, beta2, gm2;
  double P3zL, P4zL, P4tL, P5zL, P5tL;
  double Et, M, theta1, theta2;
  double P1zL, P2zL;
  double theta3, theta4, phi3, phi4, theta5;
  TVector3 tP2L, tP1L, tPtot, tbeta, tbetadir, tTrans, tP4L, tP5L;
  TVector3 tP1zCM, tP2zCM, tP3L, tPiL, tbeta2, tbetadir2, tVect, tTrans2;

  // Library instance for reference
  PDGLibrary * pLib = PDGLibrary::Instance();

  // random number generator
  RandomGen * rnd = RandomGen::Instance();

  M1 = pLib->Find(p->Pdg())->Mass();
  M2 = pLib->Find(tcode)->Mass();
  M3 = pLib->Find(s1->Pdg())->Mass();
  M4 = pLib->Find(s2->Pdg())->Mass();
  M5 = pLib->Find(s3->Pdg())->Mass();

  // set up fermi target
  Target target(ev->TargetNucleus()->Pdg());

  // handle fermi momentum
  if(DoFermi)
    {
      target.SetHitNucPdg(tcode);
      Nuclmodel->GenerateNucleon(target);
      tP2L = FermiFac * Nuclmodel->Momentum3();
      P2L = tP2L.Mag();
      E2L = TMath::Sqrt(tP2L.Mag2() + M2*M2);
    }
  else
    {
      tP2L.SetXYZ(0.0, 0.0, 0.0);
      P2L = 0;
      E2L = M2;
    }

  // first sequence, handle 4th and 5th products as composite
  E1L = p->E();

  P1L = TMath::Sqrt(E1L*E1L - M1*M1);
  tP1L = p->P4()->Vect();
  tPtot = tP1L + tP2L;

  tbeta = tPtot * (1.0 / (E1L + E2L));
  tbetadir = tbeta.Unit();
  beta = tbeta.Mag();
  gm = 1.0 / TMath::Sqrt(1.0 - beta*beta);

  theta1 = tP1L.Angle(tbeta);
  theta2 = tP2L.Angle(tbeta);
  P1zL = P1L*TMath::Cos(theta1);
  P2zL = P2L*TMath::Cos(theta2);
  tVect.SetXYZ(1,0,0);
  if(TMath::Abs((tVect - tbetadir).Mag())<.01) tVect.SetXYZ(0,1,0);
  theta5 = tVect.Angle(tbetadir);
  tTrans = (tVect - TMath::Cos(theta5)*tbetadir).Unit();

  E1CM = gm*E1L - gm*beta*P1zL;
  tP1zCM = gm*P1zL*tbetadir - gm*tbeta*E1L;
  E2CM = gm*E2L - gm*beta*P2zL;
  tP2zCM = gm*P2zL*tbetadir - gm*tbeta*E2L;
  Et = E1CM + E2CM;
  M = (rnd->RndFsi().Rndm()*(Et - M3 - M4 - M5)) + (M4 + M5);
  E3CM = (Et*Et + M3*M3 - M*M)/(2*Et);
  EiCM = Et - E3CM;
  if(E3CM*E3CM - M3*M3<0)
  {
    LOG("INukeUtils",pNOTICE)
      << "PionProduction P3 has non-real momentum - retry kinematics";
    LOG("INukeUtils",pNOTICE) << "Energy, masses of 3 fs particales:"
      << E3CM << "  " << M3 << "  " << "  " << M4 << "  " << M5;
    exceptions::INukeException exception;
    exception.SetReason("PionProduction particle 3 has non-real momentum");
    throw exception;
    return false;
  }
  P3CM = TMath::Sqrt(E3CM*E3CM - M3*M3);

  theta3 =   kPi * rnd->RndFsi().Rndm();
  theta4 =   kPi * rnd->RndFsi().Rndm();
  phi3   = 2*kPi * rnd->RndFsi().Rndm();
  phi4   = 2*kPi * rnd->RndFsi().Rndm();

  P3zL = gm*beta*E3CM + gm*P3CM*TMath::Cos(theta3);
  P3tL =  P3CM*TMath::Sin(theta3);
  PizL = gm*beta*EiCM - gm*P3CM*TMath::Cos(theta3);
  PitL = -P3CM*TMath::Sin(theta3);

  P3L = TMath::Sqrt(P3zL*P3zL + P3tL*P3tL);
  PiL = TMath::Sqrt(PizL*PizL + PitL*PitL);
  E3L = TMath::Sqrt(P3L*P3L + M3*M3);
  EiL = TMath::Sqrt(PiL*PiL + M*M);

  // handle very low momentum particles
  if(!(TMath::Finite(P3L)) || P3L < .001)
    {
      LOG("INukeUtils",pINFO)
        << "Particle 3 " << M3 << " momentum small or non-finite: " << P3L
        << "\n" << "--> Assigning .001 as new momentum";
      P3tL = 0;
      P3zL = .001;
      P3L = .001;
      E3L = TMath::Sqrt(P3L*P3L + M3*M3);
    }

  tP3L = P3zL*tbetadir + P3tL*tTrans;
  tPiL = PizL*tbetadir + PitL*tTrans;
  tP3L.Rotate(phi3,tbetadir);
  tPiL.Rotate(phi3,tbetadir);

  // second sequence, handle formally composite particles 4 and 5
  tbeta2 = tPiL * (1.0 / EiL);
  tbetadir2 = tbeta2.Unit();
  beta2 = tbeta2.Mag();
  gm2 = 1.0 / TMath::Sqrt(1.0 - beta2*beta2);

  E4CM2 = (M*M + M4*M4 - M5*M5) / (2*M);
  E5CM2 = M - E4CM2;
  P4CM2 = TMath::Sqrt(E4CM2*E4CM2 - M4*M4);

  tVect.SetXYZ(1,0,0);
  if(TMath::Abs((tVect - tbetadir2).Mag())<.01) tVect.SetXYZ(0,1,0);
  theta5 = tVect.Angle(tbetadir2);
  tTrans2 = (tVect - TMath::Cos(theta5)*tbetadir2).Unit();

  P4zL = gm2*beta2*E4CM2 + gm2*P4CM2*TMath::Cos(theta4);
  P4tL = P4CM2*TMath::Sin(theta4);
  P5zL = gm2*beta2*E5CM2 - gm2*P4CM2*TMath::Cos(theta4);
  P5tL = - P4tL;

  P4L = TMath::Sqrt(P4zL*P4zL + P4tL*P4tL);
  P5L = TMath::Sqrt(P5zL*P5zL + P5tL*P5tL);
  E4L = TMath::Sqrt(P4L*P4L + M4*M4);
  E5L = TMath::Sqrt(P5L*P5L + M5*M5);

  // handle very low momentum particles
  if(!(TMath::Finite(P4L)) || P4L < .001)
    {
      LOG("INukeUtils",pINFO)
        << "Particle 4 " << M4 << " momentum small or non-finite: " << P4L
        << "\n" << "--> Assigning .001 as new momentum";
      P4tL = 0;
      P4zL = .001;
      P4L = .001;
      E4L = TMath::Sqrt(P4L*P4L + M4*M4);
    }
  if(!(TMath::Finite(P5L)) || P5L < .001)
    {
      LOG("INukeUtils",pINFO)
        << "Particle 5 " << M5 << " momentum small or non-finite: " << P5L
        << "\n" << "--> Assigning .001 as new momentum";
      P5tL = 0;
      P5zL = .001;
      P5L = .001;
      E5L = TMath::Sqrt(P5L*P5L + M5*M5);
    }

  tP4L = P4zL*tbetadir2 + P4tL*tTrans2;
  tP5L = P5zL*tbetadir2 + P5tL*tTrans2;
  tP4L.Rotate(phi4,tbetadir2);
  tP5L.Rotate(phi4,tbetadir2);

  // pauli blocking
  if(P3L < FermiMomentum || ( pdg::IsNeutronOrProton(s2->Pdg()) && P4L < FermiMomentum ) )
  {
    LOG("INukeUtils",pNOTICE)
      << "PionProduction fails because of Pauli blocking - retry kinematics";
    exceptions::INukeException exception;
    exception.SetReason("PionProduction final state not determined");
    throw exception;
    return false;
  }

  // create scattered particles w/ appropriate momenta, code, and status
  // Set moms to be the moms of the hadron that was cloned
  s1->SetFirstMother(p->FirstMother());
  s2->SetFirstMother(p->FirstMother());
  s3->SetFirstMother(p->FirstMother());
  s1->SetLastMother(p->LastMother());
  s2->SetLastMother(p->LastMother());
  s3->SetLastMother(p->LastMother());

  TLorentzVector(tP3L,E3L);
  TLorentzVector(tP4L,E4L);
  TLorentzVector(tP5L,E5L);

  s1->SetMomentum(TLorentzVector(tP3L,E3L));
  s2->SetMomentum(TLorentzVector(tP4L,E4L));
  s3->SetMomentum(TLorentzVector(tP5L,E5L));
  int mode = kIMdHA;
  LOG ("INukeUtils",pDEBUG) << "in Pi Prod, mode =  " << mode;
  if (mode==kIMdHN)
    {
      s1->SetStatus(kIStHadronInTheNucleus);
      s2->SetStatus(kIStHadronInTheNucleus);
      s3->SetStatus(kIStHadronInTheNucleus);
    }
  else
    {
      s1->SetStatus(kIStStableFinalState);
      s2->SetStatus(kIStStableFinalState);
      s3->SetStatus(kIStStableFinalState);
    }
  return true;
}
//___________________________________________________________________________
bool genie::utils::intranuke2018::PionProduction(
  GHepRecord* ev, GHepParticle* p, GHepParticle* s1, GHepParticle* s2, GHepParticle* s3, int &RemnA, int &RemnZ,
  TLorentzVector &RemnP4, bool DoFermi, double FermiFac, double FermiMomentum, const NuclearModelI* Nuclmodel)
{
  // Aaron Meyer (7/15/2010)
  //
  // Handles pion production reactions in both hA and hN mode
  // Calculates fundamental cross sections from fit functions
  // Uses isospin relations to determine the rest of cross sections
  //
  // p is the probe particle
  // s1, s2, and s3 are the particles produced in the reaction
  // must set the status and add particles to the event record after returning from this method
  // return value for error checking


  // random number generator
  RandomGen * rnd = RandomGen::Instance();

  // library reference
  PDGLibrary * pLib = PDGLibrary::Instance();

  bool ptarg = false;
  int pcode = p->Pdg();

  int p1code = p->Pdg();
  // Randomly determine target and 1st product baryons
  int p3code = 0, p4code = 0, p5code = 0;

  //
  // Uses a fit curve log(sigma) = a - b/(T_pi - c) for pions
  // Fit parameters determined by Roman Tacik (4/3/09)
  // pi- & p cross sections are assumed to be the same as pi+ & n
  //
  // Fit curve for nucleons:
  // sigma = a*(1+b*e^(-c*(eta-d)^2))*(1-e^(-(f*eta)^g))*(1-e^(-h/eta^2))
  // 7 parameters (a,b,c,d,f,g,h)
  // eta is maximum kinematically allowed momentum of the pion, normalized by the mass
  // Uses isotopic spin decomposition of total cross sections
  //

   if ((p1code==kPdgPi0)||(p1code==kPdgPiP)||(p1code==kPdgPiM)) {

     double kine = 1000*p->KinE();

     // Determine cross sections

     // pion
     // pi- & p
     //  -> pi0 & pi0 & n
     //  a = 8.82; b = 573.2; c = 107.3;
     double xsec2pi0n = TMath::Max(0.,TMath::Exp(8.82 - (573.2/(kine-107.3))));
     //  -> pi- & pi+ & n
     //  a = 11.06; b = 985.9; c = 88.2;
     double xsecpippimn = TMath::Max(0.,TMath::Exp(11.06 - (985.9/(kine-88.2))));
     //  -> pi- & pi0 & p
     //  a = 9.58; b = 1229.4; c = 60.5;
     double xsecpimpi0p = TMath::Max(0.,TMath::Exp(9.58 - (1229.4/(kine-60.5))));
     double totpimp = xsec2pi0n + xsecpippimn + xsecpimpi0p;


     // pi+ & p
     //  -> pi+ & pi+ & n
     //  a = 5.64; b = 222.6; c = 150.0;
     double xsec2pipn = TMath::Max(0.,TMath::Exp(5.64 - (222.6/(kine-150.))));
     //  -> pi+ & pi0 & p
     //  a = 7.95; b = 852.6; c = 77.8;
     double xsecpippi0p = TMath::Max(0.,TMath::Exp(7.95 - (852.6/(kine-77.8))));
     double totpipp = xsec2pipn + xsecpippi0p;

     if (totpimp<=0 && totpipp<=0) {
       LOG("INukeUtils",pNOTICE) << "InelasticHN called below threshold energy";
       p->SetStatus(kIStHadronInTheNucleus);
       ev->AddParticle(*p);
       return false;
        }

     double xsecp, xsecn;
     switch (p1code) {
     case kPdgPi0:  xsecp = 0.5 * (totpimp + totpipp); xsecn = xsecp; break;
     case kPdgPiP:  xsecp = totpipp; xsecn = totpimp; break;
     case kPdgPiM:  xsecp = totpimp; xsecn = totpipp; break;
     default:
       LOG("INukeUtils",pWARN) << "InelasticHN cannot handle probe: "
                               << PDGLibrary::Instance()->Find(p1code)->GetName();
       exceptions::INukeException exception;
       exception.SetReason("PionProduction final state not determined");
       throw exception;
       return false;
       break;
     }

     // Normalize cross sections by Z or A-Z

     xsecp *= RemnZ;
     xsecn *= RemnA-RemnZ;

     // determine target

     double rand = rnd->RndFsi().Rndm() * (xsecp + xsecn);
     if (rand < xsecp) // proton target
       { rand /= RemnZ; ptarg = true;}
     else // neutron target
       { rand -= xsecp; rand /= RemnA-RemnZ; ptarg = false;}


     // determine final state

     if (((ptarg==true)&&(p1code==kPdgPiP))
         || ((ptarg==false)&&(p1code==kPdgPiM)))
       {
         if (rand < xsec2pipn) // pi+ & pi+ & n final state
           {
             p3code = (ptarg ? kPdgNeutron : kPdgProton);
             p4code = p1code;
             p5code = p4code;
           }
         else {  // pi+ & pi0 & p final state
           p3code = (ptarg ? kPdgProton : kPdgNeutron);
           p4code = p1code;
           p5code = kPdgPi0;
         }
       }
     else if (((ptarg==false)&&(p1code==kPdgPiP))
              || ((ptarg==true)&&(p1code==kPdgPiM)))
       {
         if (rand < xsec2pi0n) // pi0 & pi0 & n final state
           {
             p3code = (ptarg ? kPdgNeutron : kPdgProton);
             p4code = kPdgPi0;
             p5code = p4code;
           }
         else if (rand < (xsec2pi0n + xsecpippimn)) // pi+ & pi- & n final state
           {
             p3code = (ptarg ? kPdgNeutron : kPdgProton);
             p4code = p1code;
             p5code = ((p1code==kPdgPiP) ? kPdgPiM : kPdgPiP);
           }
         else // pi0 & pi- & p final state
           {
             p3code = (ptarg ? kPdgProton : kPdgNeutron);
             p4code = p1code;
             p5code = kPdgPi0;
           }
       }
     else if (p1code==kPdgPi0)
       {
         rand = rnd->RndFsi().Rndm();
         if (rand < 191./270.)
           {  // pi+ & pi- & p final state
             p3code = (ptarg ? kPdgProton : kPdgNeutron);
             p4code = kPdgPiP;
             p5code = kPdgPiM;
           }
         else if (rand < 7./135.)
           {  // pi0 & pi0 & p final state
             p3code = (ptarg ? kPdgProton : kPdgNeutron);
             p4code = kPdgPi0;
             p5code = p4code;
           }
         else
           {  // pi+ & pi0 & n final state
             p3code = (ptarg ? kPdgNeutron : kPdgProton);
             p4code = (ptarg ? kPdgPiP : kPdgPiM);
             p5code = kPdgPi0;
           }
       }
     else // unhandled
       {
         LOG("INukeUtils",pNOTICE) << "Pi production final state unable to be determined, picode, ptarg = " <<PDGLibrary::Instance()->Find(p1code)->GetName() << "  " << PDGLibrary::Instance()->Find(ptarg)->GetName();
         exceptions::INukeException exception;
         exception.SetReason("PionProduction final state not determined");
         throw exception;
         return false;
       }

   } else if(p1code==kPdgProton||p1code==kPdgNeutron) //nucleon probes
    {

      double tote = p->Energy();
      double pMass = pLib->Find(2212)->Mass();
      double nMass = pLib->Find(2112)->Mass();
      double etapp2ppPi0 =
        utils::intranuke2018::CalculateEta(pMass,tote,pMass,pMass+pMass,pLib->Find(111)->Mass());
      double etapp2pnPip =
        utils::intranuke2018::CalculateEta(pLib->Find(p1code)->Mass(),tote,((p1code==kPdgProton)?pMass:nMass),
                                       pMass+nMass,pLib->Find(211)->Mass());
      double etapn2nnPip =
        utils::intranuke2018::CalculateEta(pMass,tote,nMass,nMass+nMass,pLib->Find(211)->Mass());
      double etapn2ppPim =
        utils::intranuke2018::CalculateEta(pMass,tote,nMass,pMass+pMass,pLib->Find(211)->Mass());

      if ((etapp2ppPi0<=0.)&&(etapp2pnPip<=0.)&&(etapn2nnPip<=0.)&&(etapn2ppPim<=0.)) { // below threshold
        LOG("INukeUtils",pNOTICE) << "PionProduction() called below threshold energy";
        exceptions::INukeException exception;
        exception.SetReason("PionProduction final state not possible - below threshold");
        throw exception;
        return false;
        }

      // calculate cross sections
      double xsecppPi0=0,xsecpnPiP=0,xsecnnPiP=0,xsecppPiM=0;
      if (etapp2ppPi0>0){
      xsecppPi0 = 4511*(1-.91*TMath::Exp(-TMath::Power((etapp2ppPi0-.705),2)));
      xsecppPi0 *= (1-TMath::Exp(-TMath::Power((.556*etapp2ppPi0),3.5)));
      xsecppPi0 *= (1-TMath::Exp(-56.897/(etapp2ppPi0*etapp2ppPi0)));
      xsecppPi0 = TMath::Max(0.,xsecppPi0);}

      if (etapp2pnPip>0){
      xsecpnPiP = 18840*(1-.808*TMath::Exp(-TMath::Power((etapp2pnPip-.371),2)));
      xsecpnPiP *= (1-TMath::Exp(-TMath::Power((.568*etapp2pnPip),3.2)));
      xsecpnPiP *= (1-TMath::Exp(-39.818/(etapp2pnPip*etapp2pnPip)));
      xsecpnPiP = TMath::Max(0.,xsecpnPiP);}

      if (etapn2nnPip>0){
      xsecnnPiP = 7670*(1-.479*TMath::Exp(-TMath::Power((etapn2nnPip-.947),2)));
      xsecnnPiP *= (1-TMath::Exp(-TMath::Power((.35*etapn2nnPip),3.2)));
      xsecnnPiP *= (1-TMath::Exp(-71.53/(etapn2nnPip*etapn2nnPip)));
      xsecnnPiP = TMath::Max(0.,xsecnnPiP);}

      if (etapn2ppPim>0){
      xsecppPiM = 7670*(1-.479*TMath::Exp(-TMath::Power((etapn2ppPim-.947),2)));
      xsecppPiM *= (1-TMath::Exp(-TMath::Power((.35*etapn2ppPim),3.2)));
      xsecppPiM *= (1-TMath::Exp(-71.53/(etapn2ppPim*etapn2ppPim)));
      xsecppPiM = TMath::Max(0.,xsecppPiM);}

      // double sigma11 = xsecppPi0;
      double sigma10 = TMath::Max(0.,xsecpnPiP - xsecppPi0); // Fundamental cross sections
      double sigma01 = TMath::Max(0.,xsecppPiM + xsecnnPiP - xsecppPi0);

      double xsecpnPi0 = .5*(sigma10 + sigma01);
      xsecpnPi0 = TMath::Max(xsecpnPi0,0.);

      LOG("INukeUtils",pDEBUG) << '\n' << "Cross section values: "<<'\n'
                               << xsecppPi0 << " PP pi0"  <<'\n'
                               << xsecpnPiP << " PN pi+"  <<'\n'
                               << xsecnnPiP << " NN pi+"  <<'\n'
                               << xsecpnPi0 << " PN pi0";

      double xsecp=0,xsecn=0;
      switch (p1code) {
      case kPdgProton:  xsecp=xsecppPi0+xsecpnPiP; xsecn=xsecppPiM+xsecnnPiP+xsecpnPi0; break;
      case kPdgNeutron: xsecp=xsecppPiM+xsecnnPiP+xsecpnPi0; xsecn=xsecppPi0+xsecpnPiP; break;
      default:
        LOG("INukeUtils",pWARN) << "InelasticHN cannot handle probe: "
                                 << PDGLibrary::Instance()->Find(p1code)->GetName();
        return false;
        break;
      }

      // Normalize cross sections by Z or (A-Z)

      xsecp *= RemnZ;
      xsecn *= RemnA-RemnZ;

      // determine target

      double rand = rnd->RndFsi().Rndm() * (xsecp + xsecn);
      if (rand < xsecp) // proton target
        { rand /= RemnZ; ptarg = true;}
      else // neutron target
        { rand -= xsecp; rand /= RemnA-RemnZ; ptarg = false;}

      if(p1code==kPdgProton) // Cross sections not explicitly given are calculated from isospin relations
        {
          if(ptarg)
            {
              if   (rand<xsecppPi0) {p3code=kPdgProton; p4code=kPdgProton;  p5code=kPdgPi0;}
              else                  {p3code=kPdgProton; p4code=kPdgNeutron; p5code=kPdgPiP;}
            }
          else
            {
              if        (rand<xsecnnPiP)           {p3code=kPdgNeutron; p4code=kPdgNeutron; p5code=kPdgPiP;}
              else if   (rand<xsecppPiM+xsecnnPiP) {p3code=kPdgProton;  p4code=kPdgProton;  p5code=kPdgPiM;}
              else                                 {p3code=kPdgProton;  p4code=kPdgNeutron; p5code=kPdgPi0;}
            }
        }
      else if(p1code==kPdgNeutron)
        {
          if(ptarg)
            {
              if        (rand<xsecnnPiP)           {p3code=kPdgNeutron; p4code=kPdgNeutron; p5code=kPdgPiP;}
              else if   (rand<xsecppPiM+xsecnnPiP) {p3code=kPdgProton;  p4code=kPdgProton;  p5code=kPdgPiM;}
              else                                 {p3code=kPdgProton;  p4code=kPdgNeutron; p5code=kPdgPi0;}
            }
          else
            {
              if   (rand<xsecpnPiP) {p3code=kPdgNeutron; p4code=kPdgProton;  p5code=kPdgPiM;}
              else                  {p3code=kPdgNeutron; p4code=kPdgNeutron; p5code=kPdgPi0;}
            }
        }
    }
  else
    {
      LOG("INukeUtils",pWARN)
                << "Unable to handle probe (=" << p1code << ") in InelasticHN()";
      return false;
    }

   // determine if reaction is allowed
   if ( RemnA < 1 )
     {
       LOG("INukeUtils",pNOTICE) << "PionProduction() failed : no nucleons to produce pions";
       return false;
     }
   else if ( RemnZ + ((pcode==kPdgProton || pcode==kPdgPiP)?1:0) - ((pcode==kPdgPiM)?1:0)
             < ((p3code==kPdgProton || p3code==kPdgPiP)?1:0) - ((p3code==kPdgPiM)?1:0)
             + ((p4code==kPdgProton || p4code==kPdgPiP)?1:0) - ((p4code==kPdgPiM)?1:0)
             + ((p5code==kPdgProton || p5code==kPdgPiP)?1:0) - ((p5code==kPdgPiM)?1:0) )
     {
       LOG("INukeUtils",pNOTICE) << "PionProduction() failed : too few protons in nucleus";
       exceptions::INukeException exception;
       exception.SetReason("PionProduction fails - too few protons available");
       throw exception;
       return false;
     }

   s1->SetPdgCode(p3code);
   s2->SetPdgCode(p4code);
   s3->SetPdgCode(p5code);

   if(genie::utils::intranuke2018::ThreeBodyKinematics(
        ev,p,(ptarg?kPdgProton:kPdgNeutron),s1,s2,s3,DoFermi,FermiFac,FermiMomentum,Nuclmodel))
     {
       // okay, handle remnants and return true
       // assumes first particle is always the nucleon,
       //   second can be either nucleon or pion
       //   last always pion
       if (pcode==kPdgProton || pcode==kPdgPiP) RemnZ++;
       if (pcode==kPdgPiM) RemnZ--;
       if (pdg::IsPion(pcode)) RemnA--;
       if (pdg::IsProton(p3code)) RemnZ--;
       if (pdg::IsNeutronOrProton(p4code)) RemnA--;
       if (p4code==kPdgPiP || p4code==kPdgProton) RemnZ--;
       if (p4code==kPdgPiM) RemnZ++;
       if (p5code==kPdgPiP) RemnZ--;
       if (p5code==kPdgPiM) RemnZ++;

       LOG("INukeUtils",pDEBUG) << "Remnant (A,Z) = (" <<RemnA<<','<<RemnZ<<')';

       RemnP4 -= *s1->P4() + *s2->P4() + *s3->P4() - *p->P4();
       return true;
     }
   else {
     exceptions::INukeException exception;
     exception.SetReason("PionProduction final state not determined");
     throw exception;
     return false;
   }
}
//___________________________________________________________________________
double genie::utils::intranuke2018::CalculateEta(double Minc, double nrg, double Mtarg,
                               double Mtwopart, double Mpi)
{
  //Aaron Meyer (1/20/2010)

  //Used to calculate the maximum kinematically allowed CM frame pion momentum
  //  ke in MeV, eta normalized by pion mass
  //  approximated by taking two ejected nucleons to be one particle of the same mass
  //For pion cross sections, in utils::intranuke2018::PionProduction

  //LOG("INukeUtils",pDEBUG) << "Input values: "<<Minc<<' '<<nrg<<' '<<Mtarg<<' '<<Mtwopart<<' '<<Mpi;
  double Scm = Minc*Minc + Mtarg*Mtarg + 2*Mtarg*nrg;
  double eta = 0;
  eta= TMath::Power(Scm,2) + TMath::Power(Mtwopart,4) + TMath::Power(Mpi,4);
  eta-= 2*TMath::Power(Mtwopart*Mpi,2);
  eta-= 2*Scm*TMath::Power(Mtwopart,2);
  eta-= 2*Scm*TMath::Power(Mpi,2);
  eta = TMath::Power(eta/Scm,1./2.);
  eta/= (2*Mpi);

  eta=TMath::Max(eta,0.);
  return eta;
}
//___________________________________________________________________________
// Generic Phase Space Decay methods
//___________________________________________________________________________
bool genie::utils::intranuke2018::PhaseSpaceDecay(
  GHepRecord* ev, GHepParticle* p, const PDGCodeList & pdgv,
  TLorentzVector &RemnP4, double NucRmvE, EINukeMode mode)
{
// General method decaying the input particle system 'pdgv' with available 4-p
// given by 'pd'. The decayed system is used to populate the input TMCParticle
// array starting from the slot 'offset'.
//
  LOG("INukeUtils", pINFO) << "*** Performing a Phase Space Decay";
  assert(pdgv.size() > 1);

  LOG("INukeUtils",pINFO) << "probe mass: M = " << p->Mass();

  // Get the decay product masses & names

  ostringstream state_sstream;
  state_sstream << "( ";
  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   mass_sum = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m  = PDGLibrary::Instance()->Find(pdgc)->Mass();
    string nm = PDGLibrary::Instance()->Find(pdgc)->GetName();
    mass[i++] = m;
    mass_sum += m;
    state_sstream << nm << " ";
  }
  state_sstream << ")";

  TLorentzVector * pd = p->GetP4(); // incident particle 4p

  bool is_nuc  = pdg::IsNeutronOrProton(p->Pdg());
  bool is_kaon = p->Pdg()==kPdgKP  || p->Pdg()==kPdgKM;
  // not used // bool is_pion = p->Pdg()==kPdgPiP || p->Pdg()==kPdgPi0 || p->Pdg()==kPdgPiM;
  // update available energy -> init (mass + kinetic) + sum of f/s masses
  // for pion only.  Probe mass not available for nucleon, kaon
  double availE = pd->Energy() + mass_sum;
  if(is_nuc||is_kaon) availE -= p->Mass();
  pd->SetE(availE);

  LOG("INukeUtils",pNOTICE)
    << "size, mass_sum, availE, pd mass, energy = " << pdgv.size() << "  "
    << mass_sum << "  " << availE << "  " << p->Mass() << "  " << p->Energy() ;

  // compute the 4p transfer to the hadronic blob
  double dE = mass_sum;
  if(is_nuc||is_kaon) dE -= p->Mass();
  TLorentzVector premnsub(0,0,0,dE);
  RemnP4 -= premnsub;

  LOG("INukeUtils", pINFO)
    << "Final state = " << state_sstream.str() << " has N = " << pdgv.size()
    << " particles / total mass = " << mass_sum;
  LOG("INukeUtils", pINFO)
    << "Composite system p4 = " << utils::print::P4AsString(pd);

  // Set the decay
  TGenPhaseSpace GenPhaseSpace;
  bool permitted = GenPhaseSpace.SetDecay(*pd, pdgv.size(), mass);
  if(!permitted) {
     LOG("INukeUtils", pERROR)
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << mass_sum << "\n"
       << " Decaying system p4 = " << utils::print::P4AsString(pd);

     // clean-up and return
     RemnP4 += premnsub;
     delete [] mass;
     delete pd;
     return false;
  }

  // The decay is permitted - add the incident particle at the event record
  // and mark is as 'Nucleon Cluster Target' (used to be confusing 'Decayed State')
  p->SetStatus(kIStNucleonClusterTarget);  //kIStDecayedState);
  p->SetPdgCode(kPdgCompNuclCluster);
  ev->AddParticle(*p);
  // Get the maximum weight
  double wmax = -1;
  for(int k=0; k<200; k++) {
     double w = GenPhaseSpace.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);

  LOG("INukeUtils", pINFO)
   << "Max phase space gen. weight @ current hadronic interaction: " << wmax;

  // Generate an unweighted decay

  RandomGen * rnd = RandomGen::Instance();
  wmax *= 1.2;

  bool accept_decay=false;
  unsigned int itry=0;

  while(!accept_decay)
  {
    itry++;

    if(itry>kMaxUnweightDecayIterations) {
       // report, clean-up and return
       LOG("INukeUtils", pNOTICE)
             << "Couldn't generate an unweighted phase space decay after "
             << itry << " attempts";
       delete [] mass;
       delete pd;
       return false;
    }

    double w  = GenPhaseSpace.Generate();
    double gw = wmax * rnd->RndFsi().Rndm();

    if(w > wmax) {
       LOG("INukeUtils", pNOTICE)
           << "Decay weight = " << w << " > max decay weight = " << wmax;
    }

    LOG("INukeUtils", pNOTICE) << "Decay weight = " << w << " / R = " << gw;
    accept_decay = (gw<=w);
  }

  // Insert final state products into the event record
  // - the particles are added as daughters of the decayed state
  // - the particles are marked as final stable state (in hA mode)
  i=0;
  int mom = ev->ParticlePosition(p);
  LOG("INukeUtils", pNOTICE) << "mother index = " << mom;
  GHepStatus_t ist = kIStStableFinalState;
  GHepStatus_t ist_pi = kIStHadronInTheNucleus;

  TLorentzVector * v4 = p->GetX4();

  double checkpx = p->Px();
  double checkpy = p->Py();
  double checkpz = p->Pz();
  double checkE = p->E();

  //LOG("INukeUtils", PNOTICE)

  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {

     //-- current PDG code
     int pdgc = *pdg_iter;
     bool isnuc = pdg::IsNeutronOrProton(pdgc);

     //-- get the 4-momentum of the i-th final state particle
     TLorentzVector * p4fin = GenPhaseSpace.GetDecay(i++);

     //-- intranuke no longer throws "bindinos" but adds all the energy
     //   not going at a simulated f/s particle at a "hadronic blob"
     //   representing the remnant system: do the binding energy subtraction
     //   here & update the remnant hadronic system 4p
     double M  = PDGLibrary::Instance()->Find(pdgc)->Mass();
     double En = p4fin->Energy();

     double KE = En-M;

     //double KE = En;
     //if(is_pion) KE -= M;

     double dE_leftover = TMath::Min(NucRmvE, KE);
     KE -= dE_leftover;
     En  = KE+M;
     double pmag_old = p4fin->P();
     double pmag_new = TMath::Sqrt(TMath::Max(0.,En*En-M*M));
     double scale    = pmag_new / pmag_old;
     double pxn      = scale * p4fin->Px();
     double pyn      = scale * p4fin->Py();
     double pzn      = scale * p4fin->Pz();

     TLorentzVector p4n(pxn,pyn,pzn,En);
     //     LOG("INukeUtils", pNOTICE) << "Px = " << pxn << " Py = " << pyn
     //					<< " Pz = " << pzn << " E = " << KE;
     checkpx -= pxn;
     checkpy -= pyn;
     checkpz -= pzn;
     checkE -= KE;

     if (mode==kIMdHA &&
         (pdgc==kPdgPiP || pdgc==kPdgPi0 || pdgc==kPdgPiM) )
       {
         if (p4n.Vect().Mag()>=0.001)
           {
             GHepParticle new_particle(pdgc, ist_pi, mom,-1,-1,-1, p4n, *v4);
             ev->AddParticle(new_particle);
           }
         else
           {
             // Momentum too small, assign a non-zero momentum to the particle
             // Conserve momentum with the remnant nucleus

             LOG("INukeUtils", pINFO)<<"Momentum too small; assigning 0.001 as new momentum";

             double phi = 2*kPi*rnd->RndFsi().Rndm();
             double omega = 2*rnd->RndFsi().Rndm();
               // throw number against solid angle for uniform distribution

             double E4n = TMath::Sqrt(0.001*0.001+M*M);
             p4n.SetPxPyPzE(0.001,0,0,E4n);
             p4n.Rotate(TMath::ACos(1-omega),TVector3(0,0,1));
             p4n.Rotate(phi,TVector3(1,0,0));

             RemnP4 -= (p4n - TLorentzVector(0,0,0,M));

             GHepParticle new_particle(pdgc, ist, mom,-1,-1,-1, p4n, *v4);
             ev->AddParticle(new_particle);
           }
       }
     else
       {
         GHepParticle new_particle(pdgc, ist, mom,-1,-1,-1, p4n, *v4);

         if(isnuc) new_particle.SetRemovalEnergy(0.);
         ev->AddParticle(new_particle);
       }

     double dpx = (1-scale)*p4fin->Px();
     double dpy = (1-scale)*p4fin->Py();
     double dpz = (1-scale)*p4fin->Pz();
     TLorentzVector premnadd(dpx,dpy,dpz,dE_leftover);
     RemnP4 += premnadd;
  }
  //LOG("INukeUtils", pNOTICE) << "TEST: " << p->Mass();
  LOG("INukeUtils", pNOTICE) << "check conservation: Px = " << checkpx << " Py = " << checkpy
                             << " Pz = " << checkpz << " E = " << checkE;

  // Clean-up
  delete [] mass;
  delete pd;
  delete v4;

  return true;
}

double genie::utils::intranuke2018::sigmaTotalOset (
                                    const double &pionKineticEnergy,
                                    const double &density,
                                    const int    &pionPDG,
                                    const double &protonFraction,
                                    const bool   &isTableChosen
                                    )
{
  // ------ OsetCrossSection init (only first time function is called) ------ //
  static INukeOset *iNukeOset = NULL;

  if (iNukeOset == NULL)
  {
    if (isTableChosen)
    {
      // set directory with data on first call
      static const std::string dataDir = (gSystem->Getenv("GINUKEHADRONDATA")) ?
                                   string(gSystem->Getenv("GINUKEHADRONDATA")) :
                                              string(gSystem->Getenv("GENIE")) +
                                              string("/data/evgen/intranuke/");
      // set file with Oset table on first call
      static const std::string dataFile = dataDir + "tot_xsec/"
                                          "intranuke-xsection-pi+N-Oset.dat";
      // initialize OsetCrossSection on first call
      iNukeOset = new INukeOsetTable (dataFile.c_str());
    }
    else iNukeOset = new INukeOsetFormula();
  }
  // ------ OsetCrossSection init (only first time function is called) ------ //

  // set up Oset class (assign pion Tk, nuclear density etc)
  iNukeOset->setupOset (density, pionKineticEnergy, pionPDG, protonFraction);

  return iNukeOset->getTotalCrossSection();

}
