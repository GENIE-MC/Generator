
//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
         Aaron Meyer <asm58@pitt.edu>, Pittsburgh Univ.
	 Alex Bell, Pittsburgh Univ.
         Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, Rutherford Lab.
         September 20, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 30, 2007 - SD
   Changed the hadron tracking algorithm to take into account the radial
   nuclear density dependence. Using the somewhat empirical approach of
   increasing the nuclear radius by a const (tunable) number times the tracked 
   particle's de Broglie wavelength as this helps getting the hadron+nucleus 
   cross sections right.
 @ Mar 08, 2008 - CA
   Fixed code retrieving the remnant nucleus which stopped working as soon as
   simulation of nuclear de-excitation started pushing photons in the target
   nucleus daughter list.
 @ Jun 20, 2008 - CA
   Fix a mem leak: The (clone of the) GHepParticle being re-scattered was not 
   deleted after it was added at the GHEP event record.
 @ Jul 15, 2010 - AM
   The hN mode is now implemented in Intranuke. Similar to hA mode, but particles
   produced by reactions are stepped through the nucleus like probe particles.
   Particles react with nucleons instead of the entire nucleus, and final states
   are determined after reactions are finished, not before.
 @ Dec 15, 2014 - SD, Nick Geary 
   Update fates to include Compound Nucleus final state correctly.
 @ Jan 9, 2015 - SD, NG, Tomek Golan
   Added 2014 version of INTRANUKE codes (new class) for independent development.
 @ Oct, 2015 - TG
   Added 2015 version of INTRANUKE codes (new class) for independent development.  Include Oset model for medium corrections to piA for Tpi<350 MeV.
 @ May, 2015 Flor Blaszczyk
   K+ are now handled.
 @ July, 2016 Nicholas Suarez, Josh Kleckner, SD
   fix memory leak, fix fates, improve NNCorr binning
 & Mar, 2018  Nicholas Suarez, SD
   add compound nucleus option to populate KE<30 MeV
*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/HadronTransport/Intranuke2018.h"
#include "Physics/HadronTransport/HNIntranuke2018.h"
#include "Physics/HadronTransport/INukeException.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Physics/HadronTransport/INukeUtils2018.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/HadronTransport/INukeOset.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::utils::intranuke2018;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
//___________________________________________________________________________
// Methods specific to INTRANUKE's HN-mode
//___________________________________________________________________________
//___________________________________________________________________________
HNIntranuke2018::HNIntranuke2018() :
Intranuke2018("genie::HNIntranuke2018")
{

}
//___________________________________________________________________________
HNIntranuke2018::HNIntranuke2018(string config) :
Intranuke2018("genie::HNIntranuke2018",config)
{

}
//___________________________________________________________________________
HNIntranuke2018::~HNIntranuke2018()
{

}
//___________________________________________________________________________
void HNIntranuke2018::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("HNIntranuke2018", pNOTICE) 
     << "************ Running hN2018 MODE INTRANUKE ************";
     
  /*  LOG("HNIntranuke2018", pWARN) 
     << print::PrintFramedMesg(
         "Experimental code (INTRANUKE/hN model) - Run at your own risk");
  */

  Intranuke2018::ProcessEventRecord(evrec);

  LOG("HNIntranuke2018", pINFO) << "Done with this event";
}
//___________________________________________________________________________
void HNIntranuke2018::SimulateHadronicFinalState(GHepRecord* ev, GHepParticle* p) const
{
// Simulate a hadron interaction for the input particle p in HN mode
//
  if(!p || !ev)
    {
      LOG("HNIntranuke2018", pERROR) << "** Null input!";
      return;
    }

  // check particle id
  int pdgc = p->Pdg();
  bool is_pion    = (pdgc==kPdgPiP || pdgc==kPdgPiM || pdgc==kPdgPi0);
  bool is_kaon    = (pdgc==kPdgKP);
  bool is_baryon  = (pdgc==kPdgProton || pdgc==kPdgNeutron);
  bool is_gamma   = (pdgc==kPdgGamma);										
  if(!(is_pion || is_baryon  || is_gamma || is_kaon))
    {
      LOG("HNIntranuke2018", pERROR) << "** Cannot handle particle: " << p->Name();
      return;
    }
  try
    {
      // select a fate for the input particle
      INukeFateHN_t fate = this->HadronFateHN(p);

      // store the fate
      ev->Particle(p->FirstMother())->SetRescatterCode((int)fate);

      if(fate == kIHNFtUndefined)
	{
	  LOG("HNIntranuke2018", pERROR) << "** Couldn't select a fate";
	  LOG("HNIntranuke2018", pERROR) << "** Num Protons: " << fRemnZ 
				     << ",  Num Neutrons: "<<(fRemnA-fRemnZ);
	  LOG("HNIntranuke2018", pERROR) << "** Particle: " << "\n" << (*p);
	  //LOG("HNIntranuke2018", pERROR) << "** Event Record: " << "\n" << (*ev);
	  //p->SetStatus(kIStUndefined);
	  p->SetStatus(kIStStableFinalState);
	  ev->AddParticle(*p);
	  return;
	}

      LOG("HNIntranuke2018", pNOTICE)
	<< "Selected " << p->Name() << " fate: " << INukeHadroFates::AsString(fate);

      // handle the reaction
      if(fate == kIHNFtCEx || fate == kIHNFtElas)
	{
	  this->ElasHN(ev,p,fate);
	}
      else if(fate == kIHNFtAbs)                         {this-> AbsorbHN(ev,p,fate);}
      else if(fate == kIHNFtInelas && pdgc != kPdgGamma) 
	{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
	  LOG("HNIntranuke2018", pDEBUG)
	    << "Invoking InelasticHN() for a : " << p->Name()
	    << " whose fate is : " << INukeHadroFates::AsString(fate);
#endif

	  this-> InelasticHN(ev,p);
	}
      else if(fate == kIHNFtInelas && pdgc == kPdgGamma) {this-> GammaInelasticHN(ev,p,fate);}
      else if(fate == kIHNFtCmp){utils::intranuke2018::PreEquilibrium(ev,p,fRemnA,fRemnZ,fRemnP4,fDoFermi,fFermiFac,fNuclmodel,fNucRmvE,kIMdHN);}
      else if(fate == kIHNFtNoInteraction)
	{
	  p->SetStatus(kIStStableFinalState);
	  ev->AddParticle(*p);
	  return;
	}
    }
  catch(exceptions::INukeException exception)
    {
      this->SimulateHadronicFinalState(ev,p);
       LOG("HNIntranuke2018", pNOTICE) 
         << "retry call to SimulateHadronicFinalState ";
       LOG("HNIntranuke2018", pNOTICE) << exception;

    }
}
//___________________________________________________________________________
INukeFateHN_t HNIntranuke2018::HadronFateHN(const GHepParticle * p) const
{
// Select a hadron fate in HN mode
//
  RandomGen * rnd = RandomGen::Instance();

  // get pdgc code & kinetic energy in MeV
  int    pdgc = p->Pdg();
  double ke   = p->KinE() / units::MeV;

  bool isPion = (pdgc == kPdgPiP or pdgc == kPdgPi0 or pdgc == kPdgPiM);

  if (isPion and fUseOset and ke < 350.0) return HadronFateOset ();
 
  LOG("HNIntranuke2018", pNOTICE) 
   << "Selecting hN fate for " << p->Name() << " with KE = " << ke << " MeV";

   // try to generate a hadron fate
  unsigned int iter = 0;
  while(iter++ < kRjMaxIterations) {

    // handle pions
    //
    if (pdgc==kPdgPiP || pdgc==kPdgPiM || pdgc==kPdgPi0) {

       double frac_cex      = this->FateWeight(pdgc, kIHNFtCEx)
	                            * fHadroData2018->Frac(pdgc, kIHNFtCEx,     ke, fRemnA, fRemnZ);
       double frac_elas     = this->FateWeight(pdgc, kIHNFtElas)
	                            * fHadroData2018->Frac(pdgc, kIHNFtElas,    ke, fRemnA, fRemnZ);
       double frac_inel     = this->FateWeight(pdgc, kIHNFtInelas)
	                            * fHadroData2018->Frac(pdgc, kIHNFtInelas,  ke, fRemnA, fRemnZ);
       double frac_abs      = this->FateWeight(pdgc, kIHNFtAbs)
	                            * fHadroData2018->Frac(pdgc, kIHNFtAbs,     ke, fRemnA, fRemnZ);

       frac_cex     *= fNucCEXFac;    // scaling factors
       frac_abs     *= fNucAbsFac;
       frac_elas    *= fNucQEFac;
       if(pdgc==kPdgPi0) frac_abs*= 0.665;  //isospin factor

       LOG("HNIntranuke2018", pNOTICE) 
	 << "\n frac{" << INukeHadroFates::AsString(kIHNFtCEx)     << "} = " << frac_cex
	 << "\n frac{" << INukeHadroFates::AsString(kIHNFtElas)    << "} = " << frac_elas
	 << "\n frac{" << INukeHadroFates::AsString(kIHNFtInelas)  << "} = " << frac_inel
	 << "\n frac{" << INukeHadroFates::AsString(kIHNFtAbs)     << "} = " << frac_abs;

       // compute total fraction (can be <1 if fates have been switched off)
       double tf = frac_cex      +
                   frac_elas     +
                   frac_inel     +  
                   frac_abs;

       double r = tf * rnd->RndFsi().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("HNIntranuke2018", pDEBUG) << "r = " << r << " (max = " << tf << ")";
#endif

       double cf=0; // current fraction

       if(r < (cf += frac_cex     )) return kIHNFtCEx;    //cex
       if(r < (cf += frac_elas    )) return kIHNFtElas;   //elas
       if(r < (cf += frac_inel    )) return kIHNFtInelas; //inelas
       if(r < (cf += frac_abs     )) return kIHNFtAbs;    //abs   

       LOG("HNIntranuke2018", pWARN) 
         << "No selection after going through all fates! " 
                     << "Total fraction = " << tf << " (r = " << r << ")";
       ////////////////////////////
       return kIHNFtUndefined;
    }

    // handle nucleons
    else if (pdgc==kPdgProton || pdgc==kPdgNeutron) {

      double frac_elas     = this->FateWeight(pdgc, kIHNFtElas)
	                           * fHadroData2018->Frac(pdgc, kIHNFtElas,   ke, fRemnA, fRemnZ);
      double frac_inel     = this->FateWeight(pdgc, kIHNFtInelas)
	                           * fHadroData2018->Frac(pdgc, kIHNFtInelas, ke, fRemnA, fRemnZ);
      double frac_cmp      = this->FateWeight(pdgc, kIHNFtCmp)
	                           * fHadroData2018->Frac(pdgc, kIHNFtCmp,    ke, fRemnA , fRemnZ);

      LOG("HNIntranuke2018", pINFO) 
	<< "\n frac{" << INukeHadroFates::AsString(kIHNFtElas)    << "} = " << frac_elas
	<< "\n frac{" << INukeHadroFates::AsString(kIHNFtInelas)  << "} = " << frac_inel;

       // compute total fraction (can be <1 if fates have been switched off)
       double tf = frac_elas     +
                   frac_inel     +
	           frac_cmp;

       double r = tf * rnd->RndFsi().Rndm();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("HNIntranuke2018", pDEBUG) << "r = " << r << " (max = " << tf << ")";
#endif

       double cf=0; // current fraction
       if(r < (cf += frac_elas    )) return kIHNFtElas;    // elas
       if(r < (cf += frac_inel    )) return kIHNFtInelas;  // inelas
       if(r < (cf += frac_cmp     )) return kIHNFtCmp;     // cmp

       LOG("HNIntranuke2018", pWARN) 
         << "No selection after going through all fates! "
                        << "Total fraction = " << tf << " (r = " << r << ")";
       //////////////////////////
       return kIHNFtUndefined;
    }

    // handle gamma -- does not currently consider the elastic case 
    else if (pdgc==kPdgGamma)  return kIHNFtInelas;
    // Handle kaon -- elastic + charge exchange
    else if (pdgc==kPdgKP){
       double frac_cex      = this->FateWeight(pdgc, kIHNFtCEx)
	                            * fHadroData2018->Frac(pdgc, kIHNFtCEx,     ke, fRemnA, fRemnZ);
       double frac_elas     = this->FateWeight(pdgc, kIHNFtElas)
	                            * fHadroData2018->Frac(pdgc, kIHNFtElas,    ke, fRemnA, fRemnZ);

       //       frac_cex     *= fNucCEXFac;    // scaling factors
       //       frac_elas    *= fNucQEFac;   // Flor - Correct scaling factors?

       LOG("HNIntranuke", pINFO) 
          << "\n frac{" << INukeHadroFates::AsString(kIHNFtCEx)     << "} = " << frac_cex
          << "\n frac{" << INukeHadroFates::AsString(kIHNFtElas)    << "} = " << frac_elas;

       // compute total fraction (can be <1 if fates have been switched off)
       double tf = frac_cex      +
                   frac_elas;

       double r = tf * rnd->RndFsi().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("HNIntranuke", pDEBUG) << "r = " << r << " (max = " << tf << ")";
#endif

       double cf=0; // current fraction

       if(r < (cf += frac_cex     )) return kIHNFtCEx;    //cex
       if(r < (cf += frac_elas    )) return kIHNFtElas;   //elas  

       LOG("HNIntranuke", pWARN) 
         << "No selection after going through all fates! " 
                     << "Total fraction = " << tf << " (r = " << r << ")";
       ////////////////////////////
       return kIHNFtUndefined;
    }//End K+

    /*else if (pdgc==kPdgKM){

       return kIHNFtElas;
    }//End K-*/

  }//iterations

  return kIHNFtUndefined;
}

//___________________________________________________________________________
double HNIntranuke2018::FateWeight(int pdgc, INukeFateHN_t fate) const
{
  // turn fates off if the remnant nucleus does not have the number of p,n
  // required

  int np = fRemnZ;
  int nn = fRemnA - fRemnZ;
 
  if (np < 1 && nn < 1)
    {
      LOG("HNIntranuke2018", pERROR) << "** Nothing left in nucleus!!! **";
      return 0;
    }
  else
    {
      if (fate == kIHNFtCEx && pdgc==kPdgPiP ) { return (nn>=1) ? 1. : 0.; }
      if (fate == kIHNFtCEx && pdgc==kPdgPiM ) { return (np>=1) ? 1. : 0.; }
      if (fate == kIHNFtCEx && pdgc==kPdgKP  ) { return (nn>=1) ? 1. : 0.; } //Added, changed np to nn
      if (fate == kIHNFtAbs)      { return ((nn>=1) && (np>=1)) ? 1. : 0.; }
      if (fate == kIHNFtCmp )     { return ((pdgc==kPdgProton||pdgc==kPdgNeutron)&&fDoCompoundNucleus&&fRemnA>5) ? 1. : 0.; }

    }
  return 1.;
}
//___________________________________________________________________________
void HNIntranuke2018::AbsorbHN(
    GHepRecord * ev, GHepParticle * p, INukeFateHN_t fate) const
{
  // handles pi+d->2p, pi-d->nn, pi0 d->pn absorbtion, all using pi+d values
  
  int pdgc = p->Pdg();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HNIntranuke2018", pDEBUG)
    << "AbsorbHN() is invoked for a : " << p->Name()
    << " whose fate is : " << INukeHadroFates::AsString(fate);
#endif

  // check fate
  if(fate!=kIHNFtAbs)
    {
      LOG("HNIntranuke2018", pWARN)
        << "AbsorbHN() cannot handle fate: " << INukeHadroFates::AsString(fate);
      return;
    }

  // random number generator
  RandomGen * rnd = RandomGen::Instance();

  // Notes on the kinematics
  // -- Simple variables are used for efficiency
  // -- Variables are numbered according to particle
  // -- -- #1 -> incoming particle
  // -- -- #2 -> target (here, 2_1 and 2_2 for individual particles)
  // -- -- #3 -> scattered incoming (Particle tracked in hA mode)
  // -- -- #4 -> other scattered particle
  // -- Suffix "L" is for lab frame, suffix "CM" is for center of mass frame
  // -- Subscript "z" is for parallel component, "t" is for transverse

  int pcode, t1code, t2code, scode, s2code; // particles
  double M1, M2_1, M2_2, M3, M4;        // rest energies, in GeV
  double E1L, P1L, E2L, P2L, E3L, P3L, E4L, P4L;
  double P1zL, P2zL;
  double beta, gm; // speed and gamma for CM frame in lab
  double Et, E2CM;
  double C3CM, S3CM;  // cos and sin of scattering angle
  double Theta1, Theta2, theta5;
  double PHI3;        // transverse scattering angle
  double E1CM, E3CM, E4CM, P3CM;
  double P3zL, P3tL, P4zL, P4tL;
  double E2_1L, E2_2L;
  TVector3 tP2_1L, tP2_2L, tP1L, tP2L, tPtot, P1zCM, P2zCM;
  TVector3 tP3L, tP4L;
  TVector3 bDir, tTrans, tbeta, tVect;

  // Library instance for reference
  PDGLibrary * pLib = PDGLibrary::Instance();
 
  // Handle fermi target
  Target target(ev->TargetNucleus()->Pdg());

  // Target should be a deuteron, but for now
  // handling it as seperate nucleons
  if(pdgc==211) // pi-plus
    {
      pcode  = 211;
      t1code = 2212; // proton
      t2code = 2112; // neutron
      scode  = 2212;
      s2code = 2212;
    }
  else if(pdgc==-211) // pi-minus
    {
      pcode  = -211;
      t1code = 2212;
      t2code = 2112;
      scode  = 2112;
      s2code = 2112;
    }
  else if(pdgc==111) // pi-zero
    {
      pcode  = 111;
      t1code = 2212;
      t2code = 2112;
      scode  = 2212;
      s2code = 2112;
    }
  else
    {
      LOG("HNIntranuke2018", pWARN)
        << "AbsorbHN() cannot handle probe: " << pdgc;
      return;
    }
 
  // assign proper masses
  M1   = pLib->Find(pcode) ->Mass();
  M2_1 = pLib->Find(t1code)->Mass();
  M2_2 = pLib->Find(t2code)->Mass();
  M3   = pLib->Find(scode) ->Mass();
  M4   = pLib->Find(s2code)->Mass();

  // handle fermi momentum 
  if(fDoFermi)
    {
      target.SetHitNucPdg(t1code);
      fNuclmodel->GenerateNucleon(target);
      tP2_1L=fFermiFac * fNuclmodel->Momentum3();
      E2_1L = TMath::Sqrt(tP2_1L.Mag2() + M2_1*M2_1);
 
      target.SetHitNucPdg(t2code);
      fNuclmodel->GenerateNucleon(target);
      tP2_2L=fFermiFac * fNuclmodel->Momentum3();
      E2_2L = TMath::Sqrt(tP2_2L.Mag2() + M2_2*M2_2);
    }
  else
    {
      tP2_1L.SetXYZ(0.0, 0.0, 0.0);
      E2_1L = M2_1;

      tP2_2L.SetXYZ(0.0, 0.0, 0.0);
      E2_2L = M2_2;
    }

  E2L = E2_1L + E2_2L;

  // adjust p to reflect scattering
  // get random scattering angle
  C3CM = fHadroData2018->IntBounce(p,t1code,scode,fate);
    if (C3CM<-1.) 
    {
      p->SetStatus(kIStStableFinalState);
      ev->AddParticle(*p);
      return;
    }
  S3CM = TMath::Sqrt(1.0 - C3CM*C3CM);

  // Get lab energy and momenta
  E1L = p->E();
  if(E1L<0.001) E1L=0.001;
  P1L = TMath::Sqrt(E1L*E1L - M1*M1);
  tP1L = p->P4()->Vect();
  tP2L = tP2_1L + tP2_2L;
  P2L = tP2L.Mag();
  tPtot = tP1L + tP2L;

  // get unit vectors and angles needed for later
  bDir = tPtot.Unit();
  Theta1 = tP1L.Angle(bDir);
  Theta2 = tP2L.Angle(bDir);

  // get parallel and transverse components
  P1zL = P1L*TMath::Cos(Theta1);
  P2zL = P2L*TMath::Cos(Theta2);
  tVect.SetXYZ(1,0,0);
  if(TMath::Abs((tVect - bDir).Mag())<.01) tVect.SetXYZ(0,1,0);
  theta5 = tVect.Angle(bDir);
  tTrans = (tVect - TMath::Cos(theta5)*bDir).Unit();

  // calculate beta and gamma
  tbeta = tPtot * (1.0 / (E1L + E2L));
  beta = tbeta.Mag();
  gm = 1.0 / TMath::Sqrt(1.0 - beta*beta);

  // boost to CM frame to get scattered particle momenta
  E1CM = gm*E1L - gm*beta*P1zL;
  P1zCM = gm*P1zL*bDir - gm*tbeta*E1L;
  E2CM = gm*E2L - gm*beta*P2zL;
  P2zCM = gm*P2zL*bDir - gm*tbeta*E2L;
  Et = E1CM + E2CM;
  E3CM = (Et*Et + (M3*M3) - (M4*M4)) / (2.0*Et);
  E4CM = Et - E3CM;
  P3CM = TMath::Sqrt(E3CM*E3CM - M3*M3);

  // boost back to lab
  P3zL = gm*beta*E3CM + gm*P3CM*C3CM;
  P3tL = P3CM*S3CM;
  P4zL = gm*beta*E4CM + gm*P3CM*(-C3CM);
  P4tL = P3CM*(-S3CM);

  P3L = TMath::Sqrt(P3zL*P3zL + P3tL*P3tL);
  P4L = TMath::Sqrt(P4zL*P4zL + P4tL*P4tL);

  // check for too small values
  // may introduce error, so warn if it occurs
  if(!(TMath::Finite(P3L))||P3L<.001)
    {
      LOG("HNIntranuke2018",pINFO)
        << "Particle 3 " << M3 << " momentum small or non-finite: " << P3L
        << "\n" << "--> Assigning .001 as new momentum";
      P3tL = 0;
      P3zL = .001;
      P3L  = .001;
      E3L  = TMath::Sqrt(P3L*P3L + M3*M3);
    }

  if(!(TMath::Finite(P4L))||P4L<.001)
    {
      LOG("HNIntranuke2018",pINFO)
        << "Particle 4 " << M4 << " momentum small or non-finite: " << P4L
        << "\n" << "--> Assigning .001 as new momentum";
      P4tL = 0;
      P4zL = .001;
      P4L  = .001;
      E4L  = TMath::Sqrt(P4L*P4L + M4*M4);
    }

  // pauli blocking (do not apply PB for Oset)
  //if(!fUseOset && (P3L < fFermiMomentum || P4L < fFermiMomentum))
  double ke   = p->KinE() / units::MeV;
  if((!fUseOset || ke > 350.0 ) && (P3L < fFermiMomentum || P4L < fFermiMomentum))
    {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("HNIntranuke2018",pINFO) << "AbsorbHN failed: Pauli blocking";
#endif
      /*
      p->SetStatus(kIStHadronInTheNucleus);
      //disable until needed
      //      utils::intranuke2018::StepParticle(p,fFreeStep,fTrackingRadius);
      ev->AddParticle(*p);   
      return;
      */
      // new attempt at error handling:
      LOG("HNIntranuke2018", pINFO) << "AbsorbHN failed: Pauli blocking";
      exceptions::INukeException exception;
      exception.SetReason("hN absorption failed");
      throw exception;
    }

  // handle remnant nucleus updates
  fRemnZ--;
  fRemnA -=2;
  fRemnP4 -= TLorentzVector(tP2_1L,E2_1L);
  fRemnP4 -= TLorentzVector(tP2_2L,E2_2L);

  // get random phi angle, distributed uniformally in 360 deg
  PHI3 = 2 * kPi * rnd->RndFsi().Rndm();
  
  tP3L = P3zL*bDir + P3tL*tTrans;
  tP4L = P4zL*bDir + P4tL*tTrans;

  tP3L.Rotate(PHI3,bDir);  // randomize transverse components
  tP4L.Rotate(PHI3,bDir); 

  E3L = TMath::Sqrt(P3L*P3L + M3*M3);
  E4L = TMath::Sqrt(P4L*P4L + M4*M4);

  // create t particle w/ appropriate momenta, code, and status
  // set target's mom to be the mom of the hadron that was cloned
  GHepParticle * t = new GHepParticle(*p);
  t->SetFirstMother(p->FirstMother());
  t->SetLastMother(p->LastMother());

  TLorentzVector t4P4L(tP4L,E4L);
  t->SetPdgCode(s2code);
  t->SetMomentum(t4P4L);
  t->SetStatus(kIStHadronInTheNucleus);

  // adjust p to reflect scattering
  TLorentzVector t4P3L(tP3L,E3L);
  p->SetPdgCode(scode);
  p->SetMomentum(t4P3L);
  p->SetStatus(kIStHadronInTheNucleus);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HNIntranuke2018", pDEBUG)
    << "|p3| = " << (P3L) << ", E3 = " << (E3L);
  LOG("HNIntranuke2018", pDEBUG)
    << "|p4| = " << (P4L) << ", E4 = " << (E4L);
#endif

  ev->AddParticle(*p);
  ev->AddParticle(*t);

  delete t; // delete particle clone
}
//___________________________________________________________________________
void HNIntranuke2018::ElasHN(
	 GHepRecord * ev, GHepParticle * p, INukeFateHN_t fate) const
{
  // scatters particle within nucleus

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HNIntranuke2018", pDEBUG)
    << "ElasHN() is invoked for a : " << p->Name()
    << " whose fate is : " << INukeHadroFates::AsString(fate);
#endif

  if(fate!=kIHNFtCEx && fate!=kIHNFtElas)
    {
      LOG("HNIntranuke2018", pWARN)
	<< "ElasHN() cannot handle fate: " << INukeHadroFates::AsString(fate);
      return;
    }

  // Random number generator
  RandomGen * rnd = RandomGen::Instance();

  // vars for incoming particle, target, and scattered pdg codes
  int pcode = p->Pdg();
  int tcode, scode, s2code;
  double ppcnt = (double) fRemnZ / (double) fRemnA; // % of protons

  // Select a target randomly, weighted to #
  // -- Unless, of course, the fate is CEx,
  // -- in which case the target may be deterministic
  // Also assign scattered particle code
  if(fate==kIHNFtCEx)
    {
      if(pcode==kPdgPiP)      {tcode = kPdgNeutron; scode = kPdgPi0; s2code = kPdgProton;}
      else if(pcode==kPdgPiM) {tcode = kPdgProton;  scode = kPdgPi0; s2code = kPdgNeutron;}
      else if(pcode==kPdgKP)  {tcode = kPdgNeutron; scode = kPdgK0; s2code = kPdgProton;}
      else
	{
	  // for pi0
	  tcode  = (rnd->RndFsi().Rndm()<=ppcnt)?(kPdgProton) :(kPdgNeutron);
	  scode  = (tcode == kPdgProton)        ?(kPdgPiP)    :(kPdgPiM);
	  s2code = (tcode == kPdgProton)        ?(kPdgNeutron):(kPdgProton);
	}
    }
  else
    {
      tcode = (rnd->RndFsi().Rndm()<=ppcnt)?(kPdgProton):(kPdgNeutron);
      scode = pcode;
      s2code = tcode;
    }

  // get random scattering angle
  double C3CM = fHadroData2018->IntBounce(p,tcode,scode,fate);
  if (C3CM<-1.) 
    {
      p->SetStatus(kIStStableFinalState);
      ev->AddParticle(*p);
      return;
    }

  // create scattered particle
  GHepParticle * t = new GHepParticle(*p);
  t->SetPdgCode(tcode);
  double Mt = t->Mass();
  //t->SetMomentum(TLorentzVector(0,0,0,Mt));

  // handle fermi momentum 
  if(fDoFermi)
    {
      // Handle fermi target
      Target target(ev->TargetNucleus()->Pdg());
      //LOG("HAIntranuke2018", pNOTICE) << "Nuclmodel= " << fNuclmodel->ModelType(target) ;
      target.SetHitNucPdg(tcode);
      fNuclmodel->GenerateNucleon(target);
      TVector3 tP3L = fFermiFac * fNuclmodel->Momentum3();
      double tE = TMath::Sqrt(tP3L.Mag2() + Mt*Mt);
      t->SetMomentum(TLorentzVector(tP3L,tE));
    }
  else
    {
      t->SetMomentum(TLorentzVector(0,0,0,Mt));
    }

  bool pass = utils::intranuke2018::TwoBodyCollision(ev,pcode,tcode,scode,s2code,C3CM,
						  p,t,fRemnA,fRemnZ,fRemnP4,kIMdHN);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HNIntranuke2018",pDEBUG)
    << "|p3| = " << P3L << ", E3 = " << E3L;
  LOG("HNIntranuke2018",pDEBUG)
    << "|p4| = " << P4L << ", E4 = " << E4L;
#endif

  if (pass==true)
  {
    ev->AddParticle(*p);
    ev->AddParticle(*t);
  } else
  {
    delete t; //fixes memory leak
    LOG("HNIntranuke2018", pINFO) << "Elastic in hN failed calling TwoBodyCollision";
    exceptions::INukeException exception;
    exception.SetReason("hN scattering kinematics through TwoBodyCollision failed");
    throw exception;
  }

  delete t;

}
//___________________________________________________________________________
void HNIntranuke2018::InelasticHN(GHepRecord* ev, GHepParticle* p) const
{
  // Aaron Meyer (Jan 2010)
  // Updated version of InelasticHN 

  GHepParticle s1(*p);  
  GHepParticle s2(*p);
  GHepParticle s3(*p);
  
  
  if (utils::intranuke2018::PionProduction(ev,p,&s1,&s2,&s3,fRemnA,fRemnZ,fRemnP4,fDoFermi,fFermiFac,fFermiMomentum,fNuclmodel))
	{
	  // set status of particles and return
	  
	  s1.SetStatus(kIStHadronInTheNucleus);
	  s2.SetStatus(kIStHadronInTheNucleus);
	  s3.SetStatus(kIStHadronInTheNucleus);
	  
	  ev->AddParticle(s1);
	  ev->AddParticle(s2);
	  ev->AddParticle(s3);
	}
  else
	{
	  LOG("HNIntranuke2018", pNOTICE) << "Error: could not create pion production final state";
	  exceptions::INukeException exception;
	  exception.SetReason("PionProduction in hN failed");
	  throw exception;
	}
  return;

}
//___________________________________________________________________________
void HNIntranuke2018::GammaInelasticHN(GHepRecord* ev, GHepParticle* p, INukeFateHN_t fate) const     
{
  // This function handles pion photoproduction reactions

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HNIntranuke2018", pDEBUG)
    << "GammaInelasticHN() is invoked for a : " << p->Name()
    << " whose fate is : " << INukeHadroFates::AsString(fate);
#endif

  if(fate!=kIHNFtInelas && p->Pdg()!=kPdgGamma)
    {
      LOG("HNIntranuke2018", pWARN)
	<< "GammaInelasticHN() cannot handle fate: " << INukeHadroFates::AsString(fate);
      return;
    }

  // random number generator
  RandomGen * rnd = RandomGen::Instance();

  // vars for incoming particle, target, and scattered reaction products
  double ppcnt = (double) fRemnZ / (double) fRemnA; // % of protons
  int pcode = p->Pdg();
  int tcode = (rnd->RndFsi().Rndm()<=ppcnt)?(kPdgProton):(kPdgNeutron);
  int scode, s2code;
  double ke   = p->KinE() / units::MeV;

  LOG("HNIntranuke2018", pNOTICE)
    << "Particle code: " << pcode << ", target: " << tcode;


  if (rnd->RndFsi().Rndm() * (fHadroData2018 -> XSecGamp_fs() -> Evaluate(ke) +
			      fHadroData2018 -> XSecGamn_fs() -> Evaluate(ke)  )
      <= fHadroData2018 -> XSecGamp_fs() -> Evaluate(ke) ) { scode = kPdgProton;  }
  else                                                 { scode = kPdgNeutron; }

  //scode=fHadroData2018->AngleAndProduct(p,tcode,C3CM,fate);
  //double C3CM = 0.0; // cos of scattering angle
  double C3CM = fHadroData2018->IntBounce(p,tcode,scode,fate);

  if      ((tcode == kPdgProton ) && (scode==kPdgProton )) {s2code=kPdgPi0;}
  else if ((tcode == kPdgProton ) && (scode==kPdgNeutron)) {s2code=kPdgPiP;}
  else if ((tcode == kPdgNeutron) && (scode==kPdgProton )) {s2code=kPdgPiM;}
  else if ((tcode == kPdgNeutron) && (scode==kPdgNeutron)) {s2code=kPdgPi0;}
  else {
    LOG("HNIntranuke2018", pERROR)
      << "Error: could not determine particle final states";
    ev->AddParticle(*p);
    return;
  }    

  LOG("HNIntranuke2018", pNOTICE)
    << "GammaInelastic fate: " << INukeHadroFates::AsString(fate);
  LOG("HNIntranuke2018", pNOTICE)
    << " final state: " << scode << " and " << s2code;
  LOG("HNIntranuke2018", pNOTICE)
    << " scattering angle: " << C3CM;

  GHepParticle * t = new GHepParticle(*p);
  t->SetPdgCode(tcode);
  double Mt = t->Mass();

  // handle fermi momentum 
  if(fDoFermi)
    {
      // Handle fermi target
      Target target(ev->TargetNucleus()->Pdg());

      target.SetHitNucPdg(tcode);
      fNuclmodel->GenerateNucleon(target);
      TVector3 tP3L = fFermiFac * fNuclmodel->Momentum3();
      double tE = TMath::Sqrt(tP3L.Mag2() + Mt*Mt);
      t->SetMomentum(TLorentzVector(tP3L,tE));
    }
  else
    {
      t->SetMomentum(TLorentzVector(0,0,0,Mt));
    }

  bool pass = utils::intranuke2018::TwoBodyCollision(ev,pcode,tcode,scode,s2code,C3CM,
						  p,t,fRemnA,fRemnZ,fRemnP4,kIMdHN);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HNIntranuke2018",pDEBUG)
    << "|p3| = " << P3L << ", E3 = " << E3L;
  LOG("HNIntranuke2018",pDEBUG)
    << "|p4| = " << P4L << ", E4 = " << E4L;
#endif

  if (pass==true)
  {
    //p->SetStatus(kIStStableFinalState);
    //t->SetStatus(kIStStableFinalState);
    ev->AddParticle(*p);
    ev->AddParticle(*t);
  } else
  {
    ev->AddParticle(*p);
  }

  delete t;

}
//___________________________________________________________________________
int HNIntranuke2018::HandleCompoundNucleus(GHepRecord* ev, GHepParticle* p, int mom) const
{

  // handle compound nucleus option
  // -- Call the PreEquilibrium function
  if( fDoCompoundNucleus && IsInNucleus(p) && pdg::IsNeutronOrProton(p->Pdg())) 
    {  // random number generator
  //unused var - quiet compiler warning//RandomGen * rnd = RandomGen::Instance();

      if((p->KinE() < fEPreEq) )
	{
	  if(fRemnA>4)  //this needs to be matched to what is in PreEq and Eq
            {
              GHepParticle * sp = new GHepParticle(*p);
              sp->SetFirstMother(mom);
	      // this was PreEquilibrium - now just used for hN
	      //same arguement lists for PreEq and Eq
	      utils::intranuke2018::Equilibrium(ev,sp,fRemnA,fRemnZ,fRemnP4,
					       fDoFermi,fFermiFac,fNuclmodel,fNucRmvE,kIMdHN);

              delete sp;
              return 2;
            }
	  else
            {
              // nothing left to interact with!
              LOG("HNIntranuke2018", pNOTICE)
                << "*** Nothing left to interact with, escaping.";
              GHepParticle * sp = new GHepParticle(*p);
              sp->SetFirstMother(mom);
              sp->SetStatus(kIStStableFinalState);
              ev->AddParticle(*sp);
              delete sp;
              return 1;
            }
	}
    }
  return 0;
}
//___________________________________________________________________________
bool HNIntranuke2018::HandleCompoundNucleusHN(GHepRecord* ev, GHepParticle* p) const
{
  return (this->HandleCompoundNucleus(ev,p,p->FirstMother())==2);
}
//___________________________________________________________________________

void HNIntranuke2018::LoadConfig(void)
{
  // load hadronic cross sections
  fHadroData2018 = INukeHadroData2018::Instance();

  // fermi momentum setup
  // this is specifically set in Intranuke2018::Configure(string)
  fNuclmodel = dynamic_cast<const NuclearModelI *>( this -> SubAlg("NuclearModel") ) ;

  // other intranuke config params
  GetParam( "NUCL-R0",             fR0 );              // fm
  GetParam( "NUCL-NR",             fNR );

  GetParam( "INUKE-NucRemovalE",   fNucRmvE );        // GeV
  GetParam( "INUKE-HadStep",       fHadStep ) ;
  GetParam( "INUKE-NucAbsFac",     fNucAbsFac ) ;
  GetParam( "INUKE-NucQEFac",      fNucQEFac ) ;
  GetParam( "INUKE-NucCEXFac",     fNucCEXFac ) ;
  GetParam( "INUKE-Energy_Pre_Eq", fEPreEq ) ;
  GetParam( "INUKE-FermiFac",      fFermiFac ) ;
  GetParam( "INUKE-FermiMomentum", fFermiMomentum ) ;

  GetParam( "INUKE-DoCompoundNucleus", fDoCompoundNucleus ) ;
  GetParam( "INUKE-DoFermi",           fDoFermi ) ;
  GetParam( "INUKE-XsecNNCorr",        fXsecNNCorr ) ;
  GetParamDef( "AltOset",              fAltOset, false ) ;

  GetParam( "HNINUKE-UseOset",     fUseOset ) ;
  GetParam( "HNINUKE-DelRPion",    fDelRPion ) ;
  GetParam( "HNINUKE-DelRNucleon", fDelRNucleon ) ;

  // report
  LOG("HNIntranuke2018", pINFO) << "Settings for Intranuke2018 mode: " << INukeMode::AsString(kIMdHN);
  LOG("HNIntranuke2018", pWARN) << "R0          = " << fR0 << " fermi";
  LOG("HNIntranuke2018", pWARN) << "NR          = " << fNR;
  LOG("HNIntranuke2018", pWARN) << "DelRPion    = " << fDelRPion;
  LOG("HNIntranuke2018", pWARN) << "DelRNucleon = " << fDelRNucleon;
  LOG("HNIntranuke2018", pWARN) << "HadStep     = " << fHadStep << " fermi";
  LOG("HNIntranuke2018", pWARN) << "NucAbsFac   = " << fNucAbsFac;
  LOG("HNIntranuke2018", pWARN) << "NucQEFac    = " << fNucQEFac;
  LOG("HNIntranuke2018", pWARN) << "NucCEXFac   = " << fNucCEXFac;
  LOG("HNIntranuke2018", pWARN) << "EPreEq      = " << fEPreEq;
  LOG("HNIntranuke2018", pWARN) << "FermiFac    = " << fFermiFac;
  LOG("HNIntranuke2018", pWARN) << "FermiMomtm  = " << fFermiMomentum;
  LOG("HNIntranuke2018", pWARN) << "DoFermi?    = " << ((fDoFermi)?(true):(false));
  LOG("HNIntranuke2018", pWARN) << "DoCmpndNuc? = " << ((fDoCompoundNucleus)?(true):(false));
  LOG("HNIntranuke2018", pWARN) << "useOset     = " << fUseOset;
  LOG("HNIntranuke2018", pWARN) << "altOset     = " << fAltOset;
  LOG("HNIntranuke2018", pWARN) << "XsecNNCorr? = " << ((fXsecNNCorr)?(true):(false));
}
//___________________________________________________________________________

INukeFateHN_t HNIntranuke2018::HadronFateOset () const
{
  //LOG("HNIntranuke2018", pWARN) << "IN HadronFateOset";

  //LOG("HNIntranuke2018", pWARN) << "{ frac abs  = " << osetUtils::currentInstance->getAbsorptionFraction();
  //LOG("HNIntranuke2018", pWARN) << "  frac cex  = " << osetUtils::currentInstance->getCexFraction() << " }";

  const double fractionAbsorption = osetUtils::currentInstance->
                                    getAbsorptionFraction();
  const double fractionCex = osetUtils::currentInstance->getCexFraction ();

  RandomGen *randomGenerator = RandomGen::Instance();
  const double randomNumber  = randomGenerator->RndFsi().Rndm();

  //LOG("HNIntranuke2018", pWARN) << "{ frac abs  = " << fractionAbsorption;
  //LOG("HNIntranuke2018", pWARN) << "  frac cex  = " << fractionCex;
  //LOG("HNIntranuke2018", pWARN) << "  frac elas = " << 1-fractionAbsorption-fractionCex << " }";

  if (randomNumber < fractionAbsorption && fRemnA > 1) return kIHNFtAbs;
  else if (randomNumber < fractionAbsorption + fractionCex) return kIHNFtCEx;
  else return kIHNFtElas;
}
