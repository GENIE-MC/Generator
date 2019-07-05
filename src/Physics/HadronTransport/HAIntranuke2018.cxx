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
   Major overhaul of the function of each interaction type. Absorption fates
   changed to allow more than 6 particles at a time (up to 85 now). PiPro fates
   now allow the pion to rescatter inside the nucleus, will be changed at a
   later date. HAIntranuke class is now defined as derived from virtual class.
   Intranuke.
 @ Oct 10, 2011 - SD
   Changes to keep reweighting alive.  Add exception handling in ElasHA, InelasticHA,
   and Inelastic.
 @ Jan 24, 2012 - SD
   Add option of doing K+.  
 @ Jan 9, 2015 - SD, NG, TG
   Added 2014 version of INTRANUKE codes (new class) for independent development.
 @ Aug 30, 2016 - SD
   Fix memory leaks - Igor. 
*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>
#include <exception>

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
#include "Physics/HadronTransport/INukeException.h"
#include "Physics/HadronTransport/Intranuke2018.h"
#include "Physics/HadronTransport/HAIntranuke2018.h"
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
#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/NuclearState/NuclearModelMap.h"
//#include "Physics/HadronTransport/INukeOset.h" 

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::utils::intranuke2018;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
//___________________________________________________________________________
// Methods specific to INTRANUKE's HA-mode
//___________________________________________________________________________
//___________________________________________________________________________
HAIntranuke2018::HAIntranuke2018() :
Intranuke2018("genie::HAIntranuke2018")
{

}
//___________________________________________________________________________
HAIntranuke2018::HAIntranuke2018(string config) :
Intranuke2018("genie::HAIntranuke2018",config)
{

}
//___________________________________________________________________________
HAIntranuke2018::~HAIntranuke2018()
{

}
//___________________________________________________________________________
void HAIntranuke2018::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("HAIntranuke2018", pNOTICE) 
     << "************ Running hA2018 MODE INTRANUKE ************";
  GHepParticle * nuclearTarget = evrec -> TargetNucleus();
  nuclA = nuclearTarget -> A();

  Intranuke2018::ProcessEventRecord(evrec);

  LOG("HAIntranuke2018", pINFO) << "Done with this event";
}
//___________________________________________________________________________
void HAIntranuke2018::SimulateHadronicFinalState(
  GHepRecord* ev, GHepParticle* p) const
{
// Simulate a hadron interaction for the input particle p in HA mode
//
  // check inputs
  if(!p || !ev) {
     LOG("HAIntranuke2018", pERROR) << "** Null input!";
     return;
  }

  // get particle code and check whether this particle can be handled
  int  pdgc = p->Pdg();
  bool is_gamma   = (pdgc==kPdgGamma);						  
  bool is_pion    = (pdgc==kPdgPiP || pdgc==kPdgPiM || pdgc==kPdgPi0);
  bool is_kaon    = (pdgc==kPdgKP || pdgc==kPdgKM);
  bool is_baryon  = (pdgc==kPdgProton || pdgc==kPdgNeutron);
  bool is_handled = (is_baryon || is_pion || is_kaon || is_gamma); 
  if(!is_handled) {
     LOG("HAIntranuke2018", pERROR) << "** Can not handle particle: " << p->Name();
     return;     
  }

  // select a fate for the input particle
  INukeFateHA_t fate = this->HadronFateHA(p);

  // store the fate
  ev->Particle(p->FirstMother())->SetRescatterCode((int)fate);

  if(fate == kIHAFtUndefined) {
     LOG("HAIntranuke2018", pERROR) << "** Couldn't select a fate";
     p->SetStatus(kIStStableFinalState);
     ev->AddParticle(*p);
     return;     
  }
   LOG("HAIntranuke2018", pNOTICE)
     << "Selected "<< p->Name() << " fate: "<< INukeHadroFates::AsString(fate);

  // try to generate kinematics - repeat till is done (should seldom need >2)
   fNumIterations = 0;
   this->SimulateHadronicFinalStateKinematics(ev,p);
}
//___________________________________________________________________________
void HAIntranuke2018::SimulateHadronicFinalStateKinematics(
  GHepRecord* ev, GHepParticle* p) const
{
  // get stored fate
  INukeFateHA_t fate = (INukeFateHA_t) 
      ev->Particle(p->FirstMother())->RescatterCode();

   LOG("HAIntranuke2018", pINFO)
     << "Generating kinematics for " << p->Name() 
     << " fate: "<< INukeHadroFates::AsString(fate);

  // try to generate kinematics for the selected fate 

  try
    {
      fNumIterations++;
      /*      if (fate == kIHAFtElas)
	{ 
	  this->ElasHA(ev,p,fate);
	}
	else */ 
	if (fate == kIHAFtInelas || fate == kIHAFtCEx) 
	  {
	    this->InelasticHA(ev,p,fate);
	  }
	else if (fate == kIHAFtAbs || fate == kIHAFtPiProd)
	  {
	    this->Inelastic(ev,p,fate);
	  }
	else if (fate == kIHAFtCmp) //(suarez edit, 17 July, 2017: cmp)
	  {
	    LOG("HAIntranuke2018", pWARN) << "Running PreEquilibrium for kIHAFtCmp";
	    utils::intranuke2018::PreEquilibrium(ev,p,fRemnA,fRemnZ,fRemnP4,fDoFermi,fFermiFac,fNuclmodel,fNucRmvE,kIMdHA); //should be kiMdHA or HN?
	  }
    }
  catch(exceptions::INukeException exception)
    {     
      LOG("HAIntranuke2018", pNOTICE)  
	<< exception;
    if(fNumIterations <= 100) {
      LOG("HAIntranuke2018", pNOTICE)
	<< "Failed attempt to generate kinematics for "
	<< p->Name() << " fate: " << INukeHadroFates::AsString(fate)
	<< " - After " << fNumIterations << " tries, still retrying...";
      this->SimulateHadronicFinalStateKinematics(ev,p);
    } else {
      LOG("HAIntranuke2018", pNOTICE)
	<< "Failed attempt to generate kinematics for "
	<< p->Name() << " fate: " << INukeHadroFates::AsString(fate)
	<< " after " << fNumIterations-1 
	<< " attempts. Trying a new fate...";
      this->SimulateHadronicFinalState(ev,p);
    }
    }
}
//___________________________________________________________________________
INukeFateHA_t HAIntranuke2018::HadronFateHA(const GHepParticle * p) const
{
// Select a hadron fate in HA mode
//
  RandomGen * rnd = RandomGen::Instance();

  // get pdgc code & kinetic energy in MeV
  int    pdgc = p->Pdg();
  double ke   = p->KinE() / units::MeV;

  //bool isPion = (pdgc == kPdgPiP or pdgc == kPdgPi0 or pdgc == kPdgPiM); 
  //if (isPion and fUseOset and ke < 350.0) return this->HadronFateOset(); 

  LOG("HAIntranuke2018", pINFO) 
   << "Selecting hA fate for " << p->Name() << " with KE = " << ke << " MeV";

  // try to generate a hadron fate
  unsigned int iter = 0;
  while(iter++ < kRjMaxIterations) {

    // handle pions
    //
   if (pdgc==kPdgPiP || pdgc==kPdgPiM || pdgc==kPdgPi0) {

     double frac_cex      = fHadroData2018->FracADep(pdgc, kIHAFtCEx,     ke, nuclA);
     //     double frac_elas     = fHadroData2018->FracADep(pdgc, kIHAFtElas,    ke, nuclA);
     double frac_inel     = fHadroData2018->FracADep(pdgc, kIHAFtInelas,  ke, nuclA);
     double frac_abs      = fHadroData2018->FracADep(pdgc, kIHAFtAbs,     ke, nuclA);
     double frac_piprod   = fHadroData2018->FracADep(pdgc, kIHAFtPiProd,  ke, nuclA);
     LOG("HAIntranuke2018", pDEBUG) 
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtCEx)     << "} = " << frac_cex
       //          << "\n frac{" << INukeHadroFates::AsString(kIHAFtElas)    << "} = " << frac_elas
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtInelas)  << "} = " << frac_inel
	  << "\n frac{" << INukeHadroFates::AsString(kIHAFtAbs)     << "} = " << frac_abs
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtPiProd)  << "} = " << frac_piprod;
          
     // apply external tweaks to fractions
     frac_cex    *= fPionFracCExScale;
     frac_inel   *= fPionFracInelScale;
     frac_abs    *= fPionFracAbsScale;
     frac_piprod *= fPionFracPiProdScale;

     double frac_rescale = 1./(frac_cex + frac_inel + frac_abs + frac_piprod);

     frac_cex    *= frac_rescale;
     frac_inel   *= frac_rescale;
     frac_abs    *= frac_rescale;
     frac_piprod *= frac_rescale;


       // compute total fraction (can be <1 if fates have been switched off)
       double tf = frac_cex      +
	 //                   frac_elas     +
                   frac_inel     +  
	           frac_abs      +
                   frac_piprod;

       double r = tf * rnd->RndFsi().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("HAIntranuke2018", pDEBUG) << "r = " << r << " (max = " << tf << ")";
#endif
       double cf=0; // current fraction
       if(r < (cf += frac_cex     )) return kIHAFtCEx;     // cex
       //       if(r < (cf += frac_elas    )) return kIHAFtElas;    // elas
       if(r < (cf += frac_inel    )) return kIHAFtInelas;  // inelas
       if(r < (cf += frac_abs     )) return kIHAFtAbs;     // abs
       if(r < (cf += frac_piprod  )) return kIHAFtPiProd;  // pi prod

       LOG("HAIntranuke2018", pWARN) 
         << "No selection after going through all fates! " 
                     << "Total fraction = " << tf << " (r = " << r << ")";
    }

    // handle nucleons
    else if (pdgc==kPdgProton || pdgc==kPdgNeutron) {
      double frac_cex      = fHadroData2018->FracAIndep(pdgc, kIHAFtCEx,    ke);
      //double frac_elas     = fHadroData2018->FracAIndep(pdgc, kIHAFtElas,   ke);
      double frac_inel     = fHadroData2018->FracAIndep(pdgc, kIHAFtInelas, ke);
      double frac_abs      = fHadroData2018->FracAIndep(pdgc, kIHAFtAbs,    ke);
      double frac_pipro    = fHadroData2018->FracAIndep(pdgc, kIHAFtPiProd, ke);
      double frac_cmp      = fHadroData2018->FracAIndep(pdgc, kIHAFtCmp   , ke);

      LOG("HAIntranuke2018", pINFO)
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtCEx)     << "} = " << frac_cex
	// << "\n frac{" << INukeHadroFates::AsString(kIHAFtElas)    << "} = " << frac_elas
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtInelas)  << "} = " << frac_inel
	  << "\n frac{" << INukeHadroFates::AsString(kIHAFtAbs)     << "} = " << frac_abs
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtPiProd)  << "} = " << frac_pipro
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtCmp)     << "} = " << frac_cmp; //suarez edit, cmp

      // apply external tweaks to fractions
      frac_cex    *= fNucleonFracCExScale;
      frac_inel   *= fNucleonFracInelScale;
      frac_abs    *= fNucleonFracAbsScale;
      frac_pipro  *= fNucleonFracPiProdScale;

      double frac_rescale = 1./(frac_cex + frac_inel + frac_abs + frac_pipro);

      frac_cex    *= frac_rescale;
      frac_inel   *= frac_rescale;
      frac_abs    *= frac_rescale;
      frac_pipro  *= frac_rescale;

       // compute total fraction (can be <1 if fates have been switched off)
       double tf = frac_cex      +
	 //frac_elas     +
                   frac_inel     +  
	           frac_abs      +
	           frac_pipro    +
	           frac_cmp;  //suarez edit, cmp

       double r = tf * rnd->RndFsi().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("HAIntranuke2018", pDEBUG) << "r = " << r << " (max = " << tf << ")";
#endif
       double cf=0; // current fraction
       if(r < (cf += frac_cex     )) return kIHAFtCEx;     // cex
       //if(r < (cf += frac_elas    )) return kIHAFtElas;    // elas
       if(r < (cf += frac_inel    )) return kIHAFtInelas;  // inelas
       if(r < (cf += frac_abs     )) return kIHAFtAbs;     // abs
       if(r < (cf += frac_pipro   )) return kIHAFtPiProd;  // pi prod 
       if(r < (cf += frac_cmp     )) return kIHAFtCmp;  //suarez edit, cmp

       LOG("HAIntranuke2018", pWARN) 
         << "No selection after going through all fates! "
                        << "Total fraction = " << tf << " (r = " << r << ")";
    }
    // handle kaons
    else if (pdgc==kPdgKP || pdgc==kPdgKM) {
       double frac_inel     = fHadroData2018->FracAIndep(pdgc, kIHAFtInelas,  ke);
       double frac_abs      = fHadroData2018->FracAIndep(pdgc, kIHAFtAbs,     ke);

       LOG("HAIntranuke2018", pDEBUG) 
          << "\n frac{" << INukeHadroFates::AsString(kIHAFtInelas)  << "} = " << frac_inel
	  << "\n frac{" << INukeHadroFates::AsString(kIHAFtAbs)     << "} = " << frac_abs;
       // compute total fraction (can be <1 if fates have been switched off)
       double tf =  frac_inel     +  
	 frac_abs;
      double r = tf * rnd->RndFsi().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("HAIntranuke2018", pDEBUG) << "r = " << r << " (max = " << tf << ")";
#endif
       double cf=0; // current fraction
       if(r < (cf += frac_inel    )) return kIHAFtInelas;  // inelas
       if(r < (cf += frac_abs     )) return kIHAFtAbs;     // abs
    }
  }//iterations

  return kIHAFtUndefined; 
}
//___________________________________________________________________________
double HAIntranuke2018::PiBounce(void) const
{
// [adapted from neugen3 intranuke_bounce.F]
// [is a fortran stub / difficult to understand - needs to be improved]
//
// Generates theta in radians for elastic pion-nucleus scattering/
// Lookup table is based on Fig 17 of Freedman, Miller and Henley, Nucl.Phys.
// A389, 457 (1982)
//
  const int nprob = 25;
  double dintor = 0.0174533;
  double denom  = 47979.453;
  double rprob[nprob] = {
    5000., 4200., 3000., 2600., 2100., 1800., 1200., 750., 500., 230., 120.,
    35., 9., 3., 11., 18., 29., 27., 20., 14., 10., 6., 2., 0.14, 0.19 };

  double angles[nprob];
  for(int i=0; i<nprob; i++) angles[i] = 2.5*i;

  RandomGen * rnd = RandomGen::Instance();
  double r = rnd->RndFsi().Rndm();

  double xsum = 0.;
  double theta = 0.;
  double binl  = 0.;
  double binh  = 0.;
  int tj = 0;
  for(int i=0; i<60; i++) {
   theta = i+0.5;
   for(int j=0; j < nprob-1; j++) {
     binl = angles[j];
     binh = angles[j+1];
     tj=j;
     if(binl<=theta && binh>=theta) break;
     tj=0;
   }//j
   int itj = tj;
   double tfract = (theta-binl)/2.5;
   double delp   = rprob[itj+1] - rprob[itj];
   xsum += (rprob[itj] + tfract*delp)/denom;
   if(xsum>r) break;
   theta = 0.;
  }//i

  theta *= dintor;

  LOG("HAIntranuke2018", pNOTICE)
     << "Generated pi+A elastic scattering angle = " << theta << " radians";

  return theta;
}
//___________________________________________________________________________
double HAIntranuke2018::PnBounce(void) const
{
// [adapted from neugen3 intranuke_pnbounce.F]
// [is a fortran stub / difficult to understand - needs to be improved]
//
// Generates theta in radians for elastic nucleon-nucleus scattering.
// Use 800 MeV p+O16 as template in same (highly simplified) spirit as pi+A
// from table in Adams et al., PRL 1979. Guess value at 0-2 deg based on Ni
// data.
//
  const int nprob = 20;
  double dintor = 0.0174533;
  double denom  = 11967.0;
  double rprob[nprob] = {
    2400., 2350., 2200., 2000., 1728., 1261., 713., 312., 106., 35., 
    6., 5., 10., 12., 11., 9., 6., 1., 1., 1. };

  double angles[nprob];
  for(int i=0; i<nprob; i++) angles[i] = 1.0*i;

  RandomGen * rnd = RandomGen::Instance();
  double r = rnd->RndFsi().Rndm();

  double xsum = 0.;
  double theta = 0.;
  double binl  = 0.;
  double binh  = 0.;
  int tj = 0;
  for(int i=0; i< nprob; i++) {
   theta = i+0.5;
   for(int j=0; j < nprob-1; j++) {
     binl = angles[j];
     binh = angles[j+1];
     tj=j;
     if(binl<=theta && binh>=theta) break;
     tj=0;
   }//j
   int itj = tj;
   double tfract = (theta-binl)/2.5;
   double delp   = rprob[itj+1] - rprob[itj];
   xsum += (rprob[itj] + tfract*delp)/denom;
   if(xsum>r) break;
   theta = 0.;
  }//i

  theta *= dintor;

  LOG("HAIntranuke2018", pNOTICE)
     << "Generated N+A elastic scattering angle = " << theta << " radians";

  return theta;
}
//___________________________________________________________________________
void HAIntranuke2018::ElasHA(
         GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const
{
  // scatters particle within nucleus, copy of hN code meant to run only once
  // in hA mode

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HAIntranuke2018", pDEBUG)
    << "ElasHA() is invoked for a : " << p->Name()
    << " whose fate is : " << INukeHadroFates::AsString(fate);
#endif

  /*  if(fate!=kIHAFtElas)
    {
      LOG("HAIntranuke2018", pWARN)
	<< "ElasHA() cannot handle fate: " << INukeHadroFates::AsString(fate);
      return;
      } */

  // check remnants
  if(fRemnA<0 || fRemnZ<0) // best to stop it here and not try again.
    {
      LOG("HAIntranuke2018", pWARN) << "Invalid Nucleus! : (A,Z) = ("<<fRemnA<<','<<fRemnZ<<')';
      p->SetStatus(kIStStableFinalState);
      ev->AddParticle(*p);
      return;
    }

  // vars for incoming particle, target, and scattered pdg codes
  int pcode = p->Pdg();
  double Mp = p->Mass();
  double Mt = 0.;
  if (ev->TargetNucleus()->A()==fRemnA)
    { Mt = PDGLibrary::Instance()->Find(ev->TargetNucleus()->Pdg())->Mass(); }
  else 
    {
      Mt = fRemnP4.M();
    }
  TLorentzVector t4PpL = *p->P4();
  TLorentzVector t4PtL = fRemnP4;
  double C3CM = 0.0;

  // calculate scattering angle
  if(pcode==kPdgNeutron||pcode==kPdgProton) C3CM = TMath::Cos(this->PnBounce());
  else C3CM = TMath::Cos(this->PiBounce());

  // calculate final 4 momentum of probe
  TLorentzVector t4P3L, t4P4L;

  if (!utils::intranuke2018::TwoBodyKinematics(Mp,Mt,t4PpL,t4PtL,t4P3L,t4P4L,C3CM,fRemnP4))
    {
      LOG("HAIntranuke2018", pNOTICE) << "ElasHA() failed";
      exceptions::INukeException exception;
      exception.SetReason("TwoBodyKinematics failed in ElasHA");
      throw exception;
    }

  // Update probe particle
  p->SetMomentum(t4P3L);
  p->SetStatus(kIStStableFinalState);

  // Update Remnant nucleus
  fRemnP4 = t4P4L;
  LOG("HAIntranuke2018",pINFO)
    << "C3cm = " << C3CM;
  LOG("HAIntranuke2018",pINFO)
    << "|p3| = " << t4P3L.Vect().Mag()   << ", E3 = " << t4P3L.E() << ",Mp = " << Mp;
  LOG("HAIntranuke2018",pINFO)
    << "|p4| = " << fRemnP4.Vect().Mag() << ", E4 = " << fRemnP4.E() << ",Mt = " << Mt;

  ev->AddParticle(*p);

}
//___________________________________________________________________________
void HAIntranuke2018::InelasticHA(
	GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const
{
  // charge exch and inelastic - scatters particle within nucleus, hA version
  // each are treated as quasielastic, particle scatters off single nucleon

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HAIntranuke2018", pDEBUG)
    << "InelasticHA() is invoked for a : " << p->Name()
    << " whose fate is : " << INukeHadroFates::AsString(fate);
#endif
  if(ev->Probe() ) {
    LOG("HAIntranuke2018", pINFO) << " probe KE = " << ev->Probe()->KinE();
  }
  if(fate!=kIHAFtCEx && fate!=kIHAFtInelas)
    {
      LOG("HAIntranuke2018", pWARN)
	<< "InelasticHA() cannot handle fate: " << INukeHadroFates::AsString(fate);
      return;
    }

  // Random number generator
  RandomGen * rnd = RandomGen::Instance();

  // vars for incoming particle, target, and scattered pdg codes
  int pcode = p->Pdg();
  int tcode, scode, s2code;
  double ppcnt = (double) fRemnZ / (double) fRemnA; // % of protons

  // Select a hadron fate in HN mode
  INukeFateHN_t h_fate;
  if (fate == kIHAFtCEx) h_fate = kIHNFtCEx;
  else                   h_fate = kIHNFtElas;

  // Select a target randomly, weighted to #
  // -- Unless, of course, the fate is CEx,
  // -- in which case the target may be deterministic
  // Also assign scattered particle code
  if(fate==kIHAFtCEx)
    {
      if(pcode==kPdgPiP)         {tcode = kPdgNeutron; scode = kPdgPi0; s2code = kPdgProton;}
      else if(pcode==kPdgPiM)    {tcode = kPdgProton;  scode = kPdgPi0; s2code = kPdgNeutron;}
      else if(pcode==kPdgPi0)
	{
	  // for pi0
	  tcode  = (rnd->RndFsi().Rndm()<=ppcnt)?(kPdgProton) :(kPdgNeutron);
	  scode  = (tcode == kPdgProton)        ?(kPdgPiP)    :(kPdgPiM);
	  s2code = (tcode == kPdgProton)        ?(kPdgNeutron):(kPdgProton);
	}
      else if(pcode==kPdgProton) {tcode = kPdgNeutron; scode = kPdgNeutron; s2code = kPdgProton;}
      else if(pcode==kPdgNeutron){tcode = kPdgProton; scode = kPdgProton; s2code = kPdgNeutron;}
      else
	{ LOG("HAIntranuke2018", pWARN) << "InelasticHA() cannot handle fate: " 
				    << INukeHadroFates::AsString(fate)
				    << " for particle " << p->Name();
	  return;
	}
    }
  else
    {
      tcode = (rnd->RndFsi().Rndm()<=ppcnt)?(kPdgProton):(kPdgNeutron);
      //      if(pcode == kPdgKP || pcode == kPdgKM) tcode = kPdgProton;
      scode = pcode;
      s2code = tcode;
    }

  // check remnants
  if ( fRemnA < 1 )    //we've blown nucleus apart, no need to retry anything - exit
    {
      LOG("HAIntranuke2018",pNOTICE) << "InelasticHA() stops : not enough nucleons";
      p->SetStatus(kIStStableFinalState);
      ev->AddParticle(*p);
      return;
    }
  else if ( fRemnZ + (((pcode==kPdgProton)||(pcode==kPdgPiP))?1:0) - (pcode==kPdgPiM?1:0)
	    < ((( scode==kPdgProton)||( scode==kPdgPiP)) ?1:0) - (scode ==kPdgPiM ?1:0)
	    + (((s2code==kPdgProton)||(s2code==kPdgPiP)) ?1:0) - (s2code==kPdgPiM ?1:0) )
    {
      LOG("HAIntranuke2018",pWARN) << "InelasticHA() failed : too few protons in nucleus";
      p->SetStatus(kIStStableFinalState);
      ev->AddParticle(*p);
      return;         // another extreme case, best strategy is to exit and go to next event
    }

  GHepParticle t(*p);
  t.SetPdgCode(tcode);

  // set up fermi target
  Target target(ev->TargetNucleus()->Pdg());
  double tM = t.Mass();

  // handle fermi momentum
  if(fDoFermi)
    {
      target.SetHitNucPdg(tcode);
      fNuclmodel->GenerateNucleon(target);
      TVector3 tP3 = fFermiFac * fNuclmodel->Momentum3();
      double tE = TMath::Sqrt(tP3.Mag2()+ tM*tM);
      t.SetMomentum(TLorentzVector(tP3,tE));
    }
  else
    {
      t.SetMomentum(0,0,0,tM);
    }

  GHepParticle * cl = new GHepParticle(*p); // clone particle, to run IntBounce at proper energy
                                            // calculate energy and momentum using invariant mass
  double pM  = p->Mass();
  double E_p = ((*p->P4() + *t.P4()).Mag2() - tM*tM - pM*pM)/(2.0*tM);
  double P_p = TMath::Sqrt(E_p*E_p - pM*pM);
  cl->SetMomentum(TLorentzVector(P_p,0,0,E_p)); 
                  // momentum doesn't have to be in right direction, only magnitude
  double C3CM = fHadroData2018->IntBounce(cl,tcode,scode,h_fate);
  delete cl;
  if (C3CM<-1.)   // hope this doesn't occur too often - unphysical but we just pass it on
    {
      LOG("HAIntranuke2018", pWARN) << "unphysical angle chosen in InelasicHA - put particle outside nucleus";
      p->SetStatus(kIStStableFinalState);
      ev->AddParticle(*p);
      return;
    }  
    double KE1L = p->KinE();
    double KE2L = t.KinE();
    LOG("HAIntranuke2018",pINFO)
      <<  "  KE1L = " << KE1L << "   " << KE1L << "  KE2L = " << KE2L; 
  GHepParticle cl1(*p);
  GHepParticle cl2(t);
  bool success = utils::intranuke2018::TwoBodyCollision(ev,pcode,tcode,scode,s2code,C3CM,
					       &cl1,&cl2,fRemnA,fRemnZ,fRemnP4,kIMdHA); 
  if(success)
    {
      double P3L = TMath::Sqrt(cl1.Px()*cl1.Px() + cl1.Py()*cl1.Py() + cl1.Pz()*cl1.Pz());
      double P4L = TMath::Sqrt(cl2.Px()*cl2.Px() + cl2.Py()*cl2.Py() + cl2.Pz()*cl2.Pz());
      double E3L = cl1.KinE();
      double E4L = cl2.KinE();
      LOG ("HAIntranuke2018",pINFO) << "Successful quasielastic scattering or charge exchange";
      LOG("HAIntranuke",pINFO)
	<< "C3CM = " << C3CM << "\n  P3L, E3L = " 
	<< P3L << "   " << E3L << "  P4L, E4L = "<< P4L << "   " << E4L ;
      if(ev->Probe() ) { LOG("HAIntranuke",pINFO)
	  << "P4L = " << P4L << " ;E4L=  " << E4L << "\n probe KE = " << ev->Probe()->KinE() << "\n";
	LOG("HAIntranuke2018", pINFO) << "Nucleus : (A,Z) = ("<<fRemnA<<','<<fRemnZ<<')';
	TParticlePDG * remn = 0;
	double MassRem = 0.;
	int ipdgc = pdg::IonPdgCode(fRemnA, fRemnZ);
	remn = PDGLibrary::Instance()->Find(ipdgc);
	if(!remn) 
	  {
	    LOG("HAIntranuke2018", pINFO)
	      << "NO Particle with [A = " << fRemnA << ", Z = " << fRemnZ
	      << ", pdgc = " << ipdgc << "] in PDGLibrary!";
	  }
	else
	  {
	    MassRem = remn->Mass();
	    LOG("HAIntranuke2018", pINFO)
	      << "Particle with [A = " << fRemnA << ", Z = " << fRemnZ
	      << ", pdgc = " << ipdgc << "] in PDGLibrary!";
	  }
	double ERemn = fRemnP4.E();
	double PRemn = TMath::Sqrt(fRemnP4.Px()*fRemnP4.Px() + fRemnP4.Py()*fRemnP4.Py() + fRemnP4.Pz()*fRemnP4.Pz());
	double MRemn = TMath::Sqrt(ERemn*ERemn - PRemn*PRemn);
	LOG("HAIntranuke2018",pINFO) << "PRemn = " << PRemn << " ;ERemn=  " << ERemn;
	LOG("HAIntranuke2018",pINFO) << "MRemn=  " << MRemn << "  ;true Mass=  " << MassRem << "   ; excitation energy= " << (MRemn-MassRem)*1000. << " MeV";
      }
      if (ev->Probe() && (E3L>ev->Probe()->KinE()))  //assuming E3 is most important, definitely for pion.  what about pp?
	{
	  //	  LOG("HAIntranuke",pINFO)
	  //	    << "E3Lagain = " << E3L << " ;E4L=  " << E4L << "\n probe KE = " << ev->Probe()->KinE() << "\n";
	  exceptions::INukeException exception;
	  exception.SetReason("TwoBodyCollison gives KE> probe KE in hA simulation");
	  throw exception;
	}
      ev->AddParticle(cl1);
      ev->AddParticle(cl2);

      LOG("HAIntranuke2018", pDEBUG) << "Nucleus : (A,Z) = ("<<fRemnA<<','<<fRemnZ<<')';
    } else
    {
      exceptions::INukeException exception;
      exception.SetReason("TwoBodyCollison failed in hA simulation");
      throw exception;
    }
  
}
//___________________________________________________________________________
void HAIntranuke2018::Inelastic(
          GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const
{

  // Aaron Meyer (05/25/10)
  //
  // Called to handle all absorption and pi production reactions
  //
  // Nucleons -> Reaction approximated by exponential decay in p+n (sum) space,
  //   gaussian in p-n (difference) space
  //   -fit to hN simulations p C, Fe, Pb at 200, 800 MeV
  //   -get n from isospin, np-nn smaller by 2
  // Pions    -> Reaction approximated with a modified gaussian in p+n space,
  //   normal gaussian in p-n space
  //   -based on fits to multiplicity distributions of hN model
  //   for pi+ C, Fe, Pb at 250, 500 MeV
  //   -fit sum and diff of nn, np to Gaussian
  //   -get pi0 from isospin, np-nn smaller by 2
  //   -get pi- from isospin, np-nn smaller by 4
  //   -add 2-body absorption to better match McKeown data
  // Kaons    -> no guidance, use same code as pions.
  //
  // Normally distributed random number generated using Box-Muller transformation
  //
  // Pion production reactions rescatter pions in nucleus, otherwise unchanged from
  //   older versions of GENIE
  //

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("HAIntranuke2018", pDEBUG) 
      << "Inelastic() is invoked for a : " << p->Name() 
      << " whose fate is : " << INukeHadroFates::AsString(fate);
#endif

  bool allow_dup = true;
  PDGCodeList list(allow_dup); // list of final state particles

  // only absorption/pipro fates allowed
  if (fate == kIHAFtPiProd) {

      GHepParticle s1(*p);
      GHepParticle s2(*p);
      GHepParticle s3(*p);

      bool success = utils::intranuke2018::PionProduction(
         ev,p,&s1,&s2,&s3,fRemnA,fRemnZ,fRemnP4, fDoFermi,fFermiFac,fFermiMomentum,fNuclmodel);

      if (success){
	LOG ("HAIntranuke2018",pINFO) << " successful pion production fate";
	// set status of particles and conserve charge/baryon number
	s1.SetStatus(kIStStableFinalState);  //should be redundant
	//	  if (pdg::IsPion(s2->Pdg())) s2->SetStatus(kIStHadronInTheNucleus);
	s2.SetStatus(kIStStableFinalState);
	//	  if (pdg::IsPion(s3->Pdg())) s3->SetStatus(kIStHadronInTheNucleus);
	s3.SetStatus(kIStStableFinalState);
	
	ev->AddParticle(s1);
	ev->AddParticle(s2);
	ev->AddParticle(s3);
	
	return;
      }
      else {
	LOG("HAIntranuke2018", pNOTICE) << "Error: could not create pion production final state";
	exceptions::INukeException exception;
	exception.SetReason("PionProduction kinematics failed - retry kinematics");
	throw exception;
      }
  }

  else if (fate==kIHAFtAbs)   
//  tuned for pions - mixture of 2-body and many-body
//   use same for kaons as there is no guidance
    {
      // Instances for reference
      PDGLibrary * pLib = PDGLibrary::Instance();
      RandomGen * rnd = RandomGen::Instance();

      double ke = p->KinE() / units::MeV;
      int pdgc = p->Pdg();

      if (fRemnA<2)
      {
	  LOG("HAIntranuke2018", pNOTICE) << "stop  propagation - could not create absorption final state: too few particles";
      p->SetStatus(kIStStableFinalState);
      ev->AddParticle(*p);
      return;
      }
      if (fRemnZ<1 && (pdgc==kPdgPiM || pdgc==kPdgKM))
      {
	LOG("HAIntranuke2018", pNOTICE) << "stop propagation - could not create absorption final state: Pi- or K- cannot be absorbed by only neutrons";
	p->SetStatus(kIStStableFinalState);
	ev->AddParticle(*p);
	return;
      }
      if (fRemnA-fRemnZ<1 && (pdgc==kPdgPiP || pdgc==kPdgKP))
      {
	LOG("HAIntranuke2018", pINFO) << "stop propagation - could not create absorption final state: Pi+ or K+ cannot be absorbed by only protons";
	p->SetStatus(kIStStableFinalState);
	ev->AddParticle(*p);
	return;
      }

      // for now, empirical split between multi-nucleon absorption and pi d -> N N
      //
      // added 03/21/11 - Aaron Meyer
      //
      if (pdg::IsPion(pdgc) && rnd->RndFsi().Rndm()<1.14*(.903-0.00189*fRemnA)*(1.35-0.00467*ke))
	{  // pi d -> N N, probability determined empirically with McKeown data

	  INukeFateHN_t fate_hN=kIHNFtAbs;
	  int t1code,t2code,scode,s2code;
	  double ppcnt = (double) fRemnZ / (double) fRemnA; // % of protons

	  // choose target nucleon
	  // -- fates weighted by values from Engel, Mosel...
	  if (pdgc==kPdgPiP) { 
            double Prob_pipd_pp=2.*ppcnt*(1.-ppcnt);
            double Prob_pipnn_pn=.083*(1.-ppcnt)*(1.-ppcnt);
	    if (rnd->RndFsi().Rndm()*(Prob_pipd_pp+Prob_pipnn_pn)<Prob_pipd_pp){
	                       t1code=kPdgNeutron; t2code=kPdgProton; 
	                       scode=kPdgProton;   s2code=kPdgProton;}
	    else{
	                       t1code=kPdgNeutron; t2code=kPdgNeutron; 
	                       scode=kPdgProton;   s2code=kPdgNeutron;}
	  }
	  if (pdgc==kPdgPiM) { 
            double Prob_pimd_nn=2.*ppcnt*(1.-ppcnt);
            double Prob_pimpp_pn=.083*ppcnt*ppcnt;
	    if (rnd->RndFsi().Rndm()*(Prob_pimd_nn+Prob_pimpp_pn)<Prob_pimd_nn){
                               t1code=kPdgProton;  t2code=kPdgNeutron; 
	                       scode=kPdgNeutron;  s2code=kPdgNeutron; }
	    else{
	                       t1code=kPdgProton;  t2code=kPdgProton; 
	                       scode=kPdgProton;   s2code=kPdgNeutron;}
	  }
	  else { // pi0
            double Prob_pi0d_pn=0.88*ppcnt*(1.-ppcnt); // 2 * .44
            double Prob_pi0pp_pp=.14*ppcnt*ppcnt;
            double Prob_pi0nn_nn=.14*(1.-ppcnt)*(1.-ppcnt);
	    if (rnd->RndFsi().Rndm()*(Prob_pi0d_pn+Prob_pi0pp_pp+Prob_pi0nn_nn)<Prob_pi0d_pn){
	                       t1code=kPdgNeutron;  t2code=kPdgProton; 
	                        scode=kPdgNeutron;  s2code=kPdgProton;  }
	    else if (rnd->RndFsi().Rndm()*(Prob_pi0d_pn+Prob_pi0pp_pp+Prob_pi0nn_nn)<(Prob_pi0d_pn+Prob_pi0pp_pp)){
	                       t1code=kPdgProton;   t2code=kPdgProton; 
	                       scode=kPdgProton;    s2code=kPdgProton;  }
	    else {
	                       t1code=kPdgNeutron;  t2code=kPdgNeutron; 
	                       scode=kPdgNeutron;   s2code=kPdgNeutron;  }
	  }
	  LOG("HAIntranuke2018",pINFO) << "choose 2 body absorption, probe, fs = " << pdgc <<"  "<< scode <<"  "<<s2code;
	  // assign proper masses
	  //double M1   = pLib->Find(pdgc) ->Mass();
	  double M2_1 = pLib->Find(t1code)->Mass();
	  double M2_2 = pLib->Find(t2code)->Mass();
	  //double M2   = M2_1 + M2_2;
	  double M3   = pLib->Find(scode) ->Mass();
	  double M4   = pLib->Find(s2code)->Mass();

	  // handle fermi momentum 
	  double E2_1L, E2_2L;
	  TVector3 tP2_1L, tP2_2L;
	  //TLorentzVector dNucl_P4;
	  Target target(ev->TargetNucleus()->Pdg());
	  if(fDoFermi)
	    {
	      target.SetHitNucPdg(t1code);
	      fNuclmodel->GenerateNucleon(target);
	      //LOG("HAIntranuke2018", pNOTICE) << "Nuclmodel= " << fNuclmodel->ModelType(target) ;
	      tP2_1L=fFermiFac * fNuclmodel->Momentum3();
	      E2_1L = TMath::Sqrt(tP2_1L.Mag2() + M2_1*M2_1);
 
	      target.SetHitNucPdg(t2code);
	      fNuclmodel->GenerateNucleon(target);
	      tP2_2L=fFermiFac * fNuclmodel->Momentum3();
	      E2_2L = TMath::Sqrt(tP2_2L.Mag2() + M2_2*M2_2);

	      //dNucl_P4=TLorentzVector(tP2_1L+tP2_2L,E2_1L+E2_2L);
	    }
	  else
	    {
	      tP2_1L.SetXYZ(0.0, 0.0, 0.0);
	      E2_1L = M2_1;

	      tP2_2L.SetXYZ(0.0, 0.0, 0.0);
	      E2_2L = M2_2;
	    }
	  TLorentzVector dNucl_P4=TLorentzVector(tP2_1L+tP2_2L,E2_1L+E2_2L);

	  double E2L = E2_1L + E2_2L;

	  // adjust p to reflect scattering
	  // get random scattering angle
	  double C3CM = fHadroData2018->IntBounce(p,t1code,scode,fate_hN);
	  if (C3CM<-1.) 
	    {
	      LOG("HAIntranuke2018", pWARN) << "Inelastic() failed: IntBounce returned bad angle - try again";
	      exceptions::INukeException exception;
	      exception.SetReason("unphysical angle for hN scattering");
	      throw exception;
	      return;
	    }

	  TLorentzVector t4P1L,t4P2L,t4P3L,t4P4L;
	  t4P1L=*p->P4();
	  t4P2L=TLorentzVector(TVector3(tP2_1L+tP2_2L),E2L);
	  double bindE=0.050; // set to fit McKeown data, updated aug 18
	  //double bindE=0.0; 
	  if (utils::intranuke2018::TwoBodyKinematics(M3,M4,t4P1L,t4P2L,t4P3L,t4P4L,C3CM,fRemnP4,bindE))
	    {
	      //construct remnant nucleus and its mass

	      if (pdgc==kPdgPiP || pdgc==kPdgKP) fRemnZ++;
	      if (pdgc==kPdgPiM || pdgc==kPdgKM) fRemnZ--;
	      if (t1code==kPdgProton) fRemnZ--;
	      if (t2code==kPdgProton) fRemnZ--;
	      fRemnA-=2;

	      fRemnP4-=dNucl_P4;

	      TParticlePDG * remn = 0;
	      double MassRem = 0.;
	      int ipdgc = pdg::IonPdgCode(fRemnA, fRemnZ);
	      remn = PDGLibrary::Instance()->Find(ipdgc);
	      if(!remn) 
		{
		  LOG("HAIntranuke2018", pINFO)
		    << "NO Particle with [A = " << fRemnA << ", Z = " << fRemnZ
		    << ", pdgc = " << ipdgc << "] in PDGLibrary!";
		}
	      else
		{
		  MassRem = remn->Mass();
		  LOG("HAIntranuke2018", pINFO)
		    << "Particle with [A = " << fRemnA << ", Z = " << fRemnZ
		    << ", pdgc = " << ipdgc << "] in PDGLibrary!";
		}
	      double ERemn = fRemnP4.E();
	      double PRemn = TMath::Sqrt(fRemnP4.Px()*fRemnP4.Px() + fRemnP4.Py()*fRemnP4.Py() + fRemnP4.Pz()*fRemnP4.Pz());
	      double MRemn = TMath::Sqrt(ERemn*ERemn - PRemn*PRemn);
	      LOG("HAIntranuke2018",pINFO) << "PRemn = " << PRemn << " ;ERemn=  " << ERemn;
	      LOG("HAIntranuke2018",pINFO) << "expt MRemn=  " << MRemn << "  ;true Mass=  " << MassRem << "   ; excitation energy (>0 good)= " << (MRemn-MassRem)*1000. << " MeV";

	      // create t particles w/ appropriate momenta, code, and status
	      // Set target's mom to be the mom of the hadron that was cloned
	      GHepParticle* t1 = new GHepParticle(*p);
	      GHepParticle* t2 = new GHepParticle(*p);
	      t1->SetFirstMother(p->FirstMother());
	      t1->SetLastMother(p->LastMother());
	      t2->SetFirstMother(p->FirstMother());
	      t2->SetLastMother(p->LastMother());

	      // adjust p to reflect scattering
	      t1->SetPdgCode(scode);
	      t1->SetMomentum(t4P3L);

	      t2->SetPdgCode(s2code);
	      t2->SetMomentum(t4P4L);

	      t1->SetStatus(kIStStableFinalState);
	      t2->SetStatus(kIStStableFinalState);

	      ev->AddParticle(*t1);
	      ev->AddParticle(*t2);
	      delete t1;
	      delete t2;

	      return;
	    }
	  else
	    {
	      LOG("HAIntranuke2018", pNOTICE) << "Inelastic in hA failed calling TwoBodyKineamtics";
	      exceptions::INukeException exception;
	      exception.SetReason("Pion absorption kinematics through TwoBodyKinematics failed");
	      throw exception;

	    }

	} // end pi d -> N N
      else // multi-nucleon
	{
	  
      // declare some parameters for double gaussian and determine values chosen
      // parameters for proton and pi+, others come from isospin transformations

	  double ns0=0; // mean - sum of nucleons
	  double nd0=0; // mean - difference of nucleons
	  double Sig_ns=0; // std dev - sum
	  double Sig_nd=0; // std dev - diff
      double gam_ns=0; // exponential decay rate (for nucleons)

      if ( pdg::IsNeutronOrProton (pdgc) ) // nucleon probe
	{
	  // antisymmetric about Z=N
	  if (fRemnA-fRemnZ > fRemnZ) 
	    nd0 =  135.227 * TMath::Exp(-7.124*(fRemnA-fRemnZ)/double(fRemnA)) - 2.762;
	  else 
	    nd0 = -135.227 * TMath::Exp(-7.124*        fRemnZ /double(fRemnA)) + 4.914; 

	  Sig_nd = 2.034 + fRemnA * 0.007846;

	  double c1 = 0.041 + ke * 0.0001525;
	  double c2 = -0.003444 - ke * 0.00002324;
//change last factor from 30 to 15 so that gam_ns always larger than 0
//add check to be certain
	  double c3 = 0.064 - ke * 0.000015;  
	  gam_ns = c1 * TMath::Exp(c2*fRemnA) + c3;
	  if(gam_ns<0.002) gam_ns = 0.002;
	  //gam_ns = 10.;
	  LOG("HAIntranuke2018", pINFO) << "nucleon absorption";
	  LOG("HAIntranuke2018", pINFO) << "--> mean diff distr = " << nd0 << ", stand dev = " << Sig_nd;
	  //	  LOG("HAIntranuke2018", pINFO) << "--> mean sum distr = " << ns0 << ", Stand dev = " << Sig_ns;
	  LOG("HAIntranuke2018", pINFO) << "--> gam_ns = " << gam_ns;
	}
      else if ( pdgc==kPdgPiP || pdgc==kPdgPi0 || pdgc==kPdgPiM) //pion probe
	{
	  ns0 = .0001*(1.+ke/250.) * (fRemnA-10)*(fRemnA-10) + 3.5;
	  nd0 = (1.+ke/250.) - ((fRemnA/200.)*(1. + 2.*ke/250.));
	  Sig_ns = (10. + 4. * ke/250.)*TMath::Power(fRemnA/250.,0.9);  //(1. - TMath::Exp(-0.02*fRemnA));
	  Sig_nd = 4*(1 - TMath::Exp(-0.03*ke));
	  LOG("HAIntranuke2018", pINFO) << "pion absorption";
	  LOG("HAIntranuke2018", pINFO) << "--> mean diff distr = " << nd0 << ", stand dev = " << Sig_nd;
	  LOG("HAIntranuke2018", pINFO) << "--> mean sum distr = " << ns0 << ", Stand dev = " << Sig_ns;
	}
      else if (pdgc==kPdgKP || pdgc==kPdgKM) // kaon probe
	{
	  ns0 = (rnd->RndFsi().Rndm()>0.5?3:2);
	  nd0 = 1.;
	  Sig_ns = 0.1;
	  Sig_nd = 0.1;
	  LOG("HAIntranuke2018", pINFO) << "kaon absorption - set ns, nd later";
	  //	  LOG("HAIntranuke2018", pINFO) << "--> mean diff distr = " << nd0 << ", stand dev = " << Sig_nd;
	  //	  LOG("HAIntranuke2018", pINFO) << "--> mean sum distr = " << ns0 << ", Stand dev = " << Sig_ns;
	}
      else
	{
	  LOG("HAIntranuke2018", pWARN) << "Inelastic() cannot handle absorption reaction for " << p->Name();
	}

      // account for different isospin
      if (pdgc==kPdgPi0 || pdgc==kPdgNeutron) nd0-=2.;
      if (pdgc==kPdgPiM)                      nd0-=4.;

      int iter=0;    // counter
      int np=0,nn=0; // # of p, # of n
      bool not_done=true;
      double u1 = 0, u2 = 0;

      while (not_done)
	{
	  // infinite loop check
	  if (iter>=100) {
	    LOG("HAIntranuke2018", pNOTICE) << "Error: could not choose absorption final state";
	    LOG("HAIntranuke2018", pNOTICE) << "--> mean diff distr = " << nd0 << ", stand dev = " << Sig_nd;
	    LOG("HAIntranuke2018", pNOTICE) << "--> mean sum distr = " << ns0 << ", Stand dev = " << Sig_ns;
	    LOG("HAIntranuke2018", pNOTICE) << "--> gam_ns = " << gam_ns;
	    LOG("HAIntranuke2018", pNOTICE) << "--> A = " << fRemnA << ", Z = " << fRemnZ << ", Energy = " << ke;
	    exceptions::INukeException exception;
	    exception.SetReason("Absorption choice of # of p,n failed");
	    throw exception;
	  }
	  //here??

	  // Box-Muller transform
	  // Takes two uniform distribution random variables on (0,1]
	  // Creates two normally distributed random variables on (0,inf)

	  u1 = rnd->RndFsi().Rndm(); // uniform random variable 1
	  u2 = rnd->RndFsi().Rndm(); // "                     " 2
	  if (u1==0) u1 = rnd->RndFsi().Rndm();
	  if (u2==0) u2 = rnd->RndFsi().Rndm(); // Just in case

	  // normally distributed random variable   
	  double x2 = TMath::Sqrt(-2*TMath::Log(u1))*TMath::Sin(2*kPi*u2);

	  double ns = 0;

	  if ( pdg::IsNeutronOrProton (pdgc) ) //nucleon probe
	    {
	      ns = -TMath::Log(rnd->RndFsi().Rndm())/gam_ns; // exponential random variable
	    }	  
	  if ( pdg::IsKaon (pdgc) ) //charged kaon probe - either 2 or 3 nucleons to stay simple
	    {
	      ns =  (rnd->RndFsi().Rndm()<0.5?2:3);
	    }
	  else if ( pdgc==kPdgPiP || pdgc==kPdgPi0 || pdgc==kPdgPiM) //pion probe
	    {
	      // Pion fit for sum takes for xs*exp((xs-x0)^2 / 2*sig_xs0)
	      // Find random variable by first finding gaussian random variable
	      //   then throwing the value against a linear P.D.F.
	      //
	      // max is the maximum value allowed for the random variable (10 std + mean)
	      // minimum allowed value is 0

	      double max = ns0 + Sig_ns * 10;
	      if(max>fRemnA) max=fRemnA;
	      double x1 = 0;
	      bool not_found = true;
	      int iter2 = 0;

	      while (not_found)
		{
		  // infinite loop check
		  if (iter2>=100)
		    {
		      LOG("HAIntranuke2018", pNOTICE) << "Error: stuck in random variable loop for ns";
		      LOG("HAIntranuke2018", pNOTICE) << "--> mean of sum parent distr = " << ns0 << ", Stand dev = " << Sig_ns;
		      LOG("HAIntranuke2018", pNOTICE) << "--> A = " << fRemnA << ", Z = " << fRemnZ << ", Energy = " << ke;

		      exceptions::INukeException exception;
		      exception.SetReason("Random number generator for choice of #p,n final state failed - unusual - redo kinematics");
		      throw exception;
		    }

		  // calculate exponential random variable
		  u1 = rnd->RndFsi().Rndm();
		  u2 = rnd->RndFsi().Rndm();
		  if (u1==0) u1 = rnd->RndFsi().Rndm();
		  if (u2==0) u2 = rnd->RndFsi().Rndm();
		  x1 = TMath::Sqrt(-2*TMath::Log(u1))*TMath::Cos(2*kPi*u2);  

		  ns = ns0 + Sig_ns * x1;
		  if ( ns>max || ns<0 )                  {iter2++; continue;}
		  else if ( rnd->RndFsi().Rndm() > (ns/max) ) {iter2++; continue;}
		  else {
		    // accept this sum value
		    not_found=false;
		  }
		} //while(not_found)
	    }//else pion

	  double nd = nd0 + Sig_nd * x2; // difference (p-n) for both pion, nucleon probe
	  if (pdgc==kPdgKP)   // special for KP
	    { if (ns==2) nd=0;
	      if (ns>2) nd=1; }

	  np = int((ns+nd)/2.+.5); // Addition of .5 for rounding correction
	  nn = int((ns-nd)/2.+.5);

	  LOG("HAIntranuke2018", pINFO) << "ns = "<<ns<<", nd = "<<nd<<", np = "<<np<<", nn = "<<nn;
	  //LOG("HAIntranuke2018", pNOTICE) << "RemA = "<<fRemnA<<", RemZ = "<<fRemnZ<<", probe = "<<pdgc;

	  /*if ((ns+nd)/2. < 0 || (ns-nd)/2. < 0)  {iter++; continue;}
	    else */ 
	  //check for problems befor executing phase space 'decay'
	       if (np < 0 || nn < 0 )                 {iter++; continue;}
          else if (np + nn < 2. )                     {iter++; continue;}
          else if ((np + nn == 2.) &&  pdg::IsNeutronOrProton (pdgc))                     {iter++; continue;}
          else if (np > fRemnZ + ((pdg::IsProton(pdgc) || pdgc==kPdgPiP || pdgc==kPdgKP)?1:0)
		   - ((pdgc==kPdgPiM || pdgc==kPdgKM)?1:0)) {iter++; continue;}
          else if (nn > fRemnA-fRemnZ + ((pdg::IsNeutron(pdgc)||pdgc==kPdgPiM||pdgc==kPdgKM)?1:0)
		   - ((pdgc==kPdgPiP||pdgc==kPdgKP)?1:0)) {iter++; continue;}
	  else { 
	    not_done=false;   //success
	    LOG("HAIntranuke2018",pINFO) << "success, iter = " << iter << "  np, nn = " << np << "  " << nn; 
	    if (np+nn>86) // too many particles, scale down
	      {
		double frac = 85./double(np+nn);
		np = int(np*frac);
		nn = int(nn*frac);
	      }

	    if (  (np==fRemnZ       +((pdg::IsProton (pdgc)||pdgc==kPdgPiP||pdgc==kPdgKP)?1:0)-(pdgc==kPdgPiM||pdgc==kPdgKM?1:0))
		&&(nn==fRemnA-fRemnZ+((pdg::IsNeutron(pdgc)||pdgc==kPdgPiM||pdgc==kPdgKM)?1:0)-(pdgc==kPdgPiP||pdgc==kPdgKP?1:0)) )
	      { // leave at least one nucleon in the nucleus to prevent excess momentum
		if (rnd->RndFsi().Rndm()<np/(double)(np+nn)) np--;
		else nn--;
	      }

	    LOG("HAIntranuke2018", pNOTICE) << "Final state chosen; # protons : " 
					<< np << ", # neutrons : " << nn;
	  }
	} //while(not_done)

      // change remnants to reflect probe
      if ( pdgc==kPdgProton || pdgc==kPdgPiP || pdgc==kPdgKP)     fRemnZ++;
      if ( pdgc==kPdgPiM || pdgc==kPdgKM)                         fRemnZ--;
      if ( pdg::IsNeutronOrProton (pdgc) )         fRemnA++;

      // PhaseSpaceDecay forbids anything over 18 particles
      //
      // If over 18 particles, split into 5 groups and perform decay on each group
      // Anything over 85 particles reduced to 85 in previous step
      //
      // 4 of the nucleons are used as "probes" as well as original probe,
      // with each one sharing 1/5 of the total incident momentum
      //
      // Note: choosing 5 groups and distributing momentum evenly was arbitrary
      //   Needs to be revised later

      if ((np+nn)>18)
	{
	  // code lists
	  PDGCodeList list0(allow_dup);
	  PDGCodeList list1(allow_dup);
	  PDGCodeList list2(allow_dup);
	  PDGCodeList list3(allow_dup);
	  PDGCodeList list4(allow_dup);
	  PDGCodeList* listar[5] = {&list0, &list1, &list2, &list3, &list4};

	  //set up HadronClusters
	  // simple for now, each (of 5) in hadron cluster has 1/5 of mom and KE

       	  double probM = pLib->Find(pdgc)   ->Mass();
	  probM -= .025;   // BE correction
	  TVector3 pP3 = p->P4()->Vect() * (1./5.);
	  double probKE = p->P4()->E() -probM;
	  double clusKE = probKE * (1./5.);
	  TLorentzVector clusP4(pP3,clusKE);   //no mass
	  LOG("HAIntranuke2018",pINFO) << "probM = " << probM << " ;clusKE=  " << clusKE;
	  TLorentzVector X4(*p->X4());
	  GHepStatus_t ist = kIStNucleonClusterTarget;

	  int mom = p->FirstMother();

	  GHepParticle * p0 = new GHepParticle(kPdgCompNuclCluster,ist, mom,-1,-1,-1,clusP4,X4);
	  GHepParticle * p1 = new GHepParticle(kPdgCompNuclCluster,ist, mom,-1,-1,-1,clusP4,X4);
	  GHepParticle * p2 = new GHepParticle(kPdgCompNuclCluster,ist, mom,-1,-1,-1,clusP4,X4);
	  GHepParticle * p3 = new GHepParticle(kPdgCompNuclCluster,ist, mom,-1,-1,-1,clusP4,X4);
	  GHepParticle * p4 = new GHepParticle(kPdgCompNuclCluster,ist, mom,-1,-1,-1,clusP4,X4);

	  // To conserve 4-momenta
	  //	  fRemnP4 -= probP4 + protP4*np_p + neutP4*(4-np_p) - *p->P4();
	  fRemnP4 -= 5.*clusP4 - *p->P4();

	  for (int i=0;i<(np+nn);i++)
	    {
	      if (i<np)
		{
		  listar[i%5]->push_back(kPdgProton);
		  fRemnZ--;
		}
	      else listar[i%5]->push_back(kPdgNeutron);
	      fRemnA--;
	    }
	  for (int i=0;i<5;i++)
	    {
	      LOG("HAIntranuke2018", pINFO) << "List" << i << " size: " << listar[i]->size();
	      if (listar[i]->size() <2)
		{
		  exceptions::INukeException exception;
		  exception.SetReason("too few particles for Phase Space decay - try again");
		  throw exception;
		}	  
	    }

	  // commented out to better fit with absorption reactions
	  // Add the fermi energy of the nucleons to the phase space
	  /*if(fDoFermi)
	    {  
	      GHepParticle* p_ar[5] = {cl, p1, p2, p3, p4};
	      for (int i=0;i<5;i++)
		{
		  Target target(ev->TargetNucleus()->Pdg());
		  TVector3 pBuf = p_ar[i]->P4()->Vect();
		  double mBuf = p_ar[i]->Mass();
		  double eBuf = TMath::Sqrt(pBuf.Mag2() + mBuf*mBuf);
		  TLorentzVector tSum(pBuf,eBuf); 
		  double mSum = 0.0;
		  vector<int>::const_iterator pdg_iter;
		  for(pdg_iter=++(listar[i]->begin());pdg_iter!=listar[i]->end();++pdg_iter)
		    {
		      target.SetHitNucPdg(*pdg_iter); 
		      fNuclmodel->GenerateNucleon(target);
		      mBuf = pLib->Find(*pdg_iter)->Mass();
		      mSum += mBuf;
		      pBuf = fFermiFac * fNuclmodel->Momentum3();
		      eBuf = TMath::Sqrt(pBuf.Mag2() + mBuf*mBuf);
		      tSum += TLorentzVector(pBuf,eBuf);
		      fRemnP4 -= TLorentzVector(pBuf,eBuf-mBuf);
		    }
		  TLorentzVector dP4 = tSum + TLorentzVector(TVector3(0,0,0),-mSum);
		  p_ar[i]->SetMomentum(dP4);
		}
		}*/

	  bool success1 = utils::intranuke2018::PhaseSpaceDecay(ev,p0,*listar[0],fRemnP4,fNucRmvE,kIMdHA);
	  bool success2 = utils::intranuke2018::PhaseSpaceDecay(ev,p1,*listar[1],fRemnP4,fNucRmvE,kIMdHA);
	  bool success3 = utils::intranuke2018::PhaseSpaceDecay(ev,p2,*listar[2],fRemnP4,fNucRmvE,kIMdHA);
	  bool success4 = utils::intranuke2018::PhaseSpaceDecay(ev,p3,*listar[3],fRemnP4,fNucRmvE,kIMdHA);
	  bool success5 = utils::intranuke2018::PhaseSpaceDecay(ev,p4,*listar[4],fRemnP4,fNucRmvE,kIMdHA);
	  if(success1 && success2 && success3 && success4 && success5)
	    {
	      LOG("HAIntranuke2018", pINFO)<<"Successful many-body absorption - n>=18";
	      LOG("HAIntranuke2018", pDEBUG) << "Nucleus : (A,Z) = ("<<fRemnA<<','<<fRemnZ<<')';
	      TParticlePDG * remn = 0;
	      double MassRem = 0.;
	      int ipdgc = pdg::IonPdgCode(fRemnA, fRemnZ);
	      remn = PDGLibrary::Instance()->Find(ipdgc);
	      if(!remn) 
		{
		  LOG("HAIntranuke2018", pINFO)
		    << "NO Particle with [A = " << fRemnA << ", Z = " << fRemnZ
		    << ", pdgc = " << ipdgc << "] in PDGLibrary!";
		}
	      else
		{
		  MassRem = remn->Mass();
		  LOG("HAIntranuke2018", pINFO)
		    << "Particle with [A = " << fRemnA << ", Z = " << fRemnZ
		    << ", pdgc = " << ipdgc << "] in PDGLibrary!";
		}
	      double ERemn = fRemnP4.E();
	      double PRemn = TMath::Sqrt(fRemnP4.Px()*fRemnP4.Px() + fRemnP4.Py()*fRemnP4.Py() + fRemnP4.Pz()*fRemnP4.Pz());
	      double MRemn = TMath::Sqrt(ERemn*ERemn - PRemn*PRemn);
	      LOG("HAIntranuke2018",pINFO) << "PRemn = " << PRemn << " ;ERemn=  " << ERemn;
	      LOG("HAIntranuke2018",pINFO) << "MRemn=  " << MRemn << "  ;true Mass=  " << MassRem << "   ; excitation energy (>0 good)= " << (MRemn-MassRem)*1000. << " MeV";

	    }
	  else 
	    {
	      // try to recover 
	      LOG("HAIntranuke2018", pWARN) << "PhaseSpace decay fails for HadrCluster- recovery likely incorrect - rethrow event";
	      p->SetStatus(kIStStableFinalState);
	      ev->AddParticle(*p);
	      fRemnA+=np+nn;
	      fRemnZ+=np;
	      if ( pdgc==kPdgProton || pdgc==kPdgPiP )     fRemnZ--;
	      if ( pdgc==kPdgPiM )                         fRemnZ++;
	      if ( pdg::IsNeutronOrProton (pdgc) )         fRemnA--;	
		/*	      exceptions::INukeException exception;
	      exception.SetReason("Phase space generation of absorption final state failed");
	      throw exception;
		*/
	    }

	  //	  delete cl;
	  delete p0;
	  delete p1;
	  delete p2;
	  delete p3;
	  delete p4;

	}
      else // less than 18 particles pion
	{
	  if (pdgc==kPdgKP)  list.push_back(kPdgKP); //normally conserve strangeness
	  if (pdgc==kPdgKM)  list.push_back(kPdgKM);
	  /*
	  TParticlePDG * remn0 = 0;
	  int ipdgc0 = pdg::IonPdgCode(fRemnA, fRemnZ);
	  remn0 = PDGLibrary::Instance()->Find(ipdgc0);
	  double Mass0 = remn0->Mass();
	  TParticlePDG * remnt = 0;
	  int ipdgct = pdg::IonPdgCode(fRemnA-(nn+np), fRemnZ-np);
	  remnt = PDGLibrary::Instance()->Find(ipdgct);
	  double MassRemt = remnt->Mass();
	  LOG("HAIntranuke2018",pINFO) << "Mass0 = " << Mass0 << " ;Masst=  " << MassRemt << "  ; diff/nucleon= "<< (Mass0-MassRemt)/(np+nn);
	  */
	  //set up HadronCluster

       	  double probM = pLib->Find(pdgc)   ->Mass();
	  double probBE = (np+nn)*.005;   // BE correction
	  TVector3 pP3 = p->P4()->Vect();
	  double probKE = p->P4()->E() - (probM - probBE);
	  double clusKE = probKE;  // + np*0.9383 + nn*.9396;
	  TLorentzVector clusP4(pP3,clusKE);   //no mass is correct
	  LOG("HAIntranuke2018",pINFO) << "probM = " << probM << " ;clusKE=  " << clusKE;
	  TLorentzVector X4(*p->X4());
	  GHepStatus_t ist = kIStNucleonClusterTarget;
	  int mom = p->FirstMother();

	  GHepParticle * p0 = new GHepParticle(kPdgCompNuclCluster,ist, mom,-1,-1,-1,clusP4,X4);

	  //set up remnant nucleus
	  fRemnP4 -= clusP4 - *p->P4();

	  for (int i=0;i<np;i++)
	    {
	      list.push_back(kPdgProton);
	      fRemnA--;
	      fRemnZ--;
	    }
	  for (int i=0;i<nn;i++)
	    {
	      list.push_back(kPdgNeutron);
	      fRemnA--;
	    }
	  
	  // Library instance for reference
	  //PDGLibrary * pLib = PDGLibrary::Instance();

	  // commented out to better fit with absorption reactions	  
	  // Add the fermi energy of the nucleons to the phase space
	  /*if(fDoFermi)
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
		  fNuclmodel->GenerateNucleon(target);
		  mBuf = pLib->Find(*pdg_iter)->Mass();
		  mSum += mBuf;
		  pBuf = fFermiFac * fNuclmodel->Momentum3();
		  eBuf = TMath::Sqrt(pBuf.Mag2() + mBuf*mBuf);
		  tSum += TLorentzVector(pBuf,eBuf);
		  fRemnP4 -= TLorentzVector(pBuf,eBuf-mBuf);
		}
	      TLorentzVector dP4 = tSum + TLorentzVector(TVector3(0,0,0),-mSum);
	      p->SetMomentum(dP4);    
	      }*/
	  
	  LOG("HAIntranuke2018", pDEBUG)
	    << "Remnant nucleus (A,Z) = (" << fRemnA << ", " << fRemnZ << ")";
	  LOG("HAIntranuke2018", pINFO) << " list size: " << np+nn;
	  if (np+nn <2)
	    {
	      exceptions::INukeException exception;
	      exception.SetReason("too few particles for Phase Space decay - try again");
	      throw exception;
	    }
	  //	  GHepParticle * cl = new GHepParticle(*p);
	  //	  cl->SetPdgCode(kPdgDecayNuclCluster);
     	  //bool success1 = utils::intranuke2018::PhaseSpaceDecay(ev,p0,*listar[0],fRemnP4,fNucRmvE,kIMdHA);
	  bool success = utils::intranuke2018::PhaseSpaceDecay(ev,p0,list,fRemnP4,fNucRmvE,kIMdHA);
	  if (success)
	    {
	      LOG ("HAIntranuke2018",pINFO) << "Successful many-body absorption, n<=18";
	      LOG("HAIntranuke2018", pDEBUG) << "Nucleus : (A,Z) = ("<<fRemnA<<','<<fRemnZ<<')';
	      TParticlePDG * remn = 0;
	      double MassRem = 0.;
	      int ipdgc = pdg::IonPdgCode(fRemnA, fRemnZ);
	      remn = PDGLibrary::Instance()->Find(ipdgc);
	      if(!remn) 
		{
		  LOG("HAIntranuke2018", pINFO)
		    << "NO Particle with [A = " << fRemnA << ", Z = " << fRemnZ
		    << ", pdgc = " << ipdgc << "] in PDGLibrary!";
		}
	      else
		{
		  MassRem = remn->Mass();
		  LOG("HAIntranuke2018", pINFO)
		    << "Particle with [A = " << fRemnA << ", Z = " << fRemnZ
		    << ", pdgc = " << ipdgc << "] in PDGLibrary!";
		}
	      double ERemn = fRemnP4.E();
	      double PRemn = TMath::Sqrt(fRemnP4.Px()*fRemnP4.Px() + fRemnP4.Py()*fRemnP4.Py() + fRemnP4.Pz()*fRemnP4.Pz());
	      double MRemn = TMath::Sqrt(ERemn*ERemn - PRemn*PRemn);
	      LOG("HAIntranuke2018",pINFO) << "PRemn = " << PRemn << " ;ERemn=  " << ERemn;
	      LOG("HAIntranuke2018",pINFO) << "expt MRemn=  " << MRemn << "  ;true Mass=  " << MassRem << "   ; excitation energy (>0 good)= " << (MRemn-MassRem)*1000. << " MeV";
	    }
	  else {
	    // recover
	    p->SetStatus(kIStStableFinalState);
	    ev->AddParticle(*p);
	    fRemnA+=np+nn;
	    fRemnZ+=np;
	    if ( pdgc==kPdgProton || pdgc==kPdgPiP )     fRemnZ--;
	    if ( pdgc==kPdgPiM )                         fRemnZ++;
	    if ( pdg::IsNeutronOrProton (pdgc) )         fRemnA--;	
	    exceptions::INukeException exception;
	    exception.SetReason("Phase space generation of absorption final state failed");
	    throw exception;
	  }
	  delete p0;
	}
	} // end multi-nucleon FS
    }
  else // not absorption/pipro
    {
      LOG("HAIntranuke2018", pWARN) 
	<< "Inelastic() can not handle fate: " << INukeHadroFates::AsString(fate);
      return;
    }
}
//___________________________________________________________________________
int HAIntranuke2018::HandleCompoundNucleus(GHepRecord* /*ev*/, GHepParticle* /*p*/, int /*mom*/) const
{

  // only relevant for hN mode     - not anymore.
  return false;
  
}
//___________________________________________________________________________
void HAIntranuke2018::LoadConfig(void)
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
  GetParam( "INUKE-NucCEXFac",     fNucCEXFac ) ;
  GetParam( "INUKE-Energy_Pre_Eq", fEPreEq ) ;
  GetParam( "INUKE-FermiFac",      fFermiFac ) ;
  GetParam( "INUKE-FermiMomentum", fFermiMomentum ) ;

  GetParam( "INUKE-DoCompoundNucleus", fDoCompoundNucleus ) ;
  GetParam( "INUKE-DoFermi",           fDoFermi ) ;
  GetParam( "INUKE-XsecNNCorr",        fXsecNNCorr ) ;
  GetParamDef( "UseOset",              fUseOset, false ) ;
  GetParamDef( "AltOset",              fAltOset, false ) ;

  GetParam( "HAINUKE-DelRPion",    fDelRPion ) ;
  GetParam( "HAINUKE-DelRNucleon", fDelRNucleon ) ;

  GetParamDef( "FSI-Pion-MFPScale",              fPionMFPScale,           1.0 ) ;
  GetParamDef( "FSI-Pion-FracCExScale",          fPionFracCExScale,       1.0 ) ;
  GetParamDef( "FSI-Pion-FracAbsScale",          fPionFracAbsScale,       1.0 ) ;
  GetParamDef( "FSI-Pion-FracInelScale",         fPionFracInelScale,      1.0 ) ;
  GetParamDef( "FSI-Pion-FracPiProdScale",       fPionFracPiProdScale,    1.0 ) ;
  GetParamDef( "FSI-Nucleon-MFPScale",           fNucleonMFPScale,        1.0 ) ;
  GetParamDef( "FSI-Nucleon-FracCExScale",       fNucleonFracCExScale ,   1.0 ) ;
  GetParamDef( "FSI-Nucleon-FracInelScale",      fNucleonFracInelScale,   1.0 ) ;
  GetParamDef( "FSI-Nucleon-FracAbsScale",       fNucleonFracAbsScale,    1.0 ) ;
  GetParamDef( "FSI-Nucleon-FracPiProdScale",    fNucleonFracPiProdScale, 1.0 ) ;

    // report
  LOG("HAIntranuke2018", pINFO) << "Settings for INTRANUKE mode: " << INukeMode::AsString(kIMdHA);
  LOG("HAIntranuke2018", pINFO) << "R0          = " << fR0 << " fermi";
  LOG("HAIntranuke2018", pINFO) << "NR          = " << fNR;
  LOG("HAIntranuke2018", pINFO) << "DelRPion    = " << fDelRPion;
  LOG("HAIntranuke2018", pINFO) << "DelRNucleon = " << fDelRNucleon;
  LOG("HAIntranuke2018", pINFO) << "HadStep     = " << fHadStep << " fermi";
  LOG("HAIntranuke2018", pINFO) << "EPreEq      = " << fHadStep << " fermi";
  LOG("HAIntranuke2018", pINFO) << "NucAbsFac   = " << fNucAbsFac;
  LOG("HAIntranuke2018", pINFO) << "NucCEXFac   = " << fNucCEXFac;
  LOG("HAIntranuke2018", pINFO) << "FermiFac    = " << fFermiFac;
  LOG("HAIntranuke2018", pINFO) << "FermiMomtm  = " << fFermiMomentum;
  LOG("HAIntranuke2018", pINFO) << "DoFermi?    = " << ((fDoFermi)?(true):(false));
  LOG("HAIntranuke2018", pINFO) << "DoCmpndNuc? = " << ((fDoCompoundNucleus)?(true):(false));
  LOG("HAIntranuke2018", pINFO) << "XsecNNCorr? = " << ((fXsecNNCorr)?(true):(false));
}
//___________________________________________________________________________
/*
INukeFateHA_t HAIntranuke2018::HadronFateOset () const
{
  const double fractionAbsorption = osetUtils::currentInstance->
                                    getAbsorptionFraction();
  const double fractionCex = osetUtils::currentInstance->getCexFraction ();

  RandomGen *randomGenerator = RandomGen::Instance();
  const double randomNumber  = randomGenerator->RndFsi().Rndm();

  LOG("HAIntranuke2018", pINFO) 
    << "\n frac{" << INukeHadroFates::AsString(kIHAFtCEx)     << "} = " << fractionCex
    << "\n frac{" << INukeHadroFates::AsString(kIHAFtInelas)  << "} = " << 1-fractionCex-fractionAbsorption
    << "\n frac{" << INukeHadroFates::AsString(kIHAFtAbs)     << "} = " << fractionAbsorption;
  if (randomNumber < fractionAbsorption && fRemnA > 1) return kIHAFtAbs;
  else if (randomNumber < fractionAbsorption + fractionCex) return kIHAFtCEx;
  else return kIHAFtInelas;
}
*/
