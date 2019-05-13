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
 @ Jan 28, 2009 - CA
   The nuclear remnant is now marked as a kIStFinalStateNuclearRemnant, not
   as a kIStStableFinalState.
 @ Sep 15, 2009 - CA
   IsFake() and IsNucleus() are no longer available in GHepParticle.
   Use pdg::IsPseudoParticle() and pdg::IsIon().  
 @ Sep 15, 2009 - CA
   Store the rescattering code (hadron fate) in the GHEP record so as to
   facilitate event reweighting.
 @ Jul 15, 2010 - AM
   Split Intranuke class into two separate classes, one for each interaction mode.
   Intranuke.cxx now only contains methods common to both classes and associated
   with the stepping of the hadrons through the nucleus and with configuration.
 @ Nov 20, 2011 - CA
   Tweaked the way TransportHadrons() looks-up the nuclear environment so that
   it works for the nucleon decay mode as well.
 @ Dec 08, 2011 - CA
   Some minor structural changes. The new GEvGenMode_t is determined at the
   start of the event processing and is used throughout. fInTestMode flag and
   special INTRANUKE configs not needed. ProcessEventRecord() was added by 
   factoring out code from HNIntranuke and HAIntranuke. Some comments added.
 @ Dec 23, 2014 - TG, SD
   New 2014 class for latest Intranuke model
 @ Apr 26, 2018 - SD 
   Change year 2015 to 2018

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
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/HadronTransport/Intranuke2018.h"
#include "Physics/HadronTransport/INukeHadroData2018.h"
#include "Physics/HadronTransport/INukeHadroFates.h"
#include "Physics/HadronTransport/INukeMode.h"
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

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
Intranuke2018::Intranuke2018() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
Intranuke2018::Intranuke2018(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
Intranuke2018::Intranuke2018(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
Intranuke2018::~Intranuke2018()
{

}
//___________________________________________________________________________
void Intranuke2018::ProcessEventRecord(GHepRecord * evrec) const
{
  // Do not continue if there is no nuclear target
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("HNIntranuke2018", pINFO) << "No nuclear target found - INTRANUKE exits";
    return;
  }

  // Decide tracking radius for the current nucleus (few * R0 * A^1/3)
  this->SetTrackingRadius(nucltgt);

  // Understand what the event generation mode is (hadron/photon-nucleus,
  // lepton-nucleus, nucleon decay) from the input event.
  // The determined mode has an effect on INTRANUKE behaviour (how to lookup
  // the residual nucleus, whether to set an intranuclear vtx etc) but it
  // does not affect the INTRANUKE physics.
  fGMode = evrec->EventGenerationMode();

  // For lepton-nucleus scattering and for nucleon decay intranuclear vtx 
  // position (in the target nucleus coord system) is set elsewhere.
  // This method only takes effect in hadron/photon-nucleus interactions.
  // In this special mode, an interaction vertex is set at the periphery 
  // of the target nucleus.
  if(fGMode == kGMdHadronNucleus ||
     fGMode == kGMdPhotonNucleus)
  {
    this->GenerateVertex(evrec);  
  }

  // Now transport all hadrons outside the tracking radius.
  // Stepping part is common for both HA and HN.
  // Once it has been estabished that an interaction takes place then
  // HA and HN specific code takes over in order to simulate the final state.
  this->TransportHadrons(evrec);
}
//___________________________________________________________________________
void Intranuke2018::GenerateVertex(GHepRecord * evrec) const
{
// Sets a vertex in the nucleus periphery
// Called onlt in hadron/photon-nucleus interactions.

  GHepParticle * nucltgt = evrec->TargetNucleus();
  assert(nucltgt);

  RandomGen * rnd = RandomGen::Instance();
  TVector3 vtx(999999.,999999.,999999.);

  // *** For h+A events (test mode): 
  // Assume a hadron beam with uniform intensity across an area, 
  // so we need to choose events uniformly within that area.
  double x=999999., y=999999., epsilon = 0.001;
  double R2  = TMath::Power(fTrackingRadius,2.);
  double rp2 = TMath::Power(x,2.) + TMath::Power(y,2.);
  while(rp2 > R2-epsilon) {
      x = (fTrackingRadius-epsilon) * rnd->RndFsi().Rndm();
      y = -fTrackingRadius + 2*fTrackingRadius * rnd->RndFsi().Rndm();
      y -= ((y>0) ? epsilon : -epsilon);
      rp2 = TMath::Power(x,2.) + TMath::Power(y,2.);
  }
  vtx.SetXYZ(x,y, -1.*TMath::Sqrt(TMath::Max(0.,R2-rp2)) + epsilon);

  // get the actual unit vector along the incoming hadron direction
  TVector3 direction = evrec->Probe()->P4()->Vect().Unit();

  // rotate the vtx position
  vtx.RotateUz(direction);
  
  LOG("Intranuke2018", pNOTICE) 
     << "Generated vtx @ R = " << vtx.Mag() << " fm / " 
     << print::Vec3AsString(&vtx);

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(pdg::IsPseudoParticle(p->Pdg())) continue;
    if(pdg::IsIon           (p->Pdg())) continue;

    p->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);
  }
}
//___________________________________________________________________________
void Intranuke2018::SetTrackingRadius(const GHepParticle * p) const
{
  assert(p && pdg::IsIon(p->Pdg()));
  double A = p->A();
  fTrackingRadius = fR0 * TMath::Power(A, 1./3.);

  // multiply that by some input factor so that hadrons are tracked
  // beyond the nuclear 'boundary' since the nuclear density distribution
  // is not zero there
  fTrackingRadius *= fNR; 

  LOG("Intranuke2018", pNOTICE) 
      << "Setting tracking radius to R = " << fTrackingRadius;
}
//___________________________________________________________________________
bool Intranuke2018::NeedsRescattering(const GHepParticle * p) const
{
// checks whether the particle should be rescattered

  assert(p);

  if(fGMode == kGMdHadronNucleus ||
     fGMode == kGMdPhotonNucleus) {
    // hadron/photon-nucleus scattering propagate the incoming particle
    return (
      (p->Status() == kIStInitialState || p->Status() == kIStHadronInTheNucleus)
       && !pdg::IsIon(p->Pdg()));
  }
  else {
   // attempt to rescatter anything marked as 'hadron in the nucleus'
   return (p->Status() == kIStHadronInTheNucleus);
  }
}
//___________________________________________________________________________
bool Intranuke2018::CanRescatter(const GHepParticle * p) const
{
// checks whether a particle that needs to be rescattered, can in fact be 
// rescattered by this cascade MC

  assert(p);
  return  ( p->Pdg() == kPdgPiP     || 
            p->Pdg() == kPdgPiM     || 
            p->Pdg() == kPdgPi0     ||
            p->Pdg() == kPdgProton  ||
            p->Pdg() == kPdgNeutron ||
	    //	    p->Pdg() == kPdgGamma   ||
	    p->Pdg() == kPdgKP      //||
	    //	    p->Pdg() == kPdgKM
          );
}
//___________________________________________________________________________
bool Intranuke2018::IsInNucleus(const GHepParticle * p) const
{
// check whether the input particle is still within the nucleus
//
  return (p->X4()->Vect().Mag() < fTrackingRadius + fHadStep);
}
//___________________________________________________________________________
void Intranuke2018::TransportHadrons(GHepRecord * evrec) const
{
// transport all hadrons outside the nucleus

  int inucl = -1;
  fRemnA = -1;
  fRemnZ = -1;

  //  Get 'nuclear environment' at the beginning of hadron transport 
  //  and keep track of the remnant nucleus A,Z  

  if(fGMode == kGMdHadronNucleus ||
     fGMode == kGMdPhotonNucleus)
  {
     inucl = evrec->TargetNucleusPosition();
  }
  else if(fGMode == kGMdLeptonNucleus ||
	  fGMode == kGMdDarkMatterNucleus ||
	  fGMode == kGMdNucleonDecay) {
    inucl = evrec->RemnantNucleusPosition();
  }

  LOG("Intranuke2018", pNOTICE) 
     << "Propagating hadrons within nucleus found in position = " << inucl;
  GHepParticle * nucl = evrec->Particle(inucl);
  if(!nucl) {
    LOG("Intranuke2018", pERROR) 
       << "No nucleus found in position = " << inucl;
    LOG("Intranuke2018", pERROR) 
       << *evrec;
    return;
  }
  
  fRemnA = nucl->A();
  fRemnZ = nucl->Z();

  LOG("Intranuke2018", pNOTICE)
      << "Nucleus (A,Z) = (" << fRemnA << ", " << fRemnZ << ")";

  const TLorentzVector & p4nucl = *(nucl->P4());
  fRemnP4 = p4nucl; 

  // Loop over GHEP and run intranuclear rescattering on handled particles
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr = -1;

  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    // Check whether the particle needs rescattering, otherwise skip it
    if( ! this->NeedsRescattering(p) ) continue;
 
    LOG("Intranuke2018", pNOTICE)
      << " >> Stepping a " << p->Name() 
                        << " with kinetic E = " << p->KinE() << " GeV";

    // Rescatter a clone, not the original particle
    GHepParticle * sp = new GHepParticle(*p); 

    // Set clone's mom to be the hadron that was cloned
    sp->SetFirstMother(icurr); 

    // Check whether the particle can be rescattered 
    if(!this->CanRescatter(sp)) {

       // if I can't rescatter it, I will just take it out of the nucleus
       LOG("Intranuke2018", pNOTICE)
              << "... Current version can't rescatter a " << sp->Name();
       sp->SetFirstMother(icurr); 
       sp->SetStatus(kIStStableFinalState);
       evrec->AddParticle(*sp);
       delete sp;
       continue; // <-- skip to next GHEP entry
    }

    // Start stepping particle out of the nucleus
    bool has_interacted = false;
    while ( this-> IsInNucleus(sp) ) 
    {
      // advance the hadron by a step
      utils::intranuke2018::StepParticle(sp, fHadStep);

      // check whether it interacts
      double d = this->GenerateStep(evrec,sp);
      has_interacted = (d<fHadStep);
      if(has_interacted) break;
    }//stepping

    //updating the position of the original particle with the position of the clone
    evrec->Particle(sp->FirstMother())->SetPosition(*(sp->X4()));
 
    if(has_interacted && fRemnA>0)  {
        // the particle interacts - simulate the hadronic interaction
      LOG("Intranuke2018", pNOTICE) 
          << "Particle has interacted at location:  " 
          << sp->X4()->Vect().Mag() << " / nucl rad= " << fTrackingRadius;
	this->SimulateHadronicFinalState(evrec,sp);
    } else if(has_interacted && fRemnA<=0) {
        // nothing left to interact with!
      LOG("Intranuke2018", pNOTICE)
          << "*** Nothing left to interact with, escaping.";
	sp->SetStatus(kIStStableFinalState);
	evrec->AddParticle(*sp);
	evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
    } else {
        // the exits the nucleus without interacting - Done with it! 
        LOG("Intranuke2018", pNOTICE) 
          << "*** Hadron escaped the nucleus! Done with it.";
	sp->SetStatus(kIStStableFinalState);
	evrec->AddParticle(*sp);
	evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
    }
    delete sp;

    // Current snapshot
    //LOG("Intranuke2018", pINFO) << "Current event record snapshot: " << *evrec;

  }// GHEP entries

  // Add remnant nucleus - that 'hadronic blob' has all the remaining hadronic 
  // 4p not  put explicitly into the simulated particles
  TLorentzVector v4(0.,0.,0.,0.);
  GHepParticle remnant_nucleus(
    kPdgHadronicBlob, kIStFinalStateNuclearRemnant, inucl,-1,-1,-1, fRemnP4, v4);
  evrec->AddParticle(remnant_nucleus);
  // Mark the initial remnant nucleus as an intermediate state 
  // Don't do that in the hadron/photon-nucleus scatterig mode since the initial 
  // remnant nucleus and the target nucleus coincide.
  if(fGMode != kGMdHadronNucleus &&
     fGMode != kGMdPhotonNucleus) {
     evrec->Particle(inucl)->SetStatus(kIStIntermediateState);
  }
}
//___________________________________________________________________________
double Intranuke2018::GenerateStep(GHepRecord*  /*evrec*/, GHepParticle* p) const //Added ev to get tgt argument//
{
// Generate a step (in fermis) for particle p in the input event.
// Computes the mean free path L and generate an 'interaction' distance d 
// from an exp(-d/L) distribution

  int pdgc = p->Pdg();

  double scale = 1.;
  if (pdgc==kPdgPiP || pdgc==kPdgPiM || pdgc==kPdgPi0) {
    scale = fPionMFPScale;
  }
  else if (pdgc==kPdgProton || pdgc==kPdgNeutron) {
    scale = fNucleonMFPScale;
  }

  RandomGen * rnd = RandomGen::Instance();

  string fINukeMode = this->GetINukeMode();
  string fINukeModeGen = this->GetGenINukeMode();

  double L = utils::intranuke2018::MeanFreePath(p->Pdg(), *p->X4(), *p->P4(), fRemnA,
						fRemnZ, fDelRPion, fDelRNucleon, fUseOset, fAltOset, fXsecNNCorr, fINukeMode);

  LOG("Intranuke2018", pDEBUG)    << "mode= " << fINukeModeGen;
  if(fINukeModeGen == "hA") L *= scale;

  double d = -1.*L * TMath::Log(rnd->RndFsi().Rndm());

  /*    LOG("Intranuke2018", pDEBUG)
    << "mode= " << fINukeMode << "; Mean free path = " << L << " fm / "
                              << "Generated path length = " << d << " fm";
  */
  return d;
}
//___________________________________________________________________________
void Intranuke2018::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void Intranuke2018::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//___________________________________________________________________________
