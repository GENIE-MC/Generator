//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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
*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/Intranuke.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeUtils.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
Intranuke::Intranuke() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
Intranuke::Intranuke(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
Intranuke::Intranuke(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
Intranuke::~Intranuke()
{

}
//___________________________________________________________________________
void Intranuke::GenerateVertex(GHepRecord * evrec) const
{
// Generate a vtx and set it to all GHEP physical particles
// This method only takes effect in intranuke's h+A test mode.

  if(!fInTestMode) return;

  GHepParticle * nucltgt = evrec->TargetNucleus();
  assert(nucltgt);

  RandomGen * rnd = RandomGen::Instance();
  TVector3 vtx(999999.,999999.,999999.);

  // *** For h+A events (test mode): 
  // Assume a hadron beam with uniform intensity across an area, 
  // so we need to choose events uniformly within that area.
  double x=999999., y=999999., epsilon = 0.001;
  double R2  = TMath::Power(fNuclRadius,2.);
  double rp2 = TMath::Power(x,2.) + TMath::Power(y,2.);
  while(rp2 > R2-epsilon) {
      x = (fNuclRadius-epsilon) * rnd->RndFsi().Rndm();
      y = -fNuclRadius + 2*fNuclRadius * rnd->RndFsi().Rndm();
      y -= ((y>0) ? epsilon : -epsilon);
      rp2 = TMath::Power(x,2.) + TMath::Power(y,2.);
  }
  vtx.SetXYZ(x,y, -1.*TMath::Sqrt(TMath::Max(0.,R2-rp2)) + epsilon);

  // get the actual unit vector along the incoming hadron direction
  TVector3 direction = evrec->Probe()->P4()->Vect().Unit();

  // rotate the vtx position
  vtx.RotateUz(direction);
  
  LOG("Intranuke", pNOTICE) 
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
void Intranuke::SetNuclearRadius(const GHepParticle * p) const
{
  assert(p && pdg::IsIon(p->Pdg()));
  double A = p->A();
  fNuclRadius = fR0 * TMath::Power(A, 1./3.);

  // multiply that by some input factor so that hadrons are tracked
  // beyond the nuclear 'boundary' since the nuclear density distribution
  // is not zero there
  fNuclRadius *= fNR; 

  LOG("Intranuke", pNOTICE) 
               << "Setting nuclear radius to R = " << fNuclRadius;
}
//___________________________________________________________________________
bool Intranuke::NeedsRescattering(const GHepParticle * p) const
{
// checks whether the particle should be rescattered

  assert(p);

  if(fInTestMode) {
    // in test mode (initial state: hadron+nucleus) rescatter the initial
    // state hadron
    return ((p->Status() == kIStInitialState || p->Status() == kIStHadronInTheNucleus)
	    && !pdg::IsIon(p->Pdg()));
  }
  else {
   // attempt to rescatter anything marked as 'hadron in the nucleus'
   return (p->Status() == kIStHadronInTheNucleus);
  }
}
//___________________________________________________________________________
bool Intranuke::CanRescatter(const GHepParticle * p) const
{
// checks whether a particle that needs to be rescattered, can in fact be 
// rescattered by this cascade MC

  assert(p);
  return  ( p->Pdg() == kPdgPiP     || 
            p->Pdg() == kPdgPiM     || 
            p->Pdg() == kPdgPi0     ||
            p->Pdg() == kPdgProton  ||
            p->Pdg() == kPdgNeutron ||
	    p->Pdg() == kPdgGamma   ||
	    p->Pdg() == kPdgKP
          );
}
//___________________________________________________________________________
bool Intranuke::IsInNucleus(const GHepParticle * p) const
{
// check whether the input particle is still within the nucleus
//
  return (p->X4()->Vect().Mag() < fNuclRadius + fHadStep);
}
//___________________________________________________________________________
void Intranuke::TransportHadrons(GHepRecord * evrec) const
{
// transport all hadrons outside the nucleus

  //  Keep track of the remnant nucleus A,Z  
  //  Get 'nuclear environment' at the begginning of hadron transport 
  int iremn=-1;
  fRemnA=-1;
  fRemnZ=-1;
  if(fInTestMode) {
     // in hadron+nucleus mode get the particle at the 2nd GHEP row 
     GHepParticle * nucltgt = evrec->Particle(1);
     assert(nucltgt);
     iremn = 1;
     fRemnA = nucltgt->A();
     fRemnZ = nucltgt->Z();
  } else {
     // in neutrino+nucleus mode get the daughter of the initial
     // state nucleus that is a nucleus itself
     GHepParticle * nucltgt = evrec->Particle(1);
     assert(nucltgt);
     for(int i=nucltgt->FirstDaughter(); i<=nucltgt->LastDaughter(); i++) {
       if(i<0) continue;
       if(pdg::IsIon(evrec->Particle(i)->Pdg())) {
          fRemnA = evrec->Particle(i)->A();
          fRemnZ = evrec->Particle(i)->Z();
          iremn  = i;
       }
     }
  }
  const TLorentzVector & premn4 = *(evrec->Particle(iremn)->P4());
  fRemnP4 = premn4; 

  LOG("Intranuke", pNOTICE)
        << "Remnant nucleus (A,Z) = (" << fRemnA << ", " << fRemnZ << ")";

  // Loop over GHEP and run intranuclear rescattering on handled particles
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr = -1;

  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    // Check whether the particle needs rescattering, otherwise skip it
    if( ! this->NeedsRescattering(p) ) continue;

    if(this->HandleCompoundNucleus(evrec,p,icurr)) continue;

    LOG("Intranuke", pNOTICE)
      << " >> Stepping a " << p->Name() 
                        << " with kinetic E = " << p->KinE() << " GeV";

    // Rescatter a clone, not the original particle
    GHepParticle * sp = new GHepParticle(*p); 

    // Set clone's mom to be the hadron that was cloned
    sp->SetFirstMother(icurr); 

    // Check whether the particle can be rescattered 
    if(!this->CanRescatter(sp)) {

       // if I can't rescatter it, I will just take it out of the nucleus
       LOG("Intranuke", pNOTICE)
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
      utils::intranuke::StepParticle(sp, fHadStep);

      // check whether it interacts
      double d = this->GenerateStep(evrec,sp);
      has_interacted = (d<fHadStep);
      if(has_interacted) break;
    }//stepping
 
    if(has_interacted && fRemnA>0)  {
        // the particle interacts - simulate the hadronic interaction
      LOG("Intranuke", pNOTICE) 
          << "Particle has interacted at location:  " 
          << sp->X4()->Vect().Mag() << " / nucl rad= " << fNuclRadius;
	this->SimulateHadronicFinalState(evrec,sp);
    } else if(has_interacted && fRemnA<=0) {
        // nothing left to interact with!
      LOG("Intranuke", pNOTICE)
          << "*** Nothing left to interact with, escaping.";
	sp->SetStatus(kIStStableFinalState);
	evrec->AddParticle(*sp);
    } else {
        // the exits the nucleus without interacting - Done with it! 
        LOG("Intranuke", pNOTICE) 
          << "*** Hadron escaped the nucleus! Done with it.";
      sp->SetStatus(kIStStableFinalState);
      evrec->AddParticle(*sp);
    }
    delete sp;

    // Current snapshot
    //LOG("Intranuke", pINFO) << "Current event record snapshot: " << *evrec;

  }// GHEP entries

  // Add remnant nucleus - that 'hadronic blob' has all the remaining hadronic 
  // 4p not  put explicitly into the simulated particles
  TLorentzVector v4(0.,0.,0.,0.);
  GHepParticle remnant_nucleus(
    kPdgHadronicBlob, kIStFinalStateNuclearRemnant, 
    iremn,-1,-1,-1, fRemnP4, v4);
//     kPdgHadronicBlob, kIStStableFinalState, iremn,-1,-1,-1, fRemnP4, v4);
  evrec->AddParticle(remnant_nucleus);
  // Mark the initial remnant nucleus as an intermediate state 
  // Don't do that in the test mode sinc ethe initial remnant nucleus and
  // the initial nucleus coincide.
  if(!fInTestMode) {
     evrec->Particle(iremn)->SetStatus(kIStIntermediateState);
  }
}
//___________________________________________________________________________
double Intranuke::GenerateStep(GHepRecord* /*evrec*/, GHepParticle* p) const
{
// Generate a step (in fermis) for particle p in the input event.
// Computes the mean free path L and generate an 'interaction' distance d 
// from an exp(-d/L) distribution

  RandomGen * rnd = RandomGen::Instance();

  double L = utils::intranuke::MeanFreePath(p->Pdg(), *p->X4(), *p->P4(), fRemnA,
     fRemnZ, fDelRPion, fDelRNucleon);
  double d = -1.*L * TMath::Log(rnd->RndFsi().Rndm());

  LOG("Intranuke", pINFO)
            << "Mean free path = " << L << " fm / "
                              << "Generated path length = " << d << " fm";
  return d;
}
//___________________________________________________________________________
//___________________________________________________________________________
// Methods for configuring the INTRANUKE module:
//___________________________________________________________________________
//___________________________________________________________________________
void Intranuke::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void Intranuke::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
  }
//___________________________________________________________________________
