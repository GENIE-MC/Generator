//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 25, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 08, 2008 - CA
   Modified the global probability scale to be the maximum amongst the maximum
   interaction probabilities for each neutrino (rather than the sum of maximum
   probabilities). The modified probability scale still gives unbiased event
   generation & reduces the 'no-interaction' probability.
 @ Feb 14, 2008 - CA
   Significant speed improvements - Most of the rejected flux neutrinos, are
   rejected before having them propagated through the detector geometry 
   (flux neutrinos are pre-selected using the maximum path-lengths / for most
   realistic fluxes -high energy tail, low energy peak- most of the selection
   inefficiency is caused not because the path-lengths are not close the max
   possible ones, but because the energy is not close to the max possible one).
   Also some speed improvement was gained by properly using the total cross
   section splines (before, han repeatedly summing-up the
   The driver does not assert that _each_ flux neutrino generation & geometry
   navigation will be succesfull - In the rare event that this may happen, it
   prints an err mesg and tries again. In next revision I will limit the 
   number of successive trials something may go wrong to prevent the driver
   from hanging in truly problematic cases.
   The driver was adapted to handle flux drivers that -at some point- may stop
   generating more flux neutrinos (eg because they read flux neutrinos by 
   looping over a beam simulation ntuple and they reached its last entry).
   Code was appropriately restructured and some methods have been renamed.
 @ Feb 29, 2008 - CA
   Modified the InteractionProbability() to calculate absolute interaction
   probabilities. Added NFluxNeutrinos() and GlobProbScale() to get the
   number of neutrinos thrown by the flux driver towards the geometry and
   the global interaction probability scale so as to be able to calculate
   event sample normalization factors.
 @ Jan 15, 2009 - CA
   Stopped GMCJDriver from initializing the unphysical event mask so as not
   to overwrite the values that each GEVGDriver obtains from the environment.
 @ Jan 16, 2009 - CA
   Added methods to return pointers to the flux and geometry drivers.
*/
//____________________________________________________________________________

#include <cassert>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GEVGPool.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GeomAnalyzerI.h"
#include "GHEP/GHepFlags.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/InitialState.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/XSecSplineList.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GMCJDriver::GMCJDriver()
{
  this->InitJob();
}
//___________________________________________________________________________
GMCJDriver::~GMCJDriver()
{
  if (fGPool) delete fGPool;

  map<int,TH1D*>::iterator pmax_iter = fPmax.begin();
  for( ; pmax_iter != fPmax.end(); ++pmax_iter) {
    TH1D * pmax = pmax_iter->second;
    if(pmax) {
      delete pmax; pmax = 0;
    }
  }
  fPmax.clear();
}
//___________________________________________________________________________
void GMCJDriver::UseFluxDriver(GFluxI * flux_driver)
{
  fFluxDriver = flux_driver;
}
//___________________________________________________________________________
void GMCJDriver::UseGeomAnalyzer(GeomAnalyzerI * geom_analyzer)
{
  fGeomAnalyzer = geom_analyzer;
}
//___________________________________________________________________________
void GMCJDriver::UseSplines(bool useLogE)
{
  fUseSplines = true;
  fUseLogE    = useLogE;
}
//___________________________________________________________________________
void GMCJDriver::UseMaxPathLengths(string xml_filename)
{
// If you supply the maximum path lengths for all materials in your geometry
// (eg for ROOT/GEANT geometries they can be computed running GENIE's gmxpl 
// application, see $GENIE/src/stdapp/gMaxPathLengths.cxx ) you can speed up 
// the driver init phase by quite a bit (especially for complex geometries).

  fMaxPlXmlFilename = xml_filename;

  bool is_accessible = !(gSystem->AccessPathName(fMaxPlXmlFilename.c_str()));

  if (!is_accessible) fUseExtMaxPl = false;
  else                fUseExtMaxPl = true;
}
//___________________________________________________________________________
void GMCJDriver::KeepOnThrowingFluxNeutrinos(bool keep_on)
{
  LOG("GMCJDriver", pNOTICE)
        << "Keep on throwing flux neutrinos till one interacts? : "
                             << utils::print::BoolAsYNString(keep_on);
  fKeepThrowingFluxNu = keep_on;
}
//___________________________________________________________________________
void GMCJDriver::FilterUnphysical(const TBits & unphysmask)
{
  fUnphysMask = unphysmask;

  //if Configure() was run first configure all GEVGDrivers now
  if(fGPool) {
    GEVGPool::const_iterator giter;
    for(giter = fGPool->begin(); giter != fGPool->end(); ++giter) {
      GEVGDriver * driver = giter->second;
      driver->FilterUnphysical(fUnphysMask);
    }
  }
}
//___________________________________________________________________________
void GMCJDriver::ForceSingleProbScale()
{
// Use a single probability scale. That generates unweighted events.
// (Note that generating unweighted event kinematics is a different thing)
//
  fGenerateUnweighted = true;

  LOG("GMCJDriver", pNOTICE)
    << "GMCJDriver will generate un-weighted events. "
    << "Note: That does not force unweighted event kinematics!";
}
//___________________________________________________________________________
void GMCJDriver::Configure(void)
{
  LOG("GMCJDriver", pNOTICE)
     << utils::print::PrintFramedMesg("Configuring GMCJDriver");

  // Get the list of neutrino types from the input flux driver and the list
  // of target materials from the input geometry driver
  this->GetParticleLists();

  // Ask the input geometry driver to compute the max. path length for each
  // material in the list of target materials (or load a precomputed list)
  this->GetMaxPathLengthList();

  // Ask the input GFluxI for the max. neutrino energy (to compute Pmax)
  this->GetMaxFluxEnergy();

  // Create all possible initial states and for each one initialize, 
  // configure & store an GEVGDriver event generation driver object.
  // Once an 'initial state' has been selected from the input flux / geom,
  // the responsibility for generating the neutrino interaction will be
  // delegated to one of these drivers.
  this->PopulateEventGenDriverPool();

  // If the user wants to use cross section splines in order to speed things
  // up, then coordinate spline creation from all GEVGDriver objects pushed 
  // into GEVGPool. This will create all xsec splines needed for all (enabled)
  // processes that can be simulated involving the particles in the input flux 
  // and geometry. 
  // Spline creation will be skipped for every spline that has been pre-loaded 
  // into the the XSecSplineList.
  // Once more it is noted that computing cross section splines is a huge 
  // overhead. The user is encouraged to generate them in advance and load
  // them into the XSecSplineList
  this->BootstrapXSecSplines();

  // Create cross section splines describing the total interaction xsec
  // for a given initial state (Create them by summing all xsec splines
  // for each possible initial state)
  this->BootstrapXSecSplineSummation();

  // Compute the max. interaction probability to scale all interaction
  // probabilities to be computed by this driver
  this->ComputeProbScales();

  LOG("GMCJDriver", pNOTICE) << "Finished configuring GMCJDriver\n\n";
}
//___________________________________________________________________________
void GMCJDriver::InitJob(void)
{
  fFluxDriver         = 0;     // <-- flux driver
  fGeomAnalyzer       = 0;     // <-- geometry driver
  fGPool              = 0;     // <-- pool of GEVGDriver event generation drivers
  fEmax               = 0;     // <-- maximum neutrino energy
  fMaxPlXmlFilename   = "";    // <-- XML file with external path lengths
  fUseExtMaxPl        = false;
  fUseSplines         = false;
  fNFluxNeutrinos     = 0;     // <-- number of flux neutrinos thrown so far

  fGlobPmax           = 0;     // <-- maximum interaction probability (global prob scale)
  fPmax.clear();               // <-- maximum interaction probability per neutrino & per energy bin

  fGenerateUnweighted = false; // <-- default opt to generate weighted events

  fSelTgtPdg          = 0;
  fCurEvt             = 0;
  fCurVtx.SetXYZT(0.,0.,0.,0.);

  // Throw as many flux neutrinos as necessary till one has interacted
  // so that GenerateEvent() never  returns NULL (except when in error)
  this->KeepOnThrowingFluxNeutrinos(true);

  // Allow the selected GEVGDriver to go into recursive mode and regenerate
  // an interaction that turns out to be unphysical.
  //TBits unphysmask(GHepFlags::NFlags());
  //unphysmask.ResetAllBits(false); 
  //this->FilterUnphysical(unphysmask);

  // Force early initialization of singleton objects that are typically
  // would be initialized at their first use later on.
  // This is purely cosmetic and I do it to send the banner and some prolific
  // initialization printout at the top.
  assert( Messenger::Instance()     );
  assert( AlgConfigPool::Instance() );

  // Autoload splines (from the XML file pointed at the $GSPLOAD env. var.,
  // if the env. var. has been set);
  XSecSplineList * xspl = XSecSplineList::Instance();
  xspl->AutoLoad();

  // Clear the target and neutrino lists
  fNuList.clear();
  fTgtList.clear();

  // Clear the maximum path length list
  fMaxPathLengths.clear();
  fCurPathLengths.clear();
}
//___________________________________________________________________________
void GMCJDriver::GetParticleLists(void)
{
  // Get the list of flux neutrinos from the flux driver
  LOG("GMCJDriver", pNOTICE)
                    << "Asking the flux driver for its list of neutrinos";
  fNuList = fFluxDriver->FluxParticles();

  LOG("GMCJDriver", pNOTICE) << "Flux particles: " << fNuList;

  // Get the list of target materials from the geometry driver
  LOG("GMCJDriver", pNOTICE)
                  << "Asking the geometry driver for its list of targets";
  fTgtList = fGeomAnalyzer->ListOfTargetNuclei();

  LOG("GMCJDriver", pNOTICE) << "Target materials: " << fTgtList;
}
//___________________________________________________________________________
void GMCJDriver::GetMaxPathLengthList(void)
{
  if(fUseExtMaxPl) {
     LOG("GMCJDriver", pNOTICE)
       << "Loading external max path-length list for input geometry";
     fMaxPathLengths.LoadFromXml(fMaxPlXmlFilename);

  } else {
     LOG("GMCJDriver", pNOTICE)
       << "Querying the geometry driver to compute the max path-length list";
     fMaxPathLengths = fGeomAnalyzer->ComputeMaxPathLengths();
  }
  // Print maximum path lengths & neutrino energy
  LOG("GMCJDriver", pNOTICE)
     << "Maximum path length list: " << fMaxPathLengths;
}
//___________________________________________________________________________
void GMCJDriver::GetMaxFluxEnergy(void)
{
  LOG("GMCJDriver", pNOTICE)
     << "Querying the flux driver for the maximum energy of flux neutrinos";
  fEmax = fFluxDriver->MaxEnergy();

  LOG("GMCJDriver", pNOTICE) 
     << "Maximum flux neutrino energy = " << fEmax << " GeV";
}
//___________________________________________________________________________
void GMCJDriver::PopulateEventGenDriverPool(void)
{
  LOG("GMCJDriver", pDEBUG)
       << "Creating GEVGPool & adding a GEVGDriver object per init-state";

  if (fGPool) delete fGPool;
  fGPool = new GEVGPool;

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;

  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter) {
   for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {

     int target_pdgc   = *tgtiter;
     int neutrino_pdgc = *nuiter;

     InitialState init_state(target_pdgc, neutrino_pdgc);

     LOG("GMCJDriver", pNOTICE)
       << "\n\n ---- Creating a GEVGDriver object configured for init-state: "
       << init_state.AsString() << " ----\n\n";

     GEVGDriver * evgdriver = new GEVGDriver;
     evgdriver->Configure(init_state);
     //evgdriver->FilterUnphysical(fUnphysMask);
     evgdriver->UseSplines(); // check if all splines needed are loaded

     LOG("GMCJDriver", pDEBUG) << "Adding new GEVGDriver object to GEVGPool";
     fGPool->insert( GEVGPool::value_type(init_state.AsString(), evgdriver) );
   } // targets
  } // neutrinos

  LOG("GMCJDriver", pNOTICE)
             << "All necessary GEVGDriver object were pushed into GEVGPool\n";
}
//___________________________________________________________________________
void GMCJDriver::BootstrapXSecSplines(void)
{
// Bootstrap cross section spline generation by the event generation drivers
// that handle each initial state.

  if(!fUseSplines) return;

  LOG("GMCJDriver", pNOTICE) 
    << "Asking event generation drivers to compute all needed xsec splines";

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;
  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter){
     for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {
       int target_pdgc   = *tgtiter;
       int neutrino_pdgc = *nuiter;
       InitialState init_state(target_pdgc, neutrino_pdgc);
       LOG("GMCJDriver", pINFO)
           << "Computing all splines needed for init-state: "
           << init_state.AsString();
       GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
       evgdriver->CreateSplines(-1,-1,fUseLogE);
     } // targets
  } // neutrinos
  LOG("GMCJDriver", pINFO) << "Finished creating cross section splines\n";
}
//___________________________________________________________________________
void GMCJDriver::BootstrapXSecSplineSummation(void)
{
// Sum-up the cross section splines for all the interaction that can be
// simulated for each initial state

  LOG("GMCJDriver", pNOTICE)
    << "Summing-up splines to get total cross section for each init state";

  GEVGPool::iterator diter;
  for(diter = fGPool->begin(); diter != fGPool->end(); ++diter) {
    string       init_state = diter->first;
    GEVGDriver * evgdriver  = diter->second;
    assert(evgdriver);
    LOG("GMCJDriver", pNOTICE)
             << "**** Summing xsec splines for init-state = " << init_state;

    Range1D_t rE = evgdriver->ValidEnergyRange();
    assert(fEmax<rE.max && fEmax>rE.min);

    // decide the energy range for the sum spline - extend the spline a little
    // bit above the maximum beam energy (but below the maximum valid energy)
    // to avoid the evaluation of the cubic spline around the viscinity of
    // knots with zero y values (although the GENIE Spline object handles it)
    double dE  = fEmax/10.;
    double min = rE.min;
    double max = (fEmax+dE < rE.max) ? fEmax+dE : rE.max;
    evgdriver->CreateXSecSumSpline(100,min,max,true);
  }
  LOG("GMCJDriver", pNOTICE)
     << "Finished summing all interaction xsec splines per initial state";
}
//___________________________________________________________________________
void GMCJDriver::ComputeProbScales(void)
{
// Computing interaction probability scales.
// To minimize the numbers or trials before choosing a neutrino+target init
// state for generating an event (note: each trial means selecting a flux 
// neutrino, navigating it through the detector to compute path lengths, 
// computing  interaction probabilities for each material and so on...) 
// a set of probability scales can be used for different neutrino species 
// and at different energy bins. 
// A global probability scale is also being constructed for keeping the correct 
// proportions between differect flux neutrino species or flux neutrinos of 
// different energies.

  LOG("GMCJDriver", pNOTICE)
    << "Computing the max. interaction probability (probability scale)";

  // clean up global probability scale and maximum probabilties per neutrino
  // type & energy bin
  fGlobPmax = 0;
  map<int,TH1D*>::iterator pmax_iter = fPmax.begin();
  for( ; pmax_iter != fPmax.end(); ++pmax_iter) {
    TH1D * pmax = pmax_iter->second;
    if(pmax) {
      delete pmax; pmax = 0;    
    }
  }
  fPmax.clear();

  // for maximum interaction probability vs E /for given geometry/ I will
  // be using ~200 MeV bins
  //
  double de   = 0.2;
  double emin = 0.0;
  double emax = fEmax + de;
  int n = 1 + (int) ((emax-emin)/de);

  PDGCodeList::const_iterator nuiter;
  PDGCodeList::const_iterator tgtiter;

  // loop over all neutrino types generated by the flux driver
  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter) {
    int neutrino_pdgc = *nuiter;
    TH1D * pmax_hst = new TH1D("pmax_hst",
             "max interaction probability vs E | geom",n,emin,emax);
    pmax_hst->SetDirectory(0);

    // loop over energy bins
    for(int ie = 1; ie <= pmax_hst->GetNbinsX(); ie++) {
       double Ev = pmax_hst->GetBinCenter(ie);

       // loop over targets in input geometry, form initial state and compute
       // the sum of maximum interaction probabilities at the current energy bin
       //
       for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {
         int target_pdgc = *tgtiter;

         InitialState init_state(target_pdgc, neutrino_pdgc);

         LOG("GMCJDriver", pDEBUG)
           << "Computing Pmax for init-state: " 
           << init_state.AsString() << " at E = " << Ev;

         // get the appropriate driver
         GEVGDriver * evgdriver = fGPool->FindDriver(init_state);

         // get xsec sum over all modelled processes for given neutrino+target)
         double sxsec = evgdriver->XSecSumSpline()->Evaluate(Ev); 
         // get max{path-length x density}
         double plmax = fMaxPathLengths.PathLength(target_pdgc);  
         // compute/store the max interaction probabiity (for given energy)
         int A = pdg::IonPdgCodeToA(target_pdgc);
         double pmax = this->InteractionProbability(sxsec, plmax, A);
         pmax_hst->SetBinContent(ie, pmax_hst->GetBinContent(ie) + pmax);

         LOG("GMCJDriver", pDEBUG)
           << "Pmax[" << init_state.AsString() << ", Ev=" << Ev << "] = " << pmax;
       } // targets

       pmax_hst->SetBinContent(ie, 1.2 * pmax_hst->GetBinContent(ie));

       LOG("GMCJDriver", pINFO)
          << "Pmax[nu=" << neutrino_pdgc << ", Ev=" << Ev << "] = " 
          <<  pmax_hst->GetBinContent(ie);
    } // E

    fPmax.insert(map<int,TH1D*>::value_type(neutrino_pdgc,pmax_hst));
  } // nu

  // Compute global probability scale
  // Sum Probabilities {
  //   all neutrinos, all targets, @  max path length, @ max energy}
  //
  for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter) {
    int neutrino_pdgc = *nuiter;
    map<int,TH1D*>::const_iterator pmax_iter = fPmax.find(neutrino_pdgc);
    assert(pmax_iter != fPmax.end());
    TH1D * pmax_hst = pmax_iter->second;
    assert(pmax_hst);
//  double pmax = pmax_hst->GetBinContent(pmax_hst->FindBin(fEmax));
    double pmax = pmax_hst->GetMaximum();
    assert(pmax>0);        
//  fGlobPmax += pmax;
    fGlobPmax = TMath::Max(pmax, fGlobPmax); // ?;
  }

  LOG("GMCJDriver", pNOTICE) << "*** Probability scale = " << fGlobPmax;
}
//___________________________________________________________________________
void GMCJDriver::InitEventGeneration(void)
{
  fCurPathLengths.clear();
  fSelTgtPdg = 0;
  fCurVtx.SetXYZT(0.,0.,0.,0.);
}
//___________________________________________________________________________
EventRecord * GMCJDriver::GenerateEvent(void)
{
  LOG("GMCJDriver", pNOTICE) << "Generating next event...";

  this->InitEventGeneration();

  while(1) {
    bool flux_end = fFluxDriver->End();
    if(flux_end) {
       LOG("GMCJDriver", pNOTICE) 
           << "No more neutrinos can be thrown by the flux driver";
       return 0;
    }

    EventRecord * event = this->GenerateEvent1Try();
    if(event) return event;

    if(fKeepThrowingFluxNu) {
         LOG("GMCJDriver", pNOTICE)
             << "Flux neutrino didn't interact - Trying the next one...";
         continue;
    }
    break;
  } // (w(1)

  LOG("GMCJDriver", pINFO) << "Returning NULL event!";
  return 0;
}
//___________________________________________________________________________
EventRecord * GMCJDriver::GenerateEvent1Try(void)
{
// attempt generating a neutrino interaction by firing a single flux neutrino
//
  RandomGen * rnd = RandomGen::Instance();
  bool preselect = true;

  double Pno=0, Psum=0;
  double R = rnd->RndEvg().Rndm();
  LOG("GMCJDriver", pDEBUG) << "Rndm [0,1] = " << R;

  // Generate a neutrino using the input GFluxI & get current pdgc/p4/x4
  bool flux_ok = this->GenerateFluxNeutrino();
  if(!flux_ok) {
     LOG("GMCJDriver", pERROR) 
        << "** Rejecting current flux neutrino (flux driver err)";
     return 0;
  }

  // Compute the interaction probabilities assuming max. path lengths
  // and decide whether the neutrino would interact -- 
  // Many flux neutrinos should be rejected here, drastically reducing 
  // the number of neutrinos that I need to propagate through the 
  // actual detector geometry
  if(preselect) {
       LOG("GMCJDriver", pNOTICE) 
          << "Computing interaction probabilities for max. path lengths";

       Psum = this->ComputeInteractionProbabilities(true /* <- max PL*/);
       Pno  = 1-Psum;
       LOG("GMCJDriver", pNOTICE)
          << "The 'no interaction' probability (max. path lengths) is: " 
          << 100*Pno << " %";
        if(R>=1-Pno) {
  	   LOG("GMCJDriver", pNOTICE)  
              << "** Rejecting current flux neutrino";
	   return 0;
        }
  } // preselect 

  // Compute (pathLength x density x weight fraction) for all materials
  // in the input geometry, for the neutrino generated by the flux driver
  bool pl_ok = this->ComputePathLengths();
  if(!pl_ok) {
     LOG("GMCJDriver", pERROR) 
        << "** Rejecting current flux neutrino (err computing path-lengths)";
     return 0;
  }
  if(fCurPathLengths.AreAllZero()) {
     LOG("GMCJDriver", pNOTICE) 
        << "** Rejecting current flux neutrino (misses generation volume)";
     return 0;
  }

  Psum = this->ComputeInteractionProbabilities(false /* <- actual PL */);
  Pno  = 1-Psum;
  LOG("GMCJDriver", pNOTICE)
     << "The actual 'no interaction' probability is: " << 100*Pno << " %";
  if(R>=1-Pno) {
     LOG("GMCJDriver", pNOTICE) 
        << "** Rejecting current flux neutrino";
     return 0;
  }

  //
  // The flux neutrino interacts! Select a target material
  //
  fSelTgtPdg = this->SelectTargetMaterial(R);
  if(fSelTgtPdg==0) {
     LOG("GMCJDriver", pERROR) 
        << "** Rejecting current flux neutrino (failed to select tgt!)";
     return 0;
  }

  // Ask the GEVGDriver object to select and generate an interaction and
  // its kinematics for the selected initial state & neutrino 4-momentum
  this->GenerateEventKinematics();

  // Generate an 'interaction position' in the selected material (in the
  // detector coord system), along the direction of nup4 & set it 
  this->GenerateVertexPosition();

  // Set the event probability (probability for this event to happen given
  // the detector setup & the selected flux neutrino)
  // Note for users: 
  // The above probability is stored at GHepRecord::Probability()
  // For normalization purposes make sure that you take into account the
  // GHepRecord::Weight() -if event generation is weighted-, and
  // GFluxI::Weight() -if beam simulation is weighted-.
  this->ComputeEventProbability();

  return fCurEvt;
}
//___________________________________________________________________________
bool GMCJDriver::GenerateFluxNeutrino(void)
{
// Ask the neutrino flux driver to generate a flux neutrino and make sure
// that things look ok...
//
  LOG("GMCJDriver", pNOTICE) << "Generating a flux neutrino";

  bool ok = fFluxDriver->GenerateNext();
  if(!ok) {
     LOG("GMCJDriver", pERROR)
         << "*** The flux driver couldn't generate a flux neutrino!!";
     return false;
  }

  fNFluxNeutrinos++;
  int                    nupdg = fFluxDriver -> PdgCode  ();
  const TLorentzVector & nup4  = fFluxDriver -> Momentum ();
  const TLorentzVector & nux4  = fFluxDriver -> Position ();

  LOG("GMCJDriver", pNOTICE)
     << "\n [-] Generated flux neutrino: "
     << "\n  |----o PDG-code   : " << nupdg
     << "\n  |----o 4-momentum : " << utils::print::P4AsString(&nup4)
     << "\n  |----o 4-position : " << utils::print::X4AsString(&nux4);

  if(nup4.Energy() > fEmax) {
     LOG("GMCJDriver", pFATAL)
       << "\n *** Flux driver error ***"
       << "\n Generated flux v with E = " << nup4.Energy() << " GeV"
       << "\n Max v energy (declared by flux driver) = " << fEmax << " GeV"
       << "\n My interaction probability scaling is invalidated!!";
     return false;
  }
  if(!fNuList.ExistsInPDGCodeList(nupdg)) {
     LOG("GMCJDriver", pFATAL)
       << "\n *** Flux driver error ***"
       << "\n Generated flux v with pdg = " << nupdg
       << "\n It does not belong to the declared list of flux neutrinos"
       << "\n I was not configured to handle this!!";
     return false;
  }
  return true;
}
//___________________________________________________________________________
bool GMCJDriver::ComputePathLengths(void)
{
// Ask the geometry driver to compute (pathLength x density x weight frac.)
// for all detector materials for the neutrino generated by the flux driver
// and make sure that things look ok...

  fCurPathLengths.clear();

  const TLorentzVector & nup4  = fFluxDriver -> Momentum ();
  const TLorentzVector & nux4  = fFluxDriver -> Position ();

  fCurPathLengths = fGeomAnalyzer->ComputePathLengths(nux4, nup4);

  LOG("GMCJDriver", pNOTICE) << fCurPathLengths;

  if(fCurPathLengths.size() == 0) {
     LOG("GMCJDriver", pFATAL)
       << "\n *** Geometry driver error ***"
       << "\n Got an empty PathLengthList - No material found in geometry?";
     return false;
  }

  if(fCurPathLengths.AreAllZero()) {
         LOG("GMCJDriver", pNOTICE)
                 << "current flux v doesn't cross any geometry material...";
  }
  return true;
}
//___________________________________________________________________________
double GMCJDriver::ComputeInteractionProbabilities(bool use_max_path_length)
{
  LOG("GMCJDriver", pNOTICE)
       << "Computing relative interaction probabilities for each material";

  // current flux neutrino code & 4-p
  int                    nupdg = fFluxDriver->PdgCode();
  const TLorentzVector & nup4  = fFluxDriver->Momentum();

  fCurCumulProbMap.clear();

  const PathLengthList & path_length_list = 
        (use_max_path_length) ? fMaxPathLengths : fCurPathLengths;

  double probsum=0;
  PathLengthList::const_iterator pliter;

  for(pliter = path_length_list.begin();
                            pliter != path_length_list.end(); ++pliter) {
     int    mpdg  = pliter->first;            // material PDG code
     double pl    = pliter->second;           // density x path-length
     int    A     = pdg::IonPdgCodeToA(mpdg);
     double xsec  = 0.;                       // sum of xsecs for all modelled processes for given init state
     double prob  = 0.;                       // interaction probability
     double probn = 0.;                       // normalized interaction probability

     // find the GEVGDriver object that is handling the current init state
     InitialState init_state(mpdg, nupdg);
     GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
     if(!evgdriver) {
       LOG("GMCJDriver", pFATAL)
        << "\n * The MC Job driver isn't properly configured!"
        << "\n * No event generation driver could be found for init state: " 
        << init_state.AsString();
       exit(1);
     }
     // compute the interaction xsec and probability (if path-length>0)
     if(pl>0.) {
        const Spline * totxsecspl = evgdriver->XSecSumSpline();
        if(!totxsecspl) {
            LOG("GMCJDriver", pFATAL)
              << "\n * The MC Job driver isn't properly configured!"
              << "\n * Couldn't retrieve total cross section spline for init state: " 
              << init_state.AsString();
            exit(1);
        } else {
            xsec = totxsecspl->Evaluate( nup4.Energy() );
        }
        prob = this->InteractionProbability(xsec,pl,A);

        // scale the interaction probability to the maximum one so as not
        // to have to throw few billions of flux neutrinos before getting
        // an interaction...
        double pmax = 0;
        if(fGenerateUnweighted) pmax = fGlobPmax;
        else {
           map<int,TH1D*>::const_iterator pmax_iter = fPmax.find(nupdg);
           assert(pmax_iter != fPmax.end());
           TH1D * pmax_hst = pmax_iter->second;
           assert(pmax_hst);
           int    ie   = pmax_hst->FindBin(nup4.Energy());
           pmax = pmax_hst->GetBinContent(ie);
        }
        assert(pmax>0);        
        probn = prob/pmax;
     }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("GMCJDriver", pNOTICE)
         << "tgt: " << mpdg << " -> TotXSec = "
         << xsec/units::cm2 << " cm^2, Norm.Prob = " << 100*probn << "%";
#endif
     probsum += probn;
     fCurCumulProbMap.insert(map<int,double>::value_type(mpdg,probsum));
  }
  return probsum;
}
//___________________________________________________________________________
int GMCJDriver::SelectTargetMaterial(double R)
{
// Pick a target material using the pre-computed interaction probabilities
// for a flux neutrino that has already been determined that interacts

  LOG("GMCJDriver", pNOTICE) << "Selecting target material";
  int tgtpdg = 0;
  map<int,double>::const_iterator probiter = fCurCumulProbMap.begin();
  for( ; probiter != fCurCumulProbMap.end(); ++probiter) {
     double prob = probiter->second;
     if(R<prob) {
        tgtpdg = probiter->first;
        LOG("GMCJDriver", pNOTICE) 
          << "Selected target material = " << tgtpdg;
        return tgtpdg;
     }
  }
  LOG("GMCJDriver", pERROR)
     << "Could not select target material for an interacting neutrino";
  return 0;
}
//___________________________________________________________________________
void GMCJDriver::GenerateEventKinematics(void)
{
  int                    nupdg = fFluxDriver->PdgCode();
  const TLorentzVector & nup4  = fFluxDriver->Momentum();

  // Find the GEVGDriver object that generates interactions for the
  // given initial state (neutrino + target)
  InitialState init_state(fSelTgtPdg, nupdg);
  GEVGDriver * evgdriver = fGPool->FindDriver(init_state);
  if(!evgdriver) {
     LOG("GMCJDriver", pFATAL)
       << "No GEVGDriver object for init state: " << init_state.AsString();
     exit(1);
  }

  // Ask the GEVGDriver object to select and generate an interaction for
  // the selected initial state & neutrino 4-momentum
  LOG("GMCJDriver", pNOTICE)
          << "Asking the selected GEVGDriver object to generate an event";
  fCurEvt = evgdriver->GenerateEvent(nup4);
}
//___________________________________________________________________________
void GMCJDriver::GenerateVertexPosition(void)
{
  // Generate an 'interaction position' in the selected material, along
  // the direction of nup4
  LOG("GMCJDriver", pNOTICE)
     << "Asking the geometry analyzer to generate a vertex";

  const TLorentzVector & p4 = fFluxDriver->Momentum ();
  const TLorentzVector & x4 = fFluxDriver->Position ();

  const TVector3 & vtx = fGeomAnalyzer->GenerateVertex(x4, p4, fSelTgtPdg);

  TVector3 origin(x4.X(), x4.Y(), x4.Z());
  origin-=vtx; // computes vector dr = origin - vtx

  double dL = origin.Mag();
  double c  = kLightSpeed /(units::meter/units::second);
  double dt = dL/c;

  LOG("GMCJDriver", pNOTICE)
     << "|vtx - origin|: dL = " << dL << " m, dt = " << dt << " sec";

  fCurVtx.SetXYZT(vtx.x(), vtx.y(), vtx.z(), x4.T() + dt);

  fCurEvt->SetVertex(fCurVtx);
}
//___________________________________________________________________________
void GMCJDriver::ComputeEventProbability(void)
{
// Compute event probability for the given flux neutrino & detector geometry

  // interaction cross section
  // (convert from physical units -> cm^2)
  double xsec = fCurEvt->XSec() / units::cm2; 

  // path length in detector along v direction for specified target material
  // (convert from kgr/m2 to gr/cm2)
  PathLengthList::const_iterator pliter = fCurPathLengths.find(fSelTgtPdg);
  double path_length = pliter->second;
  path_length *= ((units::kilogram/units::m2)/(units::gram/units::cm2));

  // target material mass number
  // (in gr)
  int A = pdg::IonPdgCodeToA(fSelTgtPdg);

  // gett info on the interacted neutrino
  GHepParticle * nu = fCurEvt->Probe();
  int    nu_pdg = nu->Pdg();
  double Ev     = nu->P4()->Energy();
 
  // interaction probability
  double P = this->InteractionProbability(xsec, path_length, A);

  // weight for selected event
  double weight = 1.0;
  if(!fGenerateUnweighted) {
     map<int,TH1D*>::const_iterator pmax_iter = fPmax.find(nu_pdg);
     assert(pmax_iter != fPmax.end());
     TH1D * pmax_hst = pmax_iter->second;
     assert(pmax_hst);
     double pmax = pmax_hst->GetBinContent(pmax_hst->FindBin(Ev));
     assert(pmax>0);
     weight = pmax/fGlobPmax;
  }

  // set probability & update weight
  fCurEvt->SetProbability(P);
  fCurEvt->SetWeight(weight * fCurEvt->Weight());
}
//___________________________________________________________________________
double GMCJDriver::InteractionProbability(double xsec, double pL, int A)
{
// P = Na   (Avogadro number,                 atoms/mole) *
//     1/A  (1/mass number,                   mole/gr)    *
//     xsec (total interaction cross section, cm^2)       *
//     pL   (density-weighted path-length,    gr/cm^2)
//
  xsec = xsec / units::cm2; 
  pL   = pL   * ((units::kilogram/units::m2)/(units::gram/units::cm2));

  return kNA*(xsec*pL)/A;
}
//___________________________________________________________________________
