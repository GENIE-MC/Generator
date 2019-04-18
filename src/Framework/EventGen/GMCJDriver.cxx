//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

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
 @ Mar 11, 2009 - CA
   In GenerateEvent1Try() handle failure to generate kinematics. Added sanity
   check on the no interaction probability.
 @ Mar 04, 2010 - CA
   Remove unused FilterUnphysical(TBits) method. Now set exclusively via the
   GUNPHYSMASK env.var.
 @ Apr 15, 2010 - CA
   Fix unit error in ComputeEventProbability() - Reported by Corey Reed.
   The probability stored at the output event was wrong but this doesn't
   affect any of the existing applications as this number wasn't actually
   used anywhere.
 @ Dec 07, 2010 - CA
   Don't use a fixed bin size in ComputeProbScales() as this was causing 
   errors for low energy applications. Addresses a problem reported by
   Joachim Kopp.
 @ Feb 22, 2011 - JD
   Added a number of new methods to allow pre-calculation of exact flux 
   interaction probabilities for a given set of flux neutrinos from the 
   flux driver. See the comments for the new LoadFluxProbabilities, 
   SaveFluxProbabilities, PreCalcFluxProbabilities and PreSelectEvents 
   methods for details. Using these methods mean that there is no need 
   to generate maximum path lengths as instead use the exact interaction 
   probabilities to pre-select. This can result in very significant speed 
   increases (between factor of 5 and ~300) for event generation over complex
   detector geometries and with realistic flux drivers. See 
   src/support/t2k/EvGen/gT2KEvGen.cxx for an example of how to use.
 @ Mar, 7, 2011 - JD
   Store sum totals of the flux interaction probabilities for various neutrino 
   type in a map relating pdg code to total interaction probability. Also add
   public getter method so that this can be used in applications to work out
   expected event rates. See gT2KEvGen.cxx for an example of how to do this. 
   Also save the PDG code for each entry in the flux interaction probabilities 
   tree. 
 @ Mar, 11, 2011 - JD 
   Set the directory of fFluxIntTree to the output file fFluxIntProbFile if
   saving it later. This is so that it is incrementally saved and fixes bug
   where getting std::bad_alloc when trying to Write large trees 
   fFluxIntProbFile.   
 @ Jan 31, 2013 - CA
   Added SetEventGeneratorList(string listname). $GEVGL var no longer in use.
 @ Feb 01, 2013 - CA
   The GUNPHYSMASK env. var is no longer used. Added SetUnphysEventMask(const 
   TBits &). Input is propagated accordingly.
 @ Feb 06, 2013 - CA
   Fix small problem introduced with recent changes.
   In PopulateEventGenDriverPool() calls to GEVGDriver::SetEventGeneratorList()
   and GEVGDriver::Configure() were reversed. Problem reported by W.Huelsnitz.
 @ July 15, 2014 - HG
   Incorporated code provided by Jason Koskinen - IceCube
   Modified ComputeProbScales to evalulate the cross sections at both the high
   and low edges of the energy bin when calculating the max interaction 
   probability.
*/
//____________________________________________________________________________

#include <cassert>

#include <TVector3.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GMCJDriver.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/GEVGPool.h"
#include "Framework/EventGen/GFluxI.h"
#include "Framework/EventGen/GeomAnalyzerI.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Conventions/Constants.h"

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
  if(fUnphysEventMask) delete fUnphysEventMask;
  if (fGPool) delete fGPool;

  map<int,TH1D*>::iterator pmax_iter = fPmax.begin();
  for( ; pmax_iter != fPmax.end(); ++pmax_iter) {
    TH1D * pmax = pmax_iter->second;
    if(pmax) {
      delete pmax; pmax = 0;
    }
  }
  fPmax.clear();

  if(fFluxIntTree) delete fFluxIntTree;
  if(fFluxIntProbFile) delete fFluxIntProbFile;
}
//___________________________________________________________________________
void GMCJDriver::SetEventGeneratorList(string listname)
{
  LOG("GMCJDriver", pNOTICE)
       << "Setting event generator list: " << listname;

  fEventGenList = listname;
}
//___________________________________________________________________________
void GMCJDriver::SetUnphysEventMask(const TBits & mask)
{
  *fUnphysEventMask = mask;

  LOG("GMCJDriver", pNOTICE)
    << "Setting unphysical event mask (bits: " << GHepFlags::NFlags() - 1
    << " -> 0) : " << *fUnphysEventMask;
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
bool GMCJDriver::UseMaxPathLengths(string xml_filename)
{
// If you supply the maximum path lengths for all materials in your geometry
// (eg for ROOT/GEANT geometries they can be computed running GENIE's gmxpl 
// application, see $GENIE/src/stdapp/gMaxPathLengths.cxx ) you can speed up 
// the driver init phase by quite a bit (especially for complex geometries).

  fMaxPlXmlFilename = xml_filename;

  bool is_accessible = !(gSystem->AccessPathName(fMaxPlXmlFilename.c_str()));

  if ( is_accessible ) fUseExtMaxPl = true;
  else {
    fUseExtMaxPl = false;
    LOG("GMCJDriver", pWARN)
      << "UseMaxPathLengths could not find file: \"" << xml_filename << "\"";
  }
  return fUseExtMaxPl;

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
void GMCJDriver::PreSelectEvents(bool preselect)
{
// Set whether to pre-select events based on a max-path lengths file. This
// should be turned off if using pre-generated interaction probabilities 
// calculated from a given flux file.
  fPreSelect = preselect;
}
//___________________________________________________________________________
bool GMCJDriver::PreCalcFluxProbabilities(void)
{
// Loop over complete set of flux entries satisfying input config options 
// (such as neutrino type) and save the interaction probability in a tree 
// relating flux index (entry number in input flux tree) to interaction 
// probability. If a pre-generated flux interaction probability tree has 
// already been loaded then just returns true. Also save tree to a TFile
// for use in later jobs if flag is set 
//
  bool success = true;
 
  bool save_to_file = fFluxIntProbFile == 0 && fFluxIntFileName.size()>0;

  // Clear map storing sum(fBrFluxWeight*fBrFluxIntProb) for each neutrino pdg
  fSumFluxIntProbs.clear();

  // check if already loaded flux interaction probs using LoadFluxProbTree
  if(fFluxIntTree){
    LOG("GMCJDriver", pNOTICE) << 
         "Skipping pre-generation of flux interaction probabilities - "<<
         "using pre-generated file";
    success = true;
  }
  // otherwise create them on the fly now 
  else {

    if(save_to_file){
      fFluxIntProbFile = new TFile(fFluxIntFileName.c_str(), "CREATE");
      if(fFluxIntProbFile->IsZombie()){
        LOG("GMCJDriver", pFATAL) << "Cannot overwrite an existing file. Exiting!";
        exit(1);
      } 
    } 
  
    // Create the tree to store flux probs
    fFluxIntTree = new TTree(fFluxIntTreeName.c_str(), 
                         "Tree storing pre-calculated flux interaction probs"); 
    fFluxIntTree->Branch("FluxIndex", &fBrFluxIndex, "FluxIndex/I");
    fFluxIntTree->Branch("FluxIntProb", &fBrFluxIntProb, "FluxIntProb/D");
    fFluxIntTree->Branch("FluxEnu", &fBrFluxEnu, "FluxEnu/D"); 
    fFluxIntTree->Branch("FluxWeight", &fBrFluxWeight, "FluxWeight/D"); 
    fFluxIntTree->Branch("FluxPDG", &fBrFluxPDG, "FluxPDG/I"); 
    // Associate to file otherwise get std::bad_alloc when writing large trees 
    if(save_to_file) fFluxIntTree->SetDirectory(fFluxIntProbFile); 
 
    fFluxDriver->GenerateWeighted(true);
  
    fGlobPmax = 1.0; // Force ComputeInteractionProbabilities to return absolute value
  
    // Loop over flux entries and calculate interaction probabilities
    TStopwatch stopwatch; 
    stopwatch.Start();
    long int first_index = -1;
    bool first_loop = true;
    // loop until at end of flux ntuple
    while(fFluxDriver->End() == false){ 

      // get the next flux neutrino
      bool gotnext = fFluxDriver->GenerateNext(); 
      if(!gotnext){
        LOG("GMCJDriver", pWARN) << "*** Couldn't generate next flux ray! ";
        continue;
      }

      // stop if completed a full cycle (this check is necessary as fluxdriver
      // may be set to loop over more than one cycle before reaching end) 
      bool already_been_here = first_loop ? false : first_index == fFluxDriver->Index();
      if(already_been_here) break; 
   
      // compute the path lengths for current flux neutrino 
      if(this->ComputePathLengths() == false){ success = false; break;}
  
      // compute and store the interaction probability 
      double psum = this->ComputeInteractionProbabilities(false /*Based on actual PLs*/);
      assert(psum+controls::kASmallNum > 0.);
      fBrFluxIntProb = psum;
      fBrFluxIndex   = fFluxDriver->Index();
      fBrFluxEnu     = fFluxDriver->Momentum().E();
      fBrFluxWeight  = fFluxDriver->Weight();
      fBrFluxPDG     = fFluxDriver->PdgCode();
      fFluxIntTree->Fill();

      // store the first index so know when have cycled exactly once
      if(first_loop){
        first_index = fFluxDriver->Index();
        first_loop = false;
      }
    } // flux loop
    stopwatch.Stop();            
    LOG("GMCJDriver", pNOTICE)
                    << "Finished pre-calculating flux interaction probabilities. "
                    << "Total CPU time to process "<< fFluxIntTree->GetEntries()
                    << " entries: "<< stopwatch.CpuTime();

    // reset the flux driver so can be used at next stage. N.B. This 
    // should also reset flux driver to throw de-weighted flux neutrinos
    fFluxDriver->Clear("CycleHistory");
  }

  // If successfully calculated/loaded interaction probabilities then set global
  // probability scale and, if requested, save tree to output file
  if(success){
    fGlobPmax = 0.0;
    double safety_factor = 1.01;
    for(int i = 0; i< fFluxIntTree->GetEntries(); i++){
      fFluxIntTree->GetEntry(i);
      // Check have non-negative probabilities
      assert(fBrFluxIntProb+controls::kASmallNum > 0.0);
      assert(fBrFluxWeight+controls::kASmallNum > 0.0);
      // Update the global maximum
      fGlobPmax = TMath::Max(fGlobPmax, fBrFluxIntProb*safety_factor); 
      // Update the sum of fBrFluxIntProb*fBrFluxWeight for different species
      if(fSumFluxIntProbs.find(fBrFluxPDG) == fSumFluxIntProbs.end()){
        fSumFluxIntProbs[fBrFluxPDG] = 0.0;
      }
      fSumFluxIntProbs[fBrFluxPDG] += fBrFluxIntProb * fBrFluxWeight;
    }
    LOG("GMCJDriver", pNOTICE) <<
        "Updated global probability scale to fGlobPmax = "<< fGlobPmax; 

    if(save_to_file){
      LOG("GMCJDriver", pNOTICE) <<
          "Saving pre-generated interaction probabilities to file: "<<
          fFluxIntProbFile->GetName();
      fFluxIntProbFile->cd();
      fFluxIntTree->Write();
    }

    // Also build index for use later
    if(fFluxIntTree->BuildIndex("FluxIndex") != fFluxIntTree->GetEntries()){
      LOG("GMCJDriver", pFATAL) << 
          "Cannot build index using branch \"FluxIndex\" for flux prob tree!"; 
      exit(1);
    } 
 
    // Now that have pre-generated flux probabilities need to trun off event 
    // preselection as this is only advantages when using max path lengths
    this->PreSelectEvents(false);

    LOG("GMCJDriver", pNOTICE) << "Successfully generated/loaded pre-calculate flux interaction probabilities";
  }
  // Otherwise clean up
  else if(fFluxIntTree){ 
    delete fFluxIntTree; 
    fFluxIntTree = 0;
  }
  
  // Return whether have successfully pre-calculated flux interaction probabilities
  return success;
}
//___________________________________________________________________________
bool GMCJDriver::LoadFluxProbabilities(string filename)
{
// Load a pre-generated set of flux interaction probabilities from an external
// file. This is recommended when using large flux files (>100k entries) as  
// for these the time to calculate the interaction probabilities can exceed 
// ~20 minutes. After loading the input tree we call PreCalcFluxProbabilities
// to check that has successfully loaded
//
  if(fFluxIntProbFile){
    LOG("GMCJDriver", pWARN) 
     << "Can't load flux interaction prob file as one is already loaded"; 
    return false;
  }

  fFluxIntProbFile = new TFile(filename.c_str(), "OPEN");

  if(fFluxIntProbFile){
    fFluxIntTree = dynamic_cast<TTree*>(fFluxIntProbFile->Get(fFluxIntTreeName.c_str())); 
    if(fFluxIntTree){
      bool set_addresses = 
        fFluxIntTree->SetBranchAddress("FluxIntProb", &fBrFluxIntProb) >= 0 &&
        fFluxIntTree->SetBranchAddress("FluxIndex", &fBrFluxIndex) >= 0 &&
        fFluxIntTree->SetBranchAddress("FluxPDG", &fBrFluxPDG) >= 0 &&
        fFluxIntTree->SetBranchAddress("FluxWeight", &fBrFluxWeight) >= 0 &&
        fFluxIntTree->SetBranchAddress("FluxEnu", &fBrFluxEnu) >= 0; 
      if(set_addresses){ 
        // Finally check that can use them
        if(this->PreCalcFluxProbabilities()) {
          LOG("GMCJDriver", pNOTICE) 
           << "Successfully loaded pre-generated flux interaction probabilities";
          return true;
        }
      }
      // If cannot load then delete tree 
      LOG("GMCJDriver", pERROR) << 
          "Cannot find expected branches in input flux probability tree!"; 
      delete fFluxIntTree; fFluxIntTree = 0; 
    }
    else LOG("GMCJDriver", pERROR) 
          << "Cannot find tree: "<< fFluxIntTreeName.c_str();
  }
     
  LOG("GMCJDriver", pWARN)
     << "Unable to load flux interaction probabilities file";
  return false;
}
//___________________________________________________________________________
void GMCJDriver::SaveFluxProbabilities(string outfilename)
{
// Configue the flux driver to save the calculated flux interaction
// probabilities to the specified output file name for use in later jobs. See
// the LoadFluxProbTree method for how they are fed into a later job. 
//
  fFluxIntFileName = outfilename;
}
//___________________________________________________________________________
void GMCJDriver::Configure(bool calc_prob_scales)
{
  LOG("GMCJDriver", pNOTICE)
     << utils::print::PrintFramedMesg("Configuring GMCJDriver");

  // Get the list of neutrino types from the input flux driver and the list
  // of target materials from the input geometry driver
  this->GetParticleLists();

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

  if(calc_prob_scales){
    // Ask the input geometry driver to compute the max. path length for each
    // material in the list of target materials (or load a precomputed list)
    this->GetMaxPathLengthList();

    // Compute the max. interaction probability to scale all interaction
    // probabilities to be computed by this driver
    this->ComputeProbScales();
  }
  LOG("GMCJDriver", pNOTICE) << "Finished configuring GMCJDriver\n\n";
}
//___________________________________________________________________________
void GMCJDriver::InitJob(void)
{
  fEventGenList       = "Default";  // <-- set of event generators to be loaded by this driver

  fUnphysEventMask = new TBits(GHepFlags::NFlags()); //<-- unphysical event mask
  //fUnphysEventMask->ResetAllBits(true);
  for(unsigned int i = 0; i < GHepFlags::NFlags(); i++) {
   fUnphysEventMask->SetBitNumber(i, true);
  }

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
  fPreSelect          = true;  // <-- default to use pre-selection based on maximum path lengths 

  fSelTgtPdg          = 0;
  fCurEvt             = 0;
  fCurVtx.SetXYZT(0.,0.,0.,0.);

  fFluxIntProbFile    = 0; 
  fFluxIntTreeName    = "gFlxIntProb";
  fFluxIntFileName    = "";
  fFluxIntTree        = 0;
  fBrFluxIntProb      = -1.;
  fBrFluxIndex        = -1;
  fBrFluxEnu          = -1.;       
  fBrFluxWeight       = -1.;
  fBrFluxPDG          = 0;
  fSumFluxIntProbs.clear();

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
       << "Loading external max path-length list for input geometry from "
       << fMaxPlXmlFilename;
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
     evgdriver->SetEventGeneratorList(fEventGenList); // specify list of generators
     evgdriver->Configure(init_state);
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
    if (fEmax>rE.max || fEmax<rE.min)
      LOG("GMCJDriver",pFATAL)
        << " rE (validEnergyRange) [" << rE.min << "," << rE.max << "] "
        << " fEmax " << fEmax;
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
  {
    fGlobPmax = 0;
    map<int,TH1D*>::iterator pmax_iter = fPmax.begin();
    for( ; pmax_iter != fPmax.end(); ++pmax_iter) {
      TH1D * pmax = pmax_iter->second;
      if(pmax) {
        delete pmax; pmax = 0;    
      }
    }
    fPmax.clear();
  }

  // for maximum interaction probability vs E /for given geometry/ I will
  // be using 300 bins up to the maximum energy for the input flux
  // double de   = fEmax/300.;//djk june 5, 2013
  double de   = fEmax/300.;//djk june 5, 2013
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
      double EvLow  = pmax_hst->GetBinCenter(ie) - 0.5*pmax_hst->GetBinWidth(ie); 
      double EvHigh = pmax_hst->GetBinCenter(ie) + 0.5*pmax_hst->GetBinWidth(ie); 
      //double Ev = pmax_hst->GetBinCenter(ie);

       // loop over targets in input geometry, form initial state and compute
       // the sum of maximum interaction probabilities at the current energy bin
       //
       for(tgtiter = fTgtList.begin(); tgtiter != fTgtList.end(); ++tgtiter) {
         int target_pdgc = *tgtiter;

         InitialState init_state(target_pdgc, neutrino_pdgc);

         LOG("GMCJDriver", pDEBUG)
           << "Computing Pmax for init-state: " << init_state.AsString() << " E from " << EvLow << "-" << EvHigh;

         // get the appropriate driver
         GEVGDriver * evgdriver = fGPool->FindDriver(init_state);

         // get xsec sum over all modelled processes for given neutrino+target)
         double sxsecLow  = evgdriver->XSecSumSpline()->Evaluate(EvLow);
	 double sxsecHigh = evgdriver->XSecSumSpline()->Evaluate(EvHigh);

         // get max{path-length x density}
         double plmax = fMaxPathLengths.PathLength(target_pdgc);

         // compute/store the max interaction probabiity (for given energy)
         int A = pdg::IonPdgCodeToA(target_pdgc);
         double pmaxLow  = this->InteractionProbability(sxsecLow, plmax, A);
         double pmaxHigh = this->InteractionProbability(sxsecHigh, plmax, A);

	 double pmax = pmaxHigh;
	 if ( pmaxLow > pmaxHigh){
	   pmax = pmaxLow;
	   LOG("GMCJDriver", pWARN)
	     << "Lower energy neutrinos have a higher probability of interacting than those at higher energy."
	     << " pmaxLow(E=" << EvLow << ")=" << pmaxLow << " and " << " pmaxHigh(E=" << EvHigh << ")=" << pmaxHigh;
	 }

         pmax_hst->SetBinContent(ie, pmax_hst->GetBinContent(ie) + pmax);

         LOG("GMCJDriver", pDEBUG)
           << "Pmax[" << init_state.AsString() << ", Ev from " << EvLow << "-" << EvHigh << "] = " << pmax;
       } // targets

       pmax_hst->SetBinContent(ie, 1.2 * pmax_hst->GetBinContent(ie));

       LOG("GMCJDriver", pINFO)
	 << "Pmax[nu=" << neutrino_pdgc << ", Ev from " << EvLow << "-" << EvHigh << "] = "
          <<  pmax_hst->GetBinContent(ie);
    } // E

    fPmax.insert(map<int,TH1D*>::value_type(neutrino_pdgc,pmax_hst));
  } // nu

  // Compute global probability scale
  // Sum Probabilities {
  //   all neutrinos, all targets, @  max path length, @ max energy}
  //
  {
    for(nuiter = fNuList.begin(); nuiter != fNuList.end(); ++nuiter) {
      int neutrino_pdgc = *nuiter;
      map<int,TH1D*>::const_iterator pmax_iter = fPmax.find(neutrino_pdgc);
      assert(pmax_iter != fPmax.end());
      TH1D * pmax_hst = pmax_iter->second;
      assert(pmax_hst);
//    double pmax = pmax_hst->GetBinContent(pmax_hst->FindBin(fEmax));
      double pmax = pmax_hst->GetMaximum();
      assert(pmax>0);        
//    fGlobPmax += pmax;
      fGlobPmax = TMath::Max(pmax, fGlobPmax); // ?;
    }
    LOG("GMCJDriver", pNOTICE) << "*** Probability scale = " << fGlobPmax;
  }
}
//___________________________________________________________________________
void GMCJDriver::InitEventGeneration(void)
{
  fCurPathLengths.clear();
  fCurEvt    = 0;
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
  // actual detector geometry (this is skipped when using 
  // pre-calculated flux interaction probabilities)
  if(fPreSelect) {
       LOG("GMCJDriver", pNOTICE) 
          << "Computing interaction probabilities for max. path lengths";

       Psum = this->ComputeInteractionProbabilities(true /* <- max PL*/);
       Pno  = 1-Psum;
       LOG("GMCJDriver", pNOTICE)
          << "The no-interaction probability (max. path lengths) is: " 
          << 100*Pno << " %";
       if(Pno<0.) {
           LOG("GMCJDriver", pFATAL) 
             << "Negative no-interaction probability! (P = " << 100*Pno << " %)"
             << " Particle E=" << fFluxDriver->Momentum().E() << " type=" << fFluxDriver->PdgCode() << "Psum=" << Psum;
           gAbortingInErr=true;
           exit(1);
       }
       if(R>=1-Pno) {
  	   LOG("GMCJDriver", pNOTICE)  
              << "** Rejecting current flux neutrino";
	   return 0;
       }
  } // preselect 

  bool pl_ok = false;


  // If possible use pre-generated flux neutrino interaction probabilities 
  if(fFluxIntTree){
    Psum = this->PreGenFluxInteractionProbability(); 
  }         
  // Else compute them in the usual manner
  else {
    // Compute (pathLength x density x weight fraction) for all materials
    // in the input geometry, for the neutrino generated by the flux driver
    pl_ok = this->ComputePathLengths();
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
  }


  if(TMath::Abs(Psum) < controls::kASmallNum){
    LOG("GMCJDriver", pNOTICE)
       << "** Rejecting current flux neutrino (has null interaction probability)";
    return 0;
  } 

  // Now decide whether the current neutrino interacts
  Pno  = 1-Psum;
  LOG("GMCJDriver", pNOTICE)
     << "The actual 'no interaction' probability is: " << 100*Pno << " %";
  if(Pno<0.) {
      LOG("GMCJDriver", pFATAL) 
         << "Negative no interactin probability! (P = " << 100*Pno << " %)";

      // print info about what caused the problem
      int                    nupdg = fFluxDriver -> PdgCode  ();
      const TLorentzVector & nup4  = fFluxDriver -> Momentum ();
      const TLorentzVector & nux4  = fFluxDriver -> Position ();

      LOG("GMCJDriver", pWARN)
        << "\n [-] Problematic neutrino: "
        << "\n  |----o PDG-code   : " << nupdg
        << "\n  |----o 4-momentum : " << utils::print::P4AsString(&nup4)
        << "\n  |----o 4-position : " << utils::print::X4AsString(&nux4)
        << "\n Emax : " << fEmax;

      LOG("GMCJDriver", pWARN)
        << "\n Problematic path lengths:" << fCurPathLengths;

      LOG("GMCJDriver", pWARN)
        << "\n Maximum path lengths:" << fMaxPathLengths;

      exit(1);
  }
  if(R>=1-Pno) {
     LOG("GMCJDriver", pNOTICE) 
        << "** Rejecting current flux neutrino";
     return 0;
  }

  //
  // The flux neutrino interacts! 
  //

  // Calculate path lengths for first time and check potential mismatch if 
  // used pre-generated flux interaction probabilities
  if(fFluxIntTree){
    pl_ok = this->ComputePathLengths(); 
    if(!pl_ok) { 
      LOG("GMCJDriver", pFATAL) << "** Cannot calculate path lenths!"; 
      exit(1); 
    }  
    double Psum_curr = this->ComputeInteractionProbabilities(false /* <- actual PL */);
    bool mismatch = TMath::Abs(Psum-Psum_curr) > controls::kASmallNum;    
    if(mismatch){
      LOG("GMCJDriver", pFATAL) << 
          "** Mismatch between pre-calculated and current interaction "<<
          "probabilities!";
      exit(1);
    }
  }

  // Select a target material
  fSelTgtPdg = this->SelectTargetMaterial(R);
  if(fSelTgtPdg==0) {
     LOG("GMCJDriver", pERROR) 
        << "** Rejecting current flux neutrino (failed to select tgt!)";
     return 0;
  }

  // Ask the GEVGDriver object to select and generate an interaction and
  // its kinematics for the selected initial state & neutrino 4-momentum
  this->GenerateEventKinematics();
  if(!fCurEvt) {
     LOG("GMCJDriver", pWARN) 
        << "** Couldn't generate kinematics for selected interaction";
     return 0;
  }

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
        LOG("GMCJDriver", pDEBUG)
          << " (xsec, pl, A)=(" << xsec << "," << pl << "," << A << ")";

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
        LOG("GMCJDriver", pDEBUG)
          << "Pmax=" << pmax;
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

  // propagate current unphysical event mask 
  evgdriver->SetUnphysEventMask(*fUnphysEventMask);

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

  // get interaction cross section
  double xsec = fCurEvt->XSec();

  // get path length in detector along v direction for specified target material
  PathLengthList::const_iterator pliter = fCurPathLengths.find(fSelTgtPdg);
  double path_length = pliter->second;

  // get target material mass number
  int A = pdg::IonPdgCodeToA(fSelTgtPdg);

  // calculate interaction probability
  double P = this->InteractionProbability(xsec, path_length, A);

  //
  // get weight for selected event
  //

  GHepParticle * nu = fCurEvt->Probe();
  int    nu_pdg = nu->Pdg();
  double Ev     = nu->P4()->Energy();
 
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
double GMCJDriver::PreGenFluxInteractionProbability()
{
// Return the pre-computed interaction probability for the current flux 
// neutrino index (entry number in flux file). Exit if not possible as 
// using meaningless interaction probability leads to incorrect physics 
//
  if(!fFluxIntTree){
    LOG("GMCJDriver", pERROR) << 
         "Cannot get pre-computed flux interaction probability as no tree!";
    exit(1);
  }

  assert(fFluxDriver->Index() >= 0); // Check trying to find meaningfull index

  // Check if can find relevant entry and no mismatch in energies -->
  // using correct pre-gen interaction prob file
  bool found_entry = fFluxIntTree->GetEntryWithIndex(fFluxDriver->Index()) > 0;
  bool enu_match = false;
  if(found_entry){
    double rel_err = fBrFluxEnu-fFluxDriver->Momentum().E();
    if(fBrFluxEnu > controls::kASmallNum) rel_err /= fBrFluxEnu;
    enu_match = TMath::Abs(rel_err)<controls::kASmallNum;
    if(enu_match == false){
      LOG("GMCJDriver", pERROR) << 
           "Mismatch between: Enu_curr  = "<< fFluxDriver->Momentum().E() <<
           ", Enu_pre_gen = "<< fBrFluxEnu;
    } 
  }
  else {
    LOG("GMCJDriver", pERROR) << "Cannot find flux entry in interaction prob tree!";
  }

  // Exit if not successful
  bool success = found_entry && enu_match;
  if(!success){
    LOG("GMCJDriver", pFATAL) << 
         "Cannot find pre-generated interaction probability! Check you "<<
         "are using the correct pre-generated interaction prob file "   <<
         "generated using current flux input file with same input "     <<
         "config (same geom TopVol, neutrino species list)";
    exit(1);
  }
  assert(fGlobPmax+controls::kASmallNum>0.0);
  return fBrFluxIntProb/fGlobPmax; 
}
//___________________________________________________________________________
