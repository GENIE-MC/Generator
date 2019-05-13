//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - Feb 04, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 04, 2008 - CA
   The first implementation of this concrete flux driver was first added in 
   the development version 2.3.1
 @ Feb 19, 2008 - CA
   Extended to handle all near detector locations and super-k
 @ Feb 22, 2008 - CA
   Added method to report the actual POT.
 @ Mar 05, 2008 - CA,JD
   Added method to configure the starting z position (upstream of the detector
   face, in detector coord system). Added code to project flux neutrinos from
   z=0 to a configurable z position (somewhere upstream of the detector face)
 @ Mar 27, 2008 - CA
   Added option to recycle the flux ntuple
 @ Mar 31, 2008 - CA
   Handle the flux ntuple weights. Renamed the old implementation of GFluxI
   GenerateNext() to GenerateNext_weighted(). Coded-up a new GenerateNext()
   method using GenerateNext_weighted() + the rejection method to generate
   unweighted flux neutrinos. Added code in LoadBeamSimData() to scan for the 
   max. flux weight for the input location. The 'actual POT' norm factor is 
   updated after each generated flux neutrino taking into account the flux 
   weight variability. Added NFluxNeutrinos() and SumWeight().
 @ May 29, 2008 - CA, PG
   Protect LoadBeamSimData() against non-existent input file
 @ May 31, 2008 - CA
   Added option to keep on recyclying the flux ntuple for an 'infinite' number
   of times (SetNumOfCycles(0)) so that exiting the event generation loop can 
   be controlled by the accumulated POTs or number of events.
 @ June 1, 2008 - CA
   At LoadBeamSimData() added code to scan for number of neutrinos and sum of
   weights for a complete flux ntuple cycle, at the specified detector.
   Added POT_1cycle() to return the flux POT per flux ntuple cycle.
   Added POT_curravg() to return the current average number of POTs taking
   into account the flux ntuple recycling. At complete cycles the reported
   number is not just the average POT but the exact POT.
   Removed the older ActualPOT() method.
 @ June 4, 2008 - CA, AM
   Small modification at the POT normalization to account for the fact that
   flux neutrinos get de-weighted before passed to the event generation code.
   Now pot = file_pot/max_wght rather than file_pot*(nneutrinos/sum_wght)
 @ June 19, 2008 - CA
   Removing some LOG() mesgs speeds up GenerateNext() by a factor of 20 (!)
 @ June 18, 2009 - CA
   Demote warning mesgs if the current flux neutrino is skipped because it 
   is not in the list of neutrinos to be considered. 
   Now this can be intentional (eg. if generating nu_e only).
   In GenerateNext_weighted() Moved the code updating the number of neutrinos
   and sum of weights higher so that force-rejected flux neutrinos still
   count for normalization purposes.
 @ March 6, 2010 - JD
   Made compatible with 10a version of jnubeam flux. Maintained compatibility 
   with 07a version. This involved finding a few more branches - now look for 
   all branches and only use them if they exist. Also added check for required 
   branches. GJPARCNuFluxPassThroughInfo class now passes through flux info in 
   identical format as given in flux file. Any conversions to more usable 
   format happen at later stage (gNtpConv.cxx). In addition also store the flux
   entry number for current neutrino and the deduced flux version.
   Changed method to search for maximum flux weight to avoid seg faults when 
   have large number of flux entries in a file (~1.5E6). 
 @ March 8, 2010 - JD
   Incremented the GJPARCNuFluxPassThroughInfo class def number to reflect
   changes made in previous commit. Added a failsafe which will terminate 
   the job if go through a whole flux cycle without finding a detector location
   matching that specified by user. This avoids potential infinite number of
   cycles.
 @ Feb 4, 2011 - JD
   Made compatible with 11a version of jnubeam flux. Still compatable with older
   versions. Now set the data members of the pass-through class directly as the
   branch addresses of the input tree to reduce amount of duplicated code. 
 @ Feb 22, 2011 - JD
   Added functionality to start looping over input flux file from a random 
   offset. This is to avoid any potential biases when processing very large 
   flux files and always starting from the same position. The default is to
   apply a random offset but this can be switched off using the 
   GJPARCNuFlux::DisableOffset() method.
 @ Feb 22, 2011 - JD
   Implemented the new GFluxI::Clear, GFluxI::Index and GFluxI::GenerateWeighted
   methods needed so that can be used with the new pre-generation of flux 
   interaction probabilities methods added to GMCJDriver. 
 @ Feb 24, 2011 - JD
   Updated list of expected decay modes for the JNuBeam flux neutrinos for >10a
   flux mc. The decay mode is used to infer the neutrino pdg and previously we
   were just skipping them if we didn't recognise the mode - now the job aborts
   as this can lead to unphysical results.
 @ Feb 26, 2011 - JD
   Now check that there is at least one entry with matching flux location (idfd)
   at the LoadBeamSimData stage. Previously were only checking this after a 
   cycle of calling GenerateNext. This stops case where if only looping over a 
   single cycle the user was not warned that there were no flux location matches. 
 @ Jul 05, 2011 - TD
   Code used to match nd detector locations 1-10. Now match detector locations
   up to 51. Change made in preparation for the new sand muon flux (nd13).
 @ Feb 09, 2012 - TD
   Added ability to TChain flux files together. This is so when doing vector
   producion, we can sample all the input flux files, even if the equivalent
   hadd'ed flux file is too large.
 @ Mar 14, 2014 - TD
   Prevent an infinite loop in GenerateNext() when the flux driver has not been
   properly configured by exiting within GenerateNext_weighted().
   LoadBeamSimData() now returns bool, so that the user can catch cases when
   the flux driver has not been properly configured.
*/
//____________________________________________________________________________

#include <cstdlib>
#include <iostream>
#include <cassert>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TSystem.h>

#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Tools/Flux/GJPARCNuFlux.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/StringUtils.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GJPARCNuFlux,genie::flux::GJPARCNuFlux)

using std::endl;
using std::cout;
using std::endl;
using namespace genie;
using namespace genie::flux;

ClassImp(GJPARCNuFluxPassThroughInfo)

//____________________________________________________________________________
GJPARCNuFlux::GJPARCNuFlux()
{
  this->Initialize();
}
//___________________________________________________________________________
GJPARCNuFlux::~GJPARCNuFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
bool GJPARCNuFlux::GenerateNext(void)
{
// Get next (unweighted) flux ntuple entry on the specified detector location
//
  RandomGen * rnd = RandomGen::Instance();
  while(1) {
     // Check for end of flux ntuple
     bool end = this->End();
     if(end) return false;

     // Get next weighted flux ntuple entry
     bool nextok = this->GenerateNext_weighted();
     if(!nextok) continue;

     if(fNCycles==0) {
       LOG("Flux", pNOTICE) 
          << "Got flux entry: " << this->Index()
          << " - Cycle: "<< fICycle << "/ infinite"; 
     } else {
       LOG("Flux", pNOTICE) 
          << "Got flux entry: "<< this->Index() 
          << " - Cycle: "<< fICycle << "/"<< fNCycles; 
     }

     // If de-weighting get fractional weight & decide whether to accept curr flux neutrino
     double f = 1.0;
     if(fGenerateWeighted == false) f = this->Weight();
     LOG("Flux", pNOTICE) 
        << "Curr flux neutrino fractional weight = " << f;
     if(f > (1.+controls::kASmallNum)) {
       LOG("Flux", pERROR) 
           << "** Fractional weight = " << f << " > 1 !!";
     }
     double r = (f < 1.) ? rnd->RndFlux().Rndm() : 0;
     bool accept = (r<f);
     if(accept) {
       return true;
     }

     LOG("Flux", pNOTICE) 
       << "** Rejecting current flux neutrino based on the flux weight only";
  }
  return false;
}
//___________________________________________________________________________
bool GJPARCNuFlux::GenerateNext_weighted(void)
{
// Get next (weighted) flux ntuple entry on the specified detector location
//

  // Reset previously generated neutrino code / 4-p / 4-x
  this->ResetCurrent();

  // Check whether a jnubeam flux ntuple has been loaded
  if( (!fNuFluxTree && fNuFluxUsingTree) || (!fNuFluxChain && !fNuFluxUsingTree) ) {
     LOG("Flux", pFATAL)
          << "The flux driver has not been properly configured";
     //return false; // don't do this - creates an infinite loop!
     exit(1);
  }

  // Read next flux ntuple entry. Use fEntriesThisCycle to keep track of when
  // in new cycle as fIEntry can now have an offset
  if(fEntriesThisCycle >= fNEntries) {
     // Exit if have not found neutrino at specified location for whole cycle
     if(fNDetLocIdFound == 0){
       LOG("Flux", pFATAL)
         << "The input jnubeam flux ntuple contains no entries for detector id "
         << fDetLocId << ". Terminating job!";
       exit(1);
     }
     fNDetLocIdFound = 0; // reset the counter
     fICycle++;
     fIEntry=fOffset;
     fEntriesThisCycle = 0;
     // Run out of entries @ the current cycle.
     // Check whether more (or infinite) number of cycles is requested
     if(fICycle >= fNCycles && fNCycles != 0){
        LOG("Flux", pWARN)
            << "No more entries in input flux neutrino ntuple";
        return false;
     }
  }

  // In addition to getting info to generate event the following also
  // updates pass-through info (= info on the flux neutrino parent particle 
  // that may be stored at an extra branch of the output event tree -alongside 
  // with the generated event branch- for use further upstream in the t2k 
  // analysis chain -eg for beam reweighting etc-)
  bool found_entry;
  if (fNuFluxUsingTree)
    found_entry = fNuFluxTree->GetEntry(fIEntry) > 0;
  else
    found_entry = fNuFluxChain->GetEntry(fIEntry) > 0;
  assert(found_entry);
  fLoadedNeutrino = true;
  fEntriesThisCycle++;
  fIEntry = (fIEntry+1) % fNEntries;

  if (fNuFluxUsingTree) {
    if(fNuFluxSumTree) fNuFluxSumTree->GetEntry(0); // get entry 0 as only 1 entry in tree
  }
  else {
    // get entry corresponding to current tree number in the chain, as only 1 entry in each tree
    if(fNuFluxSumChain) fNuFluxSumChain->GetEntry(fNuFluxChain->GetTreeNumber());
  }

  // check for negative flux weights 
  if(fPassThroughInfo->norm + controls::kASmallNum < 0.0){ 
    LOG("Flux", pERROR) << "Negative flux weight! Will set weight to 0.0";
    fPassThroughInfo->norm  = 0.0;
  } 
  // remember to update fNorm as no longer set to fNuFluxTree branch address
  fNorm = (double) fPassThroughInfo->norm;

  // for 'near detector' flux ntuples make sure that the current entry
  // corresponds to a flux neutrino at the specified detector location
  if(fIsNDLoc           /* nd */  && 
     fDetLocId!=fPassThroughInfo->idfd /* doesn't match specified detector location*/) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("Flux", pNOTICE)
          << "Current flux neutrino not at specified detector location";
#endif
        return false;
  }


  //
  // Handling neutrinos at specified detector location
  //

  // count the number of times we have neutrinos at specified detector location
  fNDetLocIdFound += 1;


  // update the sum of weights & number of neutrinos
  fSumWeight += this->Weight() * fMaxWeight; // Weight returns fNorm/fMaxWeight
  fNNeutrinos++;

  // Convert the current parent paticle decay mode into a neutrino pdg code
  // See:
  // http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/nemode.h
  //  mode    description
  //  11      numu from pi+ 
  //  12      numu from K+  
  //  13      numu from mu- 
  //  21      numu_bar from pi- 
  //  22      numu_bar from K-  
  //  23      numu_bar from mu+ 
  //  31      nue from K+ (Ke3) 
  //  32      nue from K0L(Ke3) 
  //  33      nue from Mu+      
  //  41      nue_bar from K- (Ke3) 
  //  42      nue_bar from K0L(Ke3) 
  //  43      nue_bar from Mu-      
  // Since JNuBeam flux version >= 10a the following modes also expected
  //  14      numu from K+(3)
  //  15      numu from K0(3)
  //  24      numu_bar from K-(3)
  //  25      numu_bar from K0(3)
  //  34      nue from pi+     
  //  44      nue_bar from pi-  
  // In general expect more modes following the rule:
  //  11->19 --> numu
  //  21->29 --> numu_bar
  //  31->39 --> nue
  //  41->49 --> nuebar
  // This is based on example given at:
  //  http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/efill.kumac 
 
  if(fPassThroughInfo->mode >= 11 && fPassThroughInfo->mode <= 19) fgPdgC = kPdgNuMu;
  else if(fPassThroughInfo->mode >= 21 && fPassThroughInfo->mode <= 29) fgPdgC = kPdgAntiNuMu;
  else if(fPassThroughInfo->mode >= 31 && fPassThroughInfo->mode <= 39) fgPdgC = kPdgNuE;
  else if(fPassThroughInfo->mode >= 41 && fPassThroughInfo->mode <= 49) fgPdgC = kPdgAntiNuE;
  else {
    // If here then trying to process a neutrino from an unknown decay mode.
    // Rather than just skipping this flux neutrino the job is aborted to avoid
    // unphysical results. 
    LOG("Flux", pFATAL) << "Unexpected decay mode: "<< fPassThroughInfo->mode <<
                           "  --> unable to infer neutrino pdg! Aborting job!";
    exit(1);  
  }

  // Check neutrino pdg against declared list of neutrino species declared
  // by the current instance of the JPARC neutrino flux driver.
  // No undeclared neutrino species will be accepted at this point as GENIE
  // has already been configured to handle the specified list.
  // Make sure that the appropriate list of flux neutrino species was set at
  // initialization via GJPARCNuFlux::SetFluxParticles(const PDGCodeList &)

  if( ! fPdgCList->ExistsInPDGCodeList(fgPdgC) ) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("Flux", pNOTICE)
       << "Current flux neutrino "
       << "not at the list of neutrinos to be considered at this job.";
#endif
     return false;
  }

  // Check current neutrino energy against the maximum flux neutrino energy declared
  // by the current instance of the JPARC neutrino flux driver.
  // No flux neutrino exceeding that maximum energy will be accepted at this point as
  // that maximum energy has already been used for normalizing the interaction probabilities.
  // Make sure that the appropriate maximum flux neutrino energy was set at 
  // initialization via GJPARCNuFlux::SetMaxEnergy(double Ev)

  if(fPassThroughInfo->Enu > fMaxEv) {
     LOG("Flux", pWARN)
          << "Flux neutrino energy exceeds declared maximum neutrino energy";
     LOG("Flux", pWARN) 
          << "Ev = " << fPassThroughInfo->Enu << "(> Ev{max} = " << fMaxEv << ")";
  }
  
  // Set the current flux neutrino 4-momentum & 4-position

  double pxnu = fPassThroughInfo->Enu * fPassThroughInfo->nnu[0];
  double pynu = fPassThroughInfo->Enu * fPassThroughInfo->nnu[1];
  double pznu = fPassThroughInfo->Enu * fPassThroughInfo->nnu[2];
  double Enu  = fPassThroughInfo->Enu;
  fgP4.SetPxPyPzE (pxnu, pynu, pznu, Enu);

  if(fIsNDLoc) {
    double cm2m = units::cm / units::m;
    double xnu  = cm2m * fPassThroughInfo->xnu;
    double ynu  = cm2m * fPassThroughInfo->ynu;
    double znu  = 0;

    // projected 4-position (from z=0) back to a configurable plane 
    // position (fZ0) upstream of the detector face.
    xnu += (fZ0/fPassThroughInfo->nnu[2])*fPassThroughInfo->nnu[0]; 
    ynu += (fZ0/fPassThroughInfo->nnu[2])*fPassThroughInfo->nnu[1];
    znu = fZ0;

    fgX4.SetXYZT (xnu,  ynu,  znu,  0.);
  } else {
    fgX4.SetXYZT (0.,0.,0.,0.);
  
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pINFO)
	<< "Generated neutrino: "
	<< "\n pdg-code: " << fgPdgC
        << "\n p4: " << utils::print::P4AsShortString(&fgP4)
        << "\n x4: " << utils::print::X4AsString(&fgX4);
#endif
  // Update flux pass through info not set as branch addresses of flux ntuples
  if(fNuFluxUsingTree)
    fPassThroughInfo->fluxentry = this->Index();
  else
    fPassThroughInfo->fluxentry = fNuFluxChain->GetTree()->GetReadEntry();

  std::string filename;
  if(fNuFluxUsingTree)
    filename = fNuFluxFile->GetName();
  else
    filename = fNuFluxChain->GetFile()->GetName();
  std::string::size_type start_pos = filename.rfind("/");
  if (start_pos == std::string::npos) start_pos = 0; else ++start_pos;
  std::string basename(filename,start_pos);
  if(fNuFluxUsingTree)
    fPassThroughInfo->fluxfilename = basename + ":" + fNuFluxTree->GetName();
  else
    fPassThroughInfo->fluxfilename = basename + ":" + fNuFluxChain->GetName();
  return true;
}
//___________________________________________________________________________
double GJPARCNuFlux::POT_1cycle(void)
{
// Compute number of flux POTs / flux ntuple cycle
//
  if( (!fNuFluxTree && fNuFluxUsingTree) || (!fNuFluxChain && !fNuFluxUsingTree) ) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return 0;	
  }

// double pot = fFilePOT * (fNNeutrinosTot1c/fSumWeightTot1c);
//
// Use the max weight instead, since flux neutrinos get de-weighted
// before thrown to the event generation driver
//
  double pot = fFilePOT / fMaxWeight;
  return pot;
}
//___________________________________________________________________________
double GJPARCNuFlux::POT_curravg(void)
{
// Compute current number of flux POTs 
// On complete cycles, that POT number should be exact.
// Within cycles that is only an average number 

  if( (!fNuFluxTree && fNuFluxUsingTree) || (!fNuFluxChain && !fNuFluxUsingTree) ) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return 0;	
  }

// double pot = fNNeutrinos*fFilePOT/fSumWeightTot1c;
//
// See also comment at POT_1cycle()
//
  double cnt   = (double)fNNeutrinos;
  double cnt1c = (double)fNNeutrinosTot1c;
  double pot  = (cnt/cnt1c) * this->POT_1cycle();
  return pot;
}
//___________________________________________________________________________
long int GJPARCNuFlux::Index(void)
{ 
// Return the current flux entry index. If GenerateNext has not yet been 
// called then return -1. 
//
  if(fLoadedNeutrino){
    // subtract 1 as fIEntry was incremented since call to TTree::GetEntry
    // and deal with special case where fIEntry-1 is last entry in cycle 
    return fIEntry == 0 ? fNEntries - 1 : fIEntry-1;
  }
  // return -1 if no neutrino loaded since last call to this->ResetCurrent()
  return -1;
}
//___________________________________________________________________________
bool GJPARCNuFlux::LoadBeamSimData(string filename, string detector_location)
{
// Loads in a jnubeam beam simulation root file (converted from hbook format)
// into the GJPARCNuFlux driver.
// The detector location can be any of:
//  "sk","nd1" (<-2km),"nd5" (<-nd280),...,"nd10"

  LOG("Flux", pNOTICE) 
        << "Loading jnubeam flux tree from ROOT file: " << filename;
  LOG("Flux", pNOTICE) 
        << "Detector location: " << detector_location;

  // Check to see if its a single flux file (/dir/root_filename.0.root)
  // or a sequence of files to be tchained (e.g. /dir/root_filename@0@100)
  fNuFluxUsingTree = true;
  if (filename.find('@') != string::npos) { fNuFluxUsingTree = false; }

  vector<string> filenamev = utils::str::Split(filename,"@");
  string fileroot = "";
  int firstfile = -1, lastfile = -1;

  if (!fNuFluxUsingTree) {
    if (filenamev.size() != 3) {
      LOG("Flux", pFATAL)
	<< "Flux filename should be specfied as either:\n"
	<< "\t For a single input file:  /dir/root_filename.#.root\n"
	<< "\t For multiple input files: /dir/root_filename@#@#";
      exit(1);
    }
    fileroot  = filenamev[0];
    firstfile = atoi(filenamev[1].c_str());
    lastfile  = atoi(filenamev[2].c_str());
    LOG("Flux", pNOTICE) 
      << "Chaining beam simulation output files with stem: " << fileroot
      << " and run numbers in the range: [" << firstfile << ", " << firstfile << "]";
  }

  if (fNuFluxUsingTree) {
    bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
    if (!is_accessible) {
      LOG("Flux", pFATAL)
	<< "The input jnubeam flux file doesn't exist! Initialization failed!";
      exit(1);
    }
  }
  else {
    bool please_exit = false;
    for (int i = firstfile; i < lastfile+1; i++) {
      bool is_accessible = ! (gSystem->AccessPathName( Form("%s.%i.root",fileroot.c_str(),i) ));
      if (!is_accessible) {
	LOG("Flux", pFATAL)
	  << "The input jnubeam flux file " << Form("%s.%i.root",fileroot.c_str(),i) 
	  << "doesn't exist! Initialization failed!";
	please_exit = true;
      }
    }
    if(please_exit)
      exit(1);
  }

  fDetLoc   = detector_location;   
  fDetLocId = this->DLocName2Id(fDetLoc);  

  if(fDetLocId == 0) {
    LOG("Flux", pERROR) 
         << " ** Unknown input detector location: " << fDetLoc;
    return false;
  }

  fIsFDLoc = (fDetLocId==-1);
  fIsNDLoc = (fDetLocId>0);

  if (!fNuFluxUsingTree) {
    string ntuple_name = (fIsNDLoc) ? "h3002" : "h2000";
    fNuFluxChain = new TChain(ntuple_name.c_str());
    int result = fNuFluxChain->Add( Form("%s.%i.root",fileroot.c_str(),firstfile), 0);
    if (result != 1 && fIsNDLoc) {
      LOG("Flux", pINFO)
	<< "Could not find tree h3002 in file " << Form("%s.%i.root",fileroot.c_str(),firstfile)
	<< " Trying tree h3001";
      delete fNuFluxChain;
      ntuple_name = "h3001";
      fNuFluxChain = new TChain(ntuple_name.c_str());
      result = fNuFluxChain->Add( Form("%s.%i.root",fileroot.c_str(),firstfile), 0);
    }
    if (result != 1) {
      LOG("Flux", pERROR)
	<< "** Couldn't get flux tree: " << ntuple_name;
      return false;
    }
    
    for (int i = firstfile+1; i < lastfile+1; i++) {
      result = fNuFluxChain->Add( Form("%s.%i.root",fileroot.c_str(),i), 0 );
      if (result == 0)
	LOG("Flux", pERROR)
	  << "** Couldn't get flux tree " << ntuple_name << " in file " << Form("%s.%i.root",fileroot.c_str(),i);
    fNEntries = fNuFluxChain->GetEntries();
    }
  }
  
  else {
    fNuFluxFile = new TFile(filename.c_str(), "read");
    if(fNuFluxFile) {
      // nd treename can be h3002 or h3001 depending on fluxfile version
      string ntuple_name = (fIsNDLoc) ? "h3002" : "h2000";
      fNuFluxTree = (TTree*) fNuFluxFile->Get(ntuple_name.c_str());
      if(!fNuFluxTree && fIsNDLoc){
	ntuple_name = "h3001";
	fNuFluxTree = (TTree*) fNuFluxFile->Get(ntuple_name.c_str());
      }
      LOG("Flux", pINFO)   
	<< "Getting flux tree: " << ntuple_name;
      if(!fNuFluxTree) {
	LOG("Flux", pERROR) 
	  << "** Couldn't get flux tree: " << ntuple_name;
	return false;
      }
    } else {
      LOG("Flux", pERROR) << "** Couldn't open: " << filename;
      return false;
    }
    
    fNEntries = fNuFluxTree->GetEntries();
  }

  LOG("Flux", pNOTICE) 
    << "Loaded flux tree contains " <<  fNEntries << " entries";
  
  LOG("Flux", pDEBUG) 
    << "Getting tree branches & setting leaf addresses";

  // try to get all the branches that we know about and only set address if
  // they exist
  bool missing_critical = false;
  TBranch * fBr = 0;
  GJPARCNuFluxPassThroughInfo * info = fPassThroughInfo;

  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("norm")) ) fBr->SetAddress(&info->norm);
  else if( (fBr = fNuFluxChain->GetBranch("norm")) ) fNuFluxChain->SetBranchAddress("norm",&info->norm);
  else { 
    LOG("Flux", pFATAL) <<"Cannot find flux branch: norm";
    missing_critical = true;
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("Enu")) ) fBr->SetAddress(&info->Enu);
  else if( (fBr = fNuFluxChain->GetBranch("Enu")) ) fNuFluxChain->SetBranchAddress("Enu",&info->Enu);
  else {
    LOG("Flux", pFATAL) <<"Cannot find flux branch: Enu";  
    missing_critical = true;
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("ppid")) ) fBr->SetAddress(&info->ppid);
  else if( (fBr = fNuFluxChain->GetBranch("ppid")) ) fNuFluxChain->SetBranchAddress("ppid",&info->ppid);
  else {
    LOG("Flux", pFATAL) <<"Cannot find flux branch: ppid"; 
    missing_critical = true;
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("mode")) ) fBr->SetAddress(&info->mode);
  else if( (fBr = fNuFluxChain->GetBranch("mode")) ) fNuFluxChain->SetBranchAddress("mode",&info->mode);
  else { 
    LOG("Flux", pFATAL) <<"Cannot find flux branch: mode";
    missing_critical = true;
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("rnu")) ) fBr->SetAddress(&info->rnu);
  else if( (fBr = fNuFluxChain->GetBranch("rnu")) ) fNuFluxChain->SetBranchAddress("rnu",&info->rnu);
  else if(fIsNDLoc){ // Only required for ND location
    LOG("Flux", pFATAL) <<"Cannot find flux branch: rnu";
    missing_critical = true;
  } 
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("xnu")) ) fBr->SetAddress(&info->xnu);
  else if( (fBr = fNuFluxChain->GetBranch("xnu")) ) fNuFluxChain->SetBranchAddress("xnu",&info->xnu);
  else if(fIsNDLoc){ // Only required for ND location
    LOG("Flux", pFATAL) <<"Cannot find flux branch: xnu"; 
    missing_critical = true;
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("ynu")) ) fBr->SetAddress(&info->ynu);
  else if( (fBr = fNuFluxChain->GetBranch("ynu")) ) fNuFluxChain->SetBranchAddress("ynu",&info->ynu);
  else if(fIsNDLoc) { // Only required for ND location
    LOG("Flux", pFATAL) <<"Cannot find flux branch: ynu"; 
    missing_critical = true;
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("nnu")) ) fBr->SetAddress(info->nnu);
  else if( (fBr = fNuFluxChain->GetBranch("nnu")) ) fNuFluxChain->SetBranchAddress("nnu",info->nnu);
  else if(fIsNDLoc){ // Only required for ND location
    LOG("Flux", pFATAL) <<"Cannot find flux branch: nnu";
    missing_critical = true;
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("idfd")) ) fBr->SetAddress(&info->idfd);
  else if( (fBr = fNuFluxChain->GetBranch("idfd")) ) fNuFluxChain->SetBranchAddress("idfd",&info->idfd);
  else if(fIsNDLoc){ // Only required for ND location
    LOG("Flux", pFATAL) <<"Cannot find flux branch: idfd"; 
    missing_critical = true;
  }
  // check that have found essential branches
  if(missing_critical){
    LOG("Flux", pFATAL)
     << "Unable to find critical information in the flux ntuple! Initialization failed!";
    exit(1);
  }
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("ppi"))     ) fBr->SetAddress(&info->ppi);
  else if( (fBr = fNuFluxChain->GetBranch("ppi"))     ) fNuFluxChain->SetBranchAddress("ppi",&info->ppi);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("xpi"))     ) fBr->SetAddress(info->xpi);
  else if( (fBr = fNuFluxChain->GetBranch("xpi"))     ) fNuFluxChain->SetBranchAddress("xpi",info->xpi);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("npi"))     ) fBr->SetAddress(info->npi);
  else if( (fBr = fNuFluxChain->GetBranch("npi"))     ) fNuFluxChain->SetBranchAddress("npi",info->npi);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("ppi0"))    ) fBr->SetAddress(&info->ppi0);
  else if( (fBr = fNuFluxChain->GetBranch("ppi0"))    ) fNuFluxChain->SetBranchAddress("ppi0",&info->ppi0);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("xpi0"))    ) fBr->SetAddress(info->xpi0);
  else if( (fBr = fNuFluxChain->GetBranch("xpi0"))    ) fNuFluxChain->SetBranchAddress("xpi0",info->xpi0);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("npi0"))    ) fBr->SetAddress(info->npi0);
  else if( (fBr = fNuFluxChain->GetBranch("npi0"))    ) fNuFluxChain->SetBranchAddress("npi0",info->npi0);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("nvtx0"))   ) fBr->SetAddress(&info->nvtx0);
  else if( (fBr = fNuFluxChain->GetBranch("nvtx0"))   ) fNuFluxChain->SetBranchAddress("nvtx0",&info->nvtx0);
  // Following branches only present since flux version 10a
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("cospibm")) ) fBr->SetAddress(&info->cospibm);
  else if( (fBr = fNuFluxChain->GetBranch("cospibm")) ) fNuFluxChain->SetBranchAddress("cospibm",&info->cospibm);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("cospi0bm"))) fBr->SetAddress(&info->cospi0bm);
  else if( (fBr = fNuFluxChain->GetBranch("cospi0bm"))) fNuFluxChain->SetBranchAddress("cospi0bm",&info->cospi0bm);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gamom0"))  ) fBr->SetAddress(&info->gamom0);
  else if( (fBr = fNuFluxChain->GetBranch("gamom0"))  ) fNuFluxChain->SetBranchAddress("gamom0",&info->gamom0);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gipart"))  ) fBr->SetAddress(&info->gipart);
  else if( (fBr = fNuFluxChain->GetBranch("gipart"))  ) fNuFluxChain->SetBranchAddress("gipart",&info->gipart);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gvec0"))   ) fBr->SetAddress(info->gvec0);
  else if( (fBr = fNuFluxChain->GetBranch("gvec0"))   ) fNuFluxChain->SetBranchAddress("gvec0",info->gvec0);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gpos0"))   ) fBr->SetAddress(info->gpos0);
  else if( (fBr = fNuFluxChain->GetBranch("gpos0"))   ) fNuFluxChain->SetBranchAddress("gpos0",info->gpos0);
  // Following branches only present since flux vesion 10d
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("ng"))      ) fBr->SetAddress(&info->ng);
  else if( (fBr = fNuFluxChain->GetBranch("ng"))      ) fNuFluxChain->SetBranchAddress("ng",&info->ng);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gpid"))    ) fBr->SetAddress(info->gpid);
  else if( (fBr = fNuFluxChain->GetBranch("gpid"))    ) fNuFluxChain->SetBranchAddress("gpid",info->gpid);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gmec"))    ) fBr->SetAddress(info->gmec);
  else if( (fBr = fNuFluxChain->GetBranch("gmec"))    ) fNuFluxChain->SetBranchAddress("gmec",info->gmec);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gvx"))     ) fBr->SetAddress(info->gvx);
  else if( (fBr = fNuFluxChain->GetBranch("gvx"))     ) fNuFluxChain->SetBranchAddress("gvx",info->gvx);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gvy"))     ) fBr->SetAddress(info->gvy);
  else if( (fBr = fNuFluxChain->GetBranch("gvy"))     ) fNuFluxChain->SetBranchAddress("gvy",info->gvy);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gvz"))     ) fBr->SetAddress(info->gvz);
  else if( (fBr = fNuFluxChain->GetBranch("gvz"))     ) fNuFluxChain->SetBranchAddress("gvz",info->gvz);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gpx"))     ) fBr->SetAddress(info->gpx);
  else if( (fBr = fNuFluxChain->GetBranch("gpx"))     ) fNuFluxChain->SetBranchAddress("gpx",info->gpx);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gpy"))     ) fBr->SetAddress(info->gpy);
  else if( (fBr = fNuFluxChain->GetBranch("gpy"))     ) fNuFluxChain->SetBranchAddress("gpy",info->gpy);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gpz"))     ) fBr->SetAddress(info->gpz);
  else if( (fBr = fNuFluxChain->GetBranch("gpz"))     ) fNuFluxChain->SetBranchAddress("gpz",info->gpz);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gmat"))    ) fBr->SetAddress(info->gmat);
  else if( (fBr = fNuFluxChain->GetBranch("gmat"))    ) fNuFluxChain->SetBranchAddress("gmat",info->gmat);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gdistc"))  ) fBr->SetAddress(info->gdistc);
  else if( (fBr = fNuFluxChain->GetBranch("gdistc"))  ) fNuFluxChain->SetBranchAddress("gdistc",info->gdistc);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gdistal")) ) fBr->SetAddress(&info->gdistal);
  else if( (fBr = fNuFluxChain->GetBranch("gdistal")) ) fNuFluxChain->SetBranchAddress("gdistal",&info->gdistal);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gdistti")) ) fBr->SetAddress(&info->gdistti);
  else if( (fBr = fNuFluxChain->GetBranch("gdistti")) ) fNuFluxChain->SetBranchAddress("gdistti",&info->gdistti);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gdistfe")) ) fBr->SetAddress(&info->gdistfe);
  else if( (fBr = fNuFluxChain->GetBranch("gdistfe")) ) fNuFluxChain->SetBranchAddress("gdistfe",&info->gdistfe);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("gcosbm"))  ) fBr->SetAddress(info->gcosbm);
  else if( (fBr = fNuFluxChain->GetBranch("gcosbm"))  ) fNuFluxChain->SetBranchAddress("gcosbm",info->gcosbm);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("Enusk"))   ) fBr->SetAddress(&info->Enusk);
  else if( (fBr = fNuFluxChain->GetBranch("Enusk"))   ) fNuFluxChain->SetBranchAddress("Enusk",&info->Enusk);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("normsk"))  ) fBr->SetAddress(&info->normsk);
  else if( (fBr = fNuFluxChain->GetBranch("normsk"))  ) fNuFluxChain->SetBranchAddress("normsk",&info->normsk);
  if( fNuFluxUsingTree && (fBr = fNuFluxTree->GetBranch("anorm"))   ) fBr->SetAddress(&info->anorm);
  else if( (fBr = fNuFluxChain->GetBranch("anorm"))   ) fNuFluxChain->SetBranchAddress("anorm",&info->anorm);

  // Look for the flux file summary info tree (only expected for > 10a flux versions) 
  if( fNuFluxUsingTree && (fNuFluxSumTree = (TTree*) fNuFluxFile->Get("h1000")) ){
    if( (fBr = fNuFluxSumTree->GetBranch("version"))) fBr->SetAddress(&info->version);
    if( (fBr = fNuFluxSumTree->GetBranch("ntrig"))  ) fBr->SetAddress(&info->ntrig);
    if( (fBr = fNuFluxSumTree->GetBranch("tuneid")) ) fBr->SetAddress(&info->tuneid);
    if( (fBr = fNuFluxSumTree->GetBranch("pint"))   ) fBr->SetAddress(&info->pint);
    if( (fBr = fNuFluxSumTree->GetBranch("bpos"))   ) fBr->SetAddress(info->bpos);
    if( (fBr = fNuFluxSumTree->GetBranch("btilt"))  ) fBr->SetAddress(info->btilt);
    if( (fBr = fNuFluxSumTree->GetBranch("brms"))   ) fBr->SetAddress(info->brms);
    if( (fBr = fNuFluxSumTree->GetBranch("emit"))   ) fBr->SetAddress(info->emit);
    if( (fBr = fNuFluxSumTree->GetBranch("alpha"))  ) fBr->SetAddress(info->alpha);
    if( (fBr = fNuFluxSumTree->GetBranch("hcur"))   ) fBr->SetAddress(info->hcur);
    if( (fBr = fNuFluxSumTree->GetBranch("rand"))   ) fBr->SetAddress(&info->rand);
    if( (fBr = fNuFluxSumTree->GetBranch("rseed"))  ) fBr->SetAddress(info->rseed);
  }

  // Look for the flux file summary info tree (only expected for > 10a flux versions) 
  if ( !fNuFluxUsingTree ) {
    fNuFluxSumChain = new TChain("h1000");
    int result = fNuFluxSumChain->Add( Form("%s.%i.root",fileroot.c_str(),firstfile), 0 );
    if (result==1) {
      for (int i = firstfile+1; i < lastfile+1; i++) {
	result = fNuFluxSumChain->Add( Form("%s.%i.root",fileroot.c_str(),i), 0 );
      }
      if( (fBr = fNuFluxSumChain->GetBranch("version"))) fNuFluxSumChain->SetBranchAddress("version",&info->version);
      if( (fBr = fNuFluxSumChain->GetBranch("ntrig"))  ) fNuFluxSumChain->SetBranchAddress("ntrig",&info->ntrig);
      if( (fBr = fNuFluxSumChain->GetBranch("tuneid")) ) fNuFluxSumChain->SetBranchAddress("tuneid",&info->tuneid);
      if( (fBr = fNuFluxSumChain->GetBranch("pint"))   ) fNuFluxSumChain->SetBranchAddress("pint",&info->pint);
      if( (fBr = fNuFluxSumChain->GetBranch("bpos"))   ) fNuFluxSumChain->SetBranchAddress("bpos",info->bpos);
      if( (fBr = fNuFluxSumChain->GetBranch("btilt"))  ) fNuFluxSumChain->SetBranchAddress("btilt",info->btilt);
      if( (fBr = fNuFluxSumChain->GetBranch("brms"))   ) fNuFluxSumChain->SetBranchAddress("brms",info->brms);
      if( (fBr = fNuFluxSumChain->GetBranch("emit"))   ) fNuFluxSumChain->SetBranchAddress("emit",info->emit);
      if( (fBr = fNuFluxSumChain->GetBranch("alpha"))  ) fNuFluxSumChain->SetBranchAddress("alpha",info->alpha);
      if( (fBr = fNuFluxSumChain->GetBranch("hcur"))   ) fNuFluxSumChain->SetBranchAddress("hcur",info->hcur);
      if( (fBr = fNuFluxSumChain->GetBranch("rand"))   ) fNuFluxSumChain->SetBranchAddress("rand",&info->rand);
      if( (fBr = fNuFluxSumChain->GetBranch("rseed"))  ) fNuFluxSumChain->SetBranchAddress("rseed",info->rseed);
    }
  }

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle = 1;

  // sum-up weights & number of neutrinos for the specified location
  // over a complete cycle. Also record the maximum weight as previous
  // method using TTree::GetV1() seg faulted for more than ~1.5E6 entries
  fSumWeightTot1c  = 0;
  fNNeutrinosTot1c = 0;
  fNDetLocIdFound = 0;
  for(int ientry = 0; ientry < fNEntries; ientry++) {
     if (fNuFluxUsingTree)
       fNuFluxTree->GetEntry(ientry);
     else
       fNuFluxChain->GetEntry(ientry);
     // check for negative flux weights
     if(fPassThroughInfo->norm + controls::kASmallNum < 0.0){ 
       LOG("Flux", pERROR) << "Negative flux weight! Will set weight to 0.0";
       fPassThroughInfo->norm  = 0.0;
     } 
     fNorm = (double) fPassThroughInfo->norm;
     // update maximum weight
     fMaxWeight = TMath::Max(fMaxWeight, (double) fPassThroughInfo->norm); 
     // compare detector location (see GenerateNext_weighted() for details)
     if(fIsNDLoc && fDetLocId!=fPassThroughInfo->idfd) continue;
     fSumWeightTot1c += fNorm;
     fNNeutrinosTot1c++;
     fNDetLocIdFound++;
  }
  // Exit if have not found neutrino at specified location for whole cycle
  if(fNDetLocIdFound == 0){
    LOG("Flux", pFATAL)
     << "The input jnubeam flux ntuple contains no entries for detector id "
     << fDetLocId << ". Terminating job!";
    exit(1); 
  }
  fNDetLocIdFound = 0; // reset the counter

  LOG("Flux", pNOTICE) << "Maximum flux weight = " << fMaxWeight;  
  if(fMaxWeight <=0 ) {
      LOG("Flux", pFATAL) << "Non-positive maximum flux weight!";
      exit(1);
  }

  LOG("Flux", pINFO)
    << "Totals / cycle: #neutrinos = " << fNNeutrinosTot1c 
    << ", Sum{Weights} = " << fSumWeightTot1c;

  if(fUseRandomOffset){
    this->RandomOffset();  // Random start point when looping over ntuple
  }

  return true;
}
//___________________________________________________________________________
void GJPARCNuFlux::SetFluxParticles(const PDGCodeList & particles)
{
  if(!fPdgCList) {
     fPdgCList = new PDGCodeList;
  }
  fPdgCList->Copy(particles);

  LOG("Flux", pINFO) 
    << "Declared list of neutrino species: " << *fPdgCList;
}
//___________________________________________________________________________
void GJPARCNuFlux::SetMaxEnergy(double Ev)
{
  fMaxEv = TMath::Max(0.,Ev);

  LOG("Flux", pINFO) 
    << "Declared maximum flux neutrino energy: " << fMaxEv;
}
//___________________________________________________________________________
void GJPARCNuFlux::SetFilePOT(double pot)
{
// The flux normalization is in /N POT/det [ND] or /N POT/cm^2 [FD]
// This method specifies N (typically 1E+21)

  fFilePOT = pot;
}
//___________________________________________________________________________
void GJPARCNuFlux::SetUpstreamZ(double z0)
{
// The flux neutrino position (x,y) is given at the detector coord system
// at z=0. This method sets the preferred starting z position upstream of 
// the upstream detector face. Each flux neutrino will be backtracked from
// z=0 to the input z0.

  fZ0 = z0;
}
//___________________________________________________________________________
void GJPARCNuFlux::SetNumOfCycles(int n)
{
// The flux ntuples can be recycled for a number of times to boost generated
// event statistics without requiring enormous beam simulation statistics.
// That option determines how many times the driver is going to cycle through
// the input flux ntuple. 
// With n=0 the flux ntuple will be recycled an infinite amount of times so
// that the event generation loop can exit only on a POT or event num check.

  fNCycles = TMath::Max(0, n);
}
//___________________________________________________________________________
void GJPARCNuFlux::GenerateWeighted(bool gen_weighted)
{
// Ignore the flux weights when getting next flux neutrino. Always set to
// true for unweighted event generation but need option to switch off when
// using pre-calculated interaction probabilities for each flux ray to speed
// up event generation.
  fGenerateWeighted = gen_weighted;
}
//___________________________________________________________________________
void GJPARCNuFlux::RandomOffset()
{
// Choose a random number between 0-->fNEntries to set as start point for 
// looping over flux ntuple. May be necessary when looping over very large 
// flux files as always starting fromthe same point may introduce biases 
// (inversely proportional to number of cycles). This method resets the 
// starting value for fIEntry so must be called before any call to GenerateNext
// is made.
//
  double ran_frac = RandomGen::Instance()->RndFlux().Rndm();
  long int offset = (long int) floor(ran_frac * fNEntries); 
  LOG("Flux", pERROR) << "Setting flux driver to start looping over entries "
                      << "with offset of "<< offset;
  fIEntry = fOffset = offset;
}
//___________________________________________________________________________
void GJPARCNuFlux::Clear(Option_t * opt)
{
// If opt = "CycleHistory" then:
// Reset all counters and state variables from any previous cycles. This 
// should be called if, for instance, before event generation this flux
// driver had been used for pre-calculating interaction probabilities. Does
// not reset initial state variables such as flux particle types, POTs etc... 
//
  if(std::strcmp(opt, "CycleHistory") == 0){
    // Reset so that generate de-weighted events
    this->GenerateWeighted(false);
    // Reset cycle counters
    fICycle          = 0;
    fSumWeight       = 0;
    fNNeutrinos      = 0;
    fIEntry          = fOffset;
    fEntriesThisCycle = 0;
  }
}
//___________________________________________________________________________
void GJPARCNuFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GJPARCNuFlux driver";

  fMaxEv           = 0;
  fPdgCList        = new PDGCodeList;
  fPassThroughInfo = new GJPARCNuFluxPassThroughInfo;

  fNuFluxFile      = 0;
  fNuFluxTree      = 0;
  fNuFluxChain     = 0;
  fNuFluxSumTree   = 0;
  fNuFluxSumChain  = 0;
  fNuFluxUsingTree = true;
  fDetLoc          = "";
  fDetLocId        = 0;
  fNDetLocIdFound  = 0;
  fIsFDLoc         = false;
  fIsNDLoc         = false;

  fNEntries        = 0;
  fIEntry          = 0;
  fEntriesThisCycle= 0;
  fOffset          = 0;
  fNorm            = 0.;
  fMaxWeight       = 0;
  fFilePOT         = 0;
  fZ0              = 0;
  fNCycles         = 0;
  fICycle          = 0;        
  fSumWeight       = 0;
  fNNeutrinos      = 0;
  fSumWeightTot1c  = 0;
  fNNeutrinosTot1c = 0;
  fGenerateWeighted= false;
  fUseRandomOffset = true;
  fLoadedNeutrino  = false;

  this->SetDefaults();
  this->ResetCurrent();
}
//___________________________________________________________________________
void GJPARCNuFlux::SetDefaults(void)
{
// - Set default neutrino species list (nue, nuebar, numu, numubar) and 
//   maximum energy (25 GeV).
//   These defaults can be overwritten by user calls (at the driver init) to 
//   GJPARCNuFlux::SetMaxEnergy(double Ev) and
//   GJPARCNuFlux::SetFluxParticles(const PDGCodeList & particles)
// - Set the default file normalization to 1E+21 POT
// - Set the default flux neutrino start z position at -5m (z=0 is the
//   detector centre).
// - Set number of cycles to 1

  LOG("Flux", pNOTICE) << "Setting default GJPARCNuFlux driver options";

  PDGCodeList particles;
  particles.push_back(kPdgNuMu);
  particles.push_back(kPdgAntiNuMu);
  particles.push_back(kPdgNuE);
  particles.push_back(kPdgAntiNuE);

  this->SetFluxParticles(particles);
  this->SetMaxEnergy(25./*GeV*/);
  this->SetFilePOT(1E+21);
  this->SetUpstreamZ(-5.0); 
  this->SetNumOfCycles(1);
}
//___________________________________________________________________________
void GJPARCNuFlux::ResetCurrent(void)
{
// reset running values of neutrino pdg-code, 4-position & 4-momentum
// and the input ntuple leaves

  fLoadedNeutrino = false;

  fgPdgC = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);

  fPassThroughInfo->Reset();
}
//___________________________________________________________________________
void GJPARCNuFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fPdgCList)        delete fPdgCList;
  if (fPassThroughInfo) delete fPassThroughInfo;

  if (fNuFluxFile) {
	fNuFluxFile->Close();
	delete fNuFluxFile;
  }
}
//___________________________________________________________________________
int GJPARCNuFlux::DLocName2Id(string name)
{
// detector location: name -> int id
// sk  -> -1
// nd1 -> +1
// nd2 -> +2
// ...

  if(name == "sk"  ) return  -1;

  TString temp;
  for (int i=1; i<51; i++) {
    temp.Form("nd%d",i);
    if(name == temp.Data()) return i;
  }

  return 0;
}
//___________________________________________________________________________
GJPARCNuFluxPassThroughInfo::GJPARCNuFluxPassThroughInfo() :
TObject()
{
  this->Reset();
}
//___________________________________________________________________________
GJPARCNuFluxPassThroughInfo::GJPARCNuFluxPassThroughInfo(
                        const GJPARCNuFluxPassThroughInfo & info) :
TObject()
{
  fluxentry  = info.fluxentry;
  fluxfilename  = info.fluxfilename;
  Enu        = info.Enu;
  ppid       = info.ppid;
  mode       = info.mode;
  ppi        = info.ppi;
  ppi0       = info.ppi0; 
  nvtx0      = info.nvtx0;
  cospibm    = info.cospibm;
  cospi0bm   = info.cospi0bm;
  idfd       = info.idfd;
  gamom0     = info.gamom0;
  gipart     = info.gipart;
  xnu        = info.xnu;
  ynu        = info.ynu;
  rnu        = info.rnu;
  for(int i = 0; i<3; i++){
    nnu[i] = info.nnu[i];
    xpi[i] = info.xpi[i];
    xpi0[i] = info.xpi0[i];
    gpos0[i] = info.gpos0[i];
    npi[i]  = info.npi[i];
    npi0[i] = info.npi0[i];
    gvec0[i] = info.gvec0[i];
  }
  ng = info.ng;
  for(int ip = 0; ip<fNgmax; ip++){
    gpid[ip] = info.gpid[ip];
    gmec[ip] = info.gmec[ip];
    gcosbm[ip] = info.gcosbm[ip];
    gvx[ip] = info.gvx[ip];
    gvy[ip] = info.gvy[ip];
    gvz[ip] = info.gvz[ip];
    gpx[ip] = info.gpx[ip]; 
    gpy[ip] = info.gpy[ip]; 
    gpz[ip] = info.gpz[ip]; 
    gmat[ip] = info.gmat[ip];
    gdistc[ip] = info.gdistc[ip];
    gdistal[ip] = info.gdistal[ip];
    gdistti[ip] = info.gdistti[ip];
    gdistfe[ip] = info.gdistfe[ip]; 
  }
  norm = info.norm;
  Enusk = info.Enusk;
  normsk = info.normsk;
  anorm = info.anorm;
  version = info.version;
  ntrig = info.ntrig;
  tuneid = info.tuneid;
  pint = info.pint;
  rand = info.rand;
  for(int i = 0; i < 2; i++){
    bpos[i] = info.bpos[i];
    btilt[i] = info.btilt[i];
    brms[i] = info.brms[i];
    emit[i] = info.emit[i];
    alpha[i] = info.alpha[i];
    rseed[i] = info.rseed[i];
  }
  for(int i = 0; i < 3; i++) hcur[i] = info.hcur[i]; 
}
//___________________________________________________________________________
void GJPARCNuFluxPassThroughInfo::Reset(){
  fluxentry= -1;
  fluxfilename = "Not-set";
  Enu      = 0;
  ppid     = 0;
  mode     = 0;
  ppi      = 0.;
  norm     = 0;
  cospibm  = 0.;
  nvtx0    = 0;
  ppi0     = 0.;
  idfd     = 0;
  rnu      = -999999.;
  xnu      = -999999.;
  ynu      = -999999.;
  cospi0bm = -999999.;
  gipart   = -1;
  gamom0   = 0.;
  for(int i=0; i<3; i++){
    nnu[i] = 0.;
    xpi[i] = -999999.;
    npi[i] = 0.;
    xpi0[i] = -999999.;
    npi0[i] = 0.;
    gpos0[i] = -999999.;
    gvec0[i] = 0.;
  }
  ng = -1;
  for(int ip = 0; ip<fNgmax; ip++){
    gpid[ip] = -999999; 
    gmec[ip] = -999999;
    gcosbm[ip] = -999999.;
    gvx[ip] = -999999.;
    gvy[ip] = -999999.;
    gvz[ip] = -999999.;
    gpx[ip] = -999999.; 
    gpy[ip] = -999999.; 
    gpz[ip] = -999999.;
    gmat[ip] = -999999;
    gdistc[ip] = -999999.;
    gdistal[ip] = -999999.;
    gdistti[ip] = -999999.;
    gdistfe[ip] = -999999.; 
  }
  Enusk = -999999.;
  normsk = -999999.; 
  anorm = -999999.;
  version = -999999.;
  ntrig = -999999;
  tuneid = -999999; 
  pint = -999999;
  rand = -999999;
  for(int i = 0; i < 2; i++){
    bpos[i] = -999999.;
    btilt[i] = -999999.;
    brms[i] = -999999.;
    emit[i] = -999999.;
    alpha[i] = -999999.; 
    rseed[i] = -999999;
  }
  for(int i = 0; i < 3; i++) hcur[i] = -999999.;
}
//___________________________________________________________________________
namespace genie {
namespace flux  {
  ostream & operator << (
    ostream & stream, const genie::flux::GJPARCNuFluxPassThroughInfo & info) 
    {
      stream << "\n idfd       = " << info.idfd
             << "\n norm       = " << info.norm
             << "\n flux entry = " << info.fluxentry
             << "\n flux file  = " << info.fluxfilename
             << "\n Enu        = " << info.Enu
             << "\n geant code = " << info.ppid
             << "\n (pdg code) = " << pdg::GeantToPdg(info.ppid)   
             << "\n decay mode = " << info.mode
             << "\n nvtx0      = " << info.nvtx0
             << "\n |momentum| @ decay       = " << info.ppi    
             << "\n position_vector @ decay  = (" 
                 << info.xpi[0] << ", " 
                 << info.xpi[1] << ", " 
                 << info.xpi[2]  << ")"
             << "\n direction_vector @ decay = (" 
                 << info.npi[0] << ", " 
                 << info.npi[1] << ", " 
                 << info.npi[2]  << ")"
             << "\n |momentum| @ prod.       = " << info.ppi0 
             << "\n position_vector @ prod.  = (" 
                 << info.xpi0[0] << ", " 
                 << info.xpi0[1] << ", " 
                 << info.xpi0[2]  << ")"
             << "\n direction_vector @ prod. = (" 
                 << info.npi0[0] << ", " 
                 << info.npi0[1] << ", " 
                 << info.npi0[2]  << ")"
             << "\n cospibm = " << info.cospibm
             << "\n cospi0bm = " << info.cospi0bm
             << "\n Plus additional info if flux version is later than 07a" 
             << endl;

    return stream;
  }
}//flux 
}//genie
//___________________________________________________________________________


