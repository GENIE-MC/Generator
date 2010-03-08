//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Feb 04, 2008

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
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "Conventions/Units.h"
#include "Conventions/GBuild.h"
#include "FluxDrivers/GJPARCNuFlux.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "Utils/MathUtils.h"
#include "Utils/PrintUtils.h"

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
          << "Got flux entry: " << fIEntry 
          << " - Cycle: "<< fICycle << "/ infinite"; 
     } else {
       LOG("Flux", pNOTICE) 
          << "Got flux entry: "<< fIEntry 
          << " - Cycle: "<< fICycle << "/"<< fNCycles; 
     }

     // Get fractional weight & decide whether to accept curr flux neutrino
     double f = this->Weight() / fMaxWeight;
     LOG("Flux", pNOTICE) 
        << "Curr flux neutrino fractional weight = " << f;
     if(f > 1.) {
       LOG("Flux", pERROR) 
           << "** Fractional weight = " << f << " > 1 !!";
     }
     double r = (f < 1.) ? rnd->RndFlux().Rndm() : 0;
     bool accept = (r<f);
     if(accept) return true;

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
  if(!fNuFluxTree) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return false;	
  }

  // Read next flux ntuple entry
  if(fIEntry >= fNEntries) {
     // Exit if have not found neutrino at specified location for whole cycle
     if(fNDetLocIdFound == 0){
       LOG("Flux", pFATAL)
         << "The input jnubeam flux ntuple contains no entries for detector id "
         << fDetLocId << ". Terminating job!";
       exit(1);
     }
     fNDetLocIdFound = 0; // reset the counter
     // Run out of entries @ the current cycle.
     // Check whether more (or infinite) number of cycles is requested
     if(fICycle < fNCycles || fNCycles == 0 ) {
        fICycle++;
        fIEntry=0;
     } else {
        LOG("Flux", pWARN) 
            << "No more entries in input flux neutrino ntuple";
        return false;	
     }
  }

  fNuFluxTree->GetEntry(fIEntry);
  fIEntry++;

  // For 'near detector' flux ntuples make sure that the current entry
  // corresponds to a flux neutrino at the specified detector location
  if(fIsNDLoc           /* nd */  && 
     fDetLocId!=fLfIdfd /* doesn't match specified detector location*/) {
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
  fSumWeight += this->Weight();
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

  if(fLfMode >= 11 && fLfMode <= 13) fgPdgC = kPdgNuMu;
  if(fLfMode >= 21 && fLfMode <= 23) fgPdgC = kPdgAntiNuMu;
  if(fLfMode >= 31 && fLfMode <= 33) fgPdgC = kPdgNuE;
  if(fLfMode >= 41 && fLfMode <= 43) fgPdgC = kPdgAntiNuE;

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

  if(fLfEnu > fMaxEv) {
     LOG("Flux", pWARN) 
          << "Flux neutrino energy exceeds declared maximum neutrino energy";
     LOG("Flux", pWARN) 
          << "Ev = " << fLfEnu << "(> Ev{max} = " << fMaxEv << ")";
  }
  
  // Set the current flux neutrino 4-momentum & 4-position

  double pxnu = fLfEnu * fLfNnu[0];
  double pynu = fLfEnu * fLfNnu[1];
  double pznu = fLfEnu * fLfNnu[2];
  double Enu  = fLfEnu;
  fgP4.SetPxPyPzE (pxnu, pynu, pznu, Enu);

  if(fIsNDLoc) {
    double cm2m = units::cm / units::m;
    double xnu  = cm2m * fLfXnu;
    double ynu  = cm2m * fLfYnu;
    double znu  = 0;

    // projected 4-position (from z=0) back to a configurable plane 
    // position (fZ0) upstream of the detector face.
    xnu += (fZ0/fLfNnu[2])*fLfNnu[0]; 
    ynu += (fZ0/fLfNnu[2])*fLfNnu[1];
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
  // Update pass-through info (= info on the flux neutrino parent particle 
  // that may be stored at an extra branch of the output event tree -alongside 
  // with the generated event branch- for use further upstream in the t2k 
  // analysis chain -eg for beam reweighting etc-)
  fPassThroughInfo -> fluxentry= fIEntry-1; // Offset required so that start at 0 
  fPassThroughInfo -> ppid     = fLfPpid; 
  fPassThroughInfo -> mode     = fLfMode; 
  fPassThroughInfo -> ppi      = fLfPpi; 
  fPassThroughInfo -> ppi0     = fLfPpi0;  
  fPassThroughInfo -> nvtx0    = (int) fLfNVtx0;
  fPassThroughInfo -> cospibm  = fLfCospibm;
  fPassThroughInfo -> cospi0bm = fLfCospi0bm;
  fPassThroughInfo -> idfd     = fLfIdfd;
  fPassThroughInfo -> gamom0 = fLfGamom0;
  fPassThroughInfo -> gipart = (int) fLfGipart;
  for(int i = 0; i < 3; i++){
    fPassThroughInfo -> xpi[i]   = fLfXpi[i];
    fPassThroughInfo -> npi[i]   = fLfNpi[i];
    fPassThroughInfo -> xpi0[i]  = fLfXpi0[i];
    fPassThroughInfo -> npi0[i]  = fLfNpi0[i];
    fPassThroughInfo -> gpos0[i] = fLfGpos0[i];
    fPassThroughInfo -> gvec0[i] = fLfGvec0[i];
  }

  return true;
}
//___________________________________________________________________________
double GJPARCNuFlux::POT_1cycle(void)
{
// Compute number of flux POTs / flux ntuple cycle
//
  if(!fNuFluxTree) {
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

  if(!fNuFluxTree) {
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
void GJPARCNuFlux::LoadBeamSimData(string filename, string detector_location)
{
// Loads in a jnubeam beam simulation root file (converted from hbook format)
// into the GJPARCNuFlux driver.
// The detector location can be any of:
//  "sk","nd1" (<-2km),"nd5" (<-nd280),...,"nd10"

  LOG("Flux", pNOTICE) 
        << "Loading jnubeam flux tree from ROOT file: " << filename;
  LOG("Flux", pNOTICE) 
        << "Detector location: " << detector_location;

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
    LOG("Flux", pFATAL)
     << "The input jnubeam flux file doesn't exist! Initialization failed!";
    exit(1);
  }

  fDetLoc   = detector_location;   
  fDetLocId = this->DLocName2Id(fDetLoc);  

  if(fDetLocId == 0) {
    LOG("Flux", pERROR) 
         << " ** Unknown input detector location: " << fDetLoc;
    return;
  }

  fIsFDLoc = (fDetLocId==-1);
  fIsNDLoc = (fDetLocId>0);

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
      return;
    }
  } else {
    LOG("Flux", pERROR) << "** Couldn't open: " << filename;
    return;
  }

  fNEntries = fNuFluxTree->GetEntries();
  LOG("Flux", pNOTICE) 
      << "Loaded flux tree contains " <<  fNEntries << " entries";

  LOG("Flux", pDEBUG) 
      << "Getting tree branches & setting leaf addresses";

  // try to get all the branches that we know about and only set address if
  // they exist
  if( (fBrNorm = fNuFluxTree->GetBranch("norm" )) ) fBrNorm->SetAddress(&fLfNorm);
  else LOG("Flux", pFATAL)
        <<"Cannot find flux branch: norm";
  if( (fBrEnu  = fNuFluxTree->GetBranch("Enu"  )) ) fBrEnu->SetAddress(&fLfEnu);
  else LOG("Flux", pFATAL)
        <<"Cannot find flux branch: Enu";  
  if( (fBrPpid = fNuFluxTree->GetBranch("ppid" )) ) fBrPpid->SetAddress(&fLfPpid);
  else LOG("Flux", pFATAL)
        <<"Cannot find flux branch: ppid"; 
  if( (fBrMode = fNuFluxTree->GetBranch("mode" )) ) fBrMode->SetAddress(&fLfMode);
  else LOG("Flux", pFATAL)
        <<"Cannot find flux branch: mode"; 
  if( (fBrRnu  = fNuFluxTree->GetBranch("rnu"  )) ) fBrRnu->SetAddress(&fLfRnu);
  else if(fIsNDLoc) LOG("Flux", pFATAL)
        <<"Cannot find flux branch: rnu"; 
  if( (fBrXnu  = fNuFluxTree->GetBranch("xnu"  )) ) fBrXnu->SetAddress(&fLfXnu);
  else if(fIsNDLoc) LOG("Flux", pFATAL)
        <<"Cannot find flux branch: xnu"; 
  if( (fBrYnu  = fNuFluxTree->GetBranch("ynu"  )) ) fBrYnu->SetAddress(&fLfYnu);
  else if(fIsNDLoc) LOG("Flux", pFATAL)
        <<"Cannot find flux branch: ynu"; 
  if( (fBrNnu  = fNuFluxTree->GetBranch("nnu"  )) ) fBrNnu->SetAddress(fLfNnu);
  else if(fIsNDLoc) LOG("Flux", pFATAL)
        <<"Cannot find flux branch: nnu"; 
  if( (fBrIdfd = fNuFluxTree->GetBranch("idfd" )) ) fBrIdfd->SetAddress(&fLfIdfd);
  else if(fIsNDLoc) LOG("Flux", pFATAL)
        <<"Cannot find flux branch: idfd"; 
  // check that have found essential branches
  if((fBrEnu && fBrPpid && fBrMode &&  
     ((fBrRnu && fBrXnu && fBrYnu && fBrNnu && fBrIdfd) || !fIsNDLoc))==false){
    LOG("Flux", pFATAL)
     << "Unable to find critical information in the flux ntuple! Initialization failed!";
    exit(1);
  }
  if( (fBrPpi      = fNuFluxTree->GetBranch("ppi"     )) ) fBrPpi->SetAddress(&fLfPpi);
  if( (fBrXpi      = fNuFluxTree->GetBranch("xpi"     )) ) fBrXpi->SetAddress(fLfXpi);
  if( (fBrNpi      = fNuFluxTree->GetBranch("npi"     )) ) fBrNpi->SetAddress(fLfNpi);
  if( (fBrPpi0     = fNuFluxTree->GetBranch("ppi0"    )) ) fBrPpi0->SetAddress(&fLfPpi0);
  if( (fBrXpi0     = fNuFluxTree->GetBranch("xpi0"    )) ) fBrXpi0->SetAddress(fLfXpi0);
  if( (fBrNpi0     = fNuFluxTree->GetBranch("npi0"    )) ) fBrNpi0->SetAddress(fLfNpi0);
  if( (fBrNVtx0    = fNuFluxTree->GetBranch("nvtx0"   )) ) fBrNVtx0->SetAddress(&fLfNVtx0);
  if( (fBrCospibm  = fNuFluxTree->GetBranch("cospibm" )) ) fBrCospibm->SetAddress(&fLfCospibm);
  if( (fBrCospi0bm = fNuFluxTree->GetBranch("cospi0bm")) ) fBrCospi0bm->SetAddress(&fLfCospi0bm);
  if( (fBrGamom0   = fNuFluxTree->GetBranch("gamom0"  )) ) fBrGamom0->SetAddress(&fLfGamom0);
  if( (fBrGipart   = fNuFluxTree->GetBranch("gipart"  )) ) fBrGipart->SetAddress(&fLfGipart);
  if( (fBrGvec0    = fNuFluxTree->GetBranch("gvec0"   )) ) fBrGvec0->SetAddress(fLfGvec0);
  if( (fBrGpos0    = fNuFluxTree->GetBranch("gpos0"   )) ) fBrGpos0->SetAddress(fLfGpos0);

  // deduce if have 07a or 10a flux files. Not used in any logic and is just for record
  fFluxVersion = "10a";
  if(fNuFluxTree -> GetBranch("gamom0") == NULL) fFluxVersion = "07a";
  LOG("Flux", pNOTICE)
      << "Have flux version: " <<  fFluxVersion.c_str();

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle = 1;

  // sum-up weights & number of neutrinos for the specified location
  // over a complete cycle. Also record the maximum weight as previous
  // method using TTree::GetV1() seg faulted for more than ~1.5E6 entries
  fSumWeightTot1c  = 0;
  fNNeutrinosTot1c = 0;
  for(int ientry = 0; ientry < fNEntries; ientry++) {
     fNuFluxTree->GetEntry(ientry);
     // update maximum weight
     fMaxWeight = TMath::Max(fMaxWeight, (double) fLfNorm); 
     // compare detector location (see GenerateNext_weighted() for details)
     if(fIsNDLoc && fDetLocId!=fLfIdfd) continue;
     fSumWeightTot1c += this->Weight();
     fNNeutrinosTot1c++;
  }
  LOG("Flux", pNOTICE) << "Maximum flux weight = " << fMaxWeight;  
  if(fMaxWeight <=0 ) {
      LOG("Flux", pFATAL) << "Non-positive maximum flux weight!";
      exit(1);
  }

  LOG("Flux", pINFO)
    << "Totals / cycle: #neutrinos = " << fNNeutrinosTot1c 
    << ", Sum{Weights} = " << fSumWeightTot1c;
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
void GJPARCNuFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GJPARCNuFlux driver";

  fMaxEv           = 0;
  fPdgCList        = new PDGCodeList;
  fPassThroughInfo = new GJPARCNuFluxPassThroughInfo;

  fNuFluxFile      = 0;
  fFluxVersion     = "";
  fNuFluxTree      = 0;
  fDetLoc          = "";
  fDetLocId        = 0;
  fNDetLocIdFound  = 0;
  fIsFDLoc         = false;
  fIsNDLoc         = false;

  fNEntries        = 0;
  fIEntry          = 0;
  fMaxWeight       = 0;
  fFilePOT         = 0;
  fZ0              = 0;
  fNCycles         = 0;
  fICycle          = 0;        
  fSumWeight       = 0;
  fNNeutrinos      = 0;
  fSumWeightTot1c  = 0;
  fNNeutrinosTot1c = 0;

  fBrNorm          = 0;
  fBrIdfd          = 0;
  fBrEnu           = 0;
  fBrRnu           = 0;
  fBrXnu           = 0;
  fBrYnu           = 0;
  fBrNnu           = 0;
  fBrPpid          = 0;
  fBrMode          = 0;
  fBrPpi           = 0;
  fBrXpi           = 0;
  fBrNpi           = 0;
  fBrCospibm       = 0;      
  fBrPpi0          = 0;
  fBrXpi0          = 0;
  fBrNpi0          = 0;
  fBrCospi0bm      = 0;
  fBrNVtx0         = 0;
  fBrGipart        = 0;
  fBrGpos0         = 0;
  fBrGvec0         = 0;
  fBrGamom0        = 0;


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

  fgPdgC = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);

  fLfEnu      = 0;  
  fLfPpid     = 0;          
  fLfMode     = 0;          
  fLfPpi      = 0;          
  fLfNorm     = 0;          
  fLfCospibm  = 0;       
  fLfNVtx0    = 0;       
  fLfPpi0     = 0;          
  fLfIdfd     = 0;
  fLfRnu      = 0;    
  fLfXnu      = 0;    
  fLfYnu      = 0;          
  fLfCospi0bm = 0;       
  fLfGipart   = 0;
  fLfGamom0   = 0;
  for(int i=0; i<3; i++){
    fLfNnu[i] = 0.;       
    fLfXpi[i] = -999.999;      
    fLfNpi[i] = 0.;        
    fLfXpi0[i] = -999.999;    
    fLfNpi0[i] = 0.;    
    fLfGpos0[i] = -999.999;    
    fLfGvec0[i] = 0.;
  }    
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
  if(name == "nd1" ) return   1;
  if(name == "nd2" ) return   2;
  if(name == "nd3" ) return   3;
  if(name == "nd4" ) return   4;
  if(name == "nd5" ) return   5;
  if(name == "nd6" ) return   6;
  if(name == "nd7" ) return   7;
  if(name == "nd8" ) return   8;
  if(name == "nd9" ) return   9;
  if(name == "nd10") return  10;

  return 0;
}
//___________________________________________________________________________
GJPARCNuFluxPassThroughInfo::GJPARCNuFluxPassThroughInfo() :
TObject()
{
  fluxentry = -1;
  ppid      = 0;
  mode      = 0;
  ppi       = 0.;
  ppi0      = 0.;
  nvtx0     = 0;
  cospibm   = 0.; 
  cospi0bm  = 0.; 
  idfd      = 0;
  gamom0    = 0.;
  gipart    = -1;
  for(int i = 0; i<3; i++){ 
    xpi[i] = xpi0[i] = gpos0[i] = -999.999;
    npi[i]  = npi0[i] = gvec0[i] = 0.;
  }
}
//___________________________________________________________________________
GJPARCNuFluxPassThroughInfo::GJPARCNuFluxPassThroughInfo(
                        const GJPARCNuFluxPassThroughInfo & info) :
TObject()
{
  fluxentry  = info.fluxentry;
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
  for(int i = 0; i<3; i++){
    xpi[i] = info.xpi[i];
    xpi0[i] = info.xpi0[i];
    gpos0[i] = info.gpos0[i];
    npi[i]  = info.npi[i];
    npi0[i] = info.npi0[i];
    gvec0[i] = info.gvec0[i];
  }

}
//___________________________________________________________________________
namespace genie {
namespace flux  {
  ostream & operator << (
    ostream & stream, const genie::flux::GJPARCNuFluxPassThroughInfo & info) 
    {
      stream << "\n idfd       = " << info.idfd
             << "\n flux entry = " << info.fluxentry
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
             << endl;

    return stream;
  }
}//flux 
}//genie
//___________________________________________________________________________


