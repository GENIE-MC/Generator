//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
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
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "Conventions/Units.h"
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

  if(fNCycles==0) {
    LOG("Flux", pNOTICE) 
      << "Reading flux entry: " << fIEntry 
      << " - Cycle: " << fICycle << "/ infinite"; 
  } else {
    LOG("Flux", pNOTICE) 
      << "Reading flux entry: " << fIEntry 
      << " - Cycle: " << fICycle << "/" << fNCycles; 
  }

  fNuFluxTree->GetEntry(fIEntry);
  fIEntry++;

  // For 'near detector' flux ntuples make sure that the current entry
  // corresponds to a flux neutrino at the specified detector location
  if(fIsNDLoc           /* nd */  && 
     fDetLocId!=fLfIdfd /* doesn't match specified detector location*/) {
        LOG("Flux", pNOTICE) 
          << "Current flux neutrino not at specified detector location";
        return false;		
  }

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
     LOG("Flux", pWARN) 
          << "Unknown decay mode or decay mode producing an undeclared neurtino species";
     LOG("Flux", pWARN) 
          << "Declared list of neutrino species: " << *fPdgCList;
     LOG("Flux", pWARN) 
          << "Decay mode: " << fLfMode;
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
  LOG("Flux", pINFO) 
	<< "Generated neutrino: "
	<< "\n pdg-code: " << fgPdgC
        << "\n p4: " << utils::print::P4AsShortString(&fgP4)
        << "\n x4: " << utils::print::X4AsString(&fgX4);

  // Update pass-through info (= info on the flux neutrino parent particle 
  // that may be stored at an extra branch of the output event tree -alongside 
  // with the generated event branch- for use further upstream in the t2k 
  // analysis chain -eg for beam reweighting etc-)
  fPassThroughInfo -> pdg       = pdg::GeantToPdg(fLfPpid);
  fPassThroughInfo -> decayMode = fLfMode;
  fPassThroughInfo -> decayP    = fLfPpi;
  fPassThroughInfo -> decayX    = fLfXpi[0];
  fPassThroughInfo -> decayY    = fLfXpi[1];
  fPassThroughInfo -> decayZ    = fLfXpi[2];
  fPassThroughInfo -> decayDirX = fLfNpi[0];
  fPassThroughInfo -> decayDirY = fLfNpi[1];
  fPassThroughInfo -> decayDirZ = fLfNpi[2];
  fPassThroughInfo -> prodP     = fLfPpi0;
  fPassThroughInfo -> prodX     = fLfXpi0[0];
  fPassThroughInfo -> prodY     = fLfXpi0[1];
  fPassThroughInfo -> prodZ     = fLfXpi0[2];
  fPassThroughInfo -> prodDirX  = fLfNpi0[0];
  fPassThroughInfo -> prodDirY  = fLfNpi0[1];
  fPassThroughInfo -> prodDirZ  = fLfNpi0[2];
  fPassThroughInfo -> prodNVtx  = fLfNVtx0;

  // update the sum of weights & number of neutrinos
  fSumWeight += this->Weight();
  fNNeutrinos++;

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
  double pot = fNNeutrinosTot1c*fFilePOT/fSumWeightTot1c;
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
  double pot = fNNeutrinos*fFilePOT/fSumWeightTot1c;
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
      string ntuple_name = (fIsNDLoc) ? "h3001" : "h2000";
      LOG("Flux", pINFO) 
           << "Getting flux tree: " << ntuple_name;
      fNuFluxTree = (TTree*) fNuFluxFile->Get(ntuple_name.c_str());
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

  // the following branches can be found in both 'sk' and 'nd' flux ntuples
  fBrNorm       = fNuFluxTree -> GetBranch ("norm");
  fBrEnu        = fNuFluxTree -> GetBranch ("Enu");
  fBrPpid       = fNuFluxTree -> GetBranch ("ppid");
  fBrMode       = fNuFluxTree -> GetBranch ("mode");
  fBrPpi        = fNuFluxTree -> GetBranch ("ppi");
  fBrXpi        = fNuFluxTree -> GetBranch ("xpi");
  fBrNpi        = fNuFluxTree -> GetBranch ("npi");
  fBrCospibm    = fNuFluxTree -> GetBranch ("cospibm");
  fBrPpi0       = fNuFluxTree -> GetBranch ("ppi0");
  fBrXpi0       = fNuFluxTree -> GetBranch ("xpi0");
  fBrNpi0       = fNuFluxTree -> GetBranch ("npi0");
  // the following branches can be found only in 'nd' flux ntuples
  if(fIsNDLoc) {
    fBrIdfd     = fNuFluxTree -> GetBranch ("idfd");
    fBrRnu      = fNuFluxTree -> GetBranch ("rnu");
    fBrXnu      = fNuFluxTree -> GetBranch ("xnu");
    fBrYnu      = fNuFluxTree -> GetBranch ("ynu");
    fBrNnu      = fNuFluxTree -> GetBranch ("nnu");
    fBrCospi0bm = fNuFluxTree -> GetBranch ("cospi0bm");
  }
  // the following branches can be found only in 'fd' flux ntuples
  if(fIsFDLoc) {
    fBrNVtx0    = fNuFluxTree -> GetBranch ("nvtx0");
  }

  // set the leaf addresses for the above branches
  fBrNorm      -> SetAddress (&fLfNorm);
  fBrEnu       -> SetAddress (&fLfEnu);
  fBrPpid      -> SetAddress (&fLfPpid);
  fBrMode      -> SetAddress (&fLfMode);
  fBrPpi       -> SetAddress (&fLfPpi);
  fBrXpi       -> SetAddress ( fLfXpi);
  fBrNpi       -> SetAddress ( fLfNpi);
  fBrCospibm   -> SetAddress (&fLfCospibm);
  fBrPpi0      -> SetAddress (&fLfPpi0);
  fBrXpi0      -> SetAddress ( fLfXpi0);
  fBrNpi0      -> SetAddress ( fLfNpi0);
  if(fIsNDLoc) {
    fBrIdfd    -> SetAddress (&fLfIdfd);
    fBrRnu     -> SetAddress (&fLfRnu);
    fBrXnu     -> SetAddress (&fLfXnu);
    fBrYnu     -> SetAddress (&fLfYnu);
    fBrNnu     -> SetAddress ( fLfNnu);
    fBrCospi0bm-> SetAddress (&fLfCospi0bm);
  }
  if(fIsFDLoc) {
    fBrNVtx0   -> SetAddress (&fLfNVtx0);
  }

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle = 1;

  // scan for the maximum weight  
  fNuFluxTree->Draw("norm","","goff");
  Long64_t idx = TMath::LocMax(
                     fNuFluxTree->GetSelectedRows(),
                     fNuFluxTree->GetV1());
  fMaxWeight = fNuFluxTree->GetV1()[idx];
  LOG("Flux", pNOTICE) << "Maximum flux weight = " << fMaxWeight; 
  if(fMaxWeight <=0 ) {
      LOG("Flux", pFATAL) << "Non-positive maximum flux weight!";
      exit(1);
  }

  // sum-up weights & number of neutrinos for the specified location
  // over a complete cycle
  fSumWeightTot1c  = 0;
  fNNeutrinosTot1c = 0;
  for(int ientry = 0; ientry < fNEntries; ientry++) {
     fNuFluxTree->GetEntry(ientry);
     // compare detector location (see GenerateNext_weighted() for details)
     if(fIsNDLoc && fDetLocId!=fLfIdfd) continue;
     fSumWeightTot1c += this->Weight();
     fNNeutrinosTot1c++;
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
  fNuFluxTree      = 0;
  fDetLoc          = "";
  fDetLocId        = 0;
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

  fLfNorm     = 0;          
  fLfIdfd     = 0;   
  fLfEnu      = 0;  
  fLfRnu      = 0;    
  fLfXnu      = 0;    
  fLfYnu      = 0;          
  for(int i=0; i<3; i++) fLfNnu[i] = 0;       
  fLfPpid     = 0;          
  fLfMode     = 0;          
  fLfPpi      = 0;           
  for(int i=0; i<3; i++) fLfXpi[i] = 0;      
  for(int i=0; i<3; i++) fLfNpi[i] = 0;        
  fLfCospibm  = 0;       
  fLfPpi0     = 0;          
  for(int i=0; i<3; i++) fLfXpi0[i] = 0;    
  for(int i=0; i<3; i++) fLfNpi0[i] = 0;    
  fLfCospi0bm = 0;       
  fLfNVtx0    = 0;       
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
TObject(),
pdg       (0),
decayMode (0),
decayP    (0.),
decayX    (0.),
decayY    (0.), 
decayZ    (0.), 
decayDirX (0.), 
decayDirY (0.), 
decayDirZ (0.),
prodP     (0.), 
prodX     (0.), 
prodY     (0.),  
prodZ     (0.),  
prodDirX  (0.),  
prodDirY  (0.),  
prodDirZ  (0.),
prodNVtx  (0.)
{

}
//___________________________________________________________________________
GJPARCNuFluxPassThroughInfo::GJPARCNuFluxPassThroughInfo(
                        const GJPARCNuFluxPassThroughInfo & info) :
TObject(),
pdg       (info.pdg),
decayMode (info.decayMode),
decayP    (info.decayP),
decayX    (info.decayX),
decayY    (info.decayY), 
decayZ    (info.decayZ), 
decayDirX (info.decayDirX), 
decayDirY (info.decayDirY), 
decayDirZ (info.decayDirZ),
prodP     (info.prodP), 
prodX     (info.prodX), 
prodY     (info.prodY),  
prodZ     (info.prodZ),  
prodDirX  (info.prodDirX),  
prodDirY  (info.prodDirY),  
prodDirZ  (info.prodDirZ),
prodNVtx  (info.prodNVtx)
{

}
//___________________________________________________________________________
namespace genie {
namespace flux  {
  ostream & operator << (
    ostream & stream, const genie::flux::GJPARCNuFluxPassThroughInfo & info) 
    {
      stream << "\n pdg code   = " << info.pdg   
             << "\n decay mode = " << info.decayMode 
             << "\n |momentum| @ decay       = " << info.decayP    
             << "\n position_vector @ decay  = (" 
                 << info.decayX << ", " 
                 << info.decayY << ", " 
                 << info.decayZ  << ")"
             << "\n direction_vector @ decay = (" 
                 << info.decayDirX << ", " 
                 << info.decayDirY << ", " 
                 << info.decayDirZ  << ")"
             << "\n |momentum| @ prod.       = " << info.prodP    
             << "\n position_vector @ prod.  = (" 
                 << info.prodX << ", " 
                 << info.prodY << ", " 
                 << info.prodZ  << ")"
             << "\n direction_vector @ prod. = (" 
                 << info.prodDirX << ", " 
                 << info.prodDirY << ", " 
                 << info.prodDirZ  << ")"
             << "\n parent produced in vtx num (fd only) = " << info.prodNVtx
             << endl;

     return stream;
  }
}//flux 
}//genie
//___________________________________________________________________________


