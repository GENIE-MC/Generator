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
   This concrete flux driver was first added in development version 2.3.1

*/
//____________________________________________________________________________

#include <TFile.h>
#include <TTree.h>

#include "Conventions/Units.h"
#include "FluxDrivers/GJPARCNuFlux.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::flux;

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
  // Reset previously generated neutrino code / 4-p / 4-x
  this->ResetCurrent();

  // Read next flux ntuple entry
  if(fIEntry >= fNEntries) {
     LOG("Flux", pWARN) 
          << "No more entries in input flux neutrino ntuple";
     return false;	
  }
  fNuFluxTree->GetEntry(fIEntry);
  fIEntry++;

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
          << "Declared list of neutrino species: " << *fPdgCList;;
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
     return false;	
  }
  
  // Set the current flux neutrino 4-momentum & 4-position

  double pxnu = fLfEnu * fLfNnu[0];
  double pynu = fLfEnu * fLfNnu[1];
  double pznu = fLfEnu * fLfNnu[2];
  double Enu  = fLfEnu;
  double xnu  = fLfXnu * (units::m / units::cm);
  double ynu  = fLfYnu * (units::m / units::cm);
  double znu  = 1 * (units::m / units::cm);
 
  fgP4.SetPxPyPzE (pxnu, pynu, pznu, Enu);
  fgX4.SetXYZT    (xnu,  ynu,  znu,  0.);

  // Print-out & return

  LOG("Flux", pINFO) 
	<< "Generated neutrino: "
	<< "\n pdg-code: " << fgPdgC
        << "\n p4: " << utils::print::P4AsShortString(&fgP4)
        << "\n x4: " << utils::print::X4AsString(&fgX4);

  return true;
}
//___________________________________________________________________________
void GJPARCNuFlux::LoadFile(string filename)
{
  LOG("Flux", pNOTICE) 
        << "Loading jnubeam flux tree from ROOT file: " << filename;

  fNuFluxFile = new TFile(filename.c_str(), "read");
  if(fNuFluxFile) {
      LOG("Flux", pINFO) << "Getting flux tree h3001 (nd280)";
      fNuFluxTree = (TTree*) fNuFluxFile->Get("h3001");
      if(!fNuFluxTree) {
          LOG("Flux", pERROR) << "** Couldn't get flux tree h3001";
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

  fBrNorm      = fNuFluxTree -> GetBranch ("norm");
  fBrIdfd      = fNuFluxTree -> GetBranch ("idfd");
  fBrEnu       = fNuFluxTree -> GetBranch ("Enu");
  fBrRnu       = fNuFluxTree -> GetBranch ("rnu");
  fBrXnu       = fNuFluxTree -> GetBranch ("xnu");
  fBrYnu       = fNuFluxTree -> GetBranch ("ynu");
  fBrNnu       = fNuFluxTree -> GetBranch ("nnu");
  fBrPpid      = fNuFluxTree -> GetBranch ("ppid");
  fBrMode      = fNuFluxTree -> GetBranch ("mode");
  fBrPpi       = fNuFluxTree -> GetBranch ("ppi");
  fBrXpi       = fNuFluxTree -> GetBranch ("xpi");
  fBrNpi       = fNuFluxTree -> GetBranch ("npi");
  fBrCospibm   = fNuFluxTree -> GetBranch ("cospibm");
  fBrPpi0      = fNuFluxTree -> GetBranch ("ppi0");
  fBrXpi0      = fNuFluxTree -> GetBranch ("xpi0");
  fBrNpi0      = fNuFluxTree -> GetBranch ("npi0");
  fBrCospi0bm  = fNuFluxTree -> GetBranch ("cospi0bm");

  fBrNorm      -> SetAddress (&fLfNorm);
  fBrIdfd      -> SetAddress (&fLfIdfd);
  fBrEnu       -> SetAddress (&fLfEnu);
  fBrRnu       -> SetAddress (&fLfRnu);
  fBrXnu       -> SetAddress (&fLfXnu);
  fBrYnu       -> SetAddress (&fLfYnu);
  fBrNnu       -> SetAddress ( fLfNnu);
  fBrPpid      -> SetAddress (&fLfPpid);
  fBrMode      -> SetAddress (&fLfMode);
  fBrPpi       -> SetAddress (&fLfPpi);
  fBrXpi       -> SetAddress ( fLfXpi);
  fBrNpi       -> SetAddress ( fLfNpi);
  fBrCospibm   -> SetAddress (&fLfCospibm);
  fBrPpi0      -> SetAddress (&fLfPpi0);
  fBrXpi0      -> SetAddress ( fLfXpi0);
  fBrNpi0      -> SetAddress ( fLfNpi0);
  fBrCospi0bm  -> SetAddress (&fLfCospi0bm);
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
void GJPARCNuFlux::SetDetectorId(int /*detector*/)
{

}
//___________________________________________________________________________
void GJPARCNuFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GJPARCNuFlux driver";

  fMaxEv        = 0;
  fPdgCList     = new PDGCodeList;
  fNuFluxFile   = 0;
  fNuFluxTree   = 0;
  fNEntries     = 0;
  fIEntry       = 0;
  fBrNorm       = 0;
  fBrIdfd       = 0;
  fBrEnu        = 0;
  fBrRnu        = 0;
  fBrXnu        = 0;
  fBrYnu        = 0;
  fBrNnu        = 0;
  fBrPpid       = 0;
  fBrMode       = 0;
  fBrPpi        = 0;
  fBrXpi        = 0;
  fBrNpi        = 0;
  fBrCospibm    = 0;      
  fBrPpi0       = 0;
  fBrXpi0       = 0;
  fBrNpi0       = 0;
  fBrCospi0bm   = 0;

  this->SetDefaults();
  this->ResetCurrent();
}
//___________________________________________________________________________
void GJPARCNuFlux::SetDefaults(void)
{
// Set default neutrino species list (nue, nuebar, numu, numubar) and maximum
// energy (20 GeV).
// These defaults can be overwritten by user calls (at the driver init) to 
// GJPARCNuFlux::SetMaxEnergy(double Ev) and
// GJPARCNuFlux::SetFluxParticles(const PDGCodeList & particles)
//
  LOG("Flux", pNOTICE) << "Setting default GJPARCNuFlux driver options";

  PDGCodeList particles;
  particles.push_back(kPdgNuMu);
  particles.push_back(kPdgAntiNuMu);
  particles.push_back(kPdgNuE);
  particles.push_back(kPdgAntiNuE);

  this->SetFluxParticles(particles);
  this->SetMaxEnergy(20./*GeV*/);
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
}
//___________________________________________________________________________
void GJPARCNuFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fPdgCList) 
	delete fPdgCList;

  if (fNuFluxFile) {
	fNuFluxFile->Close();
	delete fNuFluxFile;
  }
}
//___________________________________________________________________________

