//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Robert Hatcher <rhatcher@fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ June 27, 2008 - CA
   The first implementation of this concrete flux driver was first added in
   the development version 2.5.1
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "Conventions/Units.h"
#include "Conventions/GBuild.h"
#include "FluxDrivers/GNuMIFlux.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "Utils/MathUtils.h"
#include "Utils/PrintUtils.h"

using std::endl;

using namespace genie;
using namespace genie::flux;

ClassImp(GNuMIFluxPassThroughInfo)

//____________________________________________________________________________
GNuMIFlux::GNuMIFlux()
{
  this->Initialize();
}
//___________________________________________________________________________
GNuMIFlux::~GNuMIFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
bool GNuMIFlux::GenerateNext(void)
{
// Get next (unweighted) flux ntuple entry on the specified detector location
//
  if(!fDetLocIsSet) {
     LOG("Flux", pERROR)
       << "Specify a detector location before generating flux neutrinos";
     return false;	
  }
  if(fMaxWeight<=0) {
     LOG("Flux", pERROR)
       << "Run ScanForMaxWeight() before generating unweighted flux neurtinos";
     return false;	
  }


  // RWH
  LOG("Flux", pWARN)
    << "GenerateNext calling GenerateNext_weighted() ... and returning";
  return this->GenerateNext_weighted();
  // RWH


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
bool GNuMIFlux::GenerateNext_weighted(void)
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

  fNuFluxTree->GetEntry(fIEntry);
  fIEntry++;
  fCurrentEntry->pcodes = 0;  // fetched entry has geant codes
  fCurrentEntry->units  = 0;  // fetched entry has original units

  // Convert the current gnumi neutrino flavor mode into a neutrino pdg code
  // Also convert other particle codes in GNuMIFluxPassThroughInfo to PDG
  fCurrentEntry->ConvertPartCodes();
  fgPdgC = fCurrentEntry->ntype;

  // Check neutrino pdg against declared list of neutrino species declared
  // by the current instance of the NuMI neutrino flux driver.
  // No undeclared neutrino species will be accepted at this point as GENIE
  // has already been configured to handle the specified list.
  // Make sure that the appropriate list of flux neutrino species was set at
  // initialization via GNuMIFlux::SetFluxParticles(const PDGCodeList &)

  if( ! fPdgCList->ExistsInPDGCodeList(fgPdgC) ) {
     LOG("Flux", pWARN)
          << "Unknown decay mode or decay mode producing an undeclared neutrino species";
     LOG("Flux", pWARN)
          << "Declared list of neutrino species: " << *fPdgCList;
     LOG("Flux", pWARN)
       << "gnumi neutrino flavor: " << fCurrentEntry->ntype 
       << "(pcodes=" << fCurrentEntry->pcodes << ")" ;
     return false;	
  }

  // Update the curr neutrino weight and energy

  // Check current neutrino energy against the maximum flux neutrino energy declared
  // by the current instance of the NuMI neutrino flux driver.
  // No flux neutrino exceeding that maximum energy will be accepted at this point as
  // that maximum energy has already been used for normalizing the interaction probabilities.
  // Make sure that the appropriate maximum flux neutrino energy was set at
  // initialization via GNuMIFlux::SetMaxEnergy(double Ev)
  fWeight = fCurrentEntry->nimpwt;   // start with importance weight
  double Ev = 0;
  switch ( fUseFluxAtDetCenter ) {
  case -1:  // near detector
    fWeight *= fCurrentEntry->nwtnear;
    Ev       = fCurrentEntry->nenergyn;
    break;
  case +1:  // far detector
    fWeight *= fCurrentEntry->nwtfar;
    Ev       = fCurrentEntry->nenergyf;
    break;
  default:  // recalculate on x-y window
     LOG("Flux", pWARN)
          << "not yet calculating on x-y window";
     return false;
    fWeight = 0;
    Ev      = 0;
    break;
  }

  if(Ev > fMaxEv) {
     LOG("Flux", pWARN)
          << "Flux neutrino energy exceeds declared maximum neutrino energy";
     LOG("Flux", pWARN)
          << "Ev = " << Ev << "(> Ev{max} = " << fMaxEv << ")";
  }

  // Set the current flux neutrino 4-momentum
  fgP4.SetPxPyPzE (0., 0., Ev, Ev); // ---> random values for now

  // Set the current flux neutrino 4-position
  // (Read curr flux ntuple entry & trace back to requested started position z0)
  double cm2m = units::cm / units::m;
  fgX4.SetXYZT    (cm2m*1.,  cm2m*2.,  cm2m*3.,  0.);  // ---> random values for now

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pINFO)
	<< "Generated neutrino: "
	<< "\n pdg-code: " << fgPdgC
        << "\n p4: " << utils::print::P4AsShortString(&fgP4)
        << "\n x4: " << utils::print::X4AsString(&fgX4);
#endif

  // update the sum of weights & number of neutrinos
  fSumWeight += this->Weight();
  fNNeutrinos++;

  return true;
}
//___________________________________________________________________________
double GNuMIFlux::POT_curr(void)
{
// Compute current number of flux POTs

  if(!fNuFluxTree) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return 0;	
  }

  return 0;
}
//___________________________________________________________________________
void GNuMIFlux::LoadBeamSimData(string filename)
{
// Loads in a gnumi beam simulation root file (converted from hbook format)
// into the GNuMIFlux driver.

  LOG("Flux", pNOTICE)
        << "Loading gnumi flux tree from ROOT file: " << filename;

  bool is_accessible = ! (gSystem->AccessPathName( filename.c_str() ));
  if (!is_accessible) {
    LOG("Flux", pFATAL)
     << "The input gnumi flux file doesn't exist! Initialization failed!";
    exit(1);
  }

  fNuFluxFile = new TFile(filename.c_str(), "read");
  if(fNuFluxFile) {
      LOG("Flux", pINFO) << "Getting flux tree: " << fNuFluxTreeName;
      fNuFluxTree = (TTree*) fNuFluxFile->Get(fNuFluxTreeName.c_str());
      if(!fNuFluxTree) {
          LOG("Flux", pERROR)
             << "** Couldn't get flux tree: " << fNuFluxTreeName;
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

  fBr_run          = fNuFluxTree -> GetBranch ( "run"      );
  fBr_evtno        = fNuFluxTree -> GetBranch ( "evtno"    );
  fBr_Ndxdz        = fNuFluxTree -> GetBranch ( "Ndxdz"    );
  fBr_Ndydz        = fNuFluxTree -> GetBranch ( "Ndydz"    );
  fBr_Npz          = fNuFluxTree -> GetBranch ( "Npz"      );
  fBr_Nenergy      = fNuFluxTree -> GetBranch ( "Nenergy"  );
  fBr_NdxdzNea     = fNuFluxTree -> GetBranch ( "Ndxdznea" );  // was NdxdzNea
  fBr_NdydzNea     = fNuFluxTree -> GetBranch ( "Ndydznea" );  // was NdydzNea
  fBr_NenergyN     = fNuFluxTree -> GetBranch ( "Nenergyn" );  // was NenergyN **
  fBr_NWtNear      = fNuFluxTree -> GetBranch ( "Nwtnear"  );  // was NWtNear
  fBr_NdxdzFar     = fNuFluxTree -> GetBranch ( "Ndxdzfar" );  // was NdxdzFar
  fBr_NdydzFar     = fNuFluxTree -> GetBranch ( "Ndydzfar" );  // was NdydzFar
  fBr_NenergyF     = fNuFluxTree -> GetBranch ( "Nenergyf" );  // was NenergyF
  fBr_NWtFar       = fNuFluxTree -> GetBranch ( "Nwtfar"   );  // was NWtFar
  fBr_Norig        = fNuFluxTree -> GetBranch ( "Norig"    );
  fBr_Ndecay       = fNuFluxTree -> GetBranch ( "Ndecay"   );
  fBr_Ntype        = fNuFluxTree -> GetBranch ( "Ntype"    );
  fBr_Vx           = fNuFluxTree -> GetBranch ( "Vx"       );
  fBr_Vy           = fNuFluxTree -> GetBranch ( "Vy"       );
  fBr_Vz           = fNuFluxTree -> GetBranch ( "Vz"       );
  fBr_pdPx         = fNuFluxTree -> GetBranch ( "pdpx"     );  // was pdPx
  fBr_pdPy         = fNuFluxTree -> GetBranch ( "pdpy"     );  // was pdPy
  fBr_pdPz         = fNuFluxTree -> GetBranch ( "pdpz"     );  // was pdPz
  fBr_ppdxdz       = fNuFluxTree -> GetBranch ( "ppdxdz"   );
  fBr_ppdydz       = fNuFluxTree -> GetBranch ( "ppdydz"   );
  fBr_pppz         = fNuFluxTree -> GetBranch ( "pppz"     );
  fBr_ppenergy     = fNuFluxTree -> GetBranch ( "ppenergy" );
  fBr_ppmedium     = fNuFluxTree -> GetBranch ( "ppmedium" );
  fBr_ptype        = fNuFluxTree -> GetBranch ( "ptype"    );
  fBr_ppvx         = fNuFluxTree -> GetBranch ( "ppvx"     );
  fBr_ppvy         = fNuFluxTree -> GetBranch ( "ppvy"     );
  fBr_ppvz         = fNuFluxTree -> GetBranch ( "ppvz"     );
  fBr_muparpx      = fNuFluxTree -> GetBranch ( "muparpx"  );
  fBr_muparpy      = fNuFluxTree -> GetBranch ( "muparpy"  );
  fBr_muparpz      = fNuFluxTree -> GetBranch ( "muparpz"  );
  fBr_mupare       = fNuFluxTree -> GetBranch ( "mupare"   );
  fBr_Necm         = fNuFluxTree -> GetBranch ( "Necm"     );
  fBr_Nimpwt       = fNuFluxTree -> GetBranch ( "Nimpwt"   );
  fBr_xpoint       = fNuFluxTree -> GetBranch ( "xpoint"   );
  fBr_ypoint       = fNuFluxTree -> GetBranch ( "ypoint"   );
  fBr_zpoint       = fNuFluxTree -> GetBranch ( "zpoint"   );
  fBr_tvx          = fNuFluxTree -> GetBranch ( "tvx"      );
  fBr_tvy          = fNuFluxTree -> GetBranch ( "tvy"      );
  fBr_tvz          = fNuFluxTree -> GetBranch ( "tvz"      );
  fBr_tpx          = fNuFluxTree -> GetBranch ( "tpx"      );
  fBr_tpy          = fNuFluxTree -> GetBranch ( "tpy"      );
  fBr_tpz          = fNuFluxTree -> GetBranch ( "tpz"      );
  fBr_tptype       = fNuFluxTree -> GetBranch ( "tptype"   );
  fBr_tgen         = fNuFluxTree -> GetBranch ( "tgen"     );
  fBr_tgptype      = fNuFluxTree -> GetBranch ( "tgptype"  );
  fBr_tgppx        = fNuFluxTree -> GetBranch ( "tgppx"    );
  fBr_tgppy        = fNuFluxTree -> GetBranch ( "tgppy"    );
  fBr_tgppz        = fNuFluxTree -> GetBranch ( "tgppz"    );
  fBr_tprivx       = fNuFluxTree -> GetBranch ( "tprivx"   );
  fBr_tprivy       = fNuFluxTree -> GetBranch ( "tprivy"   );
  fBr_tprivz       = fNuFluxTree -> GetBranch ( "tprivz"   );
  fBr_beamx        = fNuFluxTree -> GetBranch ( "beamx"    );
  fBr_beamy        = fNuFluxTree -> GetBranch ( "beamy"    );
  fBr_beamz        = fNuFluxTree -> GetBranch ( "beamz"    );
  fBr_beampx       = fNuFluxTree -> GetBranch ( "beampx"   );
  fBr_beampy       = fNuFluxTree -> GetBranch ( "beampy"   );
  fBr_beampz       = fNuFluxTree -> GetBranch ( "beampz"   );

  // set the leaf addresses for the above branches
  fBr_run          -> SetAddress ( &fCurrentEntry->run      );
  fBr_evtno        -> SetAddress ( &fCurrentEntry->evtno    );
  fBr_Ndxdz        -> SetAddress ( &fCurrentEntry->ndxdz    );
  fBr_Ndydz        -> SetAddress ( &fCurrentEntry->ndydz    );
  fBr_Npz          -> SetAddress ( &fCurrentEntry->npz      );
  fBr_Nenergy      -> SetAddress ( &fCurrentEntry->nenergy  );
  fBr_NdxdzNea     -> SetAddress ( &fCurrentEntry->ndxdznea );
  fBr_NdydzNea     -> SetAddress ( &fCurrentEntry->ndydznea );
  fBr_NenergyN     -> SetAddress ( &fCurrentEntry->nenergyn );
  fBr_NWtNear      -> SetAddress ( &fCurrentEntry->nwtnear  );
  fBr_NdxdzFar     -> SetAddress ( &fCurrentEntry->ndxdzfar );
  fBr_NdydzFar     -> SetAddress ( &fCurrentEntry->ndydzfar );
  fBr_NenergyF     -> SetAddress ( &fCurrentEntry->nenergyf );
  fBr_NWtFar       -> SetAddress ( &fCurrentEntry->nwtfar   );
  fBr_Norig        -> SetAddress ( &fCurrentEntry->norig    );
  fBr_Ndecay       -> SetAddress ( &fCurrentEntry->ndecay   );
  fBr_Ntype        -> SetAddress ( &fCurrentEntry->ntype    );
  fBr_Vx           -> SetAddress ( &fCurrentEntry->vx       );
  fBr_Vy           -> SetAddress ( &fCurrentEntry->vy       );
  fBr_Vz           -> SetAddress ( &fCurrentEntry->vz       );
  fBr_pdPx         -> SetAddress ( &fCurrentEntry->pdpx     );
  fBr_pdPy         -> SetAddress ( &fCurrentEntry->pdpy     );
  fBr_pdPz         -> SetAddress ( &fCurrentEntry->pdpz     );
  fBr_ppdxdz       -> SetAddress ( &fCurrentEntry->ppdxdz   );
  fBr_ppdydz       -> SetAddress ( &fCurrentEntry->ppdydz   );
  fBr_pppz         -> SetAddress ( &fCurrentEntry->pppz     );
  fBr_ppenergy     -> SetAddress ( &fCurrentEntry->ppenergy );
  fBr_ppmedium     -> SetAddress ( &fCurrentEntry->ppmedium );
  fBr_ptype        -> SetAddress ( &fCurrentEntry->ptype    );
  fBr_ppvx         -> SetAddress ( &fCurrentEntry->ppvx     );
  fBr_ppvy         -> SetAddress ( &fCurrentEntry->ppvy     );
  fBr_ppvz         -> SetAddress ( &fCurrentEntry->ppvz     );
  fBr_muparpx      -> SetAddress ( &fCurrentEntry->muparpx  );
  fBr_muparpy      -> SetAddress ( &fCurrentEntry->muparpy  );
  fBr_muparpz      -> SetAddress ( &fCurrentEntry->muparpz  );
  fBr_mupare       -> SetAddress ( &fCurrentEntry->mupare   );
  fBr_Necm         -> SetAddress ( &fCurrentEntry->necm     );
  fBr_Nimpwt       -> SetAddress ( &fCurrentEntry->nimpwt   );
  fBr_xpoint       -> SetAddress ( &fCurrentEntry->xpoint   );
  fBr_ypoint       -> SetAddress ( &fCurrentEntry->ypoint   );
  fBr_zpoint       -> SetAddress ( &fCurrentEntry->zpoint   );
  fBr_tvx          -> SetAddress ( &fCurrentEntry->tvx      );
  fBr_tvy          -> SetAddress ( &fCurrentEntry->tvy      );
  fBr_tvz          -> SetAddress ( &fCurrentEntry->tvz      );
  fBr_tpx          -> SetAddress ( &fCurrentEntry->tpx      );
  fBr_tpy          -> SetAddress ( &fCurrentEntry->tpy      );
  fBr_tpz          -> SetAddress ( &fCurrentEntry->tpz      );
  fBr_tptype       -> SetAddress ( &fCurrentEntry->tptype   );
  fBr_tgen         -> SetAddress ( &fCurrentEntry->tgen     );
  fBr_tgptype      -> SetAddress ( &fCurrentEntry->tgptype  );
  fBr_tgppx        -> SetAddress ( &fCurrentEntry->tgppx    );
  fBr_tgppy        -> SetAddress ( &fCurrentEntry->tgppy    );
  fBr_tgppz        -> SetAddress ( &fCurrentEntry->tgppz    );
  fBr_tprivx       -> SetAddress ( &fCurrentEntry->tprivx   );
  fBr_tprivy       -> SetAddress ( &fCurrentEntry->tprivy   );
  fBr_tprivz       -> SetAddress ( &fCurrentEntry->tprivz   );
  fBr_beamx        -> SetAddress ( &fCurrentEntry->beamx    );
  fBr_beamy        -> SetAddress ( &fCurrentEntry->beamy    );
  fBr_beamz        -> SetAddress ( &fCurrentEntry->beamz    );
  fBr_beampx       -> SetAddress ( &fCurrentEntry->beampx   );
  fBr_beampy       -> SetAddress ( &fCurrentEntry->beampy   );
  fBr_beampz       -> SetAddress ( &fCurrentEntry->beampz   );

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle = 1;
}
//___________________________________________________________________________
void GNuMIFlux::ScanForMaxWeight(void)
{
  if(!fDetLocIsSet) {
     LOG("Flux", pERROR)
       << "Specify a detector location before scanning for max weight";
     return;	
  }

  // scan for the maximum weight
  if ( fUseFluxAtDetCenter > 0 ) {
    fNuFluxTree->Draw("Nimpwt*Nwtfar","","goff");
  }
  if ( fUseFluxAtDetCenter < 0 ) {
    fNuFluxTree->Draw("Nimpwt*Nwtnear","","goff");
  }
  Long64_t idx = TMath::LocMax(
                     fNuFluxTree->GetSelectedRows(),
                     fNuFluxTree->GetV1());
  fMaxWeight = fNuFluxTree->GetV1()[idx];
  LOG("Flux", pNOTICE) << "Maximum flux weight = " << fMaxWeight;
  if(fMaxWeight <=0 ) {
      LOG("Flux", pFATAL) << "Non-positive maximum flux weight!";
      exit(1);
  }
}
//___________________________________________________________________________
void GNuMIFlux::SetFluxParticles(const PDGCodeList & particles)
{
  if(!fPdgCList) {
     fPdgCList = new PDGCodeList;
  }
  fPdgCList->Copy(particles);

  LOG("Flux", pINFO)
    << "Declared list of neutrino species: " << *fPdgCList;
}
//___________________________________________________________________________
void GNuMIFlux::SetMaxEnergy(double Ev)
{
  fMaxEv = TMath::Max(0.,Ev);

  LOG("Flux", pINFO)
    << "Declared maximum flux neutrino energy: " << fMaxEv;
}
//___________________________________________________________________________
void GNuMIFlux::SetFilePOT(double pot)
{
// POTs in input flux file

  fFilePOT = pot;
}
//___________________________________________________________________________
void GNuMIFlux::SetUpstreamZ(double z0)
{
// The flux neutrino position (x,y) is given at the detector coord system
// at z=0. This method sets the preferred starting z position upstream of
// the upstream detector face. Each flux neutrino will be backtracked from
// z=0 to the input z0.

  fZ0 = z0;
}
//___________________________________________________________________________
void GNuMIFlux::SetNumOfCycles(int n)
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
void GNuMIFlux::SetTreeName(string name)
{
  fNuFluxTreeName = name;
}
//___________________________________________________________________________
void GNuMIFlux::UseFluxAtFarDetCenter(void)
{
  fUseFluxAtDetCenter  = +1;
  fDetLocIsSet         = true;
}
//___________________________________________________________________________
void GNuMIFlux::UseFluxAtNearDetCenter(void)
{
  fUseFluxAtDetCenter  = -1;
  fDetLocIsSet         = true;
}
//___________________________________________________________________________
void GNuMIFlux::PrintCurrent(void)
{
  LOG("Flux", pNOTICE) << "CurrentEntry:" << *fCurrentEntry;
}
//___________________________________________________________________________
void GNuMIFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GNuMIFlux driver";

  fMaxEv           = 0;
  fPdgCList        = new PDGCodeList;
  fCurrentEntry    = new GNuMIFluxPassThroughInfo;

  fNuFluxFile      = 0;
  fNuFluxTree      = 0;
  fNuFluxTreeName  = "";

  fNEntries        = 0;
  fIEntry          = 0;
  fMaxWeight       =-1;
  fFilePOT         = 0;
  fZ0              = 0;
  fNCycles         = 0;
  fICycle          = 0;
  fSumWeight       = 0;
  fNNeutrinos      = 0;
  fUseFluxAtDetCenter = 0;
  fDetLocIsSet        = false;

  this->SetDefaults();
  this->ResetCurrent();
}
//___________________________________________________________________________
void GNuMIFlux::SetDefaults(void)
{
// - Set default neutrino species list (nue, nuebar, numu, numubar) and
//   maximum energy (125 GeV).
//   These defaults can be overwritten by user calls (at the driver init) to
//   GNuMIlux::SetMaxEnergy(double Ev) and
//   GNuMIFlux::SetFluxParticles(const PDGCodeList & particles)
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

  this->SetFluxParticles (particles);
  this->SetMaxEnergy     (125./*GeV*/);  // was 200, but that would be wasteful
  this->SetFilePOT       (1E+21);
  this->SetUpstreamZ     (-5.0);
  this->SetNumOfCycles   (1);
  this->SetTreeName      ("h10");
}
//___________________________________________________________________________
void GNuMIFlux::ResetCurrent(void)
{
// reset running values of neutrino pdg-code, 4-position & 4-momentum
// and the input ntuple leaves

  fgPdgC  = 0;
  fWeight = 0;
  fgP4.SetPxPyPzE (0.,0.,0.,0.);
  fgX4.SetXYZT    (0.,0.,0.,0.);

  fCurrentEntry->Reset();
}
//___________________________________________________________________________
void GNuMIFlux::CleanUp(void)
{
  LOG("Flux", pNOTICE) << "Cleaning up...";

  if (fPdgCList)        delete fPdgCList;
  if (fCurrentEntry)    delete fCurrentEntry;

  if (fNuFluxFile) {
	fNuFluxFile->Close();
	delete fNuFluxFile;
  }
}

//___________________________________________________________________________
GNuMIFluxPassThroughInfo::GNuMIFluxPassThroughInfo()
  : TObject()
{ Reset(); }
//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::Reset()
{
  pcodes = -1;
  units  = -1;
  
  run      = -1;
  evtno    = -1;
  ndxdz    = 0.;
  ndydz    = 0.;
  npz      = 0.;
  nenergy  = 0.;
  ndxdznea = 0.;
  ndydznea = 0.;
  nenergyn = 0.;
  nwtnear  = 0.;
  ndxdzfar = 0.;
  ndydzfar = 0.;
  nenergyf = 0.;
  nwtfar   = 0.;
  norig    = 0;
  ndecay   = 0;
  ntype    = 0;
  vx       = 0.;
  vy       = 0.;
  vz       = 0.;
  pdpx     = 0.;
  pdpy     = 0.;
  pdpz     = 0.;
  ppdxdz   = 0.;
  ppdydz   = 0.;
  pppz     = 0.;
  ppenergy = 0.;
  ppmedium = 0 ;
  ptype    = 0 ;
  ppvx     = 0.;
  ppvy     = 0.;
  ppvz     = 0.;
  muparpx  = 0.;
  muparpy  = 0.;
  muparpz  = 0.;
  mupare   = 0.;
  tvx      = 0.;
  tvy      = 0.;
  tvz      = 0.;
  tpx      = 0.;
  tpy      = 0.;
  tpz      = 0.;
  tptype   = 0 ;
  tgen     = 0 ;
  tgptype  = 0 ;
  tgppx    = 0.;
  tgppy    = 0.;
  tgppz    = 0.;
  tprivx   = 0.;
  tprivy   = 0.;
  tprivz   = 0.;
  beamx    = 0.;
  beamy    = 0.;
  beamz    = 0.;
  beampx   = 0.;
  beampy   = 0.;
  beampz   = 0.;

}

//___________________________________________________________________________
/*******

ALLOW DEFAULT COPY CONSTRUCTOR
GNuMIFluxPassThroughInfo::GNuMIFluxPassThroughInfo(
                        const GNuMIFluxPassThroughInfo & info) :
TObject(),
pcodes   ( info.pcodes   ),
units    ( info.units    ),
run      ( info.run      ),
evtno    ( info.evtno    ),
ndxdz    ( info.ndxdz    ),
ndydz    ( info.ndydz    ),
npz      ( info.npz      ),
nenergy  ( info.nenergy  ),
ndxdznea ( info.ndxdznea ),
ndydznea ( info.ndydznea ),
nenergyn ( info.nenergyn ),
nwtnear  ( info.nwtnear  ),
ndxdzfar ( info.ndxdzfar ),
ndydzfar ( info.ndydzfar ),
nenergyf ( info.nenergyf ),
nwtfar   ( info.nwtfar   ),
norig    ( info.norig    ),
ndecay   ( info.ndecay   ),
ntype    ( info.ntype    ),
vx       ( info.vx       ),
vy       ( info.vy       ),
vz       ( info.vz       ),
pdPx     ( info.pdPx     ),
pdPy     ( info.pdPy     ),
pdPz     ( info.pdPz     ),
ppdxdz   ( info.ppdxdz   ),
ppdydz   ( info.ppdydz   ),
pppz     ( info.pppz     ),
ppenergy ( info.ppenergy ),
ppmedium ( info.ppmedium ),
ptype    ( info.ptype    ),
ppvx     ( info.ppvx     ),
ppvy     ( info.ppvy     ),
ppvz     ( info.ppvz     ),
muparpx  ( info.muparpx  ),
muparpy  ( info.muparpy  ),
muparpz  ( info.muparpz  ),
mupare   ( info.mupare   ),
necm     ( info.necm     ),
nimpwt   ( info.nimpwt   ),
xpoint   ( info.xpoint   ),
ypoint   ( info.ypoint   ),
zpoint   ( info.zpoint   ),
tvx      ( info.tvx      ),
tvy      ( info.tvy      ),
tvz      ( info.tvz      ),
tpx      ( info.tpx      ),
tpy      ( info.tpy      ),
tpz      ( info.tpz      ),
tptype   ( info.tptype   ),
tgen     ( info.tgen     ),
tgptype  ( info.tgptype  ),
tgppx    ( info.tgppx    ),
tgppy    ( info.tgppy    ),
tgppz    ( info.tgppz    ),
tprivx   ( info.tprivx   ),
tprivy   ( info.tprivy   ),
tprivz   ( info.tprivz   ),
beamx    ( info.beamx    ),
beamy    ( info.beamy    ),
beamz    ( info.beamz    ),
beampx   ( info.beampx   ),
beampy   ( info.beampy   ),
beampz   ( info.beampz   )
{

}
*************/

//___________________________________________________________________________
void GNuMIFluxPassThroughInfo::ConvertPartCodes()
{
  if ( pcodes == 0 ) {
    pcodes = 1;  // flag that conversion has been made
    switch ( ntype ) {
    case 56: ntype = kPdgNuMu;     break;
    case 55: ntype = kPdgAntiNuMu; break;
    case 53: ntype = kPdgNuE;      break;
    case 52: ntype = kPdgAntiNuE;  break;
    default:
      LOG("Flux", pNOTICE)
        << "ConvertPartCodes saw ntype " << ntype << " -- unknown ";
    }
    if ( ptype   != 0 ) ptype   = pdg::GeantToPdg(ptype);
    if ( tptype  != 0 ) tptype  = pdg::GeantToPdg(tptype);
    if ( tgptype != 0 ) tgptype = pdg::GeantToPdg(tgptype);
  } else if ( pcodes != 1 ) {
    // flag as unknown state ...
    LOG("Flux", pNOTICE)
      << "ConvertPartCodes saw pcodes flag as " << pcodes;
  }

}

//___________________________________________________________________________
namespace genie {
namespace flux  {
  ostream & operator << (
    ostream & stream, const genie::flux::GNuMIFluxPassThroughInfo & info)
    {
      // stream << "\n ndecay   = " << info.ndecay << endl;
      stream << "\n GNuMIFlux run " << info.run << " evtno " << info.evtno
             << " pcodes " << info.pcodes << " units " << info.units << ")"
             << "\n random dk: dx/dz " << info.ndxdz 
             << " dy/dz " << info.ndydz
             << " pz " <<  info.npz << " E " << info.nenergy
             << "\n near00 dk: dx/dz " << info.ndxdznea 
             << " dy/dz " << info.ndydznea
             << " E " <<  info.nenergyn << " wgt " << info.nwtnear
             << "\n far00  dk: dx/dz " << info.ndxdzfar
             << " dy/dz " << info.ndydzfar
             << " E " <<  info.nenergyf << " wgt " << info.nwtfar
             << "\n norg " << info.norig << " ndecay " << info.ndecay
             << " ntype " << info.ntype
             << "\n had vtx " << info.vx << " " << info.vy << " " << info.vz
             << "\n parent p3 @ dk " << info.pdpx << " " << info.pdpy << " " << info.pdpz
             << "\n parent prod: dx/dz " << info.ppdxdz 
             << " dy/dz " << info.ppdydz
             << " pz " << info.pppz << " E " << info.ppenergy
             << "\n ppmedium " << info.ppmedium << " ptype " << info.ptype
             << " ppvtx " << info.ppvx << " " << info.ppvy << " " << info.ppvz
             << "\n mu parent p4 " << info.muparpx << " " << info.muparpy
             << " " << info.muparpz << " " << info.mupare
             << "\n necm " << info.necm << " nimpwt " << info.nimpwt
             << "\n point x,y,z " << info.xpoint << " " << info.ypoint
             << " " << info.zpoint
             << "\n tv x,y,z " << info.tvx << " " << info.tvy << " " << info.tvz
             << "\n tptype " << info.tptype << " tgen " << info.tgen
             << " tgptype " << info.tgptype 
             << "\n tgp px,py,pz " << info.tgppx << " " << info.tgppy
             << " " << info.tgppz
             << "\n tpriv x,y,z " << info.tprivx << " " << info.tprivy
             << " " << info.tprivz
             << "\n beam x,y,z " << info.beamx << " " << info.beamy
             << " " << info.beamz
             << "\n beam px,py,pz " << info.beampx << " " << info.beampy
             << " " << info.beampz
        ;
  /*
  //std::cout << "GNuMIFlux::PrintCurrent ....." << std::endl;
  //LOG("Flux", pINFO) 
  LOG("Flux", pNOTICE) 
    << "Current Leaf Values: "
    << " run " << fLf_run << " evtno " << fLf_evtno << "\n"
    << " NenergyN " << fLf_NenergyN << " NWtNear " << fLf_NWtNear
    << " NenergyF " << fLf_NenergyF << " NWtFar  " << fLf_NWtFar << "\n"
    << " Norig " << fLf_Norig << " Ndecay " << fLf_Ndecay << " Ntype " << fLf_Ntype << "\n"
    << " Vxyz " << fLf_Vx << " " << fLf_Vy << " " << fLf_Vz 
    << " pdPxyz " << fLf_pdPx << " " << fLf_pdPy << " " << fLf_pdPz << "\n"
    << " pp dxdz " << fLf_ppdxdz << " dydz " << fLf_ppdydz << " pz " << fLf_pppz << "\n"
    << " pp energy " << fLf_ppenergy << " medium " << fLf_ppmedium 
    << " ptype " << fLf_ptype 
    << " ppvxyz " << fLf_ppvx << " " << fLf_ppvy << " " << fLf_ppvz << "\n"
    << " muparpxyz " << fLf_muparpx << " " << fLf_muparpy << " " << fLf_muparpz
    << " mupare " << fLf_mupare << "\n"
    << ;
  */

     return stream;
  }
}//flux
}//genie
//___________________________________________________________________________

