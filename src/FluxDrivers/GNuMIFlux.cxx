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

  // Convert the current gnumi neutrino flavor mode into a neutrino pdg code

  if ( fLf_Ntype == 56 ) fgPdgC = kPdgNuMu;
  if ( fLf_Ntype == 55 ) fgPdgC = kPdgAntiNuMu;
  if ( fLf_Ntype == 53 ) fgPdgC = kPdgNuE;
  if ( fLf_Ntype == 52 ) fgPdgC = kPdgAntiNuE;

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
          << "gnumi neutrino flavor: " << fLf_Ntype;
     return false;	
  }

  // Update the curr neutrino weight and energy

  // Check current neutrino energy against the maximum flux neutrino energy declared
  // by the current instance of the NuMI neutrino flux driver.
  // No flux neutrino exceeding that maximum energy will be accepted at this point as
  // that maximum energy has already been used for normalizing the interaction probabilities.
  // Make sure that the appropriate maximum flux neutrino energy was set at
  // initialization via GNuMIFlux::SetMaxEnergy(double Ev)
  fWeight = fLf_Nimpwt;   // start with importance weight
  double Ev = 0;
  switch ( fUseFluxAtDetCenter ) {
  case -1:  // near detector
    fWeight *= fLf_NWtNear;
    Ev       = fLf_NenergyN;
    break;
  case +1:  // far detector
    fWeight *= fLf_NWtFar;
    Ev       = fLf_NenergyF;
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

  // Update pass-through info (= info on the flux neutrino parent particle
  // that may be stored at an extra branch of the output event tree -alongside
  // with the generated event branch- for use further upstream in the
  // analysis chain -eg for beam reweighting etc-)
  fPassThroughInfo -> ndecay   = fLf_Ndecay;
  fPassThroughInfo -> vx       = fLf_Vx;
  fPassThroughInfo -> vy       = fLf_Vy;
  fPassThroughInfo -> vz       = fLf_Vz;
  fPassThroughInfo -> pdPx     = fLf_pdPx;
  fPassThroughInfo -> pdPy     = fLf_pdPy;
  fPassThroughInfo -> pdPz     = fLf_pdPz;
  fPassThroughInfo -> ppdxdz   = fLf_ppdxdz;
  fPassThroughInfo -> ppdydz   = fLf_ppdydz;
  fPassThroughInfo -> pppz     = fLf_pppz;
  fPassThroughInfo -> ppenergy = fLf_ppenergy;
  fPassThroughInfo -> ppmedium = fLf_ppmedium;
  fPassThroughInfo -> ptype    = pdg::GeantToPdg( fLf_ptype );
  fPassThroughInfo -> ppvx     = fLf_ppvx;
  fPassThroughInfo -> ppvy     = fLf_ppvy;
  fPassThroughInfo -> ppvz     = fLf_ppvz;
  fPassThroughInfo -> muparpx  = fLf_muparpx;
  fPassThroughInfo -> muparpy  = fLf_muparpy;
  fPassThroughInfo -> muparpz  = fLf_muparpz;
  fPassThroughInfo -> mupare   = fLf_mupare;
  fPassThroughInfo -> tvx      = fLf_tvx;
  fPassThroughInfo -> tvy      = fLf_tvy;
  fPassThroughInfo -> tvz      = fLf_tvz;
  fPassThroughInfo -> tpx      = fLf_tpx;
  fPassThroughInfo -> tpy      = fLf_tpy;
  fPassThroughInfo -> tpz      = fLf_tpz;
  fPassThroughInfo -> tptype   = pdg::GeantToPdg( fLf_tptype );
  fPassThroughInfo -> tgen     = fLf_tgen;
  fPassThroughInfo -> tgptype  = pdg::GeantToPdg( fLf_tgptype );
  fPassThroughInfo -> tgppx    = fLf_tgppx;
  fPassThroughInfo -> tgppy    = fLf_tgppy;
  fPassThroughInfo -> tgppz    = fLf_tgppz;
  fPassThroughInfo -> tprivx   = fLf_tprivx;
  fPassThroughInfo -> tprivy   = fLf_tprivy;
  fPassThroughInfo -> tprivz   = fLf_tprivz;
  fPassThroughInfo -> beamx    = fLf_beamx;
  fPassThroughInfo -> beamy    = fLf_beamy;
  fPassThroughInfo -> beamz    = fLf_beamz;
  fPassThroughInfo -> beampx   = fLf_beampx;
  fPassThroughInfo -> beampy   = fLf_beampy;
  fPassThroughInfo -> beampz   = fLf_beampz;

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
  fBr_run          -> SetAddress ( &fLf_run      );
  fBr_evtno        -> SetAddress ( &fLf_evtno    );
  fBr_Ndxdz        -> SetAddress ( &fLf_Ndxdz    );
  fBr_Ndydz        -> SetAddress ( &fLf_Ndydz    );
  fBr_Npz          -> SetAddress ( &fLf_Npz      );
  fBr_Nenergy      -> SetAddress ( &fLf_Nenergy  );
  fBr_NdxdzNea     -> SetAddress ( &fLf_NdxdzNea );
  fBr_NdydzNea     -> SetAddress ( &fLf_NdydzNea );
  fBr_NenergyN     -> SetAddress ( &fLf_NenergyN );
  fBr_NWtNear      -> SetAddress ( &fLf_NWtNear  );
  fBr_NdxdzFar     -> SetAddress ( &fLf_NdxdzFar );
  fBr_NdydzFar     -> SetAddress ( &fLf_NdydzFar );
  fBr_NenergyF     -> SetAddress ( &fLf_NenergyF );
  fBr_NWtFar       -> SetAddress ( &fLf_NWtFar   );
  fBr_Norig        -> SetAddress ( &fLf_Norig    );
  fBr_Ndecay       -> SetAddress ( &fLf_Ndecay   );
  fBr_Ntype        -> SetAddress ( &fLf_Ntype    );
  fBr_Vx           -> SetAddress ( &fLf_Vx       );
  fBr_Vy           -> SetAddress ( &fLf_Vy       );
  fBr_Vz           -> SetAddress ( &fLf_Vz       );
  fBr_pdPx         -> SetAddress ( &fLf_pdPx     );
  fBr_pdPy         -> SetAddress ( &fLf_pdPy     );
  fBr_pdPz         -> SetAddress ( &fLf_pdPz     );
  fBr_ppdxdz       -> SetAddress ( &fLf_ppdxdz   );
  fBr_ppdydz       -> SetAddress ( &fLf_ppdydz   );
  fBr_pppz         -> SetAddress ( &fLf_pppz     );
  fBr_ppenergy     -> SetAddress ( &fLf_ppenergy );
  fBr_ppmedium     -> SetAddress ( &fLf_ppmedium );
  fBr_ptype        -> SetAddress ( &fLf_ptype    );
  fBr_ppvx         -> SetAddress ( &fLf_ppvx     );
  fBr_ppvy         -> SetAddress ( &fLf_ppvy     );
  fBr_ppvz         -> SetAddress ( &fLf_ppvz     );
  fBr_muparpx      -> SetAddress ( &fLf_muparpx  );
  fBr_muparpy      -> SetAddress ( &fLf_muparpy  );
  fBr_muparpz      -> SetAddress ( &fLf_muparpz  );
  fBr_mupare       -> SetAddress ( &fLf_mupare   );
  fBr_Necm         -> SetAddress ( &fLf_Necm     );
  fBr_Nimpwt       -> SetAddress ( &fLf_Nimpwt   );
  fBr_xpoint       -> SetAddress ( &fLf_xpoint   );
  fBr_ypoint       -> SetAddress ( &fLf_ypoint   );
  fBr_zpoint       -> SetAddress ( &fLf_zpoint   );
  fBr_tvx          -> SetAddress ( &fLf_tvx      );
  fBr_tvy          -> SetAddress ( &fLf_tvy      );
  fBr_tvz          -> SetAddress ( &fLf_tvz      );
  fBr_tpx          -> SetAddress ( &fLf_tpx      );
  fBr_tpy          -> SetAddress ( &fLf_tpy      );
  fBr_tpz          -> SetAddress ( &fLf_tpz      );
  fBr_tptype       -> SetAddress ( &fLf_tptype   );
  fBr_tgen         -> SetAddress ( &fLf_tgen     );
  fBr_tgptype      -> SetAddress ( &fLf_tgptype  );
  fBr_tgppx        -> SetAddress ( &fLf_tgppx    );
  fBr_tgppy        -> SetAddress ( &fLf_tgppy    );
  fBr_tgppz        -> SetAddress ( &fLf_tgppz    );
  fBr_tprivx       -> SetAddress ( &fLf_tprivx   );
  fBr_tprivy       -> SetAddress ( &fLf_tprivy   );
  fBr_tprivz       -> SetAddress ( &fLf_tprivz   );
  fBr_beamx        -> SetAddress ( &fLf_beamx    );
  fBr_beamy        -> SetAddress ( &fLf_beamy    );
  fBr_beamz        -> SetAddress ( &fLf_beamz    );
  fBr_beampx       -> SetAddress ( &fLf_beampx   );
  fBr_beampy       -> SetAddress ( &fLf_beampy   );
  fBr_beampz       -> SetAddress ( &fLf_beampz   );

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
  //std::cout << "GNuMIFlux::PrintCurrent ....." << std::endl;
  //LOG("Flux", pINFO) 
  LOG("Flux", pNOTICE) 
    << "Current Leaf Values: "
    << " run " << fLf_run << " evtno " << fLf_evtno << "\n"
    << " NenergyN " << fLf_NenergyN << " NWtNear " << fLf_NWtNear
    << " NenergyF " << fLf_NenergyF << " NWtFar  " << fLf_NWtFar << "\n"
    << " Ntype " << fLf_Ntype;
}
//___________________________________________________________________________
void GNuMIFlux::Initialize(void)
{
  LOG("Flux", pNOTICE) << "Initializing GNuMIFlux driver";

  fMaxEv           = 0;
  fPdgCList        = new PDGCodeList;
  fPassThroughInfo = new GNuMIFluxPassThroughInfo;

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

  fBr_run          = 0;
  fBr_evtno        = 0;
  fBr_Ndxdz        = 0;
  fBr_Ndydz        = 0;
  fBr_Npz          = 0;
  fBr_Nenergy      = 0;
  fBr_NdxdzNea     = 0;
  fBr_NdydzNea     = 0;
  fBr_NenergyN     = 0;
  fBr_NWtNear      = 0;
  fBr_NdxdzFar     = 0;
  fBr_NdydzFar     = 0;
  fBr_NenergyF     = 0;
  fBr_NWtFar       = 0;
  fBr_Norig        = 0;
  fBr_Ndecay       = 0;
  fBr_Ntype        = 0;
  fBr_Vx           = 0;
  fBr_Vy           = 0;
  fBr_Vz           = 0;
  fBr_pdPx         = 0;
  fBr_pdPy         = 0;
  fBr_pdPz         = 0;
  fBr_ppdxdz       = 0;
  fBr_ppdydz       = 0;
  fBr_pppz         = 0;
  fBr_ppenergy     = 0;
  fBr_ppmedium     = 0;
  fBr_ptype        = 0;
  fBr_ppvx         = 0;
  fBr_ppvy         = 0;
  fBr_ppvz         = 0;
  fBr_muparpx      = 0;
  fBr_muparpy      = 0;
  fBr_muparpz      = 0;
  fBr_mupare       = 0;
  fBr_Necm         = 0;
  fBr_Nimpwt       = 0;
  fBr_xpoint       = 0;
  fBr_ypoint       = 0;
  fBr_zpoint       = 0;
  fBr_tvx          = 0;
  fBr_tvy          = 0;
  fBr_tvz          = 0;
  fBr_tpx          = 0;
  fBr_tpy          = 0;
  fBr_tpz          = 0;
  fBr_tptype       = 0;
  fBr_tgen         = 0;
  fBr_tgptype      = 0;
  fBr_tgppx        = 0;
  fBr_tgppy        = 0;
  fBr_tgppz        = 0;
  fBr_tprivx       = 0;
  fBr_tprivy       = 0;
  fBr_tprivz       = 0;
  fBr_beamx        = 0;
  fBr_beamy        = 0;
  fBr_beamz        = 0;
  fBr_beampx       = 0;
  fBr_beampy       = 0;
  fBr_beampz       = 0;

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

  fLf_run          = 0;
  fLf_evtno        = 0;
  fLf_Ndxdz        = 0;
  fLf_Ndydz        = 0;
  fLf_Npz          = 0;
  fLf_Nenergy      = 0;
  fLf_NdxdzNea     = 0;
  fLf_NdydzNea     = 0;
  fLf_NenergyN     = 0;
  fLf_NWtNear      = 0;
  fLf_NdxdzFar     = 0;
  fLf_NdydzFar     = 0;
  fLf_NenergyF     = 0;
  fLf_NWtFar       = 0;
  fLf_Norig        = 0;
  fLf_Ndecay       = 0;
  fLf_Ntype        = 0;
  fLf_Vx           = 0;
  fLf_Vy           = 0;
  fLf_Vz           = 0;
  fLf_pdPx         = 0;
  fLf_pdPy         = 0;
  fLf_pdPz         = 0;
  fLf_ppdxdz       = 0;
  fLf_ppdydz       = 0;
  fLf_pppz         = 0;
  fLf_ppenergy     = 0;
  fLf_ppmedium     = 0;
  fLf_ptype        = 0;
  fLf_ppvx         = 0;
  fLf_ppvy         = 0;
  fLf_ppvz         = 0;
  fLf_muparpx      = 0;
  fLf_muparpy      = 0;
  fLf_muparpz      = 0;
  fLf_mupare       = 0;
  fLf_Necm         = 0;
  fLf_Nimpwt       = 0;
  fLf_xpoint       = 0;
  fLf_ypoint       = 0;
  fLf_zpoint       = 0;
  fLf_tvx          = 0;
  fLf_tvy          = 0;
  fLf_tvz          = 0;
  fLf_tpx          = 0;
  fLf_tpy          = 0;
  fLf_tpz          = 0;
  fLf_tptype       = 0;
  fLf_tgen         = 0;
  fLf_tgptype      = 0;
  fLf_tgppx        = 0;
  fLf_tgppy        = 0;
  fLf_tgppz        = 0;
  fLf_tprivx       = 0;
  fLf_tprivy       = 0;
  fLf_tprivz       = 0;
  fLf_beamx        = 0;
  fLf_beamy        = 0;
  fLf_beamz        = 0;
  fLf_beampx       = 0;
  fLf_beampy       = 0;
  fLf_beampz       = 0;
}
//___________________________________________________________________________
void GNuMIFlux::CleanUp(void)
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
GNuMIFluxPassThroughInfo::GNuMIFluxPassThroughInfo() :
TObject(),
ndecay   (0 ),
vx       (0.),
vy       (0.),
vz       (0.),
pdPx     (0.),
pdPy     (0.),
pdPz     (0.),
ppdxdz   (0.),
ppdydz   (0.),
pppz     (0.),
ppenergy (0.),
ppmedium (0 ),
ptype    (0 ),
ppvx     (0.),
ppvy     (0.),
ppvz     (0.),
muparpx  (0.),
muparpy  (0.),
muparpz  (0.),
mupare   (0.),
tvx      (0.),
tvy      (0.),
tvz      (0.),
tpx      (0.),
tpy      (0.),
tpz      (0.),
tptype   (0 ),
tgen     (0 ),
tgptype  (0 ),
tgppx    (0.),
tgppy    (0.),
tgppz    (0.),
tprivx   (0.),
tprivy   (0.),
tprivz   (0.),
beamx    (0.),
beamy    (0.),
beamz    (0.),
beampx   (0.),
beampy   (0.),
beampz   (0.)
{

}
//___________________________________________________________________________
GNuMIFluxPassThroughInfo::GNuMIFluxPassThroughInfo(
                        const GNuMIFluxPassThroughInfo & info) :
TObject(),
ndecay   ( info.ndecay   ),
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
//___________________________________________________________________________
namespace genie {
namespace flux  {
  ostream & operator << (
    ostream & stream, const genie::flux::GNuMIFluxPassThroughInfo & info)
    {
      stream << "\n ndecay   = " << info.ndecay
             << endl;

     return stream;
  }
}//flux
}//genie
//___________________________________________________________________________

