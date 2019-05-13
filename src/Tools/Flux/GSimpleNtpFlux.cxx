//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>
         Fermi National Accelerator Laboratory

 For the class documentation see the corresponding header file.

 @ Feb 22, 2011 - JD
   Implemented dummy versions of the new GFluxI::Clear, GFluxI::Index and 
   GFluxI::GenerateWeighted methods needed for pre-generation of flux
   interaction probabilities in GMCJDriver.
 @ Mar 14, 2014 - TD
   Prevent an infinite loop in GenerateNext() when the flux driver has not been
   properly configured by exiting within GenerateNext_weighted().

*/
//____________________________________________________________________________

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cassert>
#include <limits.h>
#include <algorithm>

#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/GBuild.h"

#include "Tools/Flux/GSimpleNtpFlux.h"

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/UnitUtils.h"

using std::endl;

#include <vector>
#include <algorithm>
#include <iomanip>
#include "TRegexp.h"
#include "TString.h"

#include "Tools/Flux/GFluxDriverFactory.h"
FLUXDRIVERREG4(genie,flux,GSimpleNtpFlux,genie::flux::GSimpleNtpFlux)

//#define __GENIE_LOW_LEVEL_MESG_ENABLED__
// next line won't work for NOvA: ROOT's Error() != DefaultErrorHandler
//#define USE_INDEX_FOR_META 

using namespace genie;
using namespace genie::flux;

ClassImp(GSimpleNtpFlux)
ClassImp(GSimpleNtpEntry)
ClassImp(GSimpleNtpNuMI)
ClassImp(GSimpleNtpAux)
ClassImp(GSimpleNtpMeta)

// static storage
UInt_t genie::flux::GSimpleNtpMeta::mxfileprint = UINT_MAX;

//____________________________________________________________________________
GSimpleNtpFlux::GSimpleNtpFlux() :
  GFluxExposureI(genie::flux::kPOTs)
{
  this->Initialize();
}
//___________________________________________________________________________
GSimpleNtpFlux::~GSimpleNtpFlux()
{
  this->CleanUp();
}
//___________________________________________________________________________
double GSimpleNtpFlux::GetTotalExposure() const
{
  // complete the GFluxExposureI interface
  return UsedPOTs();
}
//___________________________________________________________________________
long int  GSimpleNtpFlux::NFluxNeutrinos(void) const 
{ 
  ///< number of flux neutrinos looped so far
  return fNNeutrinos; 
} 
//___________________________________________________________________________
bool GSimpleNtpFlux::GenerateNext(void)
{
// Get next (unweighted) flux ntuple entry on the specified detector location
//
  RandomGen* rnd = RandomGen::Instance();
  while ( true ) {
     // Check for end of flux ntuple
     bool end = this->End();
     if ( end ) {
       LOG("Flux", pNOTICE) << "GenerateNext signaled End() ";
       return false;
     }

     // Get next weighted flux ntuple entry
     bool nextok = this->GenerateNext_weighted();
     if ( fGenWeighted ) return nextok;
     if ( ! nextok ) continue;
     if ( fAlreadyUnwgt ) return true;

     /* RWH - debug purposes
     if ( fNCycles == 0 ) {
       LOG("Flux", pNOTICE)
          << "Got flux entry: " << fIEntry
          << " - Cycle: "<< fICycle << "/ infinite";
     } else {
       LOG("Flux", pNOTICE)
          << "Got flux entry: "<< fIEntry
          << " - Cycle: "<< fICycle << "/"<< fNCycles;
     }
     */

     // Get fractional weight & decide whether to accept curr flux neutrino
     double f = this->Weight() / fMaxWeight;
     //LOG("Flux", pNOTICE)
     //   << "Curr flux neutrino fractional weight = " << f;
     if (f > 1.) {
       fMaxWeight = this->Weight() * 1.01; // bump the weight
       LOG("Flux", pERROR)
         << "** Fractional weight = " << f 
         << " > 1 !! Bump fMaxWeight estimate to " << fMaxWeight
         << GetCurrentEntry();
       std::cout << std::flush;
     }
     double r = (f < 1.) ? rnd->RndFlux().Rndm() : 0;
     bool accept = ( r < f );
     if ( accept ) {

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("Flux", pNOTICE)
         << "Generated beam neutrino: "
         << "\n pdg-code: " << fCurEntry->pdg
         << "\n p4: " << utils::print::P4AsShortString(&(fP4))
         << "\n x4: " << utils::print::X4AsString(&(fX4));
#endif

       fWeight = 1.;
       return true;
     }

     //LOG("Flux", pNOTICE)
     //  << "** Rejecting current flux neutrino based on the flux weight only";
  }
  return false;
}
//___________________________________________________________________________
bool GSimpleNtpFlux::GenerateNext_weighted(void)
{
// Get next (weighted) flux ntuple entry 
//

  // Check whether a flux ntuple has been loaded
  if ( ! fNuFluxTree ) {
     LOG("Flux", pFATAL)
          << "The flux driver has not been properly configured";
     // return false; // don't do this - creates an infinite loop!
     exit(1); 	
  }

  // Reuse an entry?
  //std::cout << " ***** iuse " << fIUse << " nuse " << fNUse
  //          << " ientry " << fIEntry << " nentry " << fNEntries
  //          << " icycle " << fICycle << " ncycle " << fNCycles << std::endl;
  if ( fIUse < fNUse && fIEntry >= 0 ) {
    // Reuse this entry
    fIUse++;
  } else {
    // Reset previously generated neutrino code / 4-p / 4-x
    this->ResetCurrent();
    // Move on, read next flux ntuple entry
    ++fIEntry;
    ++fNEntriesUsed;  // count total # used
    if ( fIEntry >= fNEntries ) {
      // Ran out of entries @ the current cycle of this flux file
      // Check whether more (or infinite) number of cycles is requested
      if (fICycle < fNCycles || fNCycles == 0 ) {
        fICycle++;
        fIEntry=0;
      } else {
        LOG("Flux", pWARN)
          << "No more entries in input flux neutrino ntuple, cycle "
          << fICycle << " of " << fNCycles;
        fEnd = true;
        //assert(0);
        return false;	
      }
    }
    
    int nbytes = fNuFluxTree->GetEntry(fIEntry);
    UInt_t metakey = fCurEntry->metakey;
    if ( fAllFilesMeta && ( fCurMeta->metakey != metakey ) ) {
      UInt_t oldkey = fCurMeta->metakey;
#ifdef USE_INDEX_FOR_META
      int nbmeta = fNuMetaTree->GetEntryWithIndex(metakey);
#else
      // unordered indices makes ROOT call Error() which might,
      // if not DefaultErrorHandler, be fatal.
      // so find the right one by a simple linear search.
      // not a large burden since it only happens infrequently and
      // the list is normally quite short.
      int nmeta = fNuMetaTree->GetEntries();
      int nbmeta = 0;
      for (int imeta = 0; imeta < nmeta; ++imeta ) {
        nbmeta = fNuMetaTree->GetEntry(imeta);
        if ( fCurMeta->metakey == metakey ) break;
      }
      // next condition should never happen
      if ( fCurMeta->metakey != metakey ) {
        fCurMeta = 0; // didn't find it!?
        LOG("Flux",pERROR) << "Failed to find right metakey=" << metakey
                           << " (was " << oldkey << ") out of " << nmeta 
                           << " entries";
      }
#endif
      LOG("Flux",pDEBUG) << "Get meta " << metakey 
                         << " (was " << oldkey << ") "
                         << fCurMeta->metakey 
                         << " nb " << nbytes << " " << nbmeta;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("Flux",pDEBUG) << "Get meta " << *fCurMeta; 
#endif
    }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    Int_t ifile = fNuFluxTree->GetFileNumber();
    LOG("Flux",pDEBUG)
      << "got " << fNNeutrinos << " nu, using fIEntry " << fIEntry 
      << " ifile " << ifile << " nbytes " << nbytes
      << *fCurEntry << *fCurMeta;
#endif

    fIUse = 1; 

    // here we might want to do flavor oscillations or simple mappings
    //RWH modify:  fPdgC = fCurEntry->pdg;

    // Check neutrino pdg against declared list of neutrino species declared
    // by the current instance of the NuMI neutrino flux driver.
    // No undeclared neutrino species will be accepted at this point as GENIE
    // has already been configured to handle the specified list.
    // Make sure that the appropriate list of flux neutrino species was set at
    // initialization via GSimpleNtpFlux::SetFluxParticles(const PDGCodeList &)

    // update the # POTs, sum of weights & number of neutrinos 
    // do this HERE (before rejecting flavors that users might be weeding out)
    // in order to keep the POT accounting correct.  This allows one to get
    // the right normalization for generating only events from the intrinsic
    // nu_e entries.
    fAccumPOTs += fEffPOTsPerNu;
    fSumWeight += this->Weight();
    fNNeutrinos++;

    if ( ! fPdgCList->ExistsInPDGCodeList(fCurEntry->pdg) ) {
      /// user might modify list via SetFluxParticles() in order to reject certain
      /// flavors, even if they're found in the file.  So don't make a big fuss.
      /// Spit out a single message and then stop reporting that flavor as problematic.
      int badpdg = fCurEntry->pdg;
      if ( ! fPdgCListRej->ExistsInPDGCodeList(badpdg) ) {
        fPdgCListRej->push_back(badpdg);
        LOG("Flux", pWARN)
          << "Encountered neutrino specie (" << badpdg 
          << ") that wasn't in SetFluxParticles() list, "
          << "\nDeclared list of neutrino species: " << *fPdgCList;
      }
      return false;	
    }

  }

  // Update the curr neutrino p4/x4 lorentz vector
  double Ev =  fCurEntry->E;
  fP4.SetPxPyPzE(fCurEntry->px,fCurEntry->py,fCurEntry->pz,Ev);
  fX4.SetXYZT(fCurEntry->vtxx,fCurEntry->vtxy,fCurEntry->vtxz,0);

  fWeight = fCurEntry->wgt;

  if (Ev > fMaxEv) {
     LOG("Flux", pWARN)
          << "Flux neutrino energy exceeds declared maximum neutrino energy"
          << "\nEv = " << Ev << "(> Ev{max} = " << fMaxEv << ")";
  }

  // if desired, move to user specified user coord z
  if ( TMath::Abs(fZ0) < 1.0e30 ) this->MoveToZ0(fZ0);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pINFO)
    << "Generated neutrino: " << fIEntry 
#ifdef NOT_YET
    << " run " << fCurNuMI->evtno
    << " evtno " << fCurNuMI->evtno 
    << " entryno " << fCurNuMI->entryno 
#endif
    << "\n pdg-code: " << fCurEntry->pdg << " wgt " << fCurEntry->wgt
    << "\n p4: " << utils::print::P4AsShortString(&fP4)
    << "\n x4: " << utils::print::X4AsString(&fX4);
#endif
  if ( Ev > fMaxEv ) {
    LOG("Flux", pFATAL)
      << "Generated neutrino had E_nu = " << Ev << " > " << fMaxEv 
      << " maximum ";
    assert(0);
  }

  return true;
}
//___________________________________________________________________________
double GSimpleNtpFlux::GetDecayDist() const
{
  // return distance (user units) between dk point and start position
  // these are in beam units
  return fCurEntry->dist;
}
//___________________________________________________________________________
void GSimpleNtpFlux::MoveToZ0(double z0usr)
{
  // move ray origin to specified user z0
  // move beam coord entry correspondingly

  double pzusr    = fP4.Pz();
  if ( TMath::Abs(pzusr) < 1.0e-30 ) {
    // neutrino is moving almost entirely in x-y plane
    LOG("Flux", pWARN)
      << "MoveToZ0(" << z0usr << ") not possible due to pz_usr (" << pzusr << ")";
    return;
  }

  double scale = (z0usr - fX4.Z()) / pzusr; 
  //LOG("Flux",pDEBUG) 
  //  << "MoveToZ0: before x4=(" << fX4.X() << "," << fX4.Y() << "," << fX4.Z()
  //  << ") z0=" << z0usr << " pzusr=" << pzusr
  //  << " p4=(" << fP4.Px() << "," << fP4.Py() << "," << fP4.Pz() << ")";
  fX4 += (scale*fP4);
  //LOG("Flux",pDEBUG)
  //  << "MoveToZ0: after (" << fX4.X() << "," << fX4.Y() << "," << fX4.Z()
  //  << ")";

  // this scaling works for distances, but not the time component
  fX4.SetT(0);

}

//___________________________________________________________________________
void GSimpleNtpFlux::CalcEffPOTsPerNu()
{
  // do this if flux window changes or # of files changes

  if (!fNuFluxTree) return;  // not yet fully configured

  fEffPOTsPerNu = ( (double)fFilePOTs / (double)fNEntries );
}

//___________________________________________________________________________
double GSimpleNtpFlux::UsedPOTs(void) const
{
// Compute current number of flux POTs

  if (!fNuFluxTree) {
     LOG("Flux", pWARN)
          << "The flux driver has not been properly configured";
     return 0;	
  }
  return fAccumPOTs;
}

//___________________________________________________________________________
void GSimpleNtpFlux::LoadBeamSimData(const std::vector<string>& patterns,
                                     const std::string&         config )
{
// Loads a beam simulation root file into the GSimpleNtpFlux driver.

  fNuFluxFilePatterns = patterns;
  std::vector<int> nfiles_from_pattern;

  // create a (sorted) set of file names
  // this also helps ensure that the same file isn't listed multiple times
  std::set<std::string> fnames;

  LOG("Flux", pINFO) << "LoadBeamSimData was passed " << patterns.size()
                       << " patterns";

  for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
    string pattern = patterns[ipatt];
    nfiles_from_pattern.push_back(0);
    LOG("Flux", pINFO)
        << "Loading flux tree from ROOT file pattern [" 
        << std::setw(3) << ipatt << "] \"" << pattern << "\"";

    // !WILDCARD only works for file name ... NOT directory
    string dirname = gSystem->UnixPathName(gSystem->WorkingDirectory());
    size_t slashpos = pattern.find_last_of("/");
    size_t fbegin;
    if ( slashpos != std::string::npos ) {
      dirname = pattern.substr(0,slashpos);
      LOG("Flux", pDEBUG) << "Look for flux using directory " << dirname;
      fbegin = slashpos + 1;
    } else { fbegin = 0; }

    void* dirp = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname.c_str()));
    if ( dirp ) {
      std::string basename = 
        pattern.substr(fbegin,pattern.size()-fbegin);
      TRegexp re(basename.c_str(),kTRUE);
      const char* onefile;
      while ( ( onefile = gSystem->GetDirEntry(dirp) ) ) {
        TString afile = onefile;
        if ( afile=="." || afile==".." ) continue;
        if ( basename!=afile && afile.Index(re) == kNPOS ) continue;
        std::string fullname = dirname + "/" + afile.Data();
        fnames.insert(fullname);
        nfiles_from_pattern[ipatt]++;
      }
      gSystem->FreeDirectory(dirp);
    } // legal directory
  } // loop over patterns

  size_t indx = 0;
  std::set<string>::const_iterator sitr = fnames.begin();
  for ( ; sitr != fnames.end(); ++sitr, ++indx) {
    string filename = *sitr;
    //std::cout << "  [" << std::setw(3) << indx << "]  \"" 
    //          << filename << "\"" << std::endl;
    bool isok = true; 
    // this next test only works for local files, so we can't do that
    // if we want to stream via xrootd
    // ! (gSystem->AccessPathName(filename.c_str()));
    if ( ! isok ) continue;
    // open the file to see what it contains
    LOG("Flux", pINFO) << "Load file " <<  filename;
    
    TFile* tf = TFile::Open(filename.c_str(),"READ");
    TTree* etree = (TTree*)tf->Get("flux");
    if ( etree ) {
      TTree* mtree = (TTree*)tf->Get("meta");
      // add the file to the chain
      LOG("Flux", pDEBUG) << "AddFile " << filename
                          << " etree " << etree << " meta " << mtree;
      this->AddFile(etree,mtree,filename);
      
    } // found a GSimpleNtpEntry "flux" tree
    tf->Close();
    delete tf;
  } // loop over sorted file names

  // this will open all files and read headers!!
  fNEntries = fNuFluxTree->GetEntries();

  if ( fNEntries == 0 ) {
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    LOG("Flux", pERROR)
      << "Loaded flux tree contains " <<  fNEntries << " entries";
    LOG("Flux", pERROR)
      << "Was passed the file patterns: ";
    for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
      string pattern = patterns[ipatt];
      LOG("Flux", pERROR)
        << "  [" << std::setw(3) << ipatt <<"] " << pattern;
    }
    LOG("Flux", pERROR)
      << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
  } else {
    LOG("Flux", pNOTICE)
      << "Loaded flux tree contains " <<  fNEntries << " entries"
      << " from " << fnames.size() << " unique files";
    for (size_t ipatt = 0; ipatt < patterns.size(); ++ipatt ) {
      string pattern = patterns[ipatt];
      LOG("Flux", pINFO)
        << " pattern: " << pattern << " yielded " 
        << nfiles_from_pattern[ipatt] << " files";
    }
  }

  int sba_status[3] = { -999, -999, -999 };
  // "entry" branch isn't optional ... contains the neutrino info
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,26,0)
  sba_status[0] = 
#endif
    fNuFluxTree->SetBranchAddress("entry",&fCurEntry);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,26,0)
  if ( sba_status[0] < 0 ) {
    LOG("Flux", pFATAL) 
      << "flux chain has no \"entry\" branch " << sba_status[0];
    assert(0);
  }
#endif
  //TBranch* bentry = fNuFluxTree->GetBranch("entry");
  //bentry->SetAutoDelete(false);

  if ( OptionalAttachBranch("numi") ) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,26,0)
    sba_status[1] = 
#endif
      fNuFluxTree->SetBranchAddress("numi",&fCurNuMI);
    //TBranch* bnumi = fNuFluxTree->GetBranch("numi");
    //bnumi->SetAutoDelete(false);
  } else { delete fCurNuMI; fCurNuMI = 0; }

  if ( OptionalAttachBranch("aux") ) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,26,0)
    sba_status[2] = 
#endif
      fNuFluxTree->SetBranchAddress("aux",&fCurAux);
    //TBranch* baux = fNuFluxTree->GetBranch("aux");
    //baux->SetAutoDelete(false);
  } else { delete fCurAux; fCurAux = 0; }

  LOG("Flux", pDEBUG)
    << " SetBranchAddress status: "
    << " \"entry\"=" << sba_status[0]
    << " \"numi\"=" << sba_status[1]
    << " \"aux\"=" << sba_status[2];

  // attach requested branches

  if (fMaxWeight<=0) {
     LOG("Flux", pDEBUG)
       << "Run ProcessMeta() as part of LoadBeamSimData";
     this->ProcessMeta();
  }

  // current ntuple cycle # (flux ntuples may be recycled)
  fICycle =  0;
  // pick a starting entry index [0:fNEntries-1]
  // pretend we just used up the the previous one
  RandomGen* rnd = RandomGen::Instance();
  fIUse   =  9999999;
  fIEntry = rnd->RndFlux().Integer(fNEntries) - 1;
  if ( config.find("no-offset-index") != string::npos ) {
    LOG("Flux",pINFO) << "Config saw \"no-offset-index\"";  
    fIEntry = -1;
  }
  LOG("Flux",pINFO) << "Start with entry fIEntry=" << fIEntry;  

  // don't count things we used to estimate max weight
  fNEntriesUsed = 0;
  fSumWeight  = 0;
  fNNeutrinos = 0;
  fAccumPOTs  = 0;

  LOG("Flux",pDEBUG) << "about to CalcEffPOTsPerNu";
  this->CalcEffPOTsPerNu();
  
}
//___________________________________________________________________________
void GSimpleNtpFlux::GetBranchInfo(std::vector<std::string>& branchNames,
                                   std::vector<std::string>& branchClassNames,
                                   std::vector<void**>&      branchObjPointers)
{
  // allow flux driver to report back current status and/or ntuple entry 
  // info for possible recording in the output file by supplying
  // the class name, and a pointer to the object that will be filled
  // as well as a suggested name for the branch.

  if ( fCurEntry ) {
    branchNames.push_back("simple");
    branchClassNames.push_back("genie::flux::GSimpleNtpEntry");
    branchObjPointers.push_back((void**)&fCurEntry);
  }

  if ( fCurNuMI ) {
    branchNames.push_back("numi");
    branchClassNames.push_back("genie::flux::GSimpleNtpNuMI");
    branchObjPointers.push_back((void**)&fCurNuMI);
  }

  if ( fCurAux ) {
    branchNames.push_back("aux");
    branchClassNames.push_back("genie::flux::GSimpleNtpAux");
    branchObjPointers.push_back((void**)&fCurAux);
  }
}
TTree* GSimpleNtpFlux::GetMetaDataTree() { return fNuMetaTree; }

//___________________________________________________________________________
void GSimpleNtpFlux::ProcessMeta(void)
{

  fAlreadyUnwgt = false;
  fFilePOTs     = 0;
  double minwgt = +1.0e10;
  double maxwgt = -1.0e10;
  double maxenu = 0.0;

  // PDGLibrary* pdglib = PDGLibrary::Instance(); // get initialized now

  if ( fAllFilesMeta ) {
    fNuMetaTree->SetBranchAddress("meta",&fCurMeta);
#ifdef USE_INDEX_FOR_META
    int nindices = fNuMetaTree->BuildIndex("metakey"); // key used to tie entries to meta data
    LOG("Flux", pDEBUG) << "ProcessMeta() BuildIndex nindices " << nindices;
#endif
    int nmeta = fNuMetaTree->GetEntries();
    for (int imeta = 0; imeta < nmeta; ++imeta ) {
      fNuMetaTree->GetEntry(imeta);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("Flux", pNOTICE) << "ProcessMeta() ifile " << imeta
                           << " (of " << fNFiles
                           << ") " << *fCurMeta;
#endif
      minwgt = TMath::Min(minwgt,fCurMeta->minWgt);
      maxwgt = TMath::Max(maxwgt,fCurMeta->maxWgt);
      maxenu = TMath::Max(maxenu,fCurMeta->maxEnergy);
      fFilePOTs += fCurMeta->protons;
      for (size_t i = 0; i < fCurMeta->pdglist.size(); ++i)
        fPdgCList->push_back(fCurMeta->pdglist[i]);
    }
    if ( minwgt == 1.0 && maxwgt == 1.0 ) fAlreadyUnwgt = true;
    fMaxWeight = maxwgt;
    this->SetMaxEnergy(maxenu);

  } else {
    //
    LOG("Flux", pFATAL) << "ProcessMeta() not all files have metadata";
    // for now PUNT ... eventually could scan all the entries
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Flux", pNOTICE) << "ProcessMeta() Maximum flux weight = " << fMaxWeight 
                       << ", energy = " << fMaxEv
                       << ", AlreadyUnwgt=" << (fAlreadyUnwgt?"true":"false");
#endif

  fCurMeta->Reset();
  fIFileNumber = -999;

}
//___________________________________________________________________________
void GSimpleNtpFlux::SetMaxEnergy(double Ev)
{
  fMaxEv = TMath::Max(0.,Ev);

  LOG("Flux", pINFO)
    << "Declared maximum flux neutrino energy: " << fMaxEv;
}
//___________________________________________________________________________
void GSimpleNtpFlux::SetEntryReuse(long int nuse)
{
// With nuse > 1 then the same entry in the file is used "nuse" times
// before moving on to the next entry in the ntuple

  fNUse    = TMath::Max(1L, nuse);
}
//___________________________________________________________________________
void GSimpleNtpFlux::GetFluxWindow(TVector3& p0, TVector3& p1, TVector3& p2) const
{
  // return flux window points
  
  p0 = TVector3(fCurMeta->windowBase[0],
                fCurMeta->windowBase[1],
                fCurMeta->windowBase[2]);
  p1 = p0 + TVector3(fCurMeta->windowDir1[0],
                     fCurMeta->windowDir1[1],
                     fCurMeta->windowDir1[2]);
  p2 = p0 + TVector3(fCurMeta->windowDir2[0],
                     fCurMeta->windowDir2[1],
                     fCurMeta->windowDir2[2]);

}

//___________________________________________________________________________
void GSimpleNtpFlux::PrintCurrent(void)
{
  LOG("Flux", pNOTICE) << "CurrentEntry:" << *fCurEntry;
}
//___________________________________________________________________________
void GSimpleNtpFlux::Clear(Option_t * opt)
{
// Clear the driver state
//
  LOG("Flux", pWARN) << "GSimpleNtpFlux::Clear(" << opt << ") called";
  // do it in all cases, but EVGDriver/GMCJDriver will pass "CycleHistory"

  fICycle     = 0;

  fSumWeight  = 0;
  fNNeutrinos = 0;
  fAccumPOTs  = 0;

}
//___________________________________________________________________________
void GSimpleNtpFlux::GenerateWeighted(bool gen_weighted)
{
// Set whether to generate weighted rays
//
  fGenWeighted = gen_weighted;
}
//___________________________________________________________________________
void GSimpleNtpFlux::Initialize(void)
{
  LOG("Flux", pINFO) << "Initializing GSimpleNtpFlux driver";

  fMaxEv           =  0;
  fEnd             =  false;

  fCurEntry        = new GSimpleNtpEntry;
  fCurNuMI         = new GSimpleNtpNuMI;
  fCurAux          = new GSimpleNtpAux;
  fCurMeta         = new GSimpleNtpMeta;

  fCurEntryCopy    = 0;
  fCurNuMICopy     = 0;
  fCurAuxCopy      = 0;

  fNuFluxTree      = new TChain("flux");
  fNuMetaTree      = new TChain("meta");

  //fNuFluxFilePatterns = "";
  fNuFluxBranchRequest = "entry,numi,aux";  // all branches

  fNFiles          =  0;

  fNEntries        =  0;
  fIEntry          = -1;
  fIFileNumber     =  0;
  fICycle          =  0;
  fNUse            =  1;
  fIUse            =  999999;

  fFilePOTs        = 0;

  fMaxWeight       = -1;

  fSumWeight       =  0;
  fNNeutrinos      =  0;
  fNEntriesUsed    =  0;
  fEffPOTsPerNu    =  0;
  fAccumPOTs       =  0;

  fGenWeighted     = false;
  fAllFilesMeta    = true;
  fAlreadyUnwgt    = false;

  this->SetDefaults();
  this->ResetCurrent();
}
//___________________________________________________________________________
void GSimpleNtpFlux::SetDefaults(void)
{

// - Set the default flux neutrino start z position at use flux window
// - Set number of cycles to 1

  LOG("Flux", pINFO) << "Setting default GSimpleNtpFlux driver options";

  this->SetUpstreamZ     (-3.4e38); // way upstream ==> use flux window
  this->SetNumOfCycles   (0);
  this->SetEntryReuse    (1);
}
//___________________________________________________________________________
void GSimpleNtpFlux::ResetCurrent(void)
{
// reset running values of neutrino pdg-code, 4-position & 4-momentum
// and the input ntuple leaves

  if (fCurEntry) fCurEntry->Reset();
  if (fCurNuMI)  fCurNuMI->Reset();
  if (fCurAux)   fCurAux->Reset();
  //allow caching//if (fCurMeta)  fCurMeta->Reset();
}
//___________________________________________________________________________
void GSimpleNtpFlux::CleanUp(void)
{
  LOG("Flux", pINFO) << "Cleaning up...";

  if (fPdgCList)    delete fPdgCList;
  if (fPdgCListRej) delete fPdgCListRej;
  if (fCurEntry)    delete fCurEntry;
  if (fCurNuMI)     delete fCurNuMI;
  if (fCurAux)      delete fCurAux;
  if (fCurMeta)     delete fCurMeta;

  if (fNuFluxTree)  delete fNuFluxTree;
  if (fNuMetaTree)  delete fNuMetaTree;

  LOG("Flux", pNOTICE)
    << " flux file cycles: " << fICycle << " of " << fNCycles 
    << ", entry " << fIEntry << " use: " << fIUse << " of " << fNUse;
}

//___________________________________________________________________________
void GSimpleNtpFlux::AddFile(TTree* fluxtree, TTree* metatree, string fname)
{
  // Add a file to the chain
  if ( ! fluxtree ) return;

  ULong64_t nentries = fluxtree->GetEntries();

  int stat = fNuFluxTree->AddFile(fname.c_str());
  if ( metatree ) fNuMetaTree->AddFile(fname.c_str());
  else            fAllFilesMeta = false;

  LOG("Flux",pINFO)
    << "flux->AddFile() of " << nentries 
    << " " << ((metatree)?"[+meta]":"[no-meta]")
    << " [status=" << stat << "]"
    << " entries in file: " << fname;

  if ( stat != 1 ) return;

  fNFiles++;

}

//___________________________________________________________________________

bool GSimpleNtpFlux::OptionalAttachBranch(std::string name)
{

  if ( fNuFluxBranchRequest.find(name) == string::npos ) {
    LOG("Flux", pINFO)
      << "no request for \"" << name <<"\" branch ";
    return false;
  }

  if ( ( fNuFluxTree->GetBranch(name.c_str()) ) ) return true;

  LOG("Flux", pINFO)
    << "no \"" << name << "\" branch in the \"flux\" tree";
  return false;
}

//___________________________________________________________________________
GSimpleNtpEntry::GSimpleNtpEntry() {  Reset(); }

void GSimpleNtpEntry::Reset()
{
  wgt      = 0.;
  vtxx     = 0.;
  vtxy     = 0.;
  vtxz     = 0.;
  dist     = 0.;
  px       = 0.;
  py       = 0.;
  pz       = 0.;
  E        = 0.;

  pdg      =  0;
  metakey  =  0;
}

void GSimpleNtpEntry::Print(const Option_t* /* opt */ ) const
{
  std::cout << *this << std::endl;
}

//___________________________________________________________________________
GSimpleNtpNuMI::GSimpleNtpNuMI() { Reset(); }

void GSimpleNtpNuMI::Reset()
{ 
  tpx  = tpy  = tpz  = 0.;
  vx   = vy   = vz   = 0.;
  pdpx = pdpy = pdpz = 0.;
  pppx = pppy = pppz = 0.;
  
  ndecay   =  0;
  ptype    =  0;
  ppmedium =  0;
  tptype   =  0;
  run      = -1;
  evtno    = -1;
  entryno  = -1;
}

void GSimpleNtpNuMI::Print(const Option_t* /* opt */ ) const
{
  std::cout << *this << std::endl;
}

//___________________________________________________________________________
GSimpleNtpAux::GSimpleNtpAux() { Reset(); }

void GSimpleNtpAux::Reset()
{ 
  auxint.clear();
  auxdbl.clear();
}

void GSimpleNtpAux::Print(const Option_t* /* opt */ ) const
{
  std::cout << *this << std::endl;
}

//___________________________________________________________________________


GSimpleNtpMeta::GSimpleNtpMeta()
  : TObject() //, nflavors(0), flavor(0)
{ 
  Reset();
}

GSimpleNtpMeta::~GSimpleNtpMeta()
{ 
  Reset();
}

void GSimpleNtpMeta::Reset()
{
  
  pdglist.clear();
  maxEnergy    = 0.;
  minWgt       = 0.;
  maxWgt       = 0.;
  protons      = 0.;
  windowBase[0]  = 0.;
  windowBase[1]  = 0.;
  windowBase[2]  = 0.;
  windowDir1[0]  = 0.;
  windowDir1[1]  = 0.;
  windowDir1[2]  = 0.;
  windowDir2[0]  = 0.;
  windowDir2[1]  = 0.;
  windowDir2[2]  = 0.;

  auxintname.clear();
  auxdblname.clear();
  infiles.clear();

  seed     = 0;
  metakey  = 0;
}

void GSimpleNtpMeta::AddFlavor(Int_t nupdg)
{
  bool found = false;
  for (size_t i=0; i < pdglist.size(); ++i) 
    if ( pdglist[i] == nupdg) found = true;
  if ( ! found ) pdglist.push_back(nupdg);

  /* // OLD fashion array
  bool found = false;
  for (int i=0; i < nflavors; ++i) if ( flavor[i] == nupdg ) found = true;
  if ( ! found ) {
    Int_t* old_list = flavor;
    flavor = new Int_t[nflavors+1];
    for (int i=0; i < nflavors; ++i) flavor[i] = old_list[i];
    flavor[nflavors] = nupdg;
    nflavors++;
    delete [] old_list;
  }
  */
}

void GSimpleNtpMeta::Print(const Option_t* /* opt */ ) const
{
  std::cout << *this << std::endl;
}

//___________________________________________________________________________


namespace genie {
namespace flux  {

  ostream & operator << (
    ostream & stream, const genie::flux::GSimpleNtpEntry & entry)
    {
      stream << "\nGSimpleNtpEntry " 
             << " PDG " << entry.pdg 
             << " wgt " << entry.wgt
             << " ( metakey " << entry.metakey << " )"
             << "\n   vtx [" << entry.vtxx << "," << entry.vtxy << ","
             << entry.vtxz << "] dist " << entry.dist
             << "\n   p4  [" << entry.px << "," << entry.py << ","
             << entry.pz << "," << entry.E << "]";
     return stream;
  }


ostream & operator << (ostream & stream, 
                       const genie::flux::GSimpleNtpNuMI & numi)
{
  stream << "\nGSimpleNtpNuMI "
         << "run " << numi.run 
         << " evtno " << numi.evtno
         << " entryno " << numi.entryno
         << "\n   ndecay " << numi.ndecay  << " ptype " << numi.ptype
         << "\n   tptype " << numi.tptype << " ppmedium " << numi.ppmedium
         << "\n  tp[xyz] [" << numi.tpx  << "," << numi.tpy  << "," << numi.tpz  << "]"
         << "\n   v[xyz] [" << numi.vx   << "," << numi.vy   << "," << numi.vz   << "]"
         << "\n  pd[xyz] [" << numi.pdpx << "," << numi.pdpy << "," << numi.pdpz << "]"
         << "\n  pp[xyz] [" << numi.pppx << "," << numi.pppy << "," << numi.pppz << "]"
    ;
  return stream;
}

  ostream & operator << (
    ostream & stream, const genie::flux::GSimpleNtpAux & aux)
    {
      stream << "\nGSimpleNtpAux ";
      size_t nInt = aux.auxint.size();
      stream << "\n   ints: ";
      for (size_t ijInt=0; ijInt < nInt; ++ijInt)
        stream << " " << aux.auxint[ijInt];
      size_t nDbl = aux.auxdbl.size();
      stream << "\n   doubles: ";
      for (size_t ijDbl=0; ijDbl < nDbl; ++ijDbl)
        stream << " " << aux.auxdbl[ijDbl];

     return stream;
  }

  ostream & operator << (
    ostream & stream, const genie::flux::GSimpleNtpMeta & meta)
    {
      size_t nf = meta.pdglist.size();
      stream << "\nGSimpleNtpMeta " << nf << " flavors: ";
      for (size_t i=0; i<nf; ++i) stream << " " << meta.pdglist[i];

      //stream << "\nGSimpleNtpMeta " << meta.nflavors
      //       << " flavors: ";
      //for (int i=0; i< meta.nflavors; ++i) stream << " " << meta.flavor[i];
             
      stream << "\n maxEnergy " << meta.maxEnergy
             << " min/maxWgt " << meta.minWgt << "/" << meta.maxWgt
             << " protons " << meta.protons
             << " metakey " << meta.metakey
             << "\n windowBase [" << meta.windowBase[0] << ","
             << meta.windowBase[1] << "," << meta.windowBase[2] << "]"
             << "\n windowDir1 [" << meta.windowDir1[0] << ","
             << meta.windowDir1[1] << "," << meta.windowDir1[2] << "]"
             << "\n windowDir2 [" << meta.windowDir2[0] << ","
             << meta.windowDir2[1] << "," << meta.windowDir2[2] << "]";

      size_t nInt = meta.auxintname.size();
      if ( nInt > 0 ) stream << "\n aux ints:    ";
      for (size_t ijInt=0; ijInt < nInt; ++ijInt)
        stream << " " << meta.auxintname[ijInt];

      size_t nDbl = meta.auxdblname.size();
      if ( nDbl > 0 ) stream << "\n aux doubles: ";
      for (size_t ijDbl=0; ijDbl < nDbl; ++ijDbl)
        stream << " " << meta.auxdblname[ijDbl];

      size_t nfiles = meta.infiles.size();
      stream << "\n " << nfiles << " input files: ";
      UInt_t nprint = TMath::Min(UInt_t(nfiles),
                                 genie::flux::GSimpleNtpMeta::mxfileprint);
      for (UInt_t ifiles=0; ifiles < nprint; ++ifiles)
        stream << "\n    " << meta.infiles[ifiles];

      stream << "\n input seed: " << meta.seed;

     return stream;
  }


}//flux
}//genie

//___________________________________________________________________________

void GSimpleNtpFlux::PrintConfig()
{
  
  std::ostringstream s;
  PDGCodeList::const_iterator itr = fPdgCList->begin();
  for ( ; itr != fPdgCList->end(); ++itr) s << (*itr) << " ";
  s << "[rejected: ";
  itr = fPdgCListRej->begin();
  for ( ; itr != fPdgCListRej->end(); ++itr) s << (*itr) << " ";
  s << " ] ";

  std::ostringstream fpattout;
  for (size_t i = 0; i < fNuFluxFilePatterns.size(); ++i)
    fpattout << "\n [" << std::setw(3) << i << "] " << fNuFluxFilePatterns[i];

  std::ostringstream flistout;
  std::vector<std::string> flist = GetFileList();
  for (size_t i = 0; i < flist.size(); ++i)
    flistout << "\n [" << std::setw(3) << i << "] " << flist[i];

  LOG("Flux", pNOTICE)
    << "GSimpleNtpFlux Config:"
    << "\n Enu_max " << fMaxEv 
    << "\n pdg-codes: " << s.str() << "\n "
    << "\"flux\" " << fNEntries << " entries, " 
    << "\"meta\" " << fNFiles << " entries" 
    << " (FilePOTs " << fFilePOTs << ") in files:"
    << flistout.str()
    <<  "\n from file patterns: "
    << fpattout.str()
    << "\n wgt max=" << fMaxWeight 
    << "\n Z0 pushback " << fZ0
    << "\n used entry " << fIEntry << " " << fIUse << "/" << fNUse
    << " times, in " << fICycle << "/" << fNCycles << " cycles"
    << "\n SumWeight " << fSumWeight << " for " << fNNeutrinos << " neutrinos"
    << " with " << fNEntriesUsed << " entries read"
    << "\n EffPOTsPerNu " << fEffPOTsPerNu << " AccumPOTs " << fAccumPOTs
    << "\n GenWeighted \"" << (fGenWeighted?"true":"false") << "\""
    << " AlreadyUnwgt \"" << (fAlreadyUnwgt?"true":"false") << "\""
    << " AllFilesMeta \"" << (fAllFilesMeta?"true":"false") << "\"";

}

//___________________________________________________________________________
std::vector<std::string> GSimpleNtpFlux::GetFileList() 
{
  std::vector<std::string> flist;
  TObjArray *fileElements=fNuFluxTree->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0;
  while (( chEl=(TChainElement*)next() )) {
    flist.push_back(chEl->GetTitle());
  }
  return flist;
}

//___________________________________________________________________________
