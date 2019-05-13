//____________________________________________________________________________
/*!

\program gevscan

\brief   A utility that reads-in a GHEP event tree and performs basic sanity 
         checks / test whether the generated events obey basic conservation laws

\syntax  gevscan
             -f ghep_event_file 
            [-o output_error_log_file]
            [-n nev1[,nev2]]
            [--add-event-printout-in-error-log]
            [--max-num-of-errors-shown n]
            [--event-record-print-level level]
            [--check-energy-momentum-conservation]
            [--check-charge-conservation]
            [--check-for-pseudoparticles-in-final-state]
            [--check-for-off-mass-shell-particles-in-final-state]
            [--check-for-num-of-final-state-nucleons-inconsistent-with-target]
            [--check-vertex-distribution]
            [--check-decayer-consistency]
            [--all]

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created August 13, 2008

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

//#define __debug__

#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLorentzVector.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"

using std::ostringstream;
using std::ofstream;
using std::string;
using std::vector;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;
using std::endl;

using namespace genie;
using namespace genie::constants;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
bool CheckRootFilename  (string filename);

// checks
void CheckEnergyMomentumConservation (void);
void CheckChargeConservation (void);
void CheckForPseudoParticlesInFinState (void);
void CheckForOffMassShellParticlesInFinState (void);
void CheckForNumFinStateNucleonsInconsistentWithTarget (void);
void CheckVertexDistribution (void);
void CheckDecayerConsistency (void);

// options
string   gOptInpFilename = "";
string   gOptOutFilename = "";
Long64_t gOptNEvtL = -1;
Long64_t gOptNEvtH = -1;
int      gOptMaxNumErrs = -1; 
bool     gOptAddEventPrintoutInErrLog = false;
bool     gOptCheckEnergyMomentumConservation = false;
bool     gOptCheckChargeConservation = false;
bool     gOptCheckForPseudoParticlesInFinState = false;
bool     gOptCheckForOffMassShellParticlesInFinState = false;
bool     gOptCheckForNumFinStateNucleonsInconsistentWithTarget = false;
bool     gOptCheckVertexDistribution = false;
bool     gOptCheckDecayerConsistency = false;

Long64_t gFirstEventNum = -1;
Long64_t gLastEventNum  = -1;

TTree *            gEventTree = 0;
NtpMCEventRecord * gMCRec = 0;
ofstream           gErrLog;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  TFile file(gOptInpFilename.c_str(),"READ");

  NtpMCTreeHeader * thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
  LOG("gevscan", pINFO) << "Input tree header: " << *thdr;
  NtpMCFormat_t format = thdr->format;
  if(format != kNFGHEP) {
      LOG("gevscan", pERROR) 
        << "*** Unsupported event-tree format : "
        << NtpMCFormat::AsString(format);
      file.Close();
      return 3;
  }

  gEventTree = dynamic_cast <TTree *> (file.Get("gtree"));
  gEventTree->SetBranchAddress("gmcrec", &gMCRec);

  Long64_t nev = gEventTree->GetEntries();
  if(gOptNEvtL == -1 && gOptNEvtH == -1) {
    // read all events
    gFirstEventNum = 0;
    gLastEventNum  = nev-1;
  }
  else {
    // read a range of events
    gFirstEventNum = TMath::Max((Long64_t)0,  gOptNEvtL);
    gLastEventNum = TMath::Min(nev-1,        gOptNEvtH);
    if(gLastEventNum - gFirstEventNum < 0) {
      LOG("gevdump", pFATAL) << "Invalid event range";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
  }

  
  if(gOptOutFilename.size() == 0) {
     ostringstream logfile;
     logfile << gOptInpFilename << ".errlog";
     gOptOutFilename = logfile.str();
  }
  if(gOptOutFilename != "none") {
     gErrLog.open(gOptOutFilename.c_str());
     gErrLog << "# ..................................................................................." << endl;
     gErrLog << "# Error log for event file " << gOptInpFilename << endl;
     gErrLog << "# ..................................................................................." << endl;
     gErrLog << "# " << endl;
  }

  if (gOptCheckEnergyMomentumConservation) {
	  CheckEnergyMomentumConservation();
  }
  if (gOptCheckChargeConservation) {
          CheckChargeConservation();
  }
  if (gOptCheckForPseudoParticlesInFinState) {
          CheckForPseudoParticlesInFinState();
  }
  if (gOptCheckForOffMassShellParticlesInFinState) {
          CheckForOffMassShellParticlesInFinState();
  }
  if (gOptCheckForNumFinStateNucleonsInconsistentWithTarget) {
          CheckForNumFinStateNucleonsInconsistentWithTarget();
  }
  if (gOptCheckVertexDistribution) {
          CheckVertexDistribution();
  }
  if (gOptCheckDecayerConsistency) {
          CheckDecayerConsistency();
  }


  if(gOptOutFilename != "none") {
     gErrLog.close();
  }

  return 0;
}
//____________________________________________________________________________
void CheckEnergyMomentumConservation (void)
{
  LOG("gevscan", pNOTICE) << "Checking energy/momentum conservation...";

  if(gErrLog.is_open()) {
    gErrLog << "# Events failing the energy-momentum conservation test:" << endl;
    gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    if(gOptMaxNumErrs != -1 && nerr >= gOptMaxNumErrs) break;

    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gevscan", pINFO) << "Checking event.... " << i;

    double E_init  = 0, E_fin  = 0; // E
    double px_init = 0, px_fin = 0; // px
    double py_init = 0, py_fin = 0; // py
    double pz_init = 0, pz_fin = 0; // pz

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      GHepStatus_t ist  = p->Status();

      if(ist == kIStInitialState) 
      {
         E_init  += p->E();
         px_init += p->Px();
         py_init += p->Py();
         pz_init += p->Pz();
       }
       if(ist == kIStStableFinalState || 
          ist == kIStFinalStateNuclearRemnant) 
       {
         E_fin   += p->E();
         px_fin  += p->Px();
         py_fin  += p->Py();
         pz_fin  += p->Pz();
       }
    }//p

    double epsilon = 1E-3; 

    bool E_conserved  = TMath::Abs(E_init  - E_fin)  < epsilon;
    bool px_conserved = TMath::Abs(px_init - px_fin) < epsilon;
    bool py_conserved = TMath::Abs(py_init - py_fin) < epsilon;
    bool pz_conserved = TMath::Abs(pz_init - pz_fin) < epsilon;

    bool ok = E_conserved  && 
              px_conserved &&
              py_conserved &&
              pz_conserved;

    if(!ok) {
       LOG("gevscan", pERROR) 
         << " ** Energy-momentum non-conservation in event: " << i 
         << "\n"
         << event;
       if(gErrLog.is_open()) {
           gErrLog << i;
           if(gOptAddEventPrintoutInErrLog) {
               gErrLog << event;
           }
       }
       nerr++;
    }
    
    gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6

  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
         gErrLog << "none" << endl;    
     }
  }

  LOG("gevscan", pNOTICE) 
     << "Found " << nerr 
     << " events failing the energy/momentum conservation test";
}
//____________________________________________________________________________
void CheckChargeConservation(void)
{
  LOG("gevscan", pNOTICE) << "Checking charge conservation...";

  if(gErrLog.is_open()) {
     gErrLog << "# Events failing the charge conservation test:" << endl;
     gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    if(gOptMaxNumErrs != -1 && nerr >= gOptMaxNumErrs) break;

    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gevscan", pINFO) << "Checking event.... " << i;

    // Can't run the test for neutrinos scattered off nuclear targets
    // because of intranuclear rescattering effects and the presence, in the event
    // record, of a charged nuclear remnant pseudo-particle whose charge is not stored.
    // To check charge conservation in the primary interaction, use a sample generated
    // for a free nucleon targets.
    GHepParticle * nucltgt = event.TargetNucleus();
    if (nucltgt) {
      LOG("gevscan", pINFO)
           << "Event in nuclear target - Skipping test...";
    }
    else {
      double Q_init  = 0;
      double Q_fin   = 0; 

      GHepParticle * p = 0;
      TIter event_iter(&event);
      while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

        GHepStatus_t ist  = p->Status();

        if(ist == kIStInitialState) 
        {
           Q_init  += p->Charge();
         }
         if(ist == kIStStableFinalState)
         {
           Q_fin  += p->Charge();
         }
      }//p

      double epsilon = 1E-3; 
      bool ok = TMath::Abs(Q_init - Q_fin)  < epsilon;
      if(!ok) {
         LOG("gevscan", pERROR) 
           << " ** Charge non-conservation in event: " << i 
           << "\n"
           << event;
         if(gErrLog.is_open()) {
            gErrLog << i << endl;    
            if(gOptAddEventPrintoutInErrLog) {
                 gErrLog << event;
            }
         }
         nerr++;
      }
      
    }
    gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6
  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
        gErrLog << "none" << endl;    
     }
  }

  LOG("gevscan", pNOTICE) 
     << "Found " << nerr 
     << " events failing the charge conservation test";
}
//____________________________________________________________________________
void CheckForPseudoParticlesInFinState(void)
{
  LOG("gevscan", pNOTICE) 
      << "Checking for pseudo-particles appearing in final state...";

  if(gErrLog.is_open()) {
     gErrLog << "# Events with pseudo-particles in final state:" << endl;
     gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    if(gOptMaxNumErrs != -1 && nerr >= gOptMaxNumErrs) break;

    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gevscan", pINFO) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    bool ok = true;
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      GHepStatus_t ist = p->Status();
      if(ist != kIStStableFinalState) continue;
      int pdgc = p->Pdg();      
      if(pdg::IsPseudoParticle(pdgc))
      {
        ok = false;
        break;
      }      
    }//p

    if(!ok) {
       LOG("gevscan", pERROR) 
         << " ** Pseudo-particle final state particle in event: " << i 
         << "\n"
         << event;
       if(gErrLog.is_open()) {
          gErrLog << i << endl;    
          if(gOptAddEventPrintoutInErrLog) {
               gErrLog << event;
          }
       }
       nerr++;
    }

    gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6

  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
        gErrLog << "none" << endl;    
     }
  }

  LOG("gevscan", pNOTICE) 
     << "Found " << nerr 
     << " events with pseudo-particles in  final state";
}
//____________________________________________________________________________
void CheckForOffMassShellParticlesInFinState(void)
{
  LOG("gevscan", pNOTICE) 
      << "Checking for off-mass-shell particles appearing in the final state...";

  if(gErrLog.is_open()) {
     gErrLog << "# Events with off-mass-shell particles in final state:" << endl;
     gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    if(gOptMaxNumErrs != -1 && nerr >= gOptMaxNumErrs) break;

    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gevscan", pINFO) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    bool ok = true;
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

      GHepStatus_t ist = p->Status();
      if(ist != kIStStableFinalState) continue;
      if(p->IsOffMassShell())
      {
        ok = false;
        break;
      }      
    }//p

    if(!ok) {
       LOG("gevscan", pERROR) 
         << " ** Off-mass-shell final state particle in event: " << i 
         << "\n"
         << event;
       if(gErrLog.is_open()) {
          gErrLog << i << endl;    
          if(gOptAddEventPrintoutInErrLog) {
               gErrLog << event;
          }
       }
       nerr++;
    }

    gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6

  }//i

  if(gErrLog.is_open()) {
     if(nerr == 0) {
        gErrLog << "none" << endl;    
     }
  }

  LOG("gevscan", pNOTICE) 
     << "Found " << nerr 
     << " events with off-mass-shell particles in final state";
}
//____________________________________________________________________________
void CheckForNumFinStateNucleonsInconsistentWithTarget(void)
{
  LOG("gevscan", pNOTICE) 
     << "Checking for number of final state nucleons inconsistent with target...";

  if(gErrLog.is_open()) {
    gErrLog << "# Events with number of final state nucleons inconsistent with target:" << endl;
    gErrLog << "# " << endl;
  }

  int nerr = 0;

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    if(gOptMaxNumErrs != -1 && nerr >= gOptMaxNumErrs) break;

    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gevscan", pINFO) << "Checking event.... " << i;

    // get target nucleus
    GHepParticle * nucltgt = event.TargetNucleus();
    if (!nucltgt) {
      LOG("gevscan", pINFO)
           << "Event not in nuclear target - Skipping test...";
    }
    else {
      GHepParticle * p = 0;

      int Z = 0;
      int N = 0;

      // get number of spectator nucleons 
      int fd = nucltgt->FirstDaughter();
      int ld = nucltgt->LastDaughter();
      for(int d = fd; d <= ld; d++) {
        p = event.Particle(d);
        if(!p) continue;
        int pdgc = p->Pdg();
        if(pdg::IsIon(pdgc)) {
          Z = p->Z();
          N = p->A() - p->Z();
        }
      }
      // add nucleons from the primary interaction
      TIter event_iter(&event);
      while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
        GHepStatus_t ist = p->Status();
        if(ist != kIStHadronInTheNucleus) continue;
        int pdgc = p->Pdg();
        if(pdg::IsProton (pdgc)) { Z++; }
        if(pdg::IsNeutron(pdgc)) { N++; }
      }//p

      LOG("gevscan", pINFO)
         << "Before intranuclear hadron transport: Z = " << Z << ", N = " << N;

      // count final state nucleons
      int Zf = 0;
      int Nf = 0;
      event_iter.Reset();
      while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
        GHepStatus_t ist = p->Status();
        if(ist != kIStStableFinalState) continue;
        int pdgc = p->Pdg();
        if(pdg::IsProton (pdgc)) { Zf++; }
        if(pdg::IsNeutron(pdgc)) { Nf++; }
      }
      LOG("gevscan", pINFO)
         << "In the final state: Z = " << Zf << ", N = " << Nf;

      bool ok = (Zf <= Z && Nf <= N);
      if(!ok) {
         LOG("gevscan", pERROR) 
           << " ** Number of final state nucleons inconsistent with target in event: " << i 
           << "\n"
           << event;
         if(gErrLog.is_open()) {
             gErrLog << i << endl;    
             if(gOptAddEventPrintoutInErrLog) {
                 gErrLog << event;
             }
         }
         nerr++;
      }
    } //nucltgt

    gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6

  }//i


  if(gErrLog.is_open()) {
     if(nerr == 0) {
         gErrLog << "none" << endl;    
     }
  }

  LOG("gevscan", pNOTICE) 
     << "Found " << nerr 
     << " events with a number of final state nucleons inconsistent with target";
}
//____________________________________________________________________________
void CheckVertexDistribution(void)
{
  LOG("gevscan", pNOTICE) 
     << "Checking intra-nuclear vertex distribution...";

  if(gErrLog.is_open()) {
    gErrLog << "# Intranuclear vertex distribution check:" << endl;
    gErrLog << "# " << endl;
  }

  TH1D * r_distr_mc       = new TH1D("r_distr_mc","",      150,0,30); //fm
  TH1D * r_distr_expected = new TH1D("r_distr_expected","",150,0,30); //fm

  int Z = -1;
  int A = -1;

  // get vertex position distribution
  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gevscan", pINFO) << "Checking event.... " << i;

    // get target nucleus
    GHepParticle * nucltgt = event.TargetNucleus();
    if (!nucltgt) {
      LOG("gevscan", pINFO)
           << "Event not in nuclear target - Skipping...";
    }
    else {
      if(Z == -1 && A == -1) {
         Z = nucltgt->Z();
         A = nucltgt->A();
      }

      // this test is run on a MC sample for a given target
      if(Z != nucltgt->Z() || A != nucltgt->A()) {
        LOG("gevscan", pINFO)
             << "Event not in nuclear target seen first - Skipping...";
      }
      else {
        GHepParticle * probe = event.Particle(0);
        double r = probe->X4()->Vect().Mag();

        r_distr_mc->Fill(r);
      }
    } //nucltgt

    gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6

  }//i

  if(A > 1) {
    // get expected vertex position distribution
    for(int ir = 1; ir <= r_distr_expected->GetNbinsX(); ir++) {
      double r = r_distr_expected->GetBinCenter(ir);
      double rho  = utils::nuclear::Density(r,A);
      double nexp = 4*kPi*r*r*rho;
      r_distr_expected->SetBinContent(ir,nexp);
    }

    // normalize 
    double N = r_distr_mc->GetEntries();
    r_distr_expected -> Scale (N / r_distr_expected -> Integral());

    // check consistency
    double pvalue = r_distr_mc->Chi2Test(r_distr_expected,"WWP");
    LOG("gevscan", pNOTICE) << "p-value {\\chi^2 test} = " << pvalue;

    if(gErrLog.is_open()) {
       if(pvalue < 0.99) {
         gErrLog << "Problem! p-value = " << pvalue << endl;    
       } else {
         gErrLog << "OK! p-value = " << pvalue << endl;    
       }
    }

#ifdef __debug__
    TFile f("./check_vtx.root","recreate");
    r_distr_mc -> Write();
    r_distr_expected -> Write();
    f.Close();
#endif

  }//A
  else {

    if(gErrLog.is_open()) {
      gErrLog << "Can not run test with current sample" << endl;   
    }

  }

}
//____________________________________________________________________________
void CheckDecayerConsistency(void)
{
// Check that particles seen in the final state in some events do not appear to 
// have decayed in other events.
// This might happen if, for example, particle decay flags which are applied to
// GENIE events do not get applied to intermediate particles appearing in the
// PYTHIA hadronization. It might also happen if the decayed particle status is
// used incorrectly in some modules (eg intranuke).
//
  LOG("gevscan", pNOTICE) 
     << "Checking decayer consistency...";

  if(gErrLog.is_open()) {
    gErrLog << "# Decayer consistency check:" << endl;
    gErrLog << "# " << endl;
  }

  bool allowdup = false;
  PDGCodeList final_state_particles(allowdup);
  PDGCodeList decayed_particles(allowdup);

  for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
  {
    gEventTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    LOG("gevscan", pINFO) << "Checking event.... " << i;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
      GHepStatus_t ist = p->Status();
      int pdgc = p->Pdg();
      if(ist == kIStStableFinalState) { final_state_particles.push_back(pdgc); }
      if(ist == kIStDecayedState    ) { decayed_particles.push_back(pdgc);     }
    }//p
    gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6
  }//i

  // find particles which appear in both lists
  PDGCodeList particles_in_both_lists(allowdup);

  PDGCodeList::const_iterator iter;
  for(iter = final_state_particles.begin(); 
      iter != final_state_particles.end(); ++iter) 
  {
     int pdgc = *iter;
     if(decayed_particles.ExistsInPDGCodeList(pdgc)) 
     {
        particles_in_both_lists.push_back(pdgc);
     }
  }

  bool ok = true;
  ostringstream mesg;
  if(particles_in_both_lists.size() == 0) {
    mesg << "OK.\n" << "No particle seen both in the final state and to have decayed.";
  } else {
    ok = false;
    mesg << "Problem!\n" << particles_in_both_lists.size() << " particles seen both final state and to have decayed.";
  }
 
  LOG("gevscan", pNOTICE) 
    << mesg.str();
  LOG("gevscan", pNOTICE) 
    << "Particles seen in final state: " << final_state_particles;
  LOG("gevscan", pNOTICE) 
    << "Particles seen to have decayed: " << decayed_particles;
  LOG("gevscan", pNOTICE) 
    << "Particles seen in both lists: " << particles_in_both_lists;

  if(gErrLog.is_open()) {
     gErrLog << mesg.str() << endl;
     gErrLog << "\nParticles seen in final state:" << final_state_particles << endl;
     gErrLog << "\nParticles seen to have decayed:" << decayed_particles << endl;
     gErrLog << "\nParticles seen in both lists:" << particles_in_both_lists << endl;
   }

   // find example events
   if(!ok) {
      if(gErrLog.is_open()) {
         gErrLog << "\nExample events: " << endl;          
      }
      for(iter  = particles_in_both_lists.begin(); 
          iter != particles_in_both_lists.end(); ++iter) 
      {
         int pdgc_bothlists = *iter;
         int iev_decay = -1;
         int iev_fs    = -1;
         bool have_example = false;
         for(Long64_t i = gFirstEventNum; i <= gLastEventNum; i++) 
         {
           have_example = (iev_decay != -1 && iev_fs != -1);
           if(have_example) break;

           gEventTree->GetEntry(i);
           EventRecord &  event = *(gMCRec->event);
           GHepParticle * p = 0;
           TIter event_iter(&event);
           while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
              int pdgc = p->Pdg();
              if(pdgc != pdgc_bothlists) continue;
              GHepStatus_t ist = p->Status();
              if(ist == kIStStableFinalState && iev_fs    == -1) { iev_fs    = i; }
              if(ist == kIStDecayedState     && iev_decay == -1) { iev_decay = i; }
           }//p
           gMCRec->Clear(); // clear out explicitly to prevent memory leak w/Root6
         }//i
         if(gErrLog.is_open()) {
            gErrLog << ">> " << PDGLibrary::Instance()->Find(pdgc_bothlists)->GetName()
                    << ": Decayed in event " << iev_decay 
                    << ". Seen in final state in event " << iev_fs << "." << endl;
            if(gOptAddEventPrintoutInErrLog) {
               gEventTree->GetEntry(iev_decay);
               EventRecord & event_dec = *(gMCRec->event);
               gErrLog << "Event " << iev_decay << ":";
               gErrLog << event_dec;
               gEventTree->GetEntry(iev_fs);
               EventRecord & event_fs = *(gMCRec->event);
               gErrLog << "Event: " << iev_fs << ":";
               gErrLog << event_fs;
            }
         }
      }//pdgc
   }//!ok

}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevscan", pNOTICE) << "*** Parsing command line arguments";

  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  CmdLnArgParser parser(argc,argv);
  
  // get input GENIE event sample
  if( parser.OptionExists('f') ) {
    LOG("gevscan", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("gevscan", pFATAL) 
        << "Unspecified input filename - Exiting";
    PrintSyntax();
    exit(1);
  }

  // get output error log
  if( parser.OptionExists('o') ) {
    LOG("gevscan", pINFO) << "Reading err log file name";
    gOptOutFilename = parser.ArgAsString('o');
  } 

  // number of events
  if ( parser.OptionExists('n') ) {
    LOG("gevdump", pINFO) << "Reading number of events to analyze";
    string nev =  parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n',",");
      if(vecn.size()!=2) {
         LOG("gevdump", pFATAL) << "Invalid syntax";
         PrintSyntax();
         gAbortingInErr = true;
         exit(1);
      }
      // read a range of events
      gOptNEvtL = vecn[0];
      gOptNEvtH = vecn[1];
    } else {
      // read single event
      gOptNEvtL = parser.ArgAsLong('n');
      gOptNEvtH = gOptNEvtL;
    }
  } else {
    LOG("gevdump", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvtL = -1;
    gOptNEvtH = -1;
  }

  gOptAddEventPrintoutInErrLog =
     parser.OptionExists("add-event-printout-in-error-log");

  if(parser.OptionExists("max-num-of-errors-shown")) {
     gOptMaxNumErrs = parser.ArgAsInt("max-num-of-errors-shown");
     gOptMaxNumErrs = TMath::Max(1,gOptMaxNumErrs);
  }
  
  bool all = parser.OptionExists("all");

  // checks
  gOptCheckEnergyMomentumConservation = all ||
     parser.OptionExists("check-energy-momentum-conservation");
  gOptCheckChargeConservation = all || 
     parser.OptionExists("check-charge-conservation");
  gOptCheckForNumFinStateNucleonsInconsistentWithTarget = all ||
     parser.OptionExists("check-for-num-of-final-state-nucleons-inconsistent-with-target");
  gOptCheckForPseudoParticlesInFinState = all ||
     parser.OptionExists("check-for-pseudoparticles-in-final-state");
  gOptCheckForOffMassShellParticlesInFinState = all ||
     parser.OptionExists("check-for-off-mass-shell-particles-in-final-state");
  gOptCheckVertexDistribution = all ||
     parser.OptionExists("check-vertex-distribution");
  gOptCheckDecayerConsistency = all ||
     parser.OptionExists("check-decayer-consistency");
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevscan", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << " gevscan -f sample.root [-n n1[,n2]] [-o errlog] [check names]\n";
}
//_________________________________________________________________________________
bool CheckRootFilename(string filename)
{
  if(filename.size() == 0) return false;
    
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gevscan", pERROR)  
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//_________________________________________________________________________________

