//____________________________________________________________________________
/*!

\program gmctest

\brief   A GENIE utility that reads a generated event sample and generates a 
         ROOT file with loads of histograms with basic MC truth quantities.
         Typically used in validating new generator releases.

         Syntax :
           gmctest -f filename -t tgt [-n events]

         Options:
           [] denotes an optional argument
           -f specifies the GENIE/ROOT file with the generated event sample
           -t target pdg code (analyzing interactions for one target at a time)
           -n specifies how many events to analyze [default: all]

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 02, 2005
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include "Conventions/Constants.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCSummary.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::constants;

// function prototypes
void   GetCommandLineArgs(int argc, char ** argv);
void   PrintSyntax       (void);
bool   CheckRootFilename (string filename);
string OutputFileName    (string input_file_name);
void   InitHistos        (void);
void   SaveHistos        (string outpfilename);
void   HistogramLimits   (Long64_t n, TTree* t, NtpMCEventRecord* r);
void   AnalyzeEvents     (Long64_t n, TTree* t, NtpMCEventRecord* r);
int    IProc             (const ProcessInfo &  pi);
int    ICurr             (const ProcessInfo &  pi);
int    INeut             (int pdgc);
string ProcName          (int i);
string CurrName          (int i);
string NeutName          (int i);
string HistoPtr          (string b, int iproc, int icurr, int ineut);
string HistoTitle        (string b, int iproc, int icurr, int ineut);

// bin sizes
const double kdE  = 0.25; // energy bin size
const double kdy  = 0.01; // inelasticity y bin size
const double kdx  = 0.01; // bjorken x bin size
const double kdQ2 = 0.10; // momentum transfer Q2 bin size
const double kdW  = 0.25; // invariant mass W bin size
const double kdv  = 0.01; // energy transfer v bin size

// output histograms
const int kNProc = 3; // QEL, RES, DIS
const int kNCurr = 2; // CC, NC
const int kNNeut = 6; // nue,numu,nutau,\bar(nue),\bar(numu),\bar(nutau)

TH1D * hEvLAB [kNProc][kNCurr][kNNeut]; // neutrino energy in LAB
TH1D * hElLAB [kNProc][kNCurr][kNNeut]; // fsl energy in LAB
TH1D * hKinX  [kNProc][kNCurr][kNNeut]; // x
TH1D * hKinY  [kNProc][kNCurr][kNNeut]; // y
TH2D * hKinXY [kNProc][kNCurr][kNNeut]; // x vs y
TH1D * hKinQ2 [kNProc][kNCurr][kNNeut]; // Q2
TH1D * hKinW  [kNProc][kNCurr][kNNeut]; // W
TH1D * hKinNu [kNProc][kNCurr][kNNeut]; // v
TH2D * hKinWQ2[kNProc][kNCurr][kNNeut]; // W vs Q2
TH2D * hMultPW   [kNCurr][kNNeut];
TH2D * hMultNW   [kNCurr][kNNeut];
TH2D * hMultPi0W [kNCurr][kNNeut];
TH2D * hMultPipW [kNCurr][kNNeut];
TH2D * hMultPimW [kNCurr][kNNeut];
TH1D * hVtxX;
TH1D * hVtxY;
TH1D * hVtxZ;
TH1D * hVtxT;
TH3D * hVtxR;

// 'detectable' histogram limits
double gEvMin   =  999999;
double gEvMax   = -999999;
double gWMin    =  999999;
double gWMax    = -999999;
double gQ2Min   =  999999;
double gQ2Max   = -999999;
double gVtxXMin =  999999;
double gVtxXMax = -999999;
double gVtxYMin =  999999;
double gVtxYMax = -999999;
double gVtxZMin =  999999;
double gVtxZMax = -999999;
double gVtxTMin =  999999;
double gVtxTMax = -999999;

// command-line arguments
Long64_t gOptN;       // process so many events, all if -1
Long64_t gOptTgtPdgC; // process events only for this nucl. target
string   gOptInpFile; // input GENIE event file

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments 
  GetCommandLineArgs(argc,argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  LOG("gmctest", pNOTICE) 
                  << "*** Opening GHEP data file: " << gOptInpFile;

  TFile inpfile(gOptInpFile.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( inpfile.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( inpfile.Get("header") );

  LOG("gmctest", pNOTICE) << "*** Input tree header: " << *thdr;

  //-- figure out the TTree format (GENIE supports multiple formats)
  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord);

  //-- set the event record branch ptr
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       tree->GetEntries() : TMath::Min( tree->GetEntries(), gOptN );

  //-- build output filename based on input filename
  string outpfilename = OutputFileName(gOptInpFile);

  //-- detect histogram limits
  HistogramLimits(nmax, tree, mcrec);

  //-- book output histograms
  InitHistos();

  //-- loop over the TTree records & fill in the histograms
  AnalyzeEvents(nmax, tree, mcrec);

  //-- save the histograms 
  SaveHistos(outpfilename);

  //-- close the input GENIE event file
  LOG("gmctest", pNOTICE) << "*** Closing input GHEP data stream";
  inpfile.Close();
  tree = 0;
  thdr = 0;

  LOG("gmctest", pINFO)  << "Done!";
  return 0;
}
//___________________________________________________________________
void AnalyzeEvents(
               Long64_t nmax, TTree * tree, NtpMCEventRecord * mcrec)
{
  LOG("gmctest", pNOTICE) << "*** Analyzing: " << nmax << " events";

  if ( nmax<0 ) return;
  if ( !tree  ) return;
  if ( !mcrec ) return;

  for(Long64_t i = 0; i < nmax; i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    // LOG("gmctest", pINFO) << rec_header;
    // LOG("gmctest", pINFO) << event;

    // go further inly if the event is physical
    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) continue;

    // nuclear target or free nucleon target
    GHepParticle * target = event.Particle(1);
    assert(target);

    // only further only if it matches the requested target
    if(target->PdgCode() != gOptTgtPdgC) continue;

    // neutrino
    GHepParticle * neutrino = event.Probe();
    assert(neutrino);
    // final state primary lepton
    GHepParticle * fsl = event.FinalStatePrimaryLepton();
    assert(fsl);
    // hit nucleon
    GHepParticle * hitnucl = event.StruckNucleon();
    if(!hitnucl) continue;

    const Interaction * interaction = event.GetInteraction();
    const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

    double weight = event.GetWeight();
    double Ev     = neutrino->Energy();
    double El     = fsl->Energy();

    const TLorentzVector & k1 = *(neutrino->P4()); // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());      // l 4-p (k2)
    const TLorentzVector & p1 = *(hitnucl->P4());  // N 4-p (p1)

    double M  = kNucleonMass;

    // see, for example, particle data booklet, page 210-211

    TLorentzVector q  = k1-k2;    // q=k1-k2, 4-p transfer
    double v  = (q*p1)/M;         // v (E transfer in hit nucleon rest frame)
    double Q2 = -1 * q.M2();      // momemtum transfer
    double x  = 0.5*Q2/(M*v);     // Bjorken x
    double y  = v*M/(k1*p1);      // Inelasticity, y = q*P1/k1*P1
    double W2 = M*M + 2*M*v - Q2; // Hadronic Invariant mass ^ 2
    double W  = TMath::Sqrt(W2);

    TLorentzVector * vtx = event.Vertex(); // Interaction vertex

    bool isDIS = proc_info.IsDeepInelastic();

    int np   = (!isDIS) ? 0 : event.NEntries(kPdgProton,  kIStStableFinalState);
    int nn   = (!isDIS) ? 0 : event.NEntries(kPdgNeutron, kIStStableFinalState);
    int npi0 = (!isDIS) ? 0 : event.NEntries(kPdgPi0,     kIStStableFinalState);
    int npip = (!isDIS) ? 0 : event.NEntries(kPdgPiPlus,  kIStStableFinalState);
    int npim = (!isDIS) ? 0 : event.NEntries(kPdgPiMinus, kIStStableFinalState);

    int iproc = IProc(proc_info);
    int icurr = ICurr(proc_info);
    int ineut = INeut(neutrino->PdgCode());

    if (iproc<0) continue;
    if (icurr<0) continue;
    if (ineut<0) continue;

    // fill histograms
    hEvLAB  [iproc][icurr][ineut] -> Fill ( Ev,  weight );
    hElLAB  [iproc][icurr][ineut] -> Fill ( El,  weight );
    hKinX   [iproc][icurr][ineut] -> Fill ( x,   weight );
    hKinY   [iproc][icurr][ineut] -> Fill ( y,   weight );
    hKinXY  [iproc][icurr][ineut] -> Fill ( x,y, weight );
    hKinNu  [iproc][icurr][ineut] -> Fill ( v,   weight );
    hKinW   [iproc][icurr][ineut] -> Fill ( W,   weight );
    hKinQ2  [iproc][icurr][ineut] -> Fill ( Q2,  weight );
    hKinWQ2 [iproc][icurr][ineut] -> Fill ( W,Q2,weight );

    if(isDIS) {
      hMultPW   [icurr][ineut] -> Fill ( np,   W, weight);
      hMultNW   [icurr][ineut] -> Fill ( nn,   W, weight);
      hMultPi0W [icurr][ineut] -> Fill ( npi0, W, weight);
      hMultPipW [icurr][ineut] -> Fill ( npip, W, weight);
      hMultPimW [icurr][ineut] -> Fill ( npim, W, weight);
    }
    hVtxX -> Fill (vtx->X());
    hVtxY -> Fill (vtx->Y());
    hVtxZ -> Fill (vtx->Z());
    hVtxT -> Fill (vtx->T());
    hVtxR -> Fill (vtx->X(), vtx->Y(), vtx->Z());

    mcrec->Clear();

  } // event loop
}
//___________________________________________________________________
void HistogramLimits(
               Long64_t nmax, TTree * tree, NtpMCEventRecord * mcrec)
{
  LOG("gmctest", pNOTICE) << "*** Computing histogram limits";

  for(Long64_t i = 0; i < nmax; i++) {
    tree->GetEntry(i);

    EventRecord &  event = *(mcrec->event);

    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) continue;
 
    double Ev = event.Particle(0)->Energy();

    TLorentzVector * vtx = event.Vertex();

    gEvMin   = TMath::Min(gEvMin,   Ev);
    gEvMax   = TMath::Max(gEvMax,   Ev);
    gVtxXMin = TMath::Min(gVtxXMin, vtx->X());
    gVtxXMax = TMath::Max(gVtxXMax, vtx->X());
    gVtxYMin = TMath::Min(gVtxYMin, vtx->Y());
    gVtxYMax = TMath::Max(gVtxYMax, vtx->Y());
    gVtxZMin = TMath::Min(gVtxZMin, vtx->Z());
    gVtxZMax = TMath::Max(gVtxZMax, vtx->Z());
    gVtxTMin = TMath::Min(gVtxTMin, vtx->T());
    gVtxTMax = TMath::Max(gVtxTMax, vtx->T());

    mcrec->Clear();
  }

  if(gEvMax-gEvMin < 0.5) { gEvMax+=1.0; gEvMin-=1.0; }
  gEvMin = TMath::Max(gEvMin,0.);

  if(gVtxXMin==gVtxXMax) gVtxXMax+=1;
  if(gVtxYMin==gVtxYMax) gVtxYMax+=1;
  if(gVtxZMin==gVtxZMax) gVtxZMax+=1;
  if(gVtxTMin==gVtxTMax) gVtxTMax+=1;

  gWMin  = 0;
  gWMax  = gEvMax;
  gQ2Min = 0;
  gQ2Max = gEvMax*gEvMax;

  LOG("gmctest", pDEBUG) << "Ev = [" << gEvMin << "," << gEvMax << "]";
  LOG("gmctest", pDEBUG) << "Vx = [" << gVtxXMin << "," << gVtxXMax << "]";
  LOG("gmctest", pDEBUG) << "Vy = [" << gVtxYMin << "," << gVtxYMax << "]";
  LOG("gmctest", pDEBUG) << "Vz = [" << gVtxZMin << "," << gVtxZMax << "]";
  LOG("gmctest", pDEBUG) << "Vt = [" << gVtxTMin << "," << gVtxTMax << "]";
}
//___________________________________________________________________
void InitHistos(void)
{
  LOG("gmctest", pNOTICE) << "*** Booking output histograms";

  int NEv = int ((gEvMax-gEvMin)/kdE);
  int Nx  = int (1./kdx);
  int Ny  = int (1./kdy);
  int NQ2 = int ((gQ2Max-gQ2Min)/kdQ2);
  int NW  = int ((gQ2Max-gQ2Min)/kdW);

  string ptrn;
  string title;

  for(int iproc=0; iproc<kNProc; iproc++) {
    for(int icurr=0; icurr<kNCurr; icurr++) {
      for(int ineut=0; ineut<kNNeut; ineut++) {

	ptrn  = HistoPtr("hEvLab",    iproc,icurr,ineut);
	title = HistoTitle("Ev(LAB)", iproc,icurr,ineut);
	hEvLAB[iproc][icurr][ineut] = new TH1D(
                      ptrn.c_str(),title.c_str(),NEv,gEvMin,gEvMax);

	ptrn  = HistoPtr("hElLab",    iproc,icurr,ineut);
	title = HistoTitle("El (LAB)",iproc,icurr,ineut);
	hElLAB[iproc][icurr][ineut] = new TH1D(
                     ptrn.c_str(), title.c_str(),NEv,gEvMin,gEvMax);

	ptrn  = HistoPtr("hKinX",iproc,icurr,ineut);
	title = HistoTitle("x",  iproc,icurr,ineut);
	hKinX[iproc][icurr][ineut] = new TH1D(
                              ptrn.c_str(), title.c_str(), Nx, 0,1);

	ptrn  = HistoPtr("hKinY",iproc,icurr,ineut);
	title = HistoTitle("y",  iproc,icurr,ineut);
	hKinY[iproc][icurr][ineut] = new TH1D(
                              ptrn.c_str(), title.c_str(), Ny, 0,1);

	ptrn  = HistoPtr("hKinXY",  iproc,icurr,ineut);
	title = HistoTitle("x vs y",iproc,icurr,ineut);
	hKinXY[iproc][icurr][ineut] = new TH2D(
		     ptrn.c_str(), title.c_str(), Nx, 0,1, Ny, 0,1);   

	ptrn  = HistoPtr("hKinQ2",iproc,icurr,ineut);
	title = HistoTitle("Q2",  iproc,icurr,ineut);
	hKinQ2[iproc][icurr][ineut] = new TH1D(
                   ptrn.c_str(), title.c_str(), NQ2, gQ2Min,gQ2Max);

	ptrn  = HistoPtr("hKinW",iproc,icurr,ineut);
	title = HistoTitle("W",  iproc,icurr,ineut);
	hKinW[iproc][icurr][ineut] = new TH1D(
                      ptrn.c_str(), title.c_str(), NW, gWMin,gWMax);

	ptrn  = HistoPtr("hKinNu",iproc,icurr,ineut);
	title = HistoTitle("v",  iproc,icurr,ineut);
	hKinNu[iproc][icurr][ineut] = new TH1D(
                   ptrn.c_str(), title.c_str(), NEv, gEvMin,gEvMax);

	ptrn  = HistoPtr("hKinWQ2",  iproc,icurr,ineut);
	title = HistoTitle("W vs Q2",iproc,icurr,ineut);
	hKinWQ2[iproc][icurr][ineut] = new TH2D(
		     ptrn.c_str(), title.c_str(), 
                               NW, gWMin,gWMax, NQ2, gQ2Min,gQ2Max);   
      }//neut
    }//curr
  }//proc

  for(int icurr=0; icurr<kNCurr; icurr++) {
    for(int ineut=0; ineut<kNNeut; ineut++) {

	ptrn  = HistoPtr("hMultPW", -1,icurr,ineut);
	title = HistoTitle("DIS Proton Multiplicity vs W", -1,icurr,ineut);
        hMultPW[icurr][ineut] = new TH2D(ptrn.c_str(), title.c_str(),
					          5,0,5, NW, gWMin, gWMax);

	ptrn  = HistoPtr("hMultNW", -1,icurr,ineut);
	title = HistoTitle("DIS Neutron Multiplicity vs W", -1,icurr,ineut);
        hMultNW[icurr][ineut] = new TH2D(ptrn.c_str(), title.c_str(),
					           5,0,5, NW, gWMin, gWMax);

	ptrn  = HistoPtr("hMultPi0W", -1,icurr,ineut);
	title = HistoTitle("DIS Pi0 Multiplicity vs W", -1,icurr,ineut);
        hMultPi0W[icurr][ineut] = new TH2D(ptrn.c_str(), title.c_str(),
					         30,0,30, NW, gWMin, gWMax);

	ptrn  = HistoPtr("hMultPipW", -1,icurr,ineut);
	title = HistoTitle("DIS Pi+ Multiplicity vs W", -1,icurr,ineut);
        hMultPipW[icurr][ineut] = new TH2D(ptrn.c_str(), title.c_str(),
					         30,0,30, NW, gWMin, gWMax);

	ptrn  = HistoPtr("hMultPimW", -1,icurr,ineut);
	title = HistoTitle("DIS Pi- Multiplicity vs W", -1,icurr,ineut);
        hMultPimW[icurr][ineut] = new TH2D(ptrn.c_str(), title.c_str(),
					         30,0,30, NW, gWMin, gWMax);
    }//neut
  }//curr

  hVtxX = new TH1D("hVtxX","interaction vtx, x", 
                          100, 0.95*gVtxXMin, 1.05*gVtxXMax);
  hVtxY = new TH1D("hVtxY","interaction vtx, y", 
                          100, 0.95*gVtxYMin, 1.05*gVtxYMax);
  hVtxZ = new TH1D("hVtxZ","interaction vtx, z", 
                          1000, 0.95*gVtxZMin, 1.05*gVtxZMax);
  hVtxT = new TH1D("hVtxT","interaction vtx, t", 
                          100, 0.95*gVtxTMin, 1.05*gVtxTMax);
  hVtxR = new TH3D("hVtxR","interaction vtx, x,y,z", 
		          100, 0.95*gVtxTMin, 1.05*gVtxTMax,
  		          100, 0.95*gVtxYMin, 1.05*gVtxYMax,
                          100, 0.95*gVtxZMin, 1.05*gVtxZMax);
}
//___________________________________________________________________
void SaveHistos(string outpfilename)
{
  LOG("gmctest", pNOTICE) 
          << "*** Saving output histograms to: " << outpfilename;

  TFile outfile(outpfilename.c_str(),"recreate");

  TDirectory * dirEvLab  = 
                   new TDirectory("Ev_Lab","Ev (LAB)");
  TDirectory * dirElLab  = 
                   new TDirectory("El_Lab","El (LAB)");
  TDirectory * dirKinX   = 
                   new TDirectory("x","Bjorken x");
  TDirectory * dirKinY   = 
                   new TDirectory("y","Inelasticity y");
  TDirectory * dirKinXY  = 
                   new TDirectory("xy","Bjorken x vs Inelasticity y");
  TDirectory * dirKinW   = 
                   new TDirectory("W","Hadronic Invariant Mass");
  TDirectory * dirKinNu  = 
                   new TDirectory("v","Energy transfer");
  TDirectory * dirKinQ2  = 
                   new TDirectory("Q2","Momentum Trasfer");
  TDirectory * dirKinWQ2 = 
                   new TDirectory("WQ2","Invariant Mass W vs Momentum Transfer Q2");
  TDirectory * dirVtx    = 
                   new TDirectory("vtx", "Interaction Vertex");
  TDirectory * dirMult   = 
                   new TDirectory("multipl","DIS Hadronic Multiplicities");

  for(int iproc=0; iproc<kNProc; iproc++) {
    for(int icurr=0; icurr<kNCurr; icurr++) {
      for(int ineut=0; ineut<kNNeut; ineut++) {

        dirEvLab -> cd();
	hEvLAB[iproc][icurr][ineut] -> Write();

        dirElLab -> cd();
	hElLAB[iproc][icurr][ineut] -> Write();

        dirKinX  -> cd();
	hKinX[iproc][icurr][ineut]  -> Write();

        dirKinY  -> cd(); 
	hKinY[iproc][icurr][ineut]  -> Write();

        dirKinXY -> cd();
	hKinXY[iproc][icurr][ineut] -> Write();

        dirKinW  -> cd();
	hKinW[iproc][icurr][ineut]  -> Write();

        dirKinNu  -> cd();
	hKinNu[iproc][icurr][ineut]  -> Write();

        dirKinQ2  -> cd(); 
	hKinQ2[iproc][icurr][ineut]  -> Write();

        dirKinWQ2 -> cd();
	hKinWQ2[iproc][icurr][ineut] -> Write();
      }
    }
  }

  dirVtx->cd();
  hVtxX -> Write();
  hVtxY -> Write();
  hVtxZ -> Write();
  hVtxT -> Write();
  hVtxR -> Write();

  dirMult->cd();
  for(int icurr=0; icurr<kNCurr; icurr++) {
    for(int ineut=0; ineut<kNNeut; ineut++) {

       hMultPW   [icurr][ineut] -> Write();
       hMultNW   [icurr][ineut] -> Write();
       hMultPi0W [icurr][ineut] -> Write();
       hMultPipW [icurr][ineut] -> Write();
       hMultPimW [icurr][ineut] -> Write();
    }
  }

  dirEvLab  -> Write();
  dirElLab  -> Write();
  dirKinX   -> Write();
  dirKinY   -> Write();
  dirKinXY  -> Write(); 
  dirKinW   -> Write();
  dirKinNu  -> Write();
  dirKinQ2  -> Write();
  dirKinWQ2 -> Write(); 
  dirVtx    -> Write();
  dirMult  -> Write();

  outfile.Close();
}
//___________________________________________________________________
int IProc(const ProcessInfo & pi)
{
  int iproc = -1;

  if      ( pi.IsQuasiElastic()  ) iproc = 0;
  else if ( pi.IsResonant()      ) iproc = 1;
  else if ( pi.IsDeepInelastic() ) iproc = 2;

  return iproc;
}
//___________________________________________________________________
int ICurr(const ProcessInfo & pi)
{
  int icurr = -1;

  if      ( pi.IsWeakCC() ) icurr = 0;
  else if ( pi.IsWeakNC() ) icurr = 1;

  return icurr;
}
//___________________________________________________________________
int INeut(int pdgc)
{
  int ineut = -1;

  if      ( pdg::IsNuE(pdgc)       ) ineut = 0;
  else if ( pdg::IsNuMu(pdgc)      ) ineut = 1;
  else if ( pdg::IsNuTau(pdgc)     ) ineut = 2;
  else if ( pdg::IsAntiNuE(pdgc)   ) ineut = 3;
  else if ( pdg::IsAntiNuMu(pdgc)  ) ineut = 4;
  else if ( pdg::IsAntiNuTau(pdgc) ) ineut = 5;

  return ineut;
}
//___________________________________________________________________
string ProcName(int i)
{
  if      (i==0) return "QEL";
  else if (i==1) return "RES";
  else if (i==2) return "DIS";
  else
    return "";
}
//___________________________________________________________________
string CurrName(int i)
{
  if      (i==0) return "CC";
  else if (i==1) return "NC";
  else
    return "";
}
//___________________________________________________________________
string NeutName(int i)
{
  if      (i==0) return "NuE";
  else if (i==1) return "NuMu";
  else if (i==2) return "NuTau";
  else if (i==3) return "NuEBar";
  else if (i==4) return "NuMuBar";
  else if (i==5) return "NuTauBar";
  else
    return "";
}
//___________________________________________________________________
string HistoPtr(string front, int iproc, int icurr, int ineut)
{
  ostringstream ptrn;
  ptrn << front << ProcName(iproc) 
       << CurrName(icurr) << NeutName(ineut);
  return ptrn.str();
}
//___________________________________________________________________
string HistoTitle(string front, int iproc, int icurr, int ineut)
{
  ostringstream title;
  title << front << ", " << ProcName(iproc) 
        << ", " << CurrName(icurr) << ", " << NeutName(ineut);
  return title.str();
}
//___________________________________________________________________
string OutputFileName(string inpname)
{
// Builds the output filename based on the name of the input filename
// Perfors the following conversion: name.root -> name.gmctest.root

  unsigned int L = inpname.length();

  // if the last 4 characters are "root" (ROOT file extension) then
  // remove them
  if(inpname.substr(L-4, L).find("root") != string::npos) {
    inpname.erase(L-4, L);
  }
  ostringstream name;
  name << inpname << "gmctest-" << gOptTgtPdgC << ".root";

  return gSystem->BaseName(name.str().c_str());
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gmctest", pNOTICE) << "*** Parsing commad line arguments";

   //number of events:
  try {
    LOG("gmctest", pINFO) << "Reading number of events to analyze";
    gOptN = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctest", pINFO)
              << "Unspecified number of events to analyze - Use all";
      gOptN = -1;
    }
  }

  //target PDG code:
  try {
    LOG("gmctest", pINFO) << "Reading target PDG code";
    gOptTgtPdgC = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctest", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //get input ROOT file (containing a GENIE ER ntuple)
  try {
    LOG("gmctest", pINFO) << "Reading input filename";
    gOptInpFile = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctst", pFATAL) << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //assert that the input ROOT file is accessible
  assert(CheckRootFilename(gOptInpFile));
}
//___________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmctest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gmctest [-n nev] -t tgtpdg -f filename\n";
}
//___________________________________________________________________
bool CheckRootFilename(string filename)
{
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gmctest", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//___________________________________________________________________
