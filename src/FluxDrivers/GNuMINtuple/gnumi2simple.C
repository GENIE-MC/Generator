#include "Numerical/RandomGen.h"
#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GSimpleNtpFlux.h"
#include "Utils/UnitUtils.h"

#include "TSystem.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TFile.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <set>

using namespace std;
using namespace genie;
using namespace genie::flux;

// worker routine
void fill_simple(string fname="fluxntgenie.root", 
                 string cfg="MINOS-Near", 
                 int flxset=0, 
                 long int nentries=0, 
                 double pots=1.0e30,
                 double enumin=0.0, bool doaux=false);


// main routine
void gnumi2simple(string fname, string det, string flux,
                  long int nentries=5000, double pots=1.0e30,
                  double enumin=0.0, bool doaux=false)
{

  cout << "output file: " << fname << endl;
  cout << "request " << nentries << " entries, or " << pots << " POTs" << endl;
  
  string cfg = "MINOS-NearDet";
  if ( det == "minos"    ) cfg = "MINOS-NearDet";
  if ( det == "novaipnd" ) cfg = "NOvA-IPNDShed";
  if ( det == "uboone"   ) cfg = "microboone-numi";

  cout << "det \"" << det << "\" ==> cfg " << cfg << endl;

  int flxset = 0;
  if ( flux == "old"         ) flxset =  0;
  if ( flux == "old1"        ) flxset = 10;
  if ( flux == "old2"        ) flxset = 20;
  if ( flux == "lowcut"      ) flxset =  1;
  if ( flux == "lowcut1"     ) flxset = 11;
  if ( flux == "lowcut2"     ) flxset = 21;
  if ( flux == "flugg"       ) flxset =  2;
  if ( flux == "flugg1"      ) flxset = 12;
  if ( flux == "flugg2"      ) flxset = 22;
  if ( flux == "flugglowth"  ) flxset =  3;
  if ( flux == "flugglowth1" ) flxset = 13;
  if ( flux == "flugglowth2" ) flxset = 23;

  cout << "flux \"" << flux << "\" ==> flxset " << flxset << endl;

  fill_simple(fname,cfg,flxset,nentries,pots,enumin,doaux);
}

void fill_simple(string fname, string cfg, int flxset, 
                 long int nentries, double pots, double enumin, bool doaux)
{
  string fluxpatt = "$FLUXPATH";

  int flxver = flxset % 10;
  switch ( flxver ) {
  case 0:
    fluxpatt += "/v19/fluka05_le010z185i/root/fluka05_le010z185i_";
    break;
  case 1:
    fluxpatt += "/v19/fluka05_le010z185i_lowcut/root/fluka05_le010z185i_";
    break;
  case 2:
    fluxpatt += "/flugg/flugg_le010z185i_run1/root/flugg_le010z185i_run1_";
    break;
  case 3:
    fluxpatt += "/flugg_lowth/flugg_le010z185i_run1_lowth/root/flugg_le010z185i_run1_lowth_";
    break;
  default:
    cout << " not a legal flux set " << endl;
    return;
  }

  if ( flxset/10 == 1 ) {
    if ( flxver == 2 || flxver == 3 ) fluxpatt += "001";
    else               fluxpatt += "1";
  } else if ( flxset/10 == 2 ) {
    if ( flxver == 2 || flxver == 3 ) fluxpatt += "00[12]";
    else               fluxpatt += "[12]";
  } else               fluxpatt += "*";
  fluxpatt += ".root";
    
  string fluxfname(gSystem->ExpandPathName(fluxpatt.c_str()));
  flux::GNuMIFlux* gnumi = new GNuMIFlux();
  gnumi->LoadBeamSimData(fluxfname,cfg);
  gnumi->SetEntryReuse(1); // don't reuse entries when reformatting

  //gnumi->PrintConfig();
  //cout << " change to cm:" << endl;
  //double cm_unit = genie::utils::units::UnitFromString("cm");
  //gnumi->SetLengthUnits(cm_unit);
  //gnumi->PrintConfig();
  //cout << " change to m:" << endl;
  ///////RWH: the following 2 lines should work ...
  //double m_unit = genie::utils::units::UnitFromString("m");
  //gnumi->SetLengthUnits(m_unit);
  gnumi->PrintConfig();
  //double scaleusr = 100.;

  gnumi->SetEntryReuse(1);     // reuse entries
  gnumi->SetUpstreamZ(-3e38);  // leave ray on flux window
  
  GFluxI* fdriver = dynamic_cast<GFluxI*>(gnumi);

  if (nentries == 0) nentries = 100000*20;

  // so as not to include scan time generate 1 nu
  fdriver->GenerateNext();

  TFile* file = TFile::Open(fname.c_str(),"RECREATE");
  TTree* fluxntp = new TTree("flux","a simple flux n-tuple");
  TTree* metantp = new TTree("meta","metadata for flux n-tuple");
  genie::flux::GSimpleNtpEntry* fentry = new genie::flux::GSimpleNtpEntry;
  genie::flux::GSimpleNtpNuMI*  fnumi  = new genie::flux::GSimpleNtpNuMI;
  genie::flux::GSimpleNtpAux*   faux   = new genie::flux::GSimpleNtpAux;
  genie::flux::GSimpleNtpMeta*  fmeta  = new genie::flux::GSimpleNtpMeta;
  fluxntp->Branch("entry",&fentry);
  fluxntp->Branch("numi",&fnumi);
  if (doaux) fluxntp->Branch("aux",&faux);
  metantp->Branch("meta",&fmeta);

  TLorentzVector p4u, x4u;
  //TLorentzVector p4b, x4b; 
  long int ngen = 0;

  set<int> pdglist;
  double maxe = 0;
  double minwgt = +1.0e10;
  double maxwgt = -1.0e10;

  UInt_t metakey = TString::Hash(fname.c_str(),strlen(fname.c_str()));
  cout << "metakey " << metakey << endl;
  // ensure that we don't get smashed by UInt_t vs Int_t
  metakey &= 0x7FFFFFFF;
  cout << "metakey " << metakey << " after 0x7FFFFFFF" << endl;
  cout << "=========================== Start " << endl;

  TStopwatch sw;
  sw.Start();
  int nok = 0, nwrite = 0;
  while ( nok < nentries && gnumi->UsedPOTs() < pots ) {
    //for (int i=0; i<nentries; ++i) {
    fdriver->GenerateNext();
    ngen++;
    fentry->Reset();
    fnumi->Reset();
    faux->Reset();

    fentry->metakey = metakey;
    fentry->pdg     = fdriver->PdgCode();
    fentry->wgt     = gnumi->Weight();
    x4u     = fdriver->Position();
    fentry->vtxx    = x4u.X();
    fentry->vtxy    = x4u.Y();
    fentry->vtxz    = x4u.Z();
    fentry->dist    = gnumi->GetDecayDist();
    p4u     = fdriver->Momentum();
    fentry->px      = p4u.Px();
    fentry->py      = p4u.Py();
    fentry->pz      = p4u.Pz();
    fentry->E       = p4u.E();

    fnumi->run      = gnumi->PassThroughInfo().run;
    fnumi->evtno    = gnumi->PassThroughInfo().evtno;
    fnumi->entryno  = gnumi->GetEntryNumber();

    fnumi->tpx      = gnumi->PassThroughInfo().tpx;
    fnumi->tpy      = gnumi->PassThroughInfo().tpy;
    fnumi->tpz      = gnumi->PassThroughInfo().tpz;
    fnumi->vx       = gnumi->PassThroughInfo().vx;
    fnumi->vy       = gnumi->PassThroughInfo().vy;
    fnumi->vz       = gnumi->PassThroughInfo().vz;

    fnumi->ndecay   = gnumi->PassThroughInfo().ndecay;
    fnumi->ppmedium = gnumi->PassThroughInfo().ppmedium;
    fnumi->tptype   = gnumi->PassThroughInfo().tptype;


    if ( doaux ) {
      faux->auxint.push_back(gnumi->PassThroughInfo().tptype);
      faux->auxdbl.push_back(gnumi->PassThroughInfo().tpx);
      faux->auxdbl.push_back(gnumi->PassThroughInfo().tpy);
      faux->auxdbl.push_back(gnumi->PassThroughInfo().tpz);
    }

    if ( fentry->E > enumin ) ++nok;
    fluxntp->Fill();
    ++nwrite;

    // accumulate meta data
    pdglist.insert(fentry->pdg);
    minwgt = TMath::Min(minwgt,fentry->wgt);
    maxwgt = TMath::Max(maxwgt,fentry->wgt);
    maxe   = TMath::Max(maxe,fentry->E);

  }
  cout << "=========================== Complete " << endl;

  //fmeta->nflavors  = pdglist.size();
  //if ( fmeta->nflavors > 0 ) {
  //  fmeta->flavor = new Int_t[fmeta->nflavors];
  //  set<int>::const_iterator itr = pdglist.begin();
  //  int indx = 0;
  //  for ( ; itr != pdglist.end(); ++itr, ++indx )  fmeta->flavor[indx] = *itr;
  //}
  fmeta->pdglist.clear();
  set<int>::const_iterator setitr = pdglist.begin();
  for ( ; setitr != pdglist.end(); ++setitr)  fmeta->pdglist.push_back(*setitr);
 
  fmeta->maxEnergy = maxe;
  fmeta->minWgt    = minwgt;
  fmeta->maxWgt    = maxwgt;
  fmeta->protons   = gnumi->UsedPOTs();
  TVector3 p0, p1, p2;
  gnumi->GetFluxWindow(p0,p1,p2);
  TVector3 d1 = p1 - p0;
  TVector3 d2 = p2 - p0;
  fmeta->windowBase[0] = p0.X();
  fmeta->windowBase[1] = p0.Y();
  fmeta->windowBase[2] = p0.Z();
  fmeta->windowDir1[0] = d1.X();
  fmeta->windowDir1[1] = d1.Y();
  fmeta->windowDir1[2] = d1.Z();
  fmeta->windowDir2[0] = d2.X();
  fmeta->windowDir2[1] = d2.Y();
  fmeta->windowDir2[2] = d2.Z();
  if ( doaux ) {
    fmeta->auxintname.push_back("tptype");
    fmeta->auxdblname.push_back("tpx");
    fmeta->auxdblname.push_back("tpy");
    fmeta->auxdblname.push_back("tpz");
  }
  string fname = "flux_pattern=" + fluxfname;
  fmeta->infiles.push_back(fluxfname); // should get this expanded from gnumi
  vector<string> flist = gnumi->GetFileList();
  for (size_t i = 0; i < flist.size(); ++i) {
    fname = "flux_infile=" + flist[i];
    fmeta->infiles.push_back(fname);
  }

  fmeta->seed    = RandomGen::Instance()->GetSeed();
  fmeta->metakey = metakey;

  metantp->Fill();

  sw.Stop();
  cout << "Generated " << nwrite << " (" << nok << " w/ E >" << enumin << " ) "
       << " ( request " << nentries << " ) "
       << endl
       << gnumi->UsedPOTs() << " POTs " << " ( request " << pots << " )"
       << endl
       << " pulled NFluxNeutrinos " << gnumi->NFluxNeutrinos()
       << endl
       << "Time to generate: " << endl;
  sw.Print();
  
  cout << "===================================================" << endl;

  cout << "Last GSimpleNtpEntry: " << endl
       << *fentry
       << endl
       << "Last GSimpleNtpMeta: " << endl
       << *fmeta
       << endl;

  file->cd();

  // write meta as simple object at top of file
  //fmeta->Write();

  // write ntuples out
  fluxntp->Write();
  metantp->Write();
  file->Close();

  cout << endl << endl;
  gnumi->PrintConfig();
}

/////////////////////// old interface /////////////////////////////
void make_simple_grid(string fname, string det, string flux,
                      long int nentries=5000, double pots=1.0e30,
                      double enumin=0.0, bool doaux=false)
{
  gnumi2simple(fname,det,flux,nentries,pots,enumin,doaux);
}
