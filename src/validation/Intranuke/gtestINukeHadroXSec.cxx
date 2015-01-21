//____________________________________________________________________________
/*!

\program gtestINukeHadroXSec

\brief   INTRANUKE test program. Reads-in a hadron-nucleus event file in GHEP
         format and prints-out hadron-nucleus cross sections.

         Syntax :
           gtestINukeHadroXSec -f input_ghep_file [-w]

         -f : Input file
         -w : If set, writes computed hadron cross sections to a text file

\authors Aaron Meyer and Steve Dytman

\created July 26, 2010

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/CmdLnArgParser.h"

using std::endl;
using std::setw;
using std::setfill;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

INukeFateHA_t FindhAFate(const GHepRecord * evrec);

// command line options
string gOptInpFilename = "";    ///< input event file
bool   gOptWriteOutput = false; ///< write out hadron cross sections
string gOptOutputFilename = "gevgen_hadron_xsection.txt";

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // parse command line arguments
  GetCommandLineArgs(argc,argv);

  //
  // initialize
  //

  const int nfates = 9;             // total number of possible fates
  int           countfate [nfates]; // total no. of events with given fate
  double        sigma     [nfates]; // cross sections
  double        sigma_err [nfates]; // cross section errors
  string        fatestr   [nfates] = "  ";
  INukeFateHA_t fatetype  [nfates];

  fatetype[0] = kIHAFtUndefined;   
  fatetype[1] = kIHAFtNoInteraction;
  fatetype[2] = kIHAFtCEx;
  fatetype[3] = kIHAFtElas;
  fatetype[4] = kIHAFtInelas;
  fatetype[5] = kIHAFtAbs;
  fatetype[6] = kIHAFtKo;
  fatetype[7] = kIHAFtPiProd;
  fatetype[8] = kIHAFtDCEx;
 
  for (int k=0; k<nfates; k++) {
    countfate[k] = 0; 
    sigma    [k] = 0.; 
    sigma_err[k] = 0.;
    fatestr  [k] = INukeHadroFates::AsString(fatetype[k]);
  }

  // event sample info (to be extracted from 1st event)
  int    nev        = 0;
  int    probe_pdg  = 0;
  int    target_pdg = 0;
  int    displayno  = 100;
  double kin_energy = 0.;

  //
  // open the input ROOT file and get the event tree
  //

  TTree *           tree   = 0;
  TTree *           ginuke = 0;
  NtpMCTreeHeader * thdr   = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree   = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  ginuke = dynamic_cast <TTree *>           ( file.Get("ginuke") );
  thdr   = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  /*if(!tree) {
    LOG("gtestINukeHadroXSec", pERROR) 
         << "No event tree found in input file!";
    return 1;
    }*/

  if (tree) {

    NtpMCEventRecord * mcrec = 0;
    tree->SetBranchAddress("gmcrec", &mcrec);

    //    int nev = (int) tree->GetEntries();
    nev = (int) tree->GetEntries();
    LOG("gtestINukeHadroXSec", pNOTICE) 
      << "Processing " << nev << " events";

    //
    // event loop
    //

    for(int ievent = 0; ievent < nev; ievent++) {

      // get next tree entry
      tree->GetEntry(ievent);

      // get the corresponding GENIE event
      EventRecord & event = *(mcrec->event);

      // extract info for the event sample
      if(ievent==0) {
	kin_energy = event.Particle(0)->KinE();
	probe_pdg  = event.Particle(0)->Pdg();
	target_pdg = event.Particle(1)->Pdg();
      }

      // analyze
      const GHepRecord * grec = dynamic_cast<const GHepRecord *> (&event);
      INukeFateHA_t fate = FindhAFate(grec);
      if(ievent<displayno) {
	LOG("gtestINukeHadroXSec", pNOTICE) 
          << "fate = " << INukeHadroFates::AsString(fate);
      }

      // We don't want the specific fate data, just the main (9) fate types
      switch (fate){
      case 0:   countfate[0]++; break;
      case 1:   countfate[1]++; break;
      case 2:   countfate[2]++; break;
      case 3:   countfate[3]++; break;
      case 4:   countfate[4]++; break;
      case 5:   countfate[5]++; break;
      case 6:   countfate[6]++; break;
      case 13:  countfate[8]++; break;
      default:  
	if (7<=fate && fate<=12) countfate[7]++;
	else {
	  LOG("gtestINukeHadroXSec", pWARN) 
             << "Undefined fate from FindhAFate() : " << fate;
        }
	break;
      }

      // clear current mc event record
      mcrec->Clear();

    } // end event loop 
  } // end if (tree)
  else if ( ginuke ) {
    // possibly a ginuke file

    LOG("gtestINukeHadroXSec", pNOTICE) 
      << "Found ginuke type file";

    nev = (int) ginuke->GetEntries();
    LOG("gtestINukeHadroXSec", pNOTICE) 
      << "Processing " << nev << " events";

    int     kmax    = 250;
    int     index   = 0;
    int     numh    = 0;
    int     numpip  = 0;
    int     numpi0  = 0;
    int     numpim  = 0;
    int     pdg_had[kmax];
    double  E_had  [kmax];
    double  energy  = 0.0;

    ginuke->SetBranchAddress("ke",   &kin_energy);
    ginuke->SetBranchAddress("probe",&probe_pdg );
    ginuke->SetBranchAddress("tgt",  &target_pdg);
    ginuke->SetBranchAddress("nh",   &numh      );
    ginuke->SetBranchAddress("npip", &numpip    );
    ginuke->SetBranchAddress("npi0", &numpi0    );
    ginuke->SetBranchAddress("npim", &numpim    );
    ginuke->SetBranchAddress("pdgh", &pdg_had   );
    ginuke->SetBranchAddress("Eh",   &E_had     );
    ginuke->SetBranchAddress("e",    &energy    );


    for(int ievent = 0; ievent < nev; ievent++) {

      // get next tree entry
      ginuke->GetEntry(ievent);

      // Determine fates (as defined in Intranuke/INukeUtils.cxx/ utils::intranuke::FindhAFate())
           if (energy==E_had[0] && numh==1) // No interaction
	{ index=1; }
      else if (energy!=E_had[0] && numh==1) // Elastic
	{ index=3; }
      else if ( pdg::IsPion(probe_pdg) && numpip+numpi0+numpim==0) // Absorption
	{ index=5; }
      else if ( (pdg::IsNucleon(probe_pdg) && numpip+numpi0+numpim==0 && numh>2 )
		|| (probe_pdg==kPdgGamma && energy!=E_had[0] && numpip+numpi0+numpim==0)) // Knock-out
	{ index=6; }
      else if ( numpip+numpi0+numpim> (pdg::IsPion(probe_pdg) ? 1 : 0) ) // Pion production
	{ index=7; }
      else if ( numh==2 ) // Inelastic or Charge Exchange
	{
	  if ( (pdg::IsPion(probe_pdg) && ( probe_pdg==pdg_had[0] || probe_pdg==pdg_had[1] ))
	       || pdg::IsNucleon(probe_pdg) ) index=4;
	  else index=2;
	}
      else //Double Charge Exchange or Undefined
	{
	  bool undef = true;
	  if ( pdg::IsPion(probe_pdg) )
	    {
	      for (int iter = 0; iter < numh; iter++)
		{
		  if      (probe_pdg==211 && pdg_had[iter]==-211) { index=8; undef=false; }
		  else if (probe_pdg==-211 && pdg_had[iter]==211) { index=8; undef=false; }
		}
	    }
	  if (undef) { index=0; }
	}
      countfate[index]++;
      if (ievent<displayno) {
	LOG("gtestINukeHadroXSec", pNOTICE) 
          << "fate = " << INukeHadroFates::AsString(fatetype[index]);
      } 

    }
  } // end if (ginuke)
  else {
    LOG("gtestINukeHadroXSec", pERROR) 
         << "Could not read input file!";
    return 1;
  }
  
  //
  // output section
  //
 
  const double fm2tomb  = units::fm2 / units::mb;
  const double dnev     = (double) nev;
  const int    NR       = 3;
  const double R0       = 1.4;

  int    A              = pdg::IonPdgCodeToA(target_pdg);
  int    Z              = pdg::IonPdgCodeToZ(target_pdg);
  double nuclear_radius = NR * R0 * TMath::Power(A, 1./3.); // fm
  double area           = TMath::Pi() * TMath::Power(nuclear_radius,2);

  PDGLibrary * pdglib = PDGLibrary::Instance();
  string probe_name  = pdglib->Find(probe_pdg)->GetName();
  string target_name = pdglib->Find(target_pdg)->GetName();

  LOG("gtestINukeHadroXSec", pNOTICE) 
     << " Probe = " << probe_name 
     << ", KE = " << kin_energy << " GeV";
  LOG("gtestINukeHadroXSec", pNOTICE) 
    << " Target = " << target_name 
    << " (Z,A) = (" << Z << ", " << A 
    << "), nuclear radius = " << nuclear_radius 
    << " fm, area = " << area << " fm**2 " << '\n';

  int    cnttot     = 0;
  int    nullint    = countfate[1]; // no interactions
  double sigtot     = 0;
  double sigtoterr  = 0;
  double sigtotScat = 0;
  double sigtotAbs  = 0;

  for(int k=0; k<nfates; k++) {
    if(k!=1) {
      cnttot += countfate[k];
      double ratio = countfate[k]/dnev;
      sigma[k]     = fm2tomb * area * ratio;
      sigma_err[k] = fm2tomb * area * TMath::Sqrt(ratio*(1-ratio)/dnev);
      if(sigma_err[k]==0) {
         sigma_err[k] = fm2tomb * area * TMath::Sqrt(countfate[k])/dnev;
      }
      if(countfate[k]>0) {
        LOG("gtestINukeHadroXSec", pNOTICE) 
            << " --> " << setw(26) << fatestr[k] 
            << ": " << setw(7) << countfate[k] << " events -> " 
            << setw(7) << sigma[k] << " +- " << sigma_err[k] << " (mb)" << '\n';
      }
      if(k==5) {
        sigtotAbs += sigma[k];
      }
      else 
      if (k!=1) { 
        sigtotScat += sigma[k];
      }
    }//k!=1
  }//k koop

  sigtot    = fm2tomb * area * cnttot/dnev;
  sigtoterr = fm2tomb * area * TMath::Sqrt(cnttot)/dnev;

  double sigtot_noelas    = fm2tomb * area * (cnttot-countfate[3])/dnev;
  double sigtoterr_noelas = fm2tomb * area * TMath::Sqrt(cnttot-countfate[3])/dnev;

  double ratio_as = (sigtotScat==0) ? 0 : sigtotAbs/(double)sigtotScat;
  
  LOG("gtestINukeHadroXSec", pNOTICE) 
    << "\n\n --------------------------------------------------- " 
    << "\n ==> " << setw(28) << " Total: " << setw(7) << cnttot 
    << " events -> " << setw(7) << sigtot << " +- " << sigtoterr << " (mb)"
    << "\n (-> " << setw(28) << " Hadrons escaped nucleus: " 
    << setw(7) << nullint << " events)"
    << "\n ==> " << setw(28) << " Ratio (abs/scat) = " 
    << setw(7) << ratio_as
    << "\n ==> " << setw(28) << " avg. num of int. = " 
    << setw(7) << cnttot/dnev
    << "\n ==> " << setw(28) << " no interaction   = " 
    << setw(7) << (dnev-cnttot)/dnev
    << "\n ------------------------------------------------------- \n";

  if(gOptWriteOutput) 
  {
    ifstream test_file;
    bool file_exists=false;
    test_file.open(gOptOutputFilename.c_str(), std::ifstream::in);
    file_exists=test_file.is_open();
    test_file.close();
    ofstream xsec_file; 
    xsec_file.open(gOptOutputFilename.c_str(), std::ios::app);
    if (!file_exists)
      {
	xsec_file << "#KE" << "\t" << "Undef" << "\t"
		  << "sig" << "\t" << "CEx"   << "\t"
		  << "sig" << "\t" << "Elas"  << "\t"
		  << "sig" << "\t" << "Inelas"<< "\t"
		  << "sig" << "\t" << "Abs"   << "\t"
		  << "sig" << "\t" << "KO"    << "\t"
		  << "sig" << "\t" << "PiPro" << "\t"
		  << "sig" << "\t" << "DCEx"  << "\t"
		  << "sig" << "\t" << "Reac"  << "\t"
		  << "sig" << "\t" << "Tot"   << "\t" << "sig" << endl;
      }
    xsec_file << kin_energy;
    for(int k=0; k<nfates; k++) {
      if (k==1) continue;
       xsec_file << "\t" << sigma[k] << "\t" << sigma_err[k];
    }
    xsec_file << "\t" << sigtot_noelas << "\t" << sigtoterr_noelas;
    xsec_file << "\t" << sigtot        << "\t" << sigtoterr << endl;
    xsec_file.close();
  }

  return 0;
}
//____________________________________________________________________________
INukeFateHA_t FindhAFate(const GHepRecord * evrec)
{
  // Determine the fate of an hA event
  // Works for ghAevgen or gntpc
  // author:        S. Dytman  -- July 30, 2007

  double p_KE  = evrec->Probe()->KinE();
  double p_pdg = evrec->Probe()->Pdg();

  // particle codes
  int numtype[] = {kPdgProton, kPdgNeutron, kPdgPiP, kPdgPiM, kPdgPi0, kPdgKP, kPdgKM, kPdgK0, kPdgGamma};
  // num of particle for numtype
  int num[]  = {0,0,0,0,0,0,0,0,0};
  int num_t  = 0;
  int num_nu = 0;
  int num_pi = 0;
  int num_k  = 0;
  // max KE for numtype
  double numKE[] = {0,0,0,0,0,0,0,0,0};

  GHepStatus_t status = kIStUndefined;

  bool hasBlob = false;
  int numFsPart = 0;

  int index = 0;
  TObjArrayIter piter(evrec);
  GHepParticle * p     = 0;
  GHepParticle * fs    = 0;
  GHepParticle * probe = evrec->Probe();
  while((p=(GHepParticle *) piter.Next()))
  {
    status=p->Status();
    if(status==kIStStableFinalState)
    {
      switch((int) p->Pdg()) 
      {
        case ((int) kPdgProton)  : index = 0; break;
        case ((int) kPdgNeutron) : index = 1; break;
        case ((int) kPdgPiP)     : index = 2; break;
        case ((int) kPdgPiM)     : index = 3; break;
        case ((int) kPdgPi0)     : index = 4; break;
        case ((int) kPdgKP)      : index = 5; break;
        case ((int) kPdgKM)      : index = 6; break;
        case ((int) kPdgK0)      : index = 7; break;
        case ((int) kPdgGamma)   : index = 8; break;
        case (2000000002)        : index = 9; hasBlob=true; break;
                          default: index = 9; break;
      }

      if(index!=9)
      {
        if(numFsPart==0) fs=p;
        numFsPart++;
        num[index]++;
        if(p->KinE() > numKE[index]) numKE[index] = p->KinE();
      }
    }
  }

  if(numFsPart==1)
  {
    double dE  = TMath::Abs( probe-> E() - fs-> E() );
    double dPz = TMath::Abs( probe->Pz() - fs->Pz() );
    double dPy = TMath::Abs( probe->Py() - fs->Py() );
    double dPx = TMath::Abs( probe->Px() - fs->Px() );

    if (dE < 1e-15 && dPz < 1e-15 && dPy < 1e-15 && dPx < 1e-15) return kIHAFtNoInteraction;
  }

  num_t  = num[0]+num[1]+num[2]+num[3]+num[4]+num[5]+num[6]+num[7];
  num_nu = num[0]+num[1];
  num_pi =               num[2]+num[3]+num[4];
  num_k  =                                    num[5]+num[6]+num[7];

  if(num_pi>((p_pdg==kPdgPiP || p_pdg==kPdgPiM || p_pdg==kPdgPi0)?(1):(0)))
  {
    /*    if(num[3]==10 && num[4]==0) return kIHAFtNPip;   //fix later
    else if(num[4]==10) return kIHAFtNPipPi0;        //fix later
    else if(num[4]>0) return kIHAFtInclPi0;
    else if(num[2]>0) return kIHAFtInclPip;
    else if(num[3]>0) return kIHAFtInclPim;
    else */
    return kIHAFtPiProd;
  }
  else if(num_pi<((p_pdg==kPdgPiP || p_pdg==kPdgPiM || p_pdg==kPdgPi0)?(1):(0)))
  {
    if     (num[0]==1 && num[1]==1) return kIHAFtAbs;
    else if(num[0]==2 && num[1]==0) return kIHAFtAbs;
    else if(num[0]==2 && num[1]==1) return kIHAFtAbs;
    else if(num[0]==1 && num[1]==2) return kIHAFtAbs;
    else if(num[0]==2 && num[1]==2) return kIHAFtAbs;
    else if(num[0]==3 && num[1]==2) return kIHAFtAbs;
    else return kIHAFtAbs;
  }
  else if(num_k<((p_pdg==kPdgKP || p_pdg==kPdgKM || p_pdg==kPdgK0)?(1):(0)))
  {
    return kIHAFtAbs;
  }  
  else
  {
    if(p_pdg==kPdgPiP || p_pdg==kPdgPiM || p_pdg==kPdgPi0
       || p_pdg==kPdgKP|| p_pdg==kPdgKM|| p_pdg==kPdgK0)
    {
      int fs_pdg, fs_ind;
      if     (num[2]==1) { fs_pdg=kPdgPiP; fs_ind=2; }
      else if(num[3]==1) { fs_pdg=kPdgPiM; fs_ind=3; }
      else if(num[4]==1) { fs_pdg=kPdgPi0; fs_ind=4; }
      else if(num[5]==1) { fs_pdg=kPdgKP; fs_ind=5; }
      else if(num[6]==1) { fs_pdg=kPdgKM; fs_ind=6; }
      else               { fs_pdg=kPdgK0; fs_ind=7; }
 
      if(p_pdg==fs_pdg)
      {
	if(num_nu==0) return kIHAFtElas;
	else return kIHAFtInelas;
      }
      else if(((p_pdg==kPdgPiP || p_pdg==kPdgPiM) && fs_ind==4) ||
              ((fs_ind==2 || fs_ind==3) && p_pdg==kPdgPi0))
      {
        return kIHAFtCEx;
      }
      else if(((p_pdg==kPdgKP || p_pdg==kPdgKM) && fs_ind==7) ||
              ((fs_ind==5 || fs_ind==6) && p_pdg==kPdgK0))
      {
        return kIHAFtCEx;
      }
      else if((p_pdg==kPdgPiP && fs_ind==3) ||
              (p_pdg==kPdgPiM &&fs_ind==2))
      {
        return kIHAFtDCEx;
      }
      else if((p_pdg==kPdgKP && fs_ind==6) ||
              (p_pdg==kPdgKM &&fs_ind==5))
      {
        return kIHAFtDCEx;
      }
    }
    else if(p_pdg==kPdgProton || p_pdg==kPdgNeutron)
    {
      int fs_ind;
      if(num[0]>=1) { fs_ind=0; }
      else          { fs_ind=1; }

      if(num_nu==1)
      {
        if(numtype[fs_ind]==p_pdg) return kIHAFtElas;
        else return kIHAFtUndefined;
      }
      else if(num_nu==2)
      {
        if(numKE[1]>numKE[0]) { fs_ind=1; }  
        
        if(numtype[fs_ind]==p_pdg)
        {
          //if(numKE[fs_ind]>=(.8*p_KE))
          //{
          //  if(num[0]==1 && num[1]==1) return kIHAFtKo;
          //  else if(num[0]==2) return kIHAFtKo;
	  //  else return kIHAFtKo;
          //}
          //else
          return kIHAFtInelas; //fix later
        }
        else
        {
          //if(numKE[fs_ind]>=(.8*p_KE)) return kIHAFtInelas;
          //else
          //{
          //  if(num[fs_ind]==2)
          //  {
          //    if(num[0]==2) return kIHAFtKo;
          //    else return kIHAFtKo;
          //  }
          //  else return kIHAFtInelas;
	  // }
	  return kIHAFtInelas; //fix later
        }
      }
      else if(num_nu>2)
      {
        if     (num[0]==2 && num[1]==1) return kIHAFtKo;
        else if(num[0]==1 && num[1]==2) return kIHAFtKo;
        else if(num[0]==2 && num[1]==2) return kIHAFtKo;
        else if(num[0]==3 && num[1]==2) return kIHAFtKo;
        else return kIHAFtKo;
      }
    }
    else if (p_pdg==kPdgKP || p_pdg==kPdgKM || p_pdg==kPdgK0)
    {
      int fs_ind;

      if (num[5]==1) fs_ind=5;
      else if (num[6]==1) fs_ind=6;
      else fs_ind=7; // num[7]==1

      if(numKE[fs_ind]>=(.8*p_KE)) return kIHAFtElas;
      else return kIHAFtInelas;
    }
    else if (p_pdg==kPdgGamma)
    {
      if     (num[0]==2 && num[1]==1) return kIHAFtKo;
      else if(num[0]==1 && num[1]==2) return kIHAFtKo;
      else if(num[0]==2 && num[1]==2) return kIHAFtKo;
      else if(num[0]==3 && num[1]==2) return kIHAFtKo;
      else if(num_nu < 1)             return kIHAFtUndefined;
      else                            return kIHAFtKo;
    }
  }

  LOG("Intranuke",pWARN) << "---> *** Undefined fate! ***" << "\n" << (*evrec);
  return kIHAFtUndefined;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gtestINukeHadroXSec", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get input ROOT file (containing a GENIE GHEP event tree)
  if( parser.OptionExists('f') ) {
    LOG("gtestINukeHadroXSec", pINFO) << "Reading input filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("gtestINukeHadroXSec", pFATAL)
       << "Unspecified input filename - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }
  if( parser.OptionExists('o') ) {
    LOG("gtestINukeHadroXSec", pINFO) << "Reading output filename";
    gOptOutputFilename = parser.ArgAsString('o');
  }

  // write-out events?
  gOptWriteOutput =  parser.OptionExists('w');
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gtestINukeHadroXSec", pNOTICE)
    << "\n\n" 
    << "Syntax:" << "\n"
    << "  gtestINukeHadroXSec -f event_file [-w]"
    << "\n";
}
//____________________________________________________________________________
