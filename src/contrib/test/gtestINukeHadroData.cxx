//____________________________________________________________________________
/*!

\program testINukeHadroData

\brief   test program used for testing the hadron cross section data used in
         INTRANUKE. 
         The program will save all the hadron cross section splines at an output 
         ROOT file (hadxs.root).
         If a hadron kinetic energy is passed as command line argument then,
         rather than saving the cross section splines at an output file, the
         program will print out the cross section values at the specified
         kinetic energy.

\syntax  gtestINukeHadroData [-e energy(MeV)]

         []  denotes an optionan argument
         -e  can be used to pass a kinetic energy for which all the hadron
             cross sections will be printed out

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 12, 2004

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE

Important Revisions:
@ Aug 17, 2010 - AM, SD
Update to print out xs for p and n separately to match what is used in
propagation.  Add photon and kaon xs ouputs.

*/
//____________________________________________________________________________

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <fstream>

#include "Conventions/Constants.h"
#include "HadronTransport/INukeHadroData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

#include <TSystem.h>
#include <TFile.h>
#include <TGraph2D.h>

using std::ostringstream;
using std::istream;
using std::ios;
using std::endl;

using namespace genie;
using namespace genie::constants;

void  SaveDataToRootFile    (bool test_intbounce);
void  SaveSplinesToRootFile (void);
void  PrintOutForInputKE    (double ke);

int main(int argc, char ** argv)
{
  double ke = -1; // input hadron kinetic energy
  bool save_data=true;

  CmdLnArgParser parser(argc,argv);

  // neutrino energy
  if( parser.OptionExists('e') ) {
    LOG("testINukeHadroData", pINFO) << "Reading Event Energy";
    ke = parser.ArgAsLong('e');
  } else {
    LOG("testINukeHadroData", pINFO) << "Unspecified energy, write to file";
  }

  if(ke<0) {
    SaveSplinesToRootFile();
    if(save_data) SaveDataToRootFile(true);
  }
  else     {
    PrintOutForInputKE(ke);
  }

  return 0;
}
//____________________________________________________________________________
void SaveDataToRootFile(bool test_intbounce)
{

  string filename = "hadxs.root";

  string data_dir = (gSystem->Getenv("GINUKEHADRONDATA")) ?
             string(gSystem->Getenv("GINUKEHADRONDATA")) :
             string(gSystem->Getenv("GENIE")) + string("/data/intranuke/");

  LOG("testINukeHadroData", pNOTICE)
	<< "Saving INTRANUKE hadron x-section data to file: " << filename;

  TFile f(filename.c_str(), "UPDATE");

  INukeHadroData * inukd = INukeHadroData::Instance();
  int    numPoints = 80;
  double cstep = 2.0 / (numPoints);

  // kIHNFtElas, kpn :
  {
    const int hN_kpNelas_nfiles = 18;
    const int hN_kpNelas_points_per_file = 37;
    const int hN_kpNelas_npoints = hN_kpNelas_points_per_file * hN_kpNelas_nfiles;

    int hN_kpNelas_energies[hN_kpNelas_nfiles] = {  
     100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800
    };

    TGraph2D * hN_kpNelas_graph = new TGraph2D(hN_kpNelas_npoints);
    TGraph2D * hN_kpNelas_int = new TGraph2D(hN_kpNelas_nfiles*numPoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_kpNelas_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_kpNelas_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/kpn/kpn" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_kpNelas_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff;
    }//loop over files

    //hN_kpNelas_graph -> Write("d_kpN_elas_data");
    //if (test_intbounce) {
    //  hN_kpNelas_int -> Write("d_kpN_elas_int");
    //}
    //delete hN_kpNelas_graph;
    //delete hN_kpNelas_int;
    //TGraph2D * hN_kpNelas_fudge = new TGraph2D(hN_kpNelas_nfiles*numPoints);

    if (test_intbounce)
	{
	  for (int ifile=0;ifile<hN_kpNelas_nfiles;ifile++) {
	    int subdiv = numPoints*ifile;
	    int ke = hN_kpNelas_energies[ifile];
	    for(int j=0;j<numPoints; j++)
	      {
		hN_kpNelas_int -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
					   inukd->XSec(kPdgKP,kPdgNeutron,kPdgKP,kIHNFtElas,double(ke)+.05,-1+j*cstep));
	      }
	  }
	}
    hN_kpNelas_graph -> Write("d_kpN_elas_data");
    if (test_intbounce) {
      hN_kpNelas_int -> Write("d_kpN_elas_int");
    }
    delete hN_kpNelas_graph;
    delete hN_kpNelas_int;
  }

  // kIHNFtElas, kpp :
  {
    const int hN_kpPelas_nfiles = 18;
    const int hN_kpPelas_points_per_file = 37;
    const int hN_kpPelas_npoints = hN_kpPelas_points_per_file * hN_kpPelas_nfiles;

    int hN_kpPelas_energies[hN_kpPelas_nfiles] = {  
     100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800
    };

    TGraph2D * hN_kpPelas_graph = new TGraph2D(hN_kpPelas_npoints);
    //TGraph2D * hN_kpPelas_fudge = new TGraph2D(hN_kpPelas_nfiles*numPoints);
    TGraph2D * hN_kpPelas_int = new TGraph2D(hN_kpPelas_nfiles*numPoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_kpPelas_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_kpPelas_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/kpp/kpp" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_kpPelas_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 

      if (test_intbounce) 
	{
	  int subdiv = numPoints*ifile;
	  //double * buffpt=0x0;
	  for(int j=0;j<numPoints; j++)
	    {
	      /*if (ifile==6||ifile==7) {/ *hN_kpPelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,0);
		else * /hN_kpPelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
		inukd->XSec(kPdgKP,kPdgProton,kPdgKP,kIHNFtElas,double(ke)+.05,-1+j*cstep));}*/
	      hN_kpPelas_int -> SetPoint(subdiv+j,double(ke),-1.+double(j)*cstep,
		inukd->XSec(kPdgKP,kPdgProton,kPdgKP,kIHNFtElas,double(ke),-1+j*cstep));
	      //if ((ifile>=3&&ifile<=5)&&(j<20||j>60)) { buffpt=hN_kpPelas_int -> GetZ();
	      //LOG("testINukeHadroData",pWARN)<<subdiv+j<<" "<<buffpt[subdiv+j]<<" "<<inukd->XSec(kPdgKP,kPdgProton,kPdgKP,kIHNFtElas,double(ke),-1+j*cstep);}
	      //hN_kpPelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,0);
	      /*if (j%3==0&&ifile%2==0) hN_kpPelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
		inukd->XSec(kPdgKP,kPdgProton,kPdgKP,kIHNFtElas,double(ke)+.05,-1+j*cstep));
		else hN_kpPelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,0);*/
	      /*if (ifile==3||ifile==4) hN_kpPelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,0);
	      else hN_kpPelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
	      inukd->XSec(kPdgKP,kPdgProton,kPdgKP,kIHNFtElas,double(ke)+.05,-1+j*cstep));*/

	    }
	}
    }//loop over files

    hN_kpPelas_graph -> Write("d_kpP_elas_data");
    if (test_intbounce)
      {
	hN_kpPelas_int -> Write("d_kpP_elas_int");
	//hN_kpPelas_fudge -> Write("d_kpP_elas_fudge");
      }
    delete hN_kpPelas_graph;
    delete hN_kpPelas_int;
    //delete hN_kpPelas_fudge;
  }

  // kIHNFtElas, pp&nn :
  {
    const int hN_ppelas_nfiles = 20;
    const int hN_ppelas_points_per_file = 21;
    const int hN_ppelas_npoints = hN_ppelas_points_per_file * hN_ppelas_nfiles;

    int hN_ppelas_energies[hN_ppelas_nfiles] = {
        50, 100, 150, 200, 250, 300, 350, 400, 450, 500,
       550, 600, 650, 700, 750, 800, 850, 900, 950, 1000
    };

    TGraph2D * hN_ppelas_graph = new TGraph2D(hN_ppelas_npoints);
    TGraph2D * hN_ppelas_int = new TGraph2D(hN_ppelas_nfiles*numPoints);
    //TGraph2D * hN_ppelas_fudge = new TGraph2D(hN_ppelas_nfiles*numPoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_ppelas_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_ppelas_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/pp/pp" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_ppelas_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 

      if (test_intbounce) 
	{
	  int subdiv = numPoints*ifile;
	  for(int j=0;j<numPoints; j++)
	    {
	      hN_ppelas_int -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
					 inukd->XSec(kPdgProton,kPdgProton,kPdgProton,kIHNFtElas,double(ke),-1+j*cstep));
	    }
	}
    }//loop over files

    hN_ppelas_graph -> Write("d_pp_elas_data");
    if (test_intbounce)
      {
	hN_ppelas_int -> Write("d_pp_elas_int");
	//hN_ppelas_fudge -> Write("d_pp_elas_fudge");
      }
    delete hN_ppelas_graph;
    delete hN_ppelas_int;
    //delete hN_ppelas_fudge;
  }

  // kIHNFtElas, pn&np :
  {
    const int hN_npelas_nfiles = 20;
    const int hN_npelas_points_per_file = 21;
    const int hN_npelas_npoints = hN_npelas_points_per_file * hN_npelas_nfiles;

    int hN_npelas_energies[hN_npelas_nfiles] = {
        50, 100, 150, 200, 250, 300, 350, 400, 450, 500,
       550, 600, 650, 700, 750, 800, 850, 900, 950, 1000
    };

    TGraph2D * hN_npelas_graph = new TGraph2D(hN_npelas_npoints);
    TGraph2D * hN_npelas_int = new TGraph2D(hN_npelas_nfiles*numPoints);
    //TGraph2D * hN_npelas_fudge = new TGraph2D(hN_npelas_nfiles*numPoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_npelas_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_npelas_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/pn/pn" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_npelas_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 

      if (test_intbounce) 
	{
	  int subdiv = numPoints*ifile;
	  for(int j=0;j<numPoints; j++)
	    {
	      hN_npelas_int -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
					 inukd->XSec(kPdgProton,kPdgNeutron,kPdgProton,kIHNFtElas,double(ke),-1+j*cstep));
	    }
	}
    }//loop over files

    hN_npelas_graph -> Write("d_np_elas_data");
    if (test_intbounce)
      {
	hN_npelas_int -> Write("d_np_elas_int");
	//hN_npelas_fudge -> Write("d_np_elas_fudge");
      }
    delete hN_npelas_graph;
    delete hN_npelas_int;
    //delete hN_npelas_fudge;
  }

  // kIHNFtElas, pipN :
  {
    const int hN_pipNelas_nfiles = 60;
    const int hN_pipNelas_points_per_file = 21;
    const int hN_pipNelas_npoints = hN_pipNelas_points_per_file * hN_pipNelas_nfiles;

    int hN_pipNelas_energies[hN_pipNelas_nfiles] = {
      10,  20,  30,  40,  50,  60,  70,  80,  90,
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660,
     700, 740, 780, 820, 860, 900, 940, 980,
     1020, 1060, 1100, 1140, 1180, 1220, 1260,
     1300, 1340, 1380, 1420, 1460, 1500
    };

    TGraph2D * hN_pipNelas_graph = new TGraph2D(hN_pipNelas_npoints);
    TGraph2D * hN_pipNelas_int = new TGraph2D(hN_pipNelas_nfiles*numPoints);
    //TGraph2D * hN_pipNelas_fudge = new TGraph2D(hN_pipNelas_nfiles*numPoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_pipNelas_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_pipNelas_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/pip/pip" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_pipNelas_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 

      if (test_intbounce) 
	{
	  int subdiv = numPoints*ifile;
	  for(int j=0;j<numPoints; j++)
	    {
	      hN_pipNelas_int -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
					 inukd->XSec(kPdgPiP,kPdgProton,kPdgPiP,kIHNFtElas,double(ke),-1+j*cstep));
	      //hN_pipNelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
	      //inukd->XSec(kPdgPiP,kPdgProton,kPdgPiP,kIHNFtElas,double(ke)+.05,-1+j*cstep));
	    }
	}
    }//loop over files

    hN_pipNelas_graph -> Write("d_pipN_elas_data");
    if (test_intbounce)
      {
	hN_pipNelas_int -> Write("d_pipN_elas_int");
	//hN_pipNelas_fudge -> Write("d_pipN_elas_fudge");
      }
    delete hN_pipNelas_graph;
    delete hN_pipNelas_int;
    //delete hN_pipNelas_fudge;
  }

  // kIHNFtElas, pi0N :
  {
    const int hN_pi0Nelas_nfiles = 60;
    const int hN_pi0Nelas_points_per_file = 21;
    const int hN_pi0Nelas_npoints = hN_pi0Nelas_points_per_file * hN_pi0Nelas_nfiles;

    int hN_pi0Nelas_energies[hN_pi0Nelas_nfiles] = {
      10,  20,  30,  40,  50,  60,  70,  80,  90,
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660,
     700, 740, 780, 820, 860, 900, 940, 980,
     1020, 1060, 1100, 1140, 1180, 1220, 1260,
     1300, 1340, 1380, 1420, 1460, 1500
    };

    TGraph2D * hN_pi0Nelas_graph = new TGraph2D(hN_pi0Nelas_npoints);
    TGraph2D * hN_pi0Nelas_int = new TGraph2D(hN_pi0Nelas_nfiles*numPoints);
    //TGraph2D * hN_pi0Nelas_fudge = new TGraph2D(hN_pi0Nelas_nfiles*numPoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_pi0Nelas_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_pi0Nelas_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/pip/pip" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_pi0Nelas_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 

      if (test_intbounce) 
	{
	  int subdiv = numPoints*ifile;
	  for(int j=0;j<numPoints; j++)
	    {
	      hN_pi0Nelas_int -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
					 inukd->XSec(kPdgPi0,kPdgProton,kPdgPi0,kIHNFtElas,double(ke),-1+j*cstep));
	      //hN_pi0Nelas_fudge -> SetPoint(subdiv+j,double(ke),-1+j*cstep,
	      //   inukd->XSec(kPdgPi0,kPdgProton,kPdgPi0,kIHNFtElas,double(ke)+.05,-1+j*cstep));
	    }
	}
    }//loop over files

    hN_pi0Nelas_graph -> Write("d_pi0N_elas_data");
    if (test_intbounce)
      {
	hN_pi0Nelas_int -> Write("d_pi0N_elas_int");
	//hN_pi0Nelas_fudge -> Write("d_pi0N_elas_fudge");
      }
    delete hN_pi0Nelas_graph;
    delete hN_pi0Nelas_int;
    //delete hN_pi0Nelas_fudge;
  }

  // kIHNFtElas, pimN :
  {
    const int hN_pimNelas_nfiles = 60;
    const int hN_pimNelas_points_per_file = 21;
    const int hN_pimNelas_npoints = hN_pimNelas_points_per_file * hN_pimNelas_nfiles;

    int hN_pimNelas_energies[hN_pimNelas_nfiles] = {      
      10,  20,  30,  40,  50,  60,  70,  80,  90,
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660,
     700, 740, 780, 820, 860, 900, 940, 980,
     1020, 1060, 1100, 1140, 1180, 1220, 1260,
     1300, 1340, 1380, 1420, 1460, 1500
    };

    TGraph2D * hN_pimNelas_graph = new TGraph2D(hN_pimNelas_npoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_pimNelas_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_pimNelas_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/pim/pim" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_pimNelas_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 
    }//loop over files

    hN_pimNelas_graph -> Write("d_pimN_elas_data");
    delete hN_pimNelas_graph;
  }

  // kIHNFtCEx, (pi+, pi0, pi-) N
  {
    const int hN_piNcex_nfiles = 60;
    const int hN_piNcex_points_per_file = 21;
    const int hN_piNcex_npoints = hN_piNcex_points_per_file * hN_piNcex_nfiles;

    int hN_piNcex_energies[hN_piNcex_nfiles] = {
      10,  20,  30,  40,  50,  60,  70,  80,  90,
     100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
     200, 210, 220, 230, 240, 250, 260, 270, 280, 290,
     300, 340, 380, 420, 460, 500, 540, 580, 620, 660,
     700, 740, 780, 820, 860, 900, 940, 980,
     1020, 1060, 1100, 1140, 1180, 1220, 1260,
     1300, 1340, 1380, 1420, 1460, 1500
    };

    TGraph2D * hN_piNcex_graph = new TGraph2D(hN_piNcex_npoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_piNcex_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_piNcex_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/pie/pie" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_piNcex_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 
    }//loop over files

    hN_piNcex_graph -> Write("d_piN_cex_data");
    delete hN_piNcex_graph;
  }

  // kIHNFtAbs, (pi+, pi0, pi-) N
  {
    const int hN_piNabs_nfiles = 19;
    const int hN_piNabs_points_per_file = 21;
    const int hN_piNabs_npoints = hN_piNabs_points_per_file * hN_piNabs_nfiles;

    int hN_piNabs_energies[hN_piNabs_nfiles] = {
      50,  75, 100, 125, 150, 175, 200, 225, 250, 275,
     300, 325, 350, 375, 400, 425, 450, 475, 500
    };

    TGraph2D * hN_piNabs_graph = new TGraph2D(hN_piNabs_npoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_piNabs_nfiles; ifile++) {
      // build filename
      ostringstream hN_datafile;
      int ke = hN_piNabs_energies[ifile];
      hN_datafile << data_dir << "/diff_ang/pid2p/pid2p" << ke << ".txt";
      // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        hN_piNabs_graph -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff; 
    }//loop over files

    hN_piNabs_graph -> Write("d_piN_abs_data");
    delete hN_piNabs_graph;
  }

  // kIHNFtInelas, gamma p -> p pi0
  {
    const int hN_gampi0pInelas_nfiles = 29;
    const int hN_gampi0pInelas_points_per_file = 37;
    const int hN_gampi0pInelas_npoints = hN_gampi0pInelas_points_per_file * hN_gampi0pInelas_nfiles;
    
    int hN_gampi0pInelas_energies[hN_gampi0pInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    TGraph2D * fhN2dXSecGamPi0P_Inelas = new TGraph2D(hN_gampi0pInelas_npoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_gampi0pInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     int ke = hN_gampi0pInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi0p/" << ke << "-pi0p.txt";
     // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        fhN2dXSecGamPi0P_Inelas -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff;
    }//loop over files

    fhN2dXSecGamPi0P_Inelas -> Write("d_gam_pi0p_inel_data");
    delete fhN2dXSecGamPi0P_Inelas;
  }

  // kIHNFtInelas, gamma n -> n pi0
  {
    const int hN_gampi0nInelas_nfiles = 29;
    const int hN_gampi0nInelas_points_per_file = 37;
    const int hN_gampi0nInelas_npoints = hN_gampi0nInelas_points_per_file * hN_gampi0nInelas_nfiles;

    int hN_gampi0nInelas_energies[hN_gampi0nInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    TGraph2D * fhN2dXSecGamPi0N_Inelas = new TGraph2D(hN_gampi0nInelas_npoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_gampi0nInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     int ke = hN_gampi0nInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi0n/" << ke << "-pi0n.txt";
     // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        fhN2dXSecGamPi0N_Inelas -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff;
    }//loop over files

    fhN2dXSecGamPi0N_Inelas -> Write("d_gam_pi0n_inel_data");
    delete fhN2dXSecGamPi0N_Inelas;
  }

  // kIHNFtInelas, gamma p -> n pi+
  {
    const int hN_gampipnInelas_nfiles = 29;
    const int hN_gampipnInelas_points_per_file = 37;
    const int hN_gampipnInelas_npoints = hN_gampipnInelas_points_per_file * hN_gampipnInelas_nfiles;

    int hN_gampipnInelas_energies[hN_gampipnInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    TGraph2D * fhN2dXSecGamPipN_Inelas = new TGraph2D(hN_gampipnInelas_npoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_gampipnInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     int ke = hN_gampipnInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi+n/" << ke << "-pi+n.txt";
     // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        fhN2dXSecGamPipN_Inelas -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff;
    }//loop over files

    fhN2dXSecGamPipN_Inelas -> Write("d_gam_pipn_inel_data");
    delete fhN2dXSecGamPipN_Inelas;
  }

  // kIHNFtInelas, gamma n -> p pi-
  {
    const int hN_gampimpInelas_nfiles = 29;
    const int hN_gampimpInelas_points_per_file = 37;
    const int hN_gampimpInelas_npoints = hN_gampimpInelas_points_per_file * hN_gampimpInelas_nfiles;

    int hN_gampimpInelas_energies[hN_gampimpInelas_nfiles] = {  
      160,  180,  200,  220,  240,  260,  280,  300,  320,  340,   
      360,  380,  400,  450,  500,  550,  600,  650,  700,  750,
      800,  850,  900,  950,  1000, 1050, 1100, 1150, 1200
    };

    TGraph2D * fhN2dXSecGamPimP_Inelas = new TGraph2D(hN_gampimpInelas_npoints);

    int ipoint=0;
    double x, y;

    for(int ifile = 0; ifile < hN_gampimpInelas_nfiles; ifile++) {
     // build filename
     ostringstream hN_datafile;
     int ke = hN_gampimpInelas_energies[ifile];
     hN_datafile << data_dir << "/diff_ang/gampi-p/" << ke << "-pi-p.txt";
     // read data
      TGraph * buff = new TGraph(hN_datafile.str().c_str());
      for(int i=0; i<buff->GetN(); i++) {
        buff -> GetPoint(i,x,y);
        fhN2dXSecGamPimP_Inelas -> SetPoint(ipoint,(double)ke,TMath::Cos(x*kPi/180.),y);
        ipoint++;
      }
      delete buff;
    }//loop over files

    fhN2dXSecGamPimP_Inelas -> Write("d_gam_pimp_inel_data");
    delete fhN2dXSecGamPimP_Inelas;
  }
 
  f.Close();
}
//____________________________________________________________________________
void SaveSplinesToRootFile(void)
{
  string filename = "hadxs.root";

  LOG("testINukeHadroData", pNOTICE) 
        << "Saving INTRANUKE hadron x-section splines to file: " << filename;

  //-- load SAID hadron cross section data
  INukeHadroData * inukd = INukeHadroData::Instance();

  //-- save splines to a ROOT file
  //
  // Note: splines are beeing saved to the ROOT files as genie::Spline
  // objects. To read them back, one in a ROOT session load the GENIE
  // libraries first using the $GENIE/src/scripts/loadlibs.C script
  // eg:
  // root [0] .x $GENIE/src/scripts/loadlibs.C
  // root [1] TFile f("./hdxs.root");
  // root [2] pN_tot->Draw();
  // root [3] pi0N_tot->Draw("same");
  //
  
  inukd -> XSecPn_Tot       () -> SaveAsROOT(filename, "pn_tot",    true);
  inukd -> XSecPn_Elas      () -> SaveAsROOT(filename, "pn_elas",   false);
  inukd -> XSecPn_Reac      () -> SaveAsROOT(filename, "pn_reac",   false);
  inukd -> XSecNn_Tot       () -> SaveAsROOT(filename, "nn_tot",    false);
  inukd -> XSecNn_Elas      () -> SaveAsROOT(filename, "nn_elas",   false);
  inukd -> XSecNn_Reac      () -> SaveAsROOT(filename, "nn_reac",   false);
  inukd -> XSecPp_Tot       () -> SaveAsROOT(filename, "pp_tot",    true);
  inukd -> XSecPp_Elas      () -> SaveAsROOT(filename, "pp_elas",   false);
  inukd -> XSecPp_Reac      () -> SaveAsROOT(filename, "pp_reac",   false);

  inukd -> XSecPipn_Tot     () -> SaveAsROOT(filename, "pipn_tot",  false);
  inukd -> XSecPipn_CEx     () -> SaveAsROOT(filename, "pipn_cex",  false);
  inukd -> XSecPipn_Elas    () -> SaveAsROOT(filename, "pipn_elas", false);
  inukd -> XSecPipn_Reac    () -> SaveAsROOT(filename, "pipn_reac", false);
  inukd -> XSecPi0n_Tot     () -> SaveAsROOT(filename, "pi0n_tot",  false);
  inukd -> XSecPi0n_CEx     () -> SaveAsROOT(filename, "pi0n_cex",  false);
  inukd -> XSecPi0n_Elas    () -> SaveAsROOT(filename, "pi0n_elas", false);
  inukd -> XSecPi0n_Reac    () -> SaveAsROOT(filename, "pi0n_reac", false);
  inukd -> XSecPipp_Tot     () -> SaveAsROOT(filename, "pipp_tot",  false);
  inukd -> XSecPipp_CEx     () -> SaveAsROOT(filename, "pipp_cex",  false);
  inukd -> XSecPipp_Elas    () -> SaveAsROOT(filename, "pipp_elas", false);
  inukd -> XSecPipp_Reac    () -> SaveAsROOT(filename, "pipp_reac", false);
  inukd -> XSecPipp_Abs     () -> SaveAsROOT(filename, "pipp_abs",  false); 
  inukd -> XSecPi0p_Tot     () -> SaveAsROOT(filename, "pi0p_tot",  false);
  inukd -> XSecPi0p_CEx     () -> SaveAsROOT(filename, "pi0p_cex",  false);
  inukd -> XSecPi0p_Elas    () -> SaveAsROOT(filename, "pi0p_elas", false);
  inukd -> XSecPi0p_Reac    () -> SaveAsROOT(filename, "pi0p_reac", false);
  inukd -> XSecPi0p_Abs     () -> SaveAsROOT(filename, "pi0p_abs",  false);

  inukd -> XSecKpn_Elas     () -> SaveAsROOT(filename, "kpn_elas",  false);
  inukd -> XSecKpp_Elas     () -> SaveAsROOT(filename, "kpp_elas",  false);
  inukd -> XSecGamp_fs      () -> SaveAsROOT(filename, "gamp_fs",   false);
  inukd -> XSecGamn_fs      () -> SaveAsROOT(filename, "gamn_fs",   false);
  inukd -> XSecGamN_Tot     () -> SaveAsROOT(filename, "gamN_tot",  false);
  
  inukd -> FracPA_Tot       () -> SaveAsROOT(filename, "pA_tot",       false);
  inukd -> FracPA_Elas      () -> SaveAsROOT(filename, "pA_elas",      false);
  inukd -> FracPA_Inel      () -> SaveAsROOT(filename, "pA_inel",      false);
  inukd -> FracPA_CEx       () -> SaveAsROOT(filename, "pA_cex",       false);
  inukd -> FracPA_Abs       () -> SaveAsROOT(filename, "pA_abs",       false);
  inukd -> FracPA_Pipro     () -> SaveAsROOT(filename, "pA_pipro",     false);
  inukd -> FracNA_Tot       () -> SaveAsROOT(filename, "nA_tot",       false);
  inukd -> FracNA_Elas      () -> SaveAsROOT(filename, "nA_elas",      false);
  inukd -> FracNA_Inel      () -> SaveAsROOT(filename, "nA_inel",      false);
  inukd -> FracNA_CEx       () -> SaveAsROOT(filename, "nA_cex",       false);
  inukd -> FracNA_Abs       () -> SaveAsROOT(filename, "nA_abs",       false);
  inukd -> FracNA_Pipro     () -> SaveAsROOT(filename, "nA_pipro",     false);
  inukd -> FracPipA_Tot     () -> SaveAsROOT(filename, "pipA_tot",     false);
  inukd -> FracPipA_Elas    () -> SaveAsROOT(filename, "pipA_elas",    false);
  inukd -> FracPipA_Inel    () -> SaveAsROOT(filename, "pipA_inel",    false);
  inukd -> FracPipA_CEx     () -> SaveAsROOT(filename, "pipA_cex",     false);
  inukd -> FracPipA_Abs     () -> SaveAsROOT(filename, "pipA_abs",     false);  
  inukd -> FracPimA_Tot     () -> SaveAsROOT(filename, "pimA_tot",     false);
  inukd -> FracPimA_Elas    () -> SaveAsROOT(filename, "pimA_elas",    false);
  inukd -> FracPimA_Inel    () -> SaveAsROOT(filename, "pimA_inel",    false);
  inukd -> FracPimA_CEx     () -> SaveAsROOT(filename, "pimA_cex",     false);
  inukd -> FracPimA_Abs     () -> SaveAsROOT(filename, "pimA_abs",     false);
  inukd -> FracPi0A_Tot     () -> SaveAsROOT(filename, "pi0A_tot",     false);
  inukd -> FracPi0A_Elas    () -> SaveAsROOT(filename, "pi0A_elas",    false);
  inukd -> FracPi0A_Inel    () -> SaveAsROOT(filename, "pi0A_inel",    false);
  inukd -> FracPi0A_CEx     () -> SaveAsROOT(filename, "pi0A_cex",     false);
  inukd -> FracPi0A_Abs     () -> SaveAsROOT(filename, "pi0A_abs",     false);
}
//____________________________________________________________________________
void PrintOutForInputKE(double ke)
{
  LOG("testINukeHadroData", pNOTICE) 
          << "Printing out INTRANUKE hadron x-sections";
 
  //-- load SAID hadron cross section data
  INukeHadroData * inukd = INukeHadroData::Instance();

  //-- Print out the hN mode hadron x-section 
  LOG("testINukeHadroData", pNOTICE)  
     << "\n hN mode x-sections:"
     << "\n XSec[pn/tot]      (K=" << ke << " MeV) = " << inukd -> XSecPn_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pn/elas]     (K=" << ke << " MeV) = " << inukd -> XSecPn_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pn/reac]     (K=" << ke << " MeV) = " << inukd -> XSecPn_Reac      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nn/tot]      (K=" << ke << " MeV) = " << inukd -> XSecNn_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nn/elas]     (K=" << ke << " MeV) = " << inukd -> XSecNn_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nn/reac]     (K=" << ke << " MeV) = " << inukd -> XSecNn_Reac      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pp/tot]      (K=" << ke << " MeV) = " << inukd -> XSecPp_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pp/elas]     (K=" << ke << " MeV) = " << inukd -> XSecPp_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pp/reac]     (K=" << ke << " MeV) = " << inukd -> XSecPp_Reac      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipn/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPipn_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipn/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPipn_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipn/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPipn_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipn/reac]   (K=" << ke << " MeV) = " << inukd -> XSecPipn_Reac    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0n/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPi0n_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0n/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPi0n_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0n/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPi0n_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0n/reac]   (K=" << ke << " MeV) = " << inukd -> XSecPi0n_Reac    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipp/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPipp_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipp/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPipp_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipp/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPipp_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipp/reac]   (K=" << ke << " MeV) = " << inukd -> XSecPipp_Reac    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipp/abs]    (K=" << ke << " MeV) = " << inukd -> XSecPipp_Abs     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0p/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPi0p_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0p/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPi0p_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0p/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPi0p_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0p/reac]   (K=" << ke << " MeV) = " << inukd -> XSecPi0p_Reac    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0p/abs]    (K=" << ke << " MeV) = " << inukd -> XSecPi0p_Abs     () -> Evaluate(ke) << " mbarn"
     << "\n"
     << "\n hA mode x-sections:"
     << "\n Frac[pA/tot]      (K=" << ke << " MeV) = " << inukd -> FracPA_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pA/elas]     (K=" << ke << " MeV) = " << inukd -> FracPA_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pA/inel]     (K=" << ke << " MeV) = " << inukd -> FracPA_Inel      () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pA/cex]      (K=" << ke << " MeV) = " << inukd -> FracPA_CEx       () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pA/abs]      (K=" << ke << " MeV) = " << inukd -> FracPA_Abs       () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pA/pipro]    (K=" << ke << " MeV) = " << inukd -> FracPA_Pipro     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[nA/tot]      (K=" << ke << " MeV) = " << inukd -> FracNA_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n Frac[nA/elas]     (K=" << ke << " MeV) = " << inukd -> FracNA_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n Frac[nA/inel]     (K=" << ke << " MeV) = " << inukd -> FracNA_Inel      () -> Evaluate(ke) << " mbarn"
     << "\n Frac[nA/reac]     (K=" << ke << " MeV) = " << inukd -> FracNA_CEx       () -> Evaluate(ke) << " mbarn"
     << "\n Frac[nA/abs]      (K=" << ke << " MeV) = " << inukd -> FracNA_Abs       () -> Evaluate(ke) << " mbarn"
     << "\n Frac[nA/pipro]    (K=" << ke << " MeV) = " << inukd -> FracNA_Pipro     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pipA/tot]    (K=" << ke << " MeV) = " << inukd -> FracPipA_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pipA/elas]   (K=" << ke << " MeV) = " << inukd -> FracPipA_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pipA/inel]   (K=" << ke << " MeV) = " << inukd -> FracPipA_Inel    () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pipA/cex]    (K=" << ke << " MeV) = " << inukd -> FracPipA_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pipA/abs]    (K=" << ke << " MeV) = " << inukd -> FracPipA_Abs     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pimA/tot]    (K=" << ke << " MeV) = " << inukd -> FracPimA_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pimA/elas]   (K=" << ke << " MeV) = " << inukd -> FracPimA_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pimA/inel]   (K=" << ke << " MeV) = " << inukd -> FracPimA_Inel    () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pimA/cex]    (K=" << ke << " MeV) = " << inukd -> FracPimA_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pimA/abs]    (K=" << ke << " MeV) = " << inukd -> FracPimA_Abs     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pi0A/tot]    (K=" << ke << " MeV) = " << inukd -> FracPi0A_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pi0A/elas]   (K=" << ke << " MeV) = " << inukd -> FracPi0A_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pi0A/inel]   (K=" << ke << " MeV) = " << inukd -> FracPi0A_Inel    () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pi0A/cex]    (K=" << ke << " MeV) = " << inukd -> FracPi0A_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n Frac[pi0A/abs]    (K=" << ke << " MeV) = " << inukd -> FracPi0A_Abs     () -> Evaluate(ke) << " mbarn";
}
//____________________________________________________________________________

