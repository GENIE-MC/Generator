//____________________________________________________________________________
/*

 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions:

*/
//____________________________________________________________________________

#include <sstream>
#include <fstream>
#include <cstdlib>

#include <TSystem.h>
#include <TMath.h>

#include "PDF/GRV98LO.h"
#include "Messenger/Messenger.h"

using namespace std;
using namespace genie;

//____________________________________________________________________________
GRV98LO::GRV98LO() :
PDFModelI("genie::GRV98LO")
{
  this->Initialize();
}
//____________________________________________________________________________
GRV98LO::GRV98LO(string config) :
PDFModelI("genie::GRV98LO", config)
{
  LOG("GRV98LO", pDEBUG) << "GRV98LO configuration:\n " << *fConfig;  

  this->Initialize();
}
//____________________________________________________________________________
GRV98LO::~GRV98LO() 
{ 

}
//____________________________________________________________________________
double GRV98LO::UpValence(double x, double Q2) const
{
  return AllPDFs(x,Q2).uval;
}
//____________________________________________________________________________
double GRV98LO::DownValence(double x, double Q2) const
{
  return AllPDFs(x,Q2).dval;
}
//____________________________________________________________________________
double GRV98LO::UpSea(double x, double Q2) const
{
  return AllPDFs(x,Q2).usea;
}
//____________________________________________________________________________
double GRV98LO::DownSea(double x, double Q2) const
{
  return AllPDFs(x,Q2).dsea;
}
//____________________________________________________________________________
double GRV98LO::Strange(double x, double Q2) const
{
  return AllPDFs(x,Q2).str;
}
//____________________________________________________________________________
double GRV98LO::Charm(double x, double Q2) const
{
  return AllPDFs(x,Q2).chm;
}
//____________________________________________________________________________
double GRV98LO::Bottom(double x, double Q2) const
{
  return AllPDFs(x,Q2).bot;
}
//____________________________________________________________________________
double GRV98LO::Top(double x, double Q2) const
{
  return AllPDFs(x,Q2).top;
}
//____________________________________________________________________________
double GRV98LO::Gluon(double x, double Q2) const
{
  return AllPDFs(x,Q2).gl;
}
//____________________________________________________________________________
PDF_t GRV98LO::AllPDFs(double x, double Q2) const
{
  PDF_t pdf;

  if(!fInitialized) {
    LOG("GRV98LO", pWARN) 
      << "GRV98LO algorithm was not initialized succesfully";
    pdf.uval = 0.;
    pdf.dval = 0.; 
    pdf.usea = 0.;
    pdf.dsea = 0.;
    pdf.str  = 0.;
    pdf.chm  = 0.;
    pdf.bot  = 0.;
    pdf.top  = 0.;
    pdf.gl   = 0.;
    return pdf;
  }

  LOG("GRV98LO", pDEBUG) 
    << "Inputs x = " << x << ", Q2 = " << Q2;

  // apply kinematical limits 
//Q2 = TMath::Max(Q2, fGridQ2[0]);
  if(Q2 <= 0.8) Q2 = 0.80001;
  Q2 = TMath::Min(Q2, fGridQ2[kNQ2-1]);
  x  = TMath::Max(x,  fGridXbj[0]);
  x  = TMath::Min(x,  fGridXbj[kNXBj-1]);

  double logx  = TMath::Log(x);
  double logQ2 = TMath::Log(Q2);
  double x1    = 1-x;
  double xv    = TMath::Power(x,  0.5);
  double xs    = TMath::Power(x, -0.2);
  double x1p3  = TMath::Power(x1, 3.);
  double x1p4  = TMath::Power(x1, 4.);
  double x1p5  = TMath::Power(x1, 5.);
  double x1p7  = TMath::Power(x1, 7.);

  double uv = fXUVF.Interpolate(logx,logQ2) * x1p3 * xv;
  double dv = fXDVF.Interpolate(logx,logQ2) * x1p4 * xv;
  double de = fXDEF.Interpolate(logx,logQ2) * x1p7 * xv;
  double ud = fXUDF.Interpolate(logx,logQ2) * x1p7 * xs;
  double us = 0.5 * (ud - de);
  double ds = 0.5 * (ud + de);
  double ss = fXSF.Interpolate(logx,logQ2)  * x1p7 * xs;
  double gl = fXGF.Interpolate(logx,logQ2)  * x1p5 * xs;
	
  pdf.uval = uv;
  pdf.dval = dv; 
  pdf.usea = us;
  pdf.dsea = ds;
  pdf.str  = ss;
  pdf.chm  = 0.;
  pdf.bot  = 0.;
  pdf.top  = 0.;
  pdf.gl   = gl;

  return pdf;                                               
}
//____________________________________________________________________________
void GRV98LO::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->Initialize();         
}
//____________________________________________________________________________
void GRV98LO::Configure(string config)
{
  Algorithm::Configure(config);
  this->Initialize();          
}
//____________________________________________________________________________
void GRV98LO::Initialize(void)
{
  fInitialized = false;

  const char * genie_dir = gSystem->Getenv("GENIE"); 
  if(!genie_dir) return;

  string grid_file_name = 
     string(gSystem->Getenv("GENIE")) + string("/data/evgen/pdfs/GRV98lo_patched.LHgrid");

  LOG("GRV98LO", pNOTICE) 
    << "Reading grid file from:\n " << grid_file_name; 

  ifstream grid_file;
  grid_file.open (grid_file_name.c_str());

  char rubbish[1000];

  const int nskip = 21;
  for(int i=0; i<nskip; i++) {
   grid_file.getline(rubbish,1000);
   LOG("GRV98LO", pDEBUG) << "Skipping: " << rubbish; 
  }

  // x's
  //

  LOG("GRV98LO", pDEBUG) << "Reading x_bj grid values";
  for(int j=0; j < kNXBj; j++) {
    double xbj = -1;
    grid_file >> xbj;
    // check against known limits
    // ...    
    fGridXbj[j] = xbj;
    fGridLogXbj[j] = TMath::Log(xbj);
  }
  ostringstream grid_values;
  grid_values << "(";
  for(int j=0; j < kNXBj; j++) {
    grid_values << fGridXbj[j];
    if(j == kNXBj - 1) { grid_values << ")";  }
    else               { grid_values << ", "; }
  }
  LOG("GRV98LO", pDEBUG) 
    << "x_bj grid values: " << grid_values.str();

  // Q^2s
  //

  LOG("GRV98LO", pDEBUG) << "Reading Q^2 grid values.";
  for(int i=0; i < kNQ2; i++) {
    double Q2 = -1;
    grid_file >> Q2;
    // check against known limits
    // ...    
    fGridQ2[i] = Q2;
    fGridLogQ2[i] = TMath::Log(Q2);
  }
  grid_values.str("");
  grid_values << "(";
  for(int i=0; i < kNQ2; i++) {
    grid_values << fGridQ2[i];
    if(i == kNQ2 - 1) { grid_values << ")";  }
    else              { grid_values << ", "; }
  }
  LOG("GRV98LO", pDEBUG) 
    << "Q^2 grid values: " << grid_values.str() << "GeV^2";

  // skip again
  grid_file.getline(rubbish,1000);
  LOG("GRV98LO", pDEBUG) << "Skipping: " << rubbish; 
  grid_file.getline(rubbish,1000);
  LOG("GRV98LO", pDEBUG) << "Skipping: " << rubbish; 

  // pdf values on grid points
  // 

  LOG("GRV98LO", pDEBUG) << "Reading PDF values on grid points";

  int k=0;
  for(int j=0; j < kNXBj-1; j++) {
    for(int i=0; i < kNQ2; i++) {
      double p0 = 0;
      double p1 = 0;
      double p2 = 0;
      double p3 = 0;
      double p4 = 0;
      double p5 = 0;
      grid_file >> p0 >> p1 >> p2 >> p3 >> p4 >> p5;
      LOG("GRV98LO", pDEBUG) 
         << "Row: " << k << ", grid point: (" << i << ", " << j << ") ->" 
         << "  p0 = " << p0
         << ", p1 = " << p1
         << ", p2 = " << p2
         << ", p3 = " << p3
         << ", p4 = " << p4
         << ", p5 = " << p5;
      fParton [0][i][j]  = p0;
      fParton [1][i][j]  = p1;
      fParton [2][i][j]  = p2;
      fParton [3][i][j]  = p3;
      fParton [4][i][j]  = p4;
      fParton [5][i][j]  = p5;
      k++;
    }
  }

  grid_file.close();

  // arrays for interpolation routines
  // 

  k=0;
  for(int i=0; i < kNQ2; i++) {
    double logQ2 = TMath::Log(fGridQ2[i]);
    for(int j=0; j < kNXBj - 1; j++) {
       double logx  = TMath::Log  (fGridXbj[j]); 
       double xb0v  = TMath::Power(fGridXbj[j],  0.5);
       double xb0s  = TMath::Power(fGridXbj[j], -0.2);
       double xb1   = 1 - fGridXbj[j];
       double xb1p3 = TMath::Power(xb1, 3.);
       double xb1p4 = TMath::Power(xb1, 4.);
       double xb1p5 = TMath::Power(xb1, 5.);
       double xb1p7 = TMath::Power(xb1, 7.);
       fXUVF.SetPoint(k, logx, logQ2, fParton[0][i][j] / (xb1p3 * xb0v) );
       fXDVF.SetPoint(k, logx, logQ2, fParton[1][i][j] / (xb1p4 * xb0v) );
       fXDEF.SetPoint(k, logx, logQ2, fParton[2][i][j] / (xb1p7 * xb0v) );
       fXUDF.SetPoint(k, logx, logQ2, fParton[3][i][j] / (xb1p7 * xb0s) );
       fXSF .SetPoint(k, logx, logQ2, fParton[4][i][j] / (xb1p7 * xb0s) );
       fXGF .SetPoint(k, logx, logQ2, fParton[5][i][j] / (xb1p5 * xb0s) );
       k++;
    }
    double logxmax = TMath::Log(fGridXbj[kNXBj-1]);
    fXUVF.SetPoint(k, logxmax, logQ2, 0.);
    fXDVF.SetPoint(k, logxmax, logQ2, 0.);
    fXDEF.SetPoint(k, logxmax, logQ2, 0.);
    fXUDF.SetPoint(k, logxmax, logQ2, 0.);
    fXSF .SetPoint(k, logxmax, logQ2, 0.);
    fXGF .SetPoint(k, logxmax, logQ2, 0.);
    k++;
  }

  fInitialized = true;
}
//____________________________________________________________________________
