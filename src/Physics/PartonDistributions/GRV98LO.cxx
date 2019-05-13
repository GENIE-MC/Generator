//____________________________________________________________________________
/*

 Copyright (c) 2003-2019, The GENIE Collaboration
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

#include "Physics/PartonDistributions/GRV98LO.h"
#include "Framework/Messenger/Messenger.h"

using namespace std;
using namespace genie;

//____________________________________________________________________________
GRV98LO::GRV98LO() :
PDFModelI("genie::GRV98LO"),
 fXUVF(NULL),
 fXDVF(NULL),
 fXDEF(NULL),
 fXUDF(NULL),
 fXSF (NULL),
 fXGF (NULL)
{
  this->Initialize();
}
//____________________________________________________________________________
GRV98LO::GRV98LO(string config) :
PDFModelI("genie::GRV98LO", config),
 fXUVF(NULL),
 fXDVF(NULL),
 fXDEF(NULL),
 fXUDF(NULL),
 fXSF (NULL),
 fXGF (NULL)
{
  LOG("GRV98LO", pDEBUG) << "GRV98LO configuration:\n " << GetConfig() ;

  this->Initialize();
}
//____________________________________________________________________________
GRV98LO::~GRV98LO() 
{ 
  if (fXUVF) {delete fXUVF; fXUVF = NULL;}
  if (fXDVF) {delete fXDVF; fXDVF = NULL;}
  if (fXDEF) {delete fXDEF; fXDEF = NULL;}
  if (fXUDF) {delete fXUDF; fXUDF = NULL;}
  if (fXSF ) {delete fXSF ; fXSF  = NULL;}
  if (fXGF ) {delete fXGF ; fXGF  = NULL;}
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
  Q2 = std::min(Q2, fGridQ2[kNQ2-1]);
  x  = std::max(x,  fGridXbj[0]);
  x  = std::min(x,  fGridXbj[kNXbj-1]);

  double logx  = std::log(x);
  double logQ2 = std::log(Q2);
  double x1    = 1-x;
  double xv    = std::sqrt(x);
  double xs    = std::pow(x, -0.2);
  double x1p3  = x1*x1*x1;
  double x1p4  = x1*x1p3;
  double x1p5  = x1*x1p4;
  double x1p7  = x1p3*x1p4;

  double uv = fXUVF->Eval(logx,logQ2) * x1p3 * xv;
  double dv = fXDVF->Eval(logx,logQ2) * x1p4 * xv;
  double de = fXDEF->Eval(logx,logQ2) * x1p7 * xv;
  double ud = fXUDF->Eval(logx,logQ2) * x1p7 * xs;
  double us = 0.5 * (ud - de);
  double ds = 0.5 * (ud + de);
  double ss = fXSF->Eval(logx,logQ2)  * x1p7 * xs;
  double gl = fXGF->Eval(logx,logQ2)  * x1p5 * xs;
  
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
  for(int j=0; j < kNXbj; j++) {
    double xbj = -1;
    grid_file >> xbj;
    // check against known limits
    // ...    
    fGridXbj[j] = xbj;
    fGridLogXbj[j] = TMath::Log(xbj);
  }
  ostringstream grid_values;
  grid_values << "(";
  for(int j=0; j < kNXbj; j++) {
    grid_values << fGridXbj[j];
    if(j == kNXbj - 1) { grid_values << ")";  }
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
  for(int j=0; j < kNXbj-1; j++) {
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
  
  vector<double> gridLogQ2 (kNQ2);
  vector<double> gridLogXbj(kNXbj);
  vector<double> knotsXUVF(kNQ2*kNXbj);
  vector<double> knotsXDVF(kNQ2*kNXbj);
  vector<double> knotsXDEF(kNQ2*kNXbj);
  vector<double> knotsXUDF(kNQ2*kNXbj);
  vector<double> knotsXSF (kNQ2*kNXbj);
  vector<double> knotsXGF (kNQ2*kNXbj);

  k=0;
  for(int i=0; i < kNQ2; i++) {
    double logQ2 = std::log(fGridQ2[i]);
    gridLogQ2[i] = logQ2;
    for(int j=0; j < kNXbj - 1; j++) {
       double logx  = std::log(fGridXbj[j]);
       gridLogXbj[j] = logx;
       double xb0v  = std::sqrt(fGridXbj[j]);
       double xb0s  = std::pow(fGridXbj[j], -0.2);
       double xb1   = 1 - fGridXbj[j];
       double xb1p3 = std::pow(xb1, 3.);
       double xb1p4 = std::pow(xb1, 4.);
       double xb1p5 = std::pow(xb1, 5.);
       double xb1p7 = std::pow(xb1, 7.);
       knotsXUVF[k] = fParton[0][i][j] / (xb1p3 * xb0v);
       knotsXDVF[k] = fParton[1][i][j] / (xb1p4 * xb0v);
       knotsXDEF[k] = fParton[2][i][j] / (xb1p7 * xb0v);
       knotsXUDF[k] = fParton[3][i][j] / (xb1p7 * xb0s);
       knotsXSF [k] = fParton[4][i][j] / (xb1p7 * xb0s);
       knotsXGF [k] = fParton[5][i][j] / (xb1p5 * xb0s);
       k++;
    }
    double logxmax = TMath::Log(fGridXbj[kNXbj-1]);
    gridLogXbj[kNXbj-1] = logxmax;
    knotsXUVF[k] = 0;
    knotsXDVF[k] = 0;
    knotsXDEF[k] = 0;
    knotsXUDF[k] = 0;
    knotsXSF [k] = 0;
    knotsXGF [k] = 0;
    k++;
  }
  
  fXUVF = new Interpolator2D(gridLogXbj.size(),&gridLogXbj[0],gridLogQ2.size(),&gridLogQ2[0],&knotsXUVF[0]);
  fXDVF = new Interpolator2D(gridLogXbj.size(),&gridLogXbj[0],gridLogQ2.size(),&gridLogQ2[0],&knotsXDVF[0]);
  fXDEF = new Interpolator2D(gridLogXbj.size(),&gridLogXbj[0],gridLogQ2.size(),&gridLogQ2[0],&knotsXDEF[0]);
  fXUDF = new Interpolator2D(gridLogXbj.size(),&gridLogXbj[0],gridLogQ2.size(),&gridLogQ2[0],&knotsXUDF[0]);
  fXSF  = new Interpolator2D(gridLogXbj.size(),&gridLogXbj[0],gridLogQ2.size(),&gridLogQ2[0],&knotsXSF [0]);
  fXGF  = new Interpolator2D(gridLogXbj.size(),&gridLogXbj[0],gridLogQ2.size(),&gridLogQ2[0],&knotsXGF [0]);
  
  fInitialized = true;
}
//____________________________________________________________________________
