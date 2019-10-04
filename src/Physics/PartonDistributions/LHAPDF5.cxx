//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TSystem.h>
#include <TMath.h>

#include "Framework/Conventions/GBuild.h"
#include "Physics/PartonDistributions/LHAPDF5.h"
#include "Framework/Messenger/Messenger.h"

#ifdef __GENIE_LHAPDF5_ENABLED__
#include "LHAPDF/LHAPDF.h"
#endif

using namespace genie;

//____________________________________________________________________________
LHAPDF5::LHAPDF5() :
PDFModelI("genie::LHAPDF5")
{
  this->Initialize();
}
//____________________________________________________________________________
LHAPDF5::LHAPDF5(string config) :
PDFModelI("genie::LHAPDF5", config)
{
  LOG("LHAPDF5", pDEBUG) << "LHAPDF5 configuration:\n " << this->GetConfig();  

  this->Initialize();
}
//____________________________________________________________________________
LHAPDF5::~LHAPDF5() 
{ 

}
//____________________________________________________________________________
void LHAPDF5::Initialize(void) const
{
#ifdef __GENIE_LHAPDF5_ENABLED__
  bool lhapath_ok = true;
  const char * lhapath = gSystem->Getenv("LHAPATH");
  if(!lhapath) lhapath_ok = false;
  else {
      void *dirp = gSystem->OpenDirectory(lhapath);
      if (dirp) gSystem->FreeDirectory(dirp);
      else lhapath_ok = false;
  }
  if(!lhapath_ok) {
   LOG("LHAPDF5", pFATAL) 
     << "\n"
     << "** LHAPDF won't be able to read-in the PDF data. \n"
     << "** The LHAPATH env. variable is not properly (or at all) defined. \n"
     << "** Please, set LHAPATH to <lhapdf_top_dir>/PDFsets/ \n"
     << "** See http://projects.hepforge.org/lhapdf/ for more details. \n\n";
   gAbortingInErr = true;
   exit(1);
  }
#endif
}
//____________________________________________________________________________
void LHAPDF5::SetPDFSetFromConfig(void) const
{
// Get PDF spec (particle type, pdf group/set) from configuration registry.
// For definitions, have a look at PDFLIB and LHAPDF manuals

#ifdef __GENIE_LHAPDF5_ENABLED__
  string name   = "";
  int    type   = 0;
  int    memset = 0;

  this->GetParam("name_lhapdf",   name);
  this->GetParam("type_lhapdf",   type);
  this->GetParam("memset_lhapdf", memset);

  LHAPDF::SetType stype = (type==0) ? LHAPDF::LHPDF :  LHAPDF::LHGRID;

  LHAPDF::initPDFByName(name, stype, memset);
  LHAPDF::extrapolate(false);
#endif
}
//____________________________________________________________________________
double LHAPDF5::UpValence(double x, double Q2) const
{
  return AllPDFs(x,Q2).uval;
}
//____________________________________________________________________________
double LHAPDF5::DownValence(double x, double Q2) const
{
  return AllPDFs(x,Q2).dval;
}
//____________________________________________________________________________
double LHAPDF5::UpSea(double x, double Q2) const
{
  return AllPDFs(x,Q2).usea;
}
//____________________________________________________________________________
double LHAPDF5::DownSea(double x, double Q2) const
{
  return AllPDFs(x,Q2).dsea;
}
//____________________________________________________________________________
double LHAPDF5::Strange(double x, double Q2) const
{
  return AllPDFs(x,Q2).str;
}
//____________________________________________________________________________
double LHAPDF5::Charm(double x, double Q2) const
{
  return AllPDFs(x,Q2).chm;
}
//____________________________________________________________________________
double LHAPDF5::Bottom(double x, double Q2) const
{
  return AllPDFs(x,Q2).bot;
}
//____________________________________________________________________________
double LHAPDF5::Top(double x, double Q2) const
{
  return AllPDFs(x,Q2).top;
}
//____________________________________________________________________________
double LHAPDF5::Gluon(double x, double Q2) const
{
  return AllPDFs(x,Q2).gl;
}
//____________________________________________________________________________
PDF_t LHAPDF5::AllPDFs(
#ifdef __GENIE_LHAPDF5_ENABLED__
  double x, double Q2
#else
  double, double
#endif
) const
{
  PDF_t pdf;

#ifdef __GENIE_LHAPDF5_ENABLED__
  // QCD scale
  double Q = TMath::Sqrt( TMath::Abs(Q2) ); 

  vector<double> pdfs = LHAPDF::xfx(x, Q);
  pdf.uval = pdfs[8] - pdfs[4];
  pdf.dval = pdfs[7] - pdfs[5];
  pdf.usea = pdfs[4];
  pdf.dsea = pdfs[5];
  pdf.str  = pdfs[9];
  pdf.chm  = pdfs[10];
  pdf.bot  = pdfs[11];
  pdf.top  = pdfs[12];
  pdf.gl   = pdfs[6];;
#endif

  return pdf;                                               
}
//____________________________________________________________________________
void LHAPDF5::Configure(const Registry & config)
{
  Algorithm::Configure(config);

  this->Initialize();         
  this->SetPDFSetFromConfig();

  fAllowReconfig=false;
}
//____________________________________________________________________________
void LHAPDF5::Configure(string config)
{
  Algorithm::Configure(config);

  this->Initialize();          
  this->SetPDFSetFromConfig(); 

  fAllowReconfig=false;
}
//____________________________________________________________________________

