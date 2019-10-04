//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 
*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>

#include <TSystem.h>
#include <TMath.h>

#include "Physics/PartonDistributions/LHAPDF6.h"
#include "Framework/Messenger/Messenger.h"

#ifdef __GENIE_LHAPDF6_ENABLED__
#include "LHAPDF/LHAPDF.h"
#endif

using namespace genie;

//____________________________________________________________________________
LHAPDF6::LHAPDF6() :
PDFModelI("genie::LHAPDF6")
{
  fSetName = "";
  fMemberID = 0;

#ifdef __GENIE_LHAPDF6_ENABLED__
  fLHAPDF = 0;
#endif
}
//____________________________________________________________________________
LHAPDF6::LHAPDF6(string config) :
PDFModelI("genie::LHAPDF6", config)
{
  fSetName = "";
  fMemberID = 0;

#ifdef __GENIE_LHAPDF6_ENABLED__
  fLHAPDF = 0;
#endif
}
//____________________________________________________________________________
LHAPDF6::~LHAPDF6() 
{ 
//#ifdef __GENIE_LHAPDF6_ENABLED__
//  if(fLHAPDF) delete fLHAPDF;
//#endif
}
//____________________________________________________________________________
double LHAPDF6::UpValence(double x, double Q2) const
{
  return AllPDFs(x,Q2).uval;
}
//____________________________________________________________________________
double LHAPDF6::DownValence(double x, double Q2) const
{
  return AllPDFs(x,Q2).dval;
}
//____________________________________________________________________________
double LHAPDF6::UpSea(double x, double Q2) const
{
  return AllPDFs(x,Q2).usea;
}
//____________________________________________________________________________
double LHAPDF6::DownSea(double x, double Q2) const
{
  return AllPDFs(x,Q2).dsea;
}
//____________________________________________________________________________
double LHAPDF6::Strange(double x, double Q2) const
{
  return AllPDFs(x,Q2).str;
}
//____________________________________________________________________________
double LHAPDF6::Charm(double x, double Q2) const
{
  return AllPDFs(x,Q2).chm;
}
//____________________________________________________________________________
double LHAPDF6::Bottom(double x, double Q2) const
{
  return AllPDFs(x,Q2).bot;
}
//____________________________________________________________________________
double LHAPDF6::Top(double x, double Q2) const
{
  return AllPDFs(x,Q2).top;
}
//____________________________________________________________________________
double LHAPDF6::Gluon(double x, double Q2) const
{
  return AllPDFs(x,Q2).gl;
}
//____________________________________________________________________________
#ifdef __GENIE_LHAPDF6_ENABLED__
PDF_t LHAPDF6::AllPDFs(double x, double Q2) const
{
  PDF_t pdf;
  vector<double> pdfvec;
  fLHAPDF->xfxQ2(x,Q2,pdfvec);
  pdf.uval = pdfvec[8] - pdfvec[4];
  pdf.dval = pdfvec[7] - pdfvec[5];
  pdf.usea = pdfvec[4];
  pdf.dsea = pdfvec[5];
  pdf.str  = pdfvec[9];
  pdf.chm  = pdfvec[10];
  pdf.bot  = pdfvec[11];
  pdf.top  = pdfvec[12];
  pdf.gl   = pdfvec[6];
  return pdf;                                               
}
#else
PDF_t LHAPDF6::AllPDFs(double, double) const
{
  LOG("LHAPDF6",pFATAL) << "LHAPDF6 not enabled.";
  exit(-1);
}
#endif
//____________________________________________________________________________
void LHAPDF6::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
  fAllowReconfig=false;
}
//____________________________________________________________________________
void LHAPDF6::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();          
  fAllowReconfig=false;
}
//____________________________________________________________________________
void LHAPDF6::LoadConfig(void) 
{
  this->GetParam("SetName",  fSetName );
  this->GetParam("MemberID", fMemberID);

#ifdef __GENIE_LHAPDF6_ENABLED__
  fLHAPDF = LHAPDF::mkPDF(fSetName, fMemberID);
  if(!fLHAPDF) {
     LOG("LHAPDF6",pFATAL)
       << "Couldn't retrieve LHADPF6 pdf set: " 
       << fSetName << ", member id: " << fMemberID;
     exit(1);
  }
#endif
}
//____________________________________________________________________________

