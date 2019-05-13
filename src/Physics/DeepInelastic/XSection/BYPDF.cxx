//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 09, 2009 - CA
   Renamed to BYPDF from BYPDFModel
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Physics/DeepInelastic/XSection/BYPDF.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
BYPDF::BYPDF() :
PDFModelI("genie::BYPDF")
{

}
//____________________________________________________________________________
BYPDF::BYPDF(string config) :
PDFModelI("genie::BYPDF", config)
{

}
//____________________________________________________________________________
BYPDF::~BYPDF()
{

}
//____________________________________________________________________________
double BYPDF::UpValence(double x, double q2) const
{
  return AllPDFs(x,q2).uval;
}
//____________________________________________________________________________
double BYPDF::DownValence(double x, double q2) const
{
  return AllPDFs(x,q2).dval;
}
//____________________________________________________________________________
double BYPDF::UpSea(double x, double q2) const
{
  return AllPDFs(x,q2).usea;
}
//____________________________________________________________________________
double BYPDF::DownSea(double x, double q2) const
{
  return AllPDFs(x,q2).dsea;
}
//____________________________________________________________________________
double BYPDF::Strange(double x, double q2) const
{
  return AllPDFs(x,q2).str;
}
//____________________________________________________________________________
double BYPDF::Charm(double x, double q2) const
{
  return AllPDFs(x,q2).chm;
}
//____________________________________________________________________________
double BYPDF::Bottom(double x, double q2) const
{
  return AllPDFs(x,q2).bot;
}
//____________________________________________________________________________
double BYPDF::Top(double x, double q2) const
{
  return AllPDFs(x,q2).top;
}
//____________________________________________________________________________
double BYPDF::Gluon(double x, double q2) const
{
  return AllPDFs(x,q2).gl;
}
//____________________________________________________________________________
PDF_t BYPDF::AllPDFs(double x, double q2) const
{
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) 
       << "Inputs: x = " << x << ", |q2| = " << TMath::Abs(q2);
#endif
  if(TMath::Abs(q2) < fQ2min) q2=fQ2min;

  // get the uncorrected PDFs
  PDF_t uncorrected_pdfs = fBasePDFModel->AllPDFs(x, q2);
  double uv = uncorrected_pdfs.uval;
  double us = uncorrected_pdfs.usea;
  double dv = uncorrected_pdfs.dval;
  double ds = uncorrected_pdfs.dsea;

  // compute correction factor delta(d/u)
  double delta = this->DeltaDU(x);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) << "delta(d/u) = " << delta;
#endif

  // compute u/(u+d) ratios for both valence & sea quarks
  double val = uv+dv;
  double sea = us+ds;
  double rv  = (val==0) ? 0. : uv/val;
  double rs  = (sea==0) ? 0. : us/sea;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG)
   << "valence[u/(u+d)] = " << rv << ", sea[u/(u+d)] = " << rs;
#endif

  // compute the corrected valence and sea quark PDFs:
  double uv_c =       uv        / ( 1 + delta*rv);
  double dv_c = (dv + uv*delta) / ( 1 + delta*rv);
  double us_c =       us        / ( 1 + delta*rs);
  double ds_c = (ds + us*delta) / ( 1 + delta*rs);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) << "Bodek-Yang PDF correction:";
  LOG("BodekYang", pDEBUG) << "uv: " << uv << " --> " << uv_c;
  LOG("BodekYang", pDEBUG) << "dv: " << dv << " --> " << dv_c;
  LOG("BodekYang", pDEBUG) << "us: " << us << " --> " << us_c;
  LOG("BodekYang", pDEBUG) << "ds: " << ds << " --> " << ds_c;
#endif

  // fill in and return the corrected PDFs:
  PDF_t corrected_pdfs;

  corrected_pdfs.uval = uv_c;
  corrected_pdfs.dval = dv_c;
  corrected_pdfs.usea = us_c;
  corrected_pdfs.dsea = ds_c;
  corrected_pdfs.str  = uncorrected_pdfs.str;
  corrected_pdfs.chm  = uncorrected_pdfs.chm;
  corrected_pdfs.bot  = uncorrected_pdfs.bot;
  corrected_pdfs.top  = uncorrected_pdfs.top;
  corrected_pdfs.gl   = uncorrected_pdfs.gl;

  return corrected_pdfs;
}
//____________________________________________________________________________
double BYPDF::DeltaDU(double x) const
{
// Computes the BY correction factor delta(d/u) 

  double d = fX0 + fX1 * x + fX2 * TMath::Power(x,2); 
  return d;
}
//____________________________________________________________________________
void BYPDF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYPDF::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYPDF::LoadConfig(void)
{

  GetParam( "BY-X0", fX0 ) ;
  GetParam( "BY-X1", fX1 ) ;
  GetParam( "BY-X2", fX2 ) ;

  GetParam( "PDF-Q2min", fQ2min ) ;

  // get the base PDF model (typically GRV9* LO)
  fBasePDFModel = 
    dynamic_cast<const PDFModelI *>(this->SubAlg("Uncorr-PDF-Set"));
}
//____________________________________________________________________________


