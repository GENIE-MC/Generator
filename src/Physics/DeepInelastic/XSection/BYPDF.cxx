//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/DeepInelastic/XSection/BYPDF.h"

using namespace genie;

//____________________________________________________________________________
BYPDF::BYPDF() : PDFModelI("genie::BYPDF") {}
//____________________________________________________________________________
BYPDF::BYPDF(string config) : PDFModelI("genie::BYPDF", config) {}
//____________________________________________________________________________
BYPDF::~BYPDF() {}
//____________________________________________________________________________
double BYPDF::UpValence(double x, double q2) const {
  return AllPDFs(x, q2).uval;
}
//____________________________________________________________________________
double BYPDF::DownValence(double x, double q2) const {
  return AllPDFs(x, q2).dval;
}
//____________________________________________________________________________
double BYPDF::UpSea(double x, double q2) const { return AllPDFs(x, q2).usea; }
//____________________________________________________________________________
double BYPDF::DownSea(double x, double q2) const { return AllPDFs(x, q2).dsea; }
//____________________________________________________________________________
double BYPDF::Strange(double x, double q2) const { return AllPDFs(x, q2).str; }
//____________________________________________________________________________
double BYPDF::Charm(double x, double q2) const { return AllPDFs(x, q2).chm; }
//____________________________________________________________________________
double BYPDF::Bottom(double x, double q2) const { return AllPDFs(x, q2).bot; }
//____________________________________________________________________________
double BYPDF::Top(double x, double q2) const { return AllPDFs(x, q2).top; }
//____________________________________________________________________________
double BYPDF::Gluon(double x, double q2) const { return AllPDFs(x, q2).gl; }
//____________________________________________________________________________
PDF_t BYPDF::AllPDFs(double x, double q2) const {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) << "Inputs: x = " << x
                           << ", |q2| = " << TMath::Abs(q2);
#endif
  if (TMath::Abs(q2) < fQ2min)
    q2 = fQ2min;

  // get the uncorrected PDFs
  PDF_t uncorrected_pdfs = fBasePDFModel->AllPDFs(x, q2);
  double uv = uncorrected_pdfs.uval;
  double us = uncorrected_pdfs.usea;
  double dv = uncorrected_pdfs.dval;
  double ds = uncorrected_pdfs.dsea;

  // The Bodek Yang model includes a parameter to scale the up and down sea
  // quark distributions by a percentage. These are added as additional degrees
  // of freedom in the model. The paper uses a 5% increase of the sea up and
  // down distribution and a decrease of the valance quarks. The value was not
  // obtained from the tune - we do not use it. We use the defalt of 0 for this
  // implementation whilst keeping the functionality it for the user to use.
  uv -= 2 * fUpScale * us ;
  dv -= 2 * fDownScale * ds ;
  us *= (1 + fUpScale);
  ds *= (1 + fDownScale);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) << "Scaling up (down), up sea (valance) contribution by = " << fUpScale;
  LOG("BodekYang", pDEBUG) << "Scaling up (down), down sea (valance) contribution by = " << fDownScale;
#endif
  
  // compute correction factor delta(d/u)
  double delta = fApplyDelta ? this->DeltaDU(x) : 0;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) << "delta(d/u) = " << delta;
#endif

  // compute u/(u+d) ratios for both valence & sea quarks
  double val = uv + dv;
  double sea = us + ds;
  double rv = (val == 0) ? 0. : uv / val;
  double rs = (sea == 0) ? 0. : us / sea;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BodekYang", pDEBUG) << "valence[u/(u+d)] = " << rv
                           << ", sea[u/(u+d)] = " << rs;
#endif

  // compute the corrected valence and sea quark PDFs:
  double uv_c = uv / (1 + delta * rv);
  double dv_c = (dv + uv * delta) / (1 + delta * rv);
  double us_c = us / (1 + delta * rs);
  double ds_c = (ds + us * delta) / (1 + delta * rs);

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
  corrected_pdfs.str = uncorrected_pdfs.str;
  corrected_pdfs.chm = uncorrected_pdfs.chm;
  corrected_pdfs.bot = uncorrected_pdfs.bot;
  corrected_pdfs.top = uncorrected_pdfs.top;
  corrected_pdfs.gl = uncorrected_pdfs.gl;

  return corrected_pdfs;
}
//____________________________________________________________________________
double BYPDF::DeltaDU(double x) const {
  // Computes the BY correction factor delta(d/u)

  double d = fX0 + fX1 * x + fX2 * TMath::Power(x, 2);
  return d;
}
//____________________________________________________________________________
void BYPDF::Configure(const Registry &config) {
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYPDF::Configure(string config) {
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYPDF::LoadConfig(void)
{

  GetParam( "BY-X0", fX0 ) ;
  GetParam( "BY-X1", fX1 ) ;
  GetParam( "BY-X2", fX2 ) ;
  GetParamDef( "BY-UpScale", fUpScale, 0. );
  GetParamDef( "BY-DownScale", fDownScale, 0. );
  GetParamDef( "BY-Delta", fApplyDelta, true );
  GetParam( "PDF-Q2min", fQ2min ) ;

  // get the base PDF model (typically GRV9* LO)
  fBasePDFModel =
      dynamic_cast<const PDFModelI *>(this->SubAlg("Uncorr-PDF-Set"));
}
//____________________________________________________________________________
