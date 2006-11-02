//____________________________________________________________________________
/*!

\class    genie::BYPDFModel

\brief    Computes corrected PDFs according to the Bodek-Yang model.

          Concrete implementation of the PDFModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  September 29, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "BodekYang/BYPDFModel.h"
#include "Messenger/Messenger.h"

using namespace genie;

//____________________________________________________________________________
BYPDFModel::BYPDFModel() :
PDFModelI("genie::BYPDFModel")
{

}
//____________________________________________________________________________
BYPDFModel::BYPDFModel(string config) :
PDFModelI("genie::BYPDFModel", config)
{

}
//____________________________________________________________________________
BYPDFModel::~BYPDFModel()
{

}
//____________________________________________________________________________
double BYPDFModel::UpValence(double x, double q2) const
{
  return AllPDFs(x,q2).uval;
}
//____________________________________________________________________________
double BYPDFModel::DownValence(double x, double q2) const
{
  return AllPDFs(x,q2).dval;
}
//____________________________________________________________________________
double BYPDFModel::UpSea(double x, double q2) const
{
  return AllPDFs(x,q2).usea;
}
//____________________________________________________________________________
double BYPDFModel::DownSea(double x, double q2) const
{
  return AllPDFs(x,q2).dsea;
}
//____________________________________________________________________________
double BYPDFModel::Strange(double x, double q2) const
{
  return AllPDFs(x,q2).str;
}
//____________________________________________________________________________
double BYPDFModel::Charm(double x, double q2) const
{
  return AllPDFs(x,q2).chm;
}
//____________________________________________________________________________
double BYPDFModel::Bottom(double x, double q2) const
{
  return AllPDFs(x,q2).bot;
}
//____________________________________________________________________________
double BYPDFModel::Top(double x, double q2) const
{
  return AllPDFs(x,q2).top;
}
//____________________________________________________________________________
double BYPDFModel::Gluon(double x, double q2) const
{
  return AllPDFs(x,q2).gl;
}
//____________________________________________________________________________
PDF_t BYPDFModel::AllPDFs(double x, double q2) const
{
  LOG("BodekYang", pDEBUG) 
       << "Inputs: x = " << x << ", |q2| = " << TMath::Abs(q2);

  if(TMath::Abs(q2) < fQ2min) q2=fQ2min;

  // get the uncorrected PDFs
  PDF_t uncorrected_pdfs = fBasePDFModel->AllPDFs(x, q2);
  double uv = uncorrected_pdfs.uval;
  double us = uncorrected_pdfs.usea;
  double dv = uncorrected_pdfs.dval;
  double ds = uncorrected_pdfs.dsea;

  // compute correction factor delta(d/u)
  double delta = this->DeltaDU(x);
  LOG("BodekYang", pDEBUG) << "delta(d/u) = " << delta;

  // compute u/(u+d) ratios for both valence & sea quarks
  double val = uv+dv;
  double sea = us+ds;
  double rv  = (val==0) ? 0. : uv/val;
  double rs  = (sea==0) ? 0. : us/sea;

  LOG("BodekYang", pDEBUG)
   << "valence[u/(u+d)] = " << rv << ", sea[u/(u+d)] = " << rs;

  // compute the corrected valence and sea quark PDFs:
  double uv_c =       uv        / ( 1 + delta*rv);
  double dv_c = (dv + uv*delta) / ( 1 + delta*rv);
  double us_c =       us        / ( 1 + delta*rs);
  double ds_c = (ds + us*delta) / ( 1 + delta*rs);

  LOG("BodekYang", pDEBUG) << "Bodek-Yang PDF correction:";
  LOG("BodekYang", pDEBUG) << "uv: " << uv << " --> " << uv_c;
  LOG("BodekYang", pDEBUG) << "dv: " << dv << " --> " << dv_c;
  LOG("BodekYang", pDEBUG) << "us: " << us << " --> " << us_c;
  LOG("BodekYang", pDEBUG) << "ds: " << ds << " --> " << ds_c;

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
double BYPDFModel::DeltaDU(double x) const
{
// Computes the BY correction factor delta(d/u) 

  double d = fX0 + fX1 * x + fX2 * TMath::Power(x,2); 
  return d;
}
//____________________________________________________________________________
void BYPDFModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYPDFModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BYPDFModel::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fX0 = fConfig->GetDoubleDef("X0", gc->GetDouble("BY-X0"));
  fX1 = fConfig->GetDoubleDef("X1", gc->GetDouble("BY-X1"));
  fX2 = fConfig->GetDoubleDef("X2", gc->GetDouble("BY-X2"));

  fQ2min = fConfig->GetDoubleDef("Q2min", gc->GetDouble("PDF-Q2min"));

  // get the base PDF model (typically GRV9* LO)
  fBasePDFModel = dynamic_cast<const PDFModelI *>(
                                      this->SubAlg("Uncorr-PDF-Set"));
}
//____________________________________________________________________________


