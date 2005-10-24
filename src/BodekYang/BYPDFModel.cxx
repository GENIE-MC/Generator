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
  LOG("BodekYang", pDEBUG) << "x = " << x << ", |q2| = " << q2;

  // get the base PDF model (typically GRV9* LO)
  const PDFModelI * base_pdf_model = dynamic_cast<const PDFModelI *>
        (this->SubAlg("base-pdf-model-alg-name","base-pdf-model-param-set"));

  // delegate the calculation request and get the uncorrected PDFs
  PDF_t uncorrected_pdfs = base_pdf_model->AllPDFs(x, q2);

  //--------- Correction to GRV9* LO up+down / valence+sea PDFs --------
  double uv  = uncorrected_pdfs.uval;
  double us  = uncorrected_pdfs.usea;
  double dv  = uncorrected_pdfs.dval;
  double ds  = uncorrected_pdfs.dsea;

  // compute correction factor delta(d/u)
  double delta = this->DeltaDU(x);
  LOG("BodekYang", pDEBUG) << "delta(d/u) = " << delta;

  // compute u/(u+d) ratios for both valence & sea quarks
  double rv = uv / (uv+dv);
  double rs = us / (us+ds);
  LOG("BodekYang", pDEBUG)
                 << "val[u/(u+d)] = " << rv << ", sea[u/(u+d)] = " << rs;

  // correct valence quark PDFs:
  double uv_c =       uv        / ( 1 + delta*rv);
  double dv_c = (dv + uv*delta) / ( 1 + delta*rv);

  // correct valence quark PDFs:
  double us_c =       us        / ( 1 + delta*rs);
  double ds_c = (ds + us*delta) / ( 1 + delta*rs);

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
// Computes correction factor delta(d/u) u+d / valence+sea PDFs
// The correction factor is computed as
//           delta(d/u) = cX0 + cX1 * x + cX2 * x^2 + cX3 * x^3
// where cX0,...cX3 are parameters that the algorithm finds from its XML
// configuration file.

  assert(
      fConfig->Exists("correction-const-X0") &&
              fConfig->Exists("correction-const-X1") &&
                      fConfig->Exists("correction-const-X2") &&
                               fConfig->Exists("correction-const-X3")
  );

  double cX0 = fConfig->GetDouble("correction-const-X0");
  double cX1 = fConfig->GetDouble("correction-const-X1");
  double cX2 = fConfig->GetDouble("correction-const-X2");
  double cX3 = fConfig->GetDouble("correction-const-X3");

  double x2 = x * x;
  double x3 = x * x2;

  double delta_d_u = cX0 + cX1*x + cX2*x2 + cX3*x3;

  return delta_d_u;
}
//____________________________________________________________________________
