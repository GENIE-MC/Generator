//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "PDF.h"

using namespace genie;

using std::endl;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const PDF & pdf_set)
  {
     pdf_set.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
PDF::PDF()
{
  this->Init();
}
//____________________________________________________________________________
PDF::PDF(const PDF & pdf_set)
{
  this->Copy(pdf_set);
}
//____________________________________________________________________________
PDF::~PDF()
{

}
//____________________________________________________________________________
void PDF::SetModel(const PDFModelI * model)
{
  this->Init();

  fModel = model;
}
//____________________________________________________________________________
void PDF::Calculate(double x, double q2)
{
  PDF_t pdfs = fModel->AllPDFs(x, q2);

  fUpValence   = pdfs.uval;
  fDownValence = pdfs.dval;
  fUpSea       = pdfs.usea;
  fDownSea     = pdfs.dsea;
  fStrange     = pdfs.str;
  fCharm       = pdfs.chm;
  fBottom      = pdfs.bot;
  fTop         = pdfs.top;
  fGluon       = pdfs.gl;
}
//____________________________________________________________________________
void PDF::ScaleValence(double kscale)
{
  fUpValence   *= kscale;
  fDownValence *= kscale;
}
//____________________________________________________________________________
void PDF::ScaleSea(double kscale)
{
  fUpSea       *= kscale;
  fDownSea     *= kscale;
  fStrange     *= kscale;
  fCharm       *= kscale;
  fBottom      *= kscale;
  fTop         *= kscale;
  fGluon       *= kscale;
}
//____________________________________________________________________________
void PDF::ScaleUpValence(double kscale)
{
  fUpValence *= kscale;
}
//____________________________________________________________________________
void PDF::ScaleDownValence(double kscale)
{
  fDownValence *= kscale;
}
//____________________________________________________________________________
void PDF::ScaleUpSea(double kscale)
{
  fUpSea *= kscale;
}
//____________________________________________________________________________
void PDF::ScaleDownSea(double kscale)
{
  fDownSea *= kscale;
}
//____________________________________________________________________________
void PDF::ScaleStrange(double kscale)
{
  fStrange *= kscale;
}
//____________________________________________________________________________
void PDF::ScaleCharm(double kscale)
{
  fCharm *= kscale;
}
//____________________________________________________________________________
void PDF::Reset(void)
{
  fUpValence   = 0.0;
  fDownValence = 0.0;
  fUpSea       = 0.0;
  fDownSea     = 0.0;
  fStrange     = 0.0;
  fCharm       = 0.0;
  fBottom      = 0.0;
  fTop         = 0.0;
  fGluon       = 0.0;
}
//____________________________________________________________________________
void PDF::Copy(const PDF & pdf_set)
{
  fModel       = pdf_set.fModel;

  fUpValence   = pdf_set.fUpValence;
  fDownValence = pdf_set.fDownValence;
  fUpSea       = pdf_set.fUpSea;
  fDownSea     = pdf_set.fDownSea;
  fStrange     = pdf_set.fStrange;
  fCharm       = pdf_set.fCharm;
  fBottom      = pdf_set.fBottom;
  fTop         = pdf_set.fTop;
  fGluon       = pdf_set.fGluon;
}
//____________________________________________________________________________
void PDF::Init(void)
{
  fModel = 0;

  fUpValence   = 0.0;
  fDownValence = 0.0;
  fUpSea       = 0.0;
  fDownSea     = 0.0;
  fStrange     = 0.0;
  fCharm       = 0.0;
  fBottom      = 0.0;
  fTop         = 0.0;
  fGluon       = 0.0;
}
//____________________________________________________________________________
void PDF::Print(ostream & stream) const
{
  stream << endl;
  stream << "UP-VAL....... " << fUpValence   << endl;
  stream << "DOWN-VAL..... " << fDownValence << endl;
  stream << "UP-SEA....... " << fUpSea       << endl;
  stream << "DOWN-SEA..... " << fDownSea     << endl;
  stream << "STRANGE...... " << fStrange     << endl;
  stream << "CHARM........ " << fCharm       << endl;
  stream << "BOTTOM....... " << fBottom      << endl;
  stream << "TOP.......... " << fTop         << endl;
  stream << "GLUON........ " << fGluon       << endl;
}
//____________________________________________________________________________

