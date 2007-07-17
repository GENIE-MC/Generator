//____________________________________________________________________________
/*!

\program testPDFLIB

\brief   test interface to PDFLIB library

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 4, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <TNtuple.h>
#include <TFile.h>
#include <TMath.h>

#include "Messenger/Messenger.h"
#include "PDF/PDFModelI.h"
#include "PDF/PDFLIB.h"
#include "PDF/PDF.h"

using namespace genie;

void fill_ntuple(TNtuple * nt, double q2, 
                                   int nptype, int ngroup, int nset);

//___________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  //-- Q^2 for which I extract PDFs

  const int nq2 = 3;
  double q2[nq2] = { 1, 2, 4 };

  //-- PDF sets / check PDFLIB definitions

  const int npdfs = 4;
                             /*CTEQ*/ /*GRV98*/
  const int nptype[npdfs] = {  1,  1,  1,  1 };
  const int ngroup[npdfs] = {  4,  4,  5,  5 };
  const int nset[npdfs]   = { 32, 41, 12, 13 };

  TNtuple * ntuple = new TNtuple("ntuple","pdfs",
                          "uv:dv:us:ds:s:g:x:Q2:nptype:ngroup:nset");

  for(int ipdf = 0; ipdf < npdfs; ipdf++) {
   for(int iq2 = 0; iq2 < nq2; iq2++) 
     fill_ntuple(ntuple, q2[iq2], nptype[ipdf], ngroup[ipdf], nset[ipdf]);
  }

  TFile f("./genie-pdflib.root","recreate");
  ntuple->Write("pdflib");

  f.Close();
  delete ntuple;

  return 0;
}
//___________________________________________________________________
void fill_ntuple(TNtuple * nt, 
                         double Q2, int nptype, int ngroup, int nset)
{
  PDFModelI * pdfmodel = new PDFLIB();
  
  Registry * conf = new Registry();

  conf->Set("Nptype", nptype);
  conf->Set("Ngroup", ngroup);
  conf->Set("Nset",   nset);

  pdfmodel->Configure(*conf);

  PDF pdf;
  pdf.SetModel(pdfmodel);

  const int nx = 300;

  const double xmin_idx = -2;
  const double xmax_idx =  0;
  const double dx_idx   = (xmax_idx-xmin_idx)/(nx-1);

  for(int ix=0; ix<nx; ix++) {

      double x = TMath::Power(10, xmin_idx + ix * dx_idx);
      pdf.Calculate(x, Q2);

      LOG("Main", pINFO) << ENDL << pdf;

      nt->Fill(pdf.UpValence(),
               pdf.DownValence(),
               pdf.UpSea(),
               pdf.DownSea(),
               pdf.Strange(),
               pdf.Gluon(),
               x, Q2, nptype, ngroup, nset);
  }

  delete pdfmodel;
  delete conf;
}
//___________________________________________________________________

