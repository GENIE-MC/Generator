//____________________________________________________________________________
/*!

\program gtestPDFLIB

\brief   Test interface to LHAPDF library

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created May 4, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________

#include <TNtuple.h>
#include <TFile.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/PartonDistributions/PDFModelI.h"
#include "Physics/PartonDistributions/LHAPDF5.h"
#include "Physics/PartonDistributions/PDF.h"

using namespace genie;

//___________________________________________________________________
int main(int /*argc*/, char ** /*argv*/)
{
  // x, Q2 range
  const int nQ2 = 3;
  const int nx  = 300;
  double Q2arr[nQ2] = { 1.0, 2.0, 4.0 };
  const double xmin_idx = -2;
  const double xmax_idx =  0;
  const double dx_idx   = (xmax_idx-xmin_idx)/(nx-1);

  // Output ntuple
  TNtuple * nt = new TNtuple("nt","pdfs","uv:dv:us:ds:s:g:x:Q2");

  // PDF model
  AlgFactory * algf = AlgFactory::Instance();
  const PDFModelI * pdfmodel = 
        dynamic_cast<const PDFModelI *> (
                    algf->GetAlgorithm("genie::LHAPDF5","GRVLO"));
  
  PDF pdf;
  pdf.SetModel(pdfmodel);

  // Extract PDFs
  for(int iq2 = 0; iq2 < nQ2; iq2++) {
    for(int ix = 0; ix < nx; ix++) {

      double Q2 = Q2arr[iq2];
      double x  = TMath::Power(10, xmin_idx + ix * dx_idx);

      pdf.Calculate(x, Q2);
      LOG("test", pINFO) << "PDFs:\n" << pdf;

      nt->Fill(pdf.UpValence(),
               pdf.DownValence(),
               pdf.UpSea(),
               pdf.DownSea(),
               pdf.Strange(),
               pdf.Gluon(),
               x, Q2);
    }
  }
  TFile f("./genie-pdfs.root","recreate");
  nt->Write("pdfs");

  f.Close();
  delete nt;

  return 0;
}
//___________________________________________________________________
