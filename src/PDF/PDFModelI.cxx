//____________________________________________________________________________
/*!

\class    genie::PDFModelI

\brief    Pure abstract base class. Defines the PDFModelI interface to be
          implemented by any algorithmic class computing Parton Density
          Functions.

          Wrapper classes to existing Parton Density Function
          libraries (PDFLIB, LHAPDF) should also adopt this interface so as
          to be integrated into the GENIE framework.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/ 
//____________________________________________________________________________

#include "PDF/PDFModelI.h"

using namespace genie;

//____________________________________________________________________________
PDFModelI::PDFModelI() :
Algorithm()
{

}
//____________________________________________________________________________
PDFModelI::PDFModelI(const char * param_set) :
Algorithm(param_set)
{

}
//____________________________________________________________________________
PDFModelI::~PDFModelI() 
{ 

}
//____________________________________________________________________________



