//____________________________________________________________________________
/*!

\class    genie::PDFLIB

\brief    PDFLIB library interface.

          Concrete implementation of the PDFModelI interface.
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 06, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "PDF/PDFLIB.h"
#include "Messenger/Messenger.h"

using namespace genie;

//the actual PDFLIB calls

extern "C" void pdfset_(const char param[20][20], double val[20]);

extern "C" void structm_(double *, double *, double *, double *, double *, 
                  double *, double *, double *, double *, double *, double *);

//____________________________________________________________________________
PDFLIB::PDFLIB() :
PDFModelI("genie::PDFLIB")
{
  this->Initialize();
}
//____________________________________________________________________________
PDFLIB::PDFLIB(string config) :
PDFModelI("genie::PDFLIB", config)
{
  LOG("PDF", pDEBUG) << "PDFLIB configuration: " << ENDL << *fConfig;  

  this->Initialize();
}
//____________________________________________________________________________
PDFLIB::~PDFLIB() 
{ 

}
//____________________________________________________________________________
void PDFLIB::Initialize(void) const
{
  char   param[20][20];
  double val[20];

  strcpy(param[0], "Init0");

  pdfset_(param, val); // call pdfset from the fortran PDFLIB library
}
//____________________________________________________________________________
void PDFLIB::SetPDFSetFromConfig(void) const
{
  // Get configuration (particle type, pdf group/set) from registry.
  // (for definitions, have a look at PDFLIB manual)

  int nptype = -1; // particle type
  int ngroup = -1; // PDF author group
  int nset   = -1; // PDF set --within PDF author group--

  fConfig->Get("Nptype", nptype);
  fConfig->Get("Ngroup", ngroup);
  fConfig->Get("Nset",   nset  );

  // Call pdfset from the fortran PDFLIB library to set PDF group/set

  char   param[20][20];
  double val[20];

  strcpy(param[0],"Nptype");
  val[0] = nptype;
  strcpy(param[1],"Ngroup");
  val[1] = ngroup;
  strcpy(param[2],"Nset");
  val[2] = nset;

  pdfset_(param, val);
}
//____________________________________________________________________________
double PDFLIB::UpValence(double x, double q2) const
{
  return AllPDFs(x,q2).uval;
}
//____________________________________________________________________________
double PDFLIB::DownValence(double x, double q2) const
{
  return AllPDFs(x,q2).dval;
}
//____________________________________________________________________________
double PDFLIB::UpSea(double x, double q2) const
{
  return AllPDFs(x,q2).usea;
}
//____________________________________________________________________________
double PDFLIB::DownSea(double x, double q2) const
{
  return AllPDFs(x,q2).dsea;
}
//____________________________________________________________________________
double PDFLIB::Strange(double x, double q2) const
{
  return AllPDFs(x,q2).str;
}
//____________________________________________________________________________
double PDFLIB::Charm(double x, double q2) const
{
  return AllPDFs(x,q2).chm;
}
//____________________________________________________________________________
double PDFLIB::Bottom(double x, double q2) const
{
  return AllPDFs(x,q2).bot;
}
//____________________________________________________________________________
double PDFLIB::Top(double x, double q2) const
{
  return AllPDFs(x,q2).top;
}
//____________________________________________________________________________
double PDFLIB::Gluon(double x, double q2) const
{
  return AllPDFs(x,q2).gl;
}
//____________________________________________________________________________
PDF_t PDFLIB::AllPDFs(double x, double q2) const
{
  double scale = TMath::Sqrt( TMath::Abs(q2) ); // QCD scale

  // call structm from the fortran PDFLIB library

  double uval, dval, usea, dsea, str, chm, bot, top, gl;

  structm_(&x, &scale, &uval, &dval, &usea, &dsea, 
                                                &str, &chm, &bot, &top, &gl);

  PDF_t pdf;

  pdf.uval = uval;
  pdf.dval = dval;
  pdf.usea = usea;
  pdf.dsea = dsea;
  pdf.str  = str;
  pdf.chm  = chm;
  pdf.bot  = bot;
  pdf.top  = top;
  pdf.gl   = gl;

  return pdf;                                               
}
//____________________________________________________________________________
void PDFLIB::Configure(const Registry & config)
{
  Algorithm::Configure(config);

  this->Initialize();          // re-initialize PDFLIB
  this->SetPDFSetFromConfig(); // call PDFLIB's pdfset with new configuration
}
//____________________________________________________________________________
void PDFLIB::Configure(string config)
{
  Algorithm::Configure(config);

  this->Initialize();          // re-initialize PDFLIB
  this->SetPDFSetFromConfig(); // call PDFLIB's pdfset with new configuration
}
//____________________________________________________________________________

