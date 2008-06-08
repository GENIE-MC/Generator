//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author Anselmo Meregaglia <anselmo.meregaglia@cern.ch>, IPHC Strasbourg
        Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, STFC, Rutherford Lab
	22-Jan-2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 08, 2008 - AM,CA
   This GENIE wrapper class for the LHAPDF library was first added in 2.3.1 
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/GBuild.h"
#include "PDF/LHAPDF.h"
#include "Messenger/Messenger.h"

using namespace genie;

#ifdef __GENIE_LHAPDF_ENABLED__
//
//the actual LHAPDF calls
//
extern "C" void initpdfsetbyname_(char name[20]);
extern "C" void initpdfmem_(int nmem);
extern "C" void evolvepdf_(double x, double Q, double f[12]);

#endif

//____________________________________________________________________________
LHAPDF::LHAPDF() :
PDFModelI("genie::LHAPDF")
{
  this->Initialize();
}
//____________________________________________________________________________
LHAPDF::LHAPDF(string config) :
PDFModelI("genie::LHAPDF", config)
{
  LOG("PDF", pDEBUG) 
      << "LHAPDF configuration:\n " << *fConfig;  

  this->Initialize();
}
//____________________________________________________________________________
LHAPDF::~LHAPDF() 
{ 

}
//____________________________________________________________________________
double LHAPDF::UpValence(double x, double q2) const
{
  return AllPDFs(x,q2).uval;
}
//____________________________________________________________________________
double LHAPDF::DownValence(double x, double q2) const
{
  return AllPDFs(x,q2).dval;
}
//____________________________________________________________________________
double LHAPDF::UpSea(double x, double q2) const
{
  return AllPDFs(x,q2).usea;
}
//____________________________________________________________________________
double LHAPDF::DownSea(double x, double q2) const
{
  return AllPDFs(x,q2).dsea;
}
//____________________________________________________________________________
double LHAPDF::Strange(double x, double q2) const
{
  return AllPDFs(x,q2).str;
}
//____________________________________________________________________________
double LHAPDF::Charm(double x, double q2) const
{
  return AllPDFs(x,q2).chm;
}
//____________________________________________________________________________
double LHAPDF::Bottom(double x, double q2) const
{
  return AllPDFs(x,q2).bot;
}
//____________________________________________________________________________
double LHAPDF::Top(double x, double q2) const
{
  return AllPDFs(x,q2).top;
}
//____________________________________________________________________________
double LHAPDF::Gluon(double x, double q2) const
{
  return AllPDFs(x,q2).gl;
}
//____________________________________________________________________________
PDF_t LHAPDF::AllPDFs(double x, double q2) const
{
  PDF_t pdf;

  LOG("PDF", pDEBUG) << "Getting LHAPDF pdfs at x = " << x << ", q2 = " << q2;

#ifdef __GENIE_LHAPDF_ENABLED__

  double f[12];
//  evolvePDF_(x,q2,f);
  evolvepdf_(x,q2,f);

  // set PDF_t
  pdf.uval = f[2]-f[-2];
  pdf.dval = f[1]-f[-1];
  pdf.usea = f[-2];
  pdf.dsea = f[-1];
  pdf.str  = f[3];
  pdf.chm  = f[4];
  pdf.bot  = f[5];
  pdf.top  = f[6];
  pdf.gl   = f[0];

#else
  LOG("PDF", pERROR) 
     << " ** LHAPDF was not enabled during the GENIE configuration";
  LOG("PDF", pERROR) 
     << " ** To use the LHAPDF library use"
     << " --enable-lhapdf --with-lhapdf-lib=/some/path/libLHAPDF.a";
#endif

  return pdf;                                               
}
//____________________________________________________________________________
void LHAPDF::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->SetPDFSetFromConfig(); 
}
//____________________________________________________________________________
void LHAPDF::Configure(string config)
{
  Algorithm::Configure(config);
  this->SetPDFSetFromConfig(); 
}
//____________________________________________________________________________
void LHAPDF::SetPDFSetFromConfig(void) const
{
// Get configuration (pdf member/set) from registry.
// (for definitions, have a look at LHAPDF manual)

#ifdef __GENIE_LHAPDF_ENABLED__

  int nmem = -1; // PDF member
  int nset   = -1; // PDF set 

  fConfig->Get("Nmember", nmem);
  fConfig->Get("Nset",   nset  );

  // Call pdfset from the fortran LHAPDF library to set PDF set/member

  char   param[20];

  if(nset==1)
    strcpy(param,"GRV98= LO");

  initpdfsetbyname_(param); // call pdfset from the fortran LHAPDF library
  initpdfmem_(nmem);
//  InitPDFsetByName_(param); // call pdfset from the fortran LHAPDF library
//  InitPDFmem_(nmem);

#else
  LOG("PDF", pERROR) 
     << " ** LHAPDF was not enabled during the GENIE configuration";
  LOG("PDF", pERROR) 
     << " ** To use the LHAPDF library use"
     << " --enable-lhapdf --with-lhapdf-lib=/some/path/libLHAPDF.a";
#endif
}
//____________________________________________________________________________

