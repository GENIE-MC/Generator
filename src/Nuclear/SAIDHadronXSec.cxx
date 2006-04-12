//____________________________________________________________________________
/*!

\class    genie::SAIDHadronXSec

\brief    Singleton class to load & serve SAID hadron cross section splines.
          Data and corrections provided by S.Dytmnan.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  April 12, 2006

*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <iostream>

#include <TSystem.h>

#include "Messenger/Messenger.h"
#include "Nuclear/SAIDHadronXSec.h"
#include "Numerical/Spline.h"

using std::cout;
using std::endl;

using namespace genie;

//____________________________________________________________________________
SAIDHadronXSec * SAIDHadronXSec::fInstance = 0;
//____________________________________________________________________________
SAIDHadronXSec::SAIDHadronXSec()
{
  this->LoadTables();
  fInstance = 0;
}
//____________________________________________________________________________
SAIDHadronXSec::~SAIDHadronXSec()
{
  cout << "SAIDHadronXSec singleton dtor: "
                      << "Deleting all hadron cross section splines" << endl;

  delete fSplPiplusPElas;
  delete fSplPiplusPInel;
  delete fSplPiminusPElas;
  delete fSplPiminusPInel;
  delete fSplPiminusPCEx; 
  delete fSplPiplusAbs;    
  delete fSplNPElas;      
  delete fSplNPInel;      
  delete fSplPPElas;      
  delete fSplPPInel;      
}
//____________________________________________________________________________
SAIDHadronXSec * SAIDHadronXSec::Instance()
{
  if(fInstance == 0) {
    LOG("FermiP", pINFO) << "SAIDHadronXSec late initialization";

    static SAIDHadronXSec::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new SAIDHadronXSec;
  }
  return fInstance;
}
//____________________________________________________________________________
void SAIDHadronXSec::LoadTables(void)
{
  //-- Get the directory with the SAID hadron cross section data (search for
  //   $GSAIDHADRONDATA or use default: $GENIE/data/hadron_xsec)

  string said_data_dir = (gSystem->Getenv("GSAIDHADRONDATA")) ?
            string(gSystem->Getenv("GSAIDHADRONDATA")) :
            string(gSystem->Getenv("GENIE")) + string("/data/hadron_xsec");

  LOG("SAID", pNOTICE)  
       << "Loading SAID hadron cross section data from: " << said_data_dir;

  //-- Build filenames

  string fdata_piplus_p_elas  = said_data_dir + "/pi+p-xselas-said.dat";
  string fdata_piplus_p_inel  = said_data_dir + "/pi+p-xsreac-said.dat";
  string fdata_piminus_p_elas = said_data_dir + "/pi-p-xselas-said.dat";
  string fdata_piminus_p_inel = said_data_dir + "/pi-p-xsreac-said.dat";
  string fdata_piminus_p_cex  = said_data_dir + "/pi-p-cex-xstot-said.dat";
  string fdata_pid_pp_tot     = said_data_dir + "/pid-pp-xstot-said.dat";
  string fdata_np_elas        = said_data_dir + "/np-xselas-said.dat";
  string fdata_np_inel        = said_data_dir + "/np-xsreac-said.dat";
  string fdata_pp_elas        = said_data_dir + "/pp-xselas-said.dat";
  string fdata_pp_inel        = said_data_dir + "/pp-xsreac-said.dat";

  //-- Make sure that all tables are available

  assert( ! gSystem->AccessPathName(fdata_piplus_p_elas.c_str())  );
  assert( ! gSystem->AccessPathName(fdata_piplus_p_inel.c_str())  );
  assert( ! gSystem->AccessPathName(fdata_piminus_p_elas.c_str()) );
  assert( ! gSystem->AccessPathName(fdata_piminus_p_inel.c_str()) );
  assert( ! gSystem->AccessPathName(fdata_piminus_p_cex.c_str())  );
  assert( ! gSystem->AccessPathName(fdata_pid_pp_tot.c_str())     );
  assert( ! gSystem->AccessPathName(fdata_np_elas.c_str())        );
  assert( ! gSystem->AccessPathName(fdata_np_inel.c_str())        );
  assert( ! gSystem->AccessPathName(fdata_pp_elas.c_str())        );
  assert( ! gSystem->AccessPathName(fdata_pp_inel.c_str())        );

  //-- Load the data into splines
  
  fSplPiplusPElas  = new Spline ( fdata_piplus_p_elas  );
  fSplPiplusPInel  = new Spline ( fdata_piplus_p_inel  );
  fSplPiminusPElas = new Spline ( fdata_piminus_p_elas );
  fSplPiminusPInel = new Spline ( fdata_piminus_p_inel );
  fSplPiminusPCEx  = new Spline ( fdata_piminus_p_cex  );
  fSplNPElas       = new Spline ( fdata_np_elas        );
  fSplNPInel       = new Spline ( fdata_np_inel        );
  fSplPPElas       = new Spline ( fdata_pp_elas        );   
  fSplPPInel       = new Spline ( fdata_pp_inel        );   

  // Especially for pion absorption cross section copy S.Dytmnan's correction.
  // His comments: 'model not valid above 600 MeV, slow linear falloff matches data'

  Spline xsabs(fdata_pid_pp_tot); 

  double   frac = 1; //?
  int      nk   = xsabs.NKnots();
  double * xsec = new double[nk];
  double * E    = new double[nk];

  for(int i=0; i<nk; i++) {
    xsec[i] = (E[i] < 600.) ? frac * xsabs.Evaluate(E[i]) : 
                              frac * xsabs.Evaluate(600.)*(1800.-E[i])/1200.;
  }
  fSplPiplusAbs = new Spline(nk,E,xsec); // corrected pi abs xsec

  delete [] xsec;
  delete [] E;

  LOG("SAID", pNOTICE) << "... Done loading SAID hadron cross section data";
};
//____________________________________________________________________________

