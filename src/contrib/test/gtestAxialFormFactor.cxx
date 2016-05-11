//____________________________________________________________________________
/*!

\program gtestAxialFormFactor

\brief   Program used for testing / debugging the axial form factor

\author   Aaron Meyer <asmeyer2012 \at uchicago.edu>
          based off testElFormFactors by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created August 20, 2013

\syntax  gtestAxialFormFactor [-l] [-o output_filename_prefix]
                              [-c min,max,inc[,min,max,inc...]]

         []  denotes an optionan argument
         -l  switch which evaluates over logarithmic Q^2
         -o  output filename prefix for root file output
         -c  allows scanning over z-expansion coefficients
             -- must give a multiple of 3 numbers as arguments: min, max, increment
             -- will scan over at most MAX_COEF coefficients and fill ntuple tree with all entries
             -- can change MAX_COEF and recompile as necessary
             -- to skip scanning over a coefficient, put max < min for that coefficient
             note that Kmax in UserPhysicsOptions should match the number of scanned coeffs

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
//#include <TNtupleD.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "GHEP/GHepParticle.h"
#include "LlewellynSmith/AxialFormFactorModelI.h"
#include "LlewellynSmith/AxialFormFactor.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"

// number of coefficient values stored in tree
#define MAX_COEF 3
// Q^2 ranges (GeV^2) and number of entries for tree
#define Q2_LIN     4
#define Q2_LOG_MIN 0.1
#define Q2_LOG_MAX 100
#define N_TREE     100

using namespace genie;
using std::string;
using std::ostringstream;

struct tdata // addresses of data which will be added to tree
{
 int*     pKmax;
 int*     pmod;
 double*  pQ2;
 double*  pz;
 double*  pFA;
 double*  pcAn; // passed as just the pointer to array
};

void GetCommandLineArgs (int argc, char ** argv);
void CalculateFormFactor(AxialFormFactor axff,            \
                         TTree * affnt, tdata params,     \
                         AlgFactory * algf, Registry * r, \
                         const Registry * gc,             \
                         Interaction* interaction, double t0);//, int Kmax);

bool IncrementCoefficients(double* coefmin, double* coefmax, double* coefinc, \
                           int kmaxinc, tdata params, AlgFactory * algf,
                           Registry* r, const Registry * gc);

string  kDefOptEvFilePrefix = "test.axialff";
string  gOptEvFilePrefix;
bool    gOptDoLogarithmicQ2;
double  gOptCoeffMin[MAX_COEF] = {0.};
double  gOptCoeffMax[MAX_COEF] = {0.};
double  gOptCoeffInc[MAX_COEF] = {0.};
int     gOptKmaxInc = 0;

//__________________________________________________________________________
int main(int argc, char ** argv)
{

  GetCommandLineArgs(argc,argv);

  if (gOptKmaxInc > 0)
  {
    LOG("testAxialFormFactor", pINFO) << "Found coefficient ranges:";
    for (int ik=0;ik<gOptKmaxInc;ik++)
    {
      LOG("testAxialFormFactor", pINFO) << "( min:"  << gOptCoeffMin[ik]  \
                                        << " , max:" << gOptCoeffMax[ik]  \
                                        << " , inc:" << gOptCoeffInc[ik] << " )";
    }
  }

  TTree * affnt = new TTree("affnt","Axial Form Factor Test Tree");
  AxialFormFactor axff;
  AlgFactory * algf = AlgFactory::Instance();
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  AlgId id("genie::ZExpAxialFormFactorModel","Default");
  const Algorithm * alg = algf->GetAlgorithm(id);
  Registry * r = confp->FindRegistry(alg);

  int Kmax;
  int mod;
  double Q2;
  double z;
  double FA;
  double cAn[MAX_COEF];
  affnt = new TTree("axff","Axial Form Factor Test");
  affnt->Branch("Kmax", &Kmax, "Kmax/I"); // number of coefficients
  affnt->Branch("mod",  &mod,  "mod/I" ); // model
  affnt->Branch("Q2",   &Q2,   "Q2/D"  ); // Q^2
  affnt->Branch("z" ,   &z ,   "z/D"   ); // z
  affnt->Branch("FA",   &FA,   "FA/D"  ); // F_A axial form factor
  affnt->Branch("cAn",   cAn,  "cAn[Kmax]/D"); // z-expansion coefficients

  // load all parameters into a struct to make passing them as arguments easier
  tdata params; //defined tdata struct above
  params.pKmax = &Kmax;
  params.pmod  = &mod;
  params.pQ2   = &Q2;
  params.pz    = &z;
  params.pFA   = &FA;
  params.pcAn  = cAn;

  // flag for single output or loop
  bool do_once = (gOptKmaxInc < 1);
  bool do_Q4limit = r->GetBoolDef("QEL-Q4limit", gc->GetBool("QEL-Q4limit"));

  if (gOptKmaxInc < 1) Kmax = r->GetIntDef("QEL-Kmax", gc->GetInt("QEL-Kmax")); 
  else                 Kmax = gOptKmaxInc;
  //if (do_Q4limit) Kmax = Kmax+4;
  if (do_Q4limit) LOG("testAxialFormFactor",pWARN) \
   << "Q4limit specified, note that coefficients will be altered";

  // initialize to zero
  for (int ip=0;ip<MAX_COEF;ip++)
  {
    params.pcAn[ip] = 0.;
  }

  // Get T0 from registry value
  double t0 = r->GetDoubleDef("QEL-T0", gc->GetDouble("QEL-T0")); 

  // make sure coefficients are getting incremented correctly
  if (MAX_COEF < gOptKmaxInc ||
     (gOptKmaxInc == 0 && MAX_COEF < Kmax) )
  {
    LOG("testAxialFormFactor",pWARN) \
      << "Too many coefficient ranges, reducing to " << MAX_COEF;
    Kmax        = MAX_COEF;
    gOptKmaxInc = MAX_COEF;
  }
  LOG("testAxialFormFactor",pWARN) << "pKmax = " << *params.pKmax;

  // set up interaction
  Interaction * interaction = 
         Interaction::QELCC(kPdgTgtFe56,kPdgProton,kPdgNuMu,0.02);

  if (do_once)
  { // calculate once and be done with it

    // pre-load parameters
    ostringstream alg_key;
    for (int ip=1;ip<Kmax+1;ip++)
    {
      alg_key.str("");
      alg_key << "QEL-Z_A" << ip;
      params.pcAn[ip-1] =
        r->GetDoubleDef(alg_key.str(), gc->GetDouble(alg_key.str())); 
      LOG("testAxialFormFactor",pINFO) << "Loading " << alg_key.str() \
       << " : " << params.pcAn[ip-1];
    }

    CalculateFormFactor(axff,affnt,params,algf,r,gc,interaction,t0);//,t0,Kmax);
  } else {
    // override and control coefficient values
    r->UnLock();
    r->Set("QEL-Kmax",*params.pKmax);

    ostringstream alg_key;
    for (int ip=1;ip<gOptKmaxInc+1;ip++)
    {
      alg_key.str("");
      alg_key << "QEL-Z_A" << ip;
      r->Set(alg_key.str(),gOptCoeffMin[ip-1]);
      params.pcAn[ip-1] = gOptCoeffMin[ip-1];
      LOG("testAxialFormFactor",pWARN) << "Set parameter: " << params.pcAn[ip-1];
    }
    algf->ForceReconfiguration();

    do
    { CalculateFormFactor(axff,affnt,params,algf,r,gc,interaction,t0); }//,t0,Kmax); }
    while (IncrementCoefficients(gOptCoeffMin,gOptCoeffMax,gOptCoeffInc,gOptKmaxInc,
                                 params,algf,r,gc));

  }

  TFile f((gOptEvFilePrefix + ".root").c_str(),"recreate");
  affnt->Write();
  f.Close();

  return 0;
}
//__________________________________________________________________________
void CalculateFormFactor(AxialFormFactor axff,
                         TTree * affnt, tdata params, \
                         AlgFactory * algf, Registry * r, \
                         const Registry * gc, \
                         Interaction* interaction, double t0) //, int Kmax)
{
  const AxialFormFactorModelI * dipole =
      dynamic_cast<const AxialFormFactorModelI *> (
        algf->GetAlgorithm("genie::DipoleAxialFormFactorModel", "Default"));
  const AxialFormFactorModelI * zexp =
      dynamic_cast<const AxialFormFactorModelI *> (
        algf->GetAlgorithm("genie::ZExpAxialFormFactorModel", "Default"));

  interaction->KinePtr()->SetQ2(0.);
  double t     = interaction->KinePtr()->q2();
  double tcut  = 9.0 * TMath::Power(constants::kPi0Mass, 2);
  double Q2    = 0.;
                                                                                            
  for(int iq=0; iq<N_TREE; iq++) { 
                                                                                            
   if (gOptDoLogarithmicQ2)
   {
     Q2 = TMath::Exp(  TMath::Log(double(Q2_LOG_MIN))
                     + (double(iq+1)/double(N_TREE))
                     * TMath::Log(double(Q2_LOG_MAX)/double(Q2_LOG_MIN)));
     //Q2 = TMath::Exp(((iq+1)*double(Q2_RANGE_LOG/N_TREE))*TMath::Log(10.)); // logarithmic
   } else { Q2 = (iq+1)*double(Q2_LIN)/double(N_TREE) ; } // linear 

   interaction->KinePtr()->SetQ2(Q2);
   *params.pQ2 = Q2;
                                                                                            
   // calculate z parameter used in expansion
   t = interaction->KinePtr()->q2();
   double znum  = TMath::Sqrt(tcut - t) - TMath::Sqrt(tcut - t0);
   double zden  = TMath::Sqrt(tcut - t) + TMath::Sqrt(tcut - t0);
   double zpar = znum/zden;
                                                                                            
   if (zpar != zpar) LOG("testAxialFormFactor", pERROR) << "Undefined expansion parameter";
   *params.pz = zpar;

   axff.SetModel(dipole);
   axff.Calculate(interaction);
   *params.pmod = 0;
   *params.pFA  = axff.FA();
   affnt->Fill();
                                                                                            
   axff.SetModel(zexp);
   axff.Calculate(interaction);
   *params.pmod = 1;
   *params.pFA = axff.FA();
   affnt->Fill();
   
  }
}
//__________________________________________________________________________
bool IncrementCoefficients(double* coefmin, double* coefmax, double* coefinc, \
          int kmaxinc, tdata params, \
          AlgFactory * algf, Registry* r, const Registry * gc)

{

  if (kmaxinc < 1)
  {
    LOG("testAxialFormFactor",pERROR) << "No coefficients to increment";
    return false;
  } else {

  ostringstream alg_key;
  int ip = -1;
  bool stopflag;
  do
  {
    if (ip > -1)
    { // a previous iteration went over max
      params.pcAn[ip] = coefmin[ip];
      r->Set(alg_key.str(),params.pcAn[ip]);
    }
    stopflag = true;
    alg_key.str("");                                      // reset sstream
    ip++;                                                 // increment index
    if (ip == kmaxinc) return false;                      // stop if gone too far
    alg_key << "QEL-Z_A" << ip+1;                         // get new name
    params.pcAn[ip] += coefinc[ip];                       // increment parameter
    r->Set(alg_key.str(),params.pcAn[ip]);
    if (params.pcAn[ip] > coefmax[ip]) stopflag=false;    

  } while (! stopflag); // loop


  } // if kmaxinc < 1

  algf->ForceReconfiguration();
  return true;

}
//__________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("testAxialFormFactor", pINFO) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("testAxialFormFactor", pINFO) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("testAxialFormFactor", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // logarithmic vs linear in q2
   if( parser.OptionExists('l') ) {                                             
     LOG("testAxialFormFactor", pINFO) << "Looping over logarithmic Q2";
     gOptDoLogarithmicQ2 = true;
   } else {                                                                     
     LOG("testAxialFormFactor", pINFO) << "Looping over linear Q2";
     gOptDoLogarithmicQ2 = false;
   } //-l                                                                      

  if( parser.OptionExists('c') ) {
    LOG("testAxialFormFactor", pINFO) << "Reading Coefficient Ranges";
    string coef = parser.ArgAsString('c');
    
    // split into sections of min,max,inc(rement)
    vector<string> coefrange = utils::str::Split(coef, ",");
    assert(coefrange.size() % 3 == 0);
    gOptKmaxInc = coefrange.size() / 3;
    LOG("testAxialFormFactor", pINFO) << "Number of ranges to implement : " << gOptKmaxInc;
    for (int ik = 0;ik<gOptKmaxInc;ik++)
    {
      gOptCoeffMin[ik] = atof(coefrange[ik*3  ].c_str());
      gOptCoeffMax[ik] = atof(coefrange[ik*3+1].c_str());
      gOptCoeffInc[ik] = atof(coefrange[ik*3+2].c_str());
    }
    
  }

}
