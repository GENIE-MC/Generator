//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenWrapper

\brief    NeuGEN's C++ Wrapper class

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#include <cassert>
#include <iostream>

#include <TSystem.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "Facades/NeuGenWrapper.h"
#include "Facades/NGInitState.h"

using std::endl;
using std::cout;

using genie::Messenger;
using genie::EventRecord;
using genie::GHepParticle;
using namespace genie::nuvld::facades;

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  namespace facades {
    ostream & operator << (ostream & stream, const NeuGenWrapper & conf)
    {
      conf.Print(stream);
      return stream;
    }
  }
 }
}    
//____________________________________________________________________________
NeuGenWrapper::NeuGenWrapper() :
_name("default")
{
  this->Init();
}
//____________________________________________________________________________
NeuGenWrapper::NeuGenWrapper(const char * name) :
_name(name)
{
  this->Init();
}
//____________________________________________________________________________
NeuGenWrapper::NeuGenWrapper(const NeuGenConfig * config) 
{
  this->Init();
  this->Reconfigure(config);
}
//____________________________________________________________________________
NeuGenWrapper::~NeuGenWrapper()
{
  if(_config) delete _config;
}
//____________________________________________________________________________
void NeuGenWrapper::Reconfigure(const NeuGenConfig * config)
{
  if(_config) delete _config;

  _config = new NeuGenConfig(config);

  int five  = 5;
  int six   = 6;
  int seven = 7;
  int eight = 8;
  int nine  = 9;
  
  float val = 0;

  //! Set Axial mass for QEL scattering  
  val = _config->MaQel();
  set_parameters_("QEL_MA",&six,&val);

  //! Set Axial mass for RES scattering
  val = _config->MaRes();
  set_parameters_("RES_MA",&six,&val);

  //! Set Axial mass for COH scattering
  val = _config->MaCoh();
  set_parameters_("COH_MA",&six,&val);

  //! Set QEL axial form factor at Q^2 = 0
  val = _config->Fa0Qel();
  set_parameters_("QEL_FA0",&seven,&val);

  //! Set Elastic scattering parameter
  val = _config->EtaQel();
  set_parameters_("QEL_ETA",&seven,&val);

  //! Set R-S Model parameter Omega
  val = _config->OmegaRes();
  set_parameters_("RES_OMEGA",&nine,&val);

  //! Set R-S Model parameter Z
  val = _config->ZRes();
  set_parameters_("RES_Z",&five,&val);

  //! Set Nuclear size scale param in COH scattering
  val = _config->R0Coh();
  set_parameters_("COH_R0",&six,&val);

  //! Set Re/Im for pion scattering amplitude  
  val = _config->REICoh();
  set_parameters_("COH_REI",&seven,&val);

  //! KNO Hadronization model parameters

  val = _config->KnoA(e_vp);
  set_parameters_("KNO_A1",&six,&val);
  val = _config->KnoA(e_vn);
  set_parameters_("KNO_A2",&six,&val);
  val = _config->KnoA(e_vbp);
  set_parameters_("KNO_A3",&six,&val);
  val = _config->KnoA(e_vbn);
  set_parameters_("KNO_A4",&six,&val);
  val = _config->KnoB();
  set_parameters_("KNO_B",&five,&val);
  
  //! DIS-RES tuning params for multiplicity = 2 states
  
  val = _config->DisRes(2, e_vp);
  set_parameters_("KNO_R112",&eight,&val);
  val = _config->DisRes(2, e_vn);
  set_parameters_("KNO_R122",&eight,&val);
  val = _config->DisRes(2, e_vbp);
  set_parameters_("KNO_R132",&eight,&val);
  val = _config->DisRes(2, e_vbn);
  set_parameters_("KNO_R142",&eight,&val);

  //! DIS-RES tuning params for multiplicity = 3 states

  val = _config->DisRes(3, e_vp);
  set_parameters_("KNO_R113",&eight,&val);
  val = _config->DisRes(3, e_vn);
  set_parameters_("KNO_R123",&eight,&val);
  val = _config->DisRes(3, e_vbp);
  set_parameters_("KNO_R133",&eight,&val);
  val = _config->DisRes(3, e_vbn);
  set_parameters_("KNO_R143",&eight,&val);

  float q2min = 0.;
  
  int pdfgroup = _config->PdfGroup();
  int pdfset   = _config->PdfSet();
  
  LOG("NeuGen", pINFO) << "Reconfiguring: pdfgroup [" << pdfgroup <<
  "], pdset [" << pdfset <<"]";
  set_pdfset_(&pdfgroup,&pdfset,&q2min);   
}  
//____________________________________________________________________________
void NeuGenWrapper::SetDefaultConfig(void)
{
  int  nchar = 0;
  int  nver = 2;
  bool ok    = false;
 
  // Get the user-prefered NeuGEN configuration set from the GNEUGENCONF env. 
  // variable and initialize NeuGEN with it.
  // If the env.variable is not set, use the "MODBYEF2" NeuGEN configuration.
  // For the names of NeuGEN configuration sets, see the NeuGEN documentation.

  string gneugenconf =
              ( gSystem->Getenv("GNEUGENCONF") ? 
                                gSystem->Getenv("GNEUGENCONF") : "MODBYRS");
  nchar = gneugenconf.size();

  string gneugenconfv =
              ( gSystem->Getenv("GNEUGENCONFV") ? 
                                gSystem->Getenv("GNEUGENCONFV") : "2");

  LOG("NeuGen", pINFO)
          << "Using NeuGEN configuration set: [" << gneugenconf << "]" << " version: ["
          << gneugenconfv << "]" ;
  char * conf = (char*)gneugenconf.c_str(); // get what NeuGEN expects
  char * confv = (char*)gneugenconfv.c_str(); // get what NeuGEN expects
  int  vrs = atoi(confv); 
  LOG("NeuGen", pINFO) << "version " << vrs;

  initialize_configuration_(conf, &nchar, &nver, &ok);
  assert(ok);
  print_configuration_();

//  set_default_parameters_();
}
//____________________________________________________________________________
float NeuGenWrapper::XSec(float e, NGInteraction * ni, NeuGenCuts * cuts)
{
  float val=0.;

  // Get Process and nucleus  
  int         process = ni->GetProcess();
  NGNucleus_t nucleus = ni->GetNucleus();

  // Set cuts  
  if(cuts)  this->SetCuts(cuts);
  else      this->NoCuts();
  
  // Get cross section
  int nuc = nucleus;
  sig_value_(&e,&nuc,&process,&val);

  return val;
}
//____________________________________________________________________________
float NeuGenWrapper::ExclusiveXSec(float e, NGInteraction * ni,
                                     NGFinalState * final, NeuGenCuts * cuts )
{
  float val=0.;

  // Get Process  and nucleus  
  int         process = ni->GetProcess(final);
  NGNucleus_t nucleus = ni->GetNucleus();

  writestate_(&process);

  // Set cuts  
  if(cuts)  this->SetCuts(cuts);
  else      this->NoCuts();

  // Get cross section
  int nuc = nucleus;
  sig_value_(&e,&nuc,&process,&val);

  return val;
}
//____________________________________________________________________________
float NeuGenWrapper::DiffXSec(float e, NGKineVar_t  kid,
                              float kval, NGInteraction *ni, NeuGenCuts * cuts)
{
  float val=0.;

  // Get Process  and nucleus  
  int         process = ni->GetProcess();
  NGNucleus_t nucleus = ni->GetNucleus();

  writestate_(&process);

  // Set cuts  
  if(cuts)  this->SetCuts(cuts);
  else      this->NoCuts();
  
  // Get cross section  
  int kvid = kid;
  int nuc  = nucleus;

  dsig_value_(&e,&kvid,&kval,&nuc,&process,&val);

  return val;
}
//____________________________________________________________________________
float NeuGenWrapper::ExclusiveDiffXSec(float e,NGKineVar_t  kid,
           float kval, NGInteraction *ni, NGFinalState * final, NeuGenCuts * cuts)
{
  float val=0.;

  // Get Process  and nucleus  
  int         process = ni->GetProcess(final);
  NGNucleus_t nucleus = ni->GetNucleus();

  writestate_(&process);

  // Set cuts  
  if(cuts)  this->SetCuts(cuts);
  else      this->NoCuts();

  // Get cross section
  int nuc  = nucleus;
  int kvid = kid;

  dsig_value_(&e,&kvid,&kval,&nuc,&process,&val);

  return val;
}
//____________________________________________________________________________
float NeuGenWrapper::eDiff2Xsec(float e, NGKineVar_t  kv1, float kval1,
         NGKineVar_t  kv2, float kval2, NGInteraction * ni, NeuGenCuts * cuts)
{
  float val=0.;

  // Get Process  and nucleus
  int         process = ni->GetProcess();
  NGNucleus_t nucleus = ni->GetNucleus();

  writestate_(&process);

  // Set cuts
  if(cuts)  this->SetCuts(cuts);
  else      this->NoCuts();

  // Get cross section
  int kvid1 = kv1;
  int kvid2 = kv2;  
  int nuc   = nucleus;

  LOG("NeuGen", pINFO)
      << "\nComputing electron d^2xsec / da*db with a = "
         << NGKineVar::AsString(kv1) << ", b = " << NGKineVar::AsString(kv2);

  if(kv1 == e_costh && kv2 == e_ep) {

     ddsig_e_value_(&e,&kvid1,&kval1,&kvid2,&kval2,&nuc,&process,&val);
     return val;
  }
  
  LOG("NeuGen", pWARN) << "Can't compute e cross section / return 0";
  return 0.;
}
//____________________________________________________________________________
void NeuGenWrapper::SetCuts(NeuGenCuts * cuts)
{
// Set Cuts

  float kvmin = 0.;
  float kvmax = 0.;
  int   kvid  = 0;

  NGKineVar_t kv;

  if(cuts->IsCutSet()) {

    kv    = cuts->CutKVId();
    kvmin = cuts->KVMin();
    kvmax = cuts->KVMax();
    kvid  = kv;

  } else {
    kvmin = -9999999.;
    kvmax =  9999999.;
  }

  int  procmask = cuts->ProcMask();
  bool sumqel   = cuts->SumQel();
  bool sumres   = cuts->SumRes();
  bool sumdis   = cuts->SumDis();

  set_kvcuts_(&kvid,&kvmin,&kvmax);

  set_masks_(&procmask,&sumqel,&sumres,&sumdis);
}
//____________________________________________________________________________
void NeuGenWrapper::NoCuts(void)
{
// Turn off cuts

  int   procmask =  0;
  int   kvid     =  0;
  float kvmin    = -9999999.;
  float kvmax    =  9999999.;
  bool  valid    =  true;

  set_kvcuts_(&kvid,&kvmin,&kvmax);

  set_masks_(&procmask,&valid,&valid,&valid);
}
//____________________________________________________________________________
float NeuGenWrapper::StrucFunc(float x, float Q2, int A,
                      NGInitState_t init, NGCcNc_t ccnc, int isf, int raw_dis)
{
// raw_dis :  1 = calculate sf from y-dependence subject to normal
//                neugen cuts and flags (i.e. turn off certain proceses)
//            2 = calculate sf purely from the DIS model with no
//                resonance contribution or DIS subtraction
// isf: 1-6, the structure function number F1-F6
                    
  float f1=0., f2=0., f3=0., f4=0., f5=0., f6=0.;

  int  im    = -1;
  int  iccnc = -1;
  bool isNu  = false;

  int   lA  = A;
  int   lrd = raw_dis;
  float lx  = x;
  float lQ2 = Q2;
  
  switch(init) {
      case e_vp:   im =  1; isNu = true;   break;
      case e_vn:   im =  2; isNu = true;   break;
      case e_vN:   im =  5; isNu = true;   break;
      case e_vbp:  im =  3; isNu = true;   break;
      case e_vbn:  im =  4; isNu = true;   break;
      case e_vbN:  im =  6; isNu = true;   break;
      case e_lp:   im =  1; isNu = false;  break;
      case e_ln:   im =  2; isNu = false;  break;
      case e_lN:   im =  3; isNu = false;  break;
      default:     im = -1; isNu = false;
  }
  switch(ccnc) {
      case e_cc: iccnc = 1; break;
      case e_nc: iccnc = 2; break;
      default:   iccnc = 0; break;
  }
   
  if(isNu) {

     SLOG("NeuGen", pINFO)
                << "IM = " << im << ", CCNC = " << iccnc << ", A = " << A
                << ", Q2 = " << Q2 << ", x = " << x;
                
     nu_structurefunctions_(&im, &iccnc, &lA, &lx, &lQ2, &lrd, &f1,&f2,&f3,&f4,&f5,&f6);

     SLOG("NeuGen", pINFO)
                <<  " F1 = " << f1 << ", F2 = " << f2 << ", F3 = " << f3
                << ", F4 = " << f4 << ", F5 = " << f5 << ", F6 = " << f6;
  } else {                                
     e_structurefunctions_(&im, &lA, &lx, &lQ2, &lrd, &f1,&f2);
  }

  switch(isf) {
      case 1:  return f1; break;
      case 2:  return f2; break;
      case 3:  return f3; break;
      case 4:  return f4; break;
      case 5:  return f5; break;
      case 6:  return f6; break;
      default: return 0.;
 }
 return 0.;     
}
//____________________________________________________________________________
Spline * NeuGenWrapper::XSecSpline (
      float emin, float emax, int nbins, NGInteraction * ni, NeuGenCuts * cuts)
{
  SLOG("NeuGen", pINFO) << "Building XSec Spline";

  assert(emin < emax && nbins > 1);

  double * E    = new double[nbins];
  double * xsec = new double[nbins];

  double de = (emax - emin) / ( nbins - 1);

  for(int ie = 0; ie < nbins; ie++) {

      E[ie]    = emin + ie * de;
      xsec[ie] = this->XSec (E[ie], ni, cuts);

      LOG("NeuGen", pINFO) << "xsec(" << E[ie] << ") = " << xsec[ie];
  }

  Spline * spline = new Spline(nbins, E, xsec);

  //delete [] E;
  //delete [] xsec;

  return spline;
}
//____________________________________________________________________________
Spline * NeuGenWrapper::ExclusiveXSecSpline(
    float emin, float emax, int nbins,
                    NGInteraction * ni, NGFinalState * final, NeuGenCuts * cuts)
{
  SLOG("NeuGen", pINFO) << "Building Exclusive XSec Spline";

  assert(emin < emax && nbins > 0);

  double * E    = new double[nbins];
  double * xsec = new double[nbins];

  double de = (emax - emin) / ( nbins - 1);

  for(int ie = 0; ie < nbins; ie++) {

      E[ie]    = emin + ie * de;
      xsec[ie] = this->ExclusiveXSec(E[ie], ni, final, cuts);

      LOG("NeuGen", pINFO) << "xsec(" << E[ie] << ") = " << xsec[ie];
  }

  Spline * spline = new Spline(nbins, E, xsec);

  //delete [] E;
  //delete [] xsec;

  return spline;
}
//____________________________________________________________________________
Spline * NeuGenWrapper::StrucFuncSpline(
   NGKineVar_t xvar, float varmin, float varmax, int nbins, float fixedvar,
             int A, NGInitState_t init, NGCcNc_t ccnc, NGSF_t sf, int raw_dis)
{
  SLOG("NeuGen", pINFO) << "Building Structure Function Spline";

  assert(varmin < varmax && nbins > 0);

  Spline * spl = 0;

  float * SF   = new float[nbins];
  float * var  = new float[nbins];
  float   dvar = (varmax - varmin) / (nbins - 1);

  // tmp - to be removed when non raw-DIS SFs can be computed in NeuGEN.
  // For now, force raw-DIS so that NeuGEN does not quit.
  if(raw_dis != 2) {
    LOG("NeuGen", pWARN)
             << "RAW-DIS = " << raw_dis
                            << " -> switched to 2: forcing raw-dis SFs";
    raw_dis = 2;
  }

  switch(xvar) {
    
  case (e_qqs) :
    //----- compute: structrure function = f(Q2) at fixed x
    
    for(int i = 0; i < nbins; i++) {
      var [i]  = varmin + i * dvar; // Q2
      float Q2 = (float) var [i];
      if(sf == e_F2) {
        SF [i] = this->StrucFunc(fixedvar, Q2, A, init, ccnc, 2, raw_dis);
      } else
      if (sf == e_xF3) {
        SF [i] = this->StrucFunc(fixedvar, Q2, A, init, ccnc, 3, raw_dis);
        SF [i] *= fixedvar;
      } else
        SF[i] = 0;
        
      LOG("NeuGen", pINFO)
           << "SF(Q2 = " << var[i] << ", x = " << fixedvar << ") = " << SF[i];
    }
    spl = new Spline(nbins, var, SF);
    return spl;
    break;

  case (e_x)   :
    //----- compute: structrure function = f(x) at fixed Q2

    for(int i = 0; i < nbins; i++) {
      var [i] = varmin + i * dvar; // x
      float x = (float) var[i];
      if(sf == e_F2) {
        SF [i] = this->StrucFunc(x, fixedvar, A, init, ccnc, 2, raw_dis);
      } else
      if (sf == e_xF3) {
        SF [i] = this->StrucFunc(x, fixedvar, A, init, ccnc, 3, raw_dis);
        SF [i] *= var[i];
      } else
        SF[i] = 0;

      LOG("NeuGen", pINFO)
           << "SF(Q2 = " << fixedvar << ", x = " << var[i] << ") = " << SF[i];
    }
    spl = new Spline(nbins, var, SF);
    return spl;
    break;

  default:
    delete [] SF;
    delete [] var;
    spl = 0;    
  }
  return spl;
}               
//____________________________________________________________________________
void NeuGenWrapper::GenControl(char * flag, int var)
{
  bool ok = false;
  
  LOG("NeuGen", pINFO) 
             << "Setting NeuGEN control flag: [" << flag << "] --> " << var;

  gen_control_(flag, &var, &ok);
  //assert(ok);
}
//____________________________________________________________________________
EventRecord * NeuGenWrapper::GenerateEvent(int nupdgc, float E, int A, int Z)
{
  EventRecord * evrec = new EventRecord;

  //-- ask NeuGEN to generate an event

  SLOG("NeuGen", pINFO)
     << "Generating neutrino event for v(" << nupdgc
         << ") + tgt (A=" << A << ", Z=" << Z << ") at E = " << E << " GeV";
  generate_nu_event_(&nupdgc,&E,&A,&Z);

  //-- get the number of STDHEP entries NeuGENs event record

  int N = 0;
  get_n_stdhep_entries_(&N);

  SLOG("NeuGen", pINFO)
                 << "NeuGEN's STDHEP record contains: " << N << " entries";
  
  //-- get particle information and build GENIE event record

  for (int i=0; i<N; i++) {
    
    //-- get i^th particle status & pdg code
    int ist=0, pdgc=0;
    get_stdhep_particle_info_(&i, &ist, &pdgc);

    SLOG("NeuGen", pINFO)
        << "Adding NeuGEN's STDHEP entry = " << i << " with PDGC = " << pdgc;

    // set all Soudan2 codes to 'Rootino'...
    if(pdgc>1e7 && pdgc<1e8) {
       SLOG("NeuGen", pWARN)
                      << "GENIE resets Soudan-2 PDGC = " << pdgc << " to 0";
       pdgc=0; 
    }
    
    //-- get i^th particle mothers
    int mom1=0, mom2=0;
    get_stdhep_particle_mothers_(&i, &mom1, &mom2);
    
    //-- get i^th particle daughters
    int dau1=0, dau2=0;
    get_stdhep_particle_daughters_(&i, &dau1, &dau2);

    //-- translate NeuGEN/fortran indices [1,N] to GENIE/C++ indices [0,N-1]
    mom1--;
    mom2--;
    dau1--;
    dau2--;

    //-- get i^th particle 4-p
    double px=0, py=0, pz=0, E=0;
    get_stdhep_particle_p4_(&i, &px, &py, &pz, &E);

    //-- get i^th particle vertex
    double x=0, y=0, z=0, t=0;
    get_stdhep_particle_v4_(&i, &x, &y, &z, &z);

    //-- add the particle at GENIE's event record
    new ( (*evrec)[i] ) GHepParticle(
                         pdgc, ist, mom1,mom2,dau1,dau2, px,py,pz,E, x,y,z,t);    
  }

  // Since the event was not created from within the GENIE framework it will
  // be missing the summary (Interaction object). 
  evrec->AttachInteraction(0);
  
  return evrec;
}
//____________________________________________________________________________
void NeuGenWrapper::PrintEvent(void)
{
// Prints NeuGEN STDHEP common block.
// Note that for printing GENIE's EventRecord you must: "cout << event_rec;"

  int one = 1;
  heplst_(&one);
}
//____________________________________________________________________________
void NeuGenWrapper::Print(ostream & stream) const
{
  stream << "Name:..........." << _name  << endl;
  stream << "Current config: " << endl;
  stream << *_config << endl;
}
//____________________________________________________________________________
void NeuGenWrapper::Init(void) 
{
  this->SetDefaultConfig();
  _config = 0;
}  
//____________________________________________________________________________
