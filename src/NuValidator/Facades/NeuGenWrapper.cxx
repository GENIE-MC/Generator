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

#include "Facades/NeuGenWrapper.h"
#include "Facades/NGInitState.h"

using std::endl;
using std::cout;

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
  
  set_pdfset_(&pdfgroup,&pdfset,&q2min);   
}  
//____________________________________________________________________________
void NeuGenWrapper::SetDefaultConfig(void)
{
  set_default_parameters_();
}
//____________________________________________________________________________
float NeuGenWrapper::XSec(float e, NGInteraction * ni, NeuGenCuts * cuts)
{
  float val=0.;

  // Get Process and nucleus
  
  int         process = ni->GetProcess();
  NGNucleus_t nucleus = ni->GetNucleus();

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
XSecVsEnergy * NeuGenWrapper::XSecSpline (
      float emin, float emax, int nbins, NGInteraction * ni, NeuGenCuts * cuts)
{
  assert(emin < emax && nbins > 1);

  double * E    = new double[nbins];
  double * xsec = new double[nbins];

  double de = (emax - emin) / ( nbins - 1);

  for(int ie = 0; ie < nbins; ie++) {

      E[ie]    = emin + ie * de;
      xsec[ie] = this->XSec (E[ie], ni, cuts);

      cout << "xsec(" << E[ie] << ") = " << xsec[ie] << endl;
  }

  XSecVsEnergy * spline = new XSecVsEnergy(nbins, E, xsec);

  //delete [] E;
  //delete [] xsec;

  return spline;
}
//____________________________________________________________________________
XSecVsEnergy * NeuGenWrapper::ExclusiveXSecSpline(
    float emin, float emax, int nbins,
                    NGInteraction * ni, NGFinalState * final, NeuGenCuts * cuts)
{
  assert(emin < emax && nbins > 0);

  double * E    = new double[nbins];
  double * xsec = new double[nbins];

  double de = (emax - emin) / ( nbins - 1);

  for(int ie = 0; ie < nbins; ie++) {

      E[ie]    = emin + ie * de;
      xsec[ie] = this->ExclusiveXSec(E[ie], ni, final, cuts);

      cout << "xsec(" << E[ie] << ") = " << xsec[ie] << endl;
  }

  XSecVsEnergy * spline = new XSecVsEnergy(nbins, E, xsec);

  //delete [] E;
  //delete [] xsec;

  return spline;
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
