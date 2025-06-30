//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/XSection/HEDISPXSec.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

#include <TMath.h>

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
HEDISPXSec::HEDISPXSec() :
XSecAlgorithmI("genie::HEDISPXSec")
{

}
//____________________________________________________________________________
HEDISPXSec::HEDISPXSec(string config) :
XSecAlgorithmI("genie::HEDISPXSec", config)
{

}
//____________________________________________________________________________
HEDISPXSec::~HEDISPXSec()
{

}
//____________________________________________________________________________
double HEDISPXSec::XSec(
     const Interaction * interaction, KinePhaseSpace_t kps) const
{

  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Load SF tables 
  HEDISStrucFunc * sf_tbl = HEDISStrucFunc::Instance(fSFinfo);

  // W limits are computed using kinematics assumption.
  // The lower limit is tuneable because hadronization might have issues at low W (as in PYTHIA6).
  // To be consistent the cross section must be computed in the W range where the events are generated.
  const Kinematics  & kinematics = interaction -> Kine();
  const KPhaseSpace & ps         = interaction -> PhaseSpace();
  double W = kinematics.W();
  Range1D_t Wl  = ps.WLim();
  Wl.min = TMath::Max(Wl.min,fWmin);
  if      ( W<Wl.min ) return 0.;
  else if ( W>Wl.max ) return 0.;

  const InitialState & init_state = interaction -> InitState();

  double y     = kinematics.y();
  double Q2    = kinematics.Q2();
  double x     = kinematics.x();
  double E     = init_state.ProbeE(kRfLab);
  double Mnuc  = init_state.Tgt().HitNucMass();
  double Mlep2 = TMath::Power(interaction->FSPrimLepton()->Mass(),2);

  // Get F1,F2,F3 for particular quark channel and compute differential xsec
  SF_xQ2 sf = sf_tbl->EvalQrkSFLO( interaction, x, Q2 );
  double xsec = (fMassTerms) ? ds_dxdy_mass( sf, x, y, E, Mnuc, Mlep2 ) : ds_dxdy( sf, x, y );

  // If NLO is enable we compute sigma_NLO/sigma_LO. Then the quark xsec 
  // is multiplied by this ratio.
  // This is done because at NLO we can only compute the nucleon xsec. But
  // for the hadronization we need the different quark contributions.
  // This could be avoid if a NLO parton showering is introduced.
  if (fSFinfo.IsNLO && xsec>0.) {
    SF_xQ2 sflo  = sf_tbl->EvalNucSFLO(interaction,x,Q2);
    SF_xQ2 sfnlo = sf_tbl->EvalNucSFNLO(interaction,x,Q2);
    double lo  = (fMassTerms) ? ds_dxdy_mass( sflo, x, y, E, Mnuc, Mlep2 ) : ds_dxdy( sflo, x, y );
    if (lo>0.) {
      double nlo = (fMassTerms) ? ds_dxdy_mass( sfnlo, x, y, E, Mnuc, Mlep2 ) : ds_dxdy( sfnlo, x, y );
      xsec *= nlo / lo;
    }
  }

  // Compute the front factor
  double propagator = 0;
  if (interaction -> ProcInfo().IsWeakCC()) propagator = TMath::Power( fSFinfo.MassW*fSFinfo.MassW/(Q2+fSFinfo.MassW*fSFinfo.MassW), 2);
  else                                      propagator = TMath::Power( fSFinfo.MassZ*fSFinfo.MassZ/(Q2+fSFinfo.MassZ*fSFinfo.MassZ)/(1.-fSFinfo.Rho), 2);

  xsec *= kGF2/(2*kPi*x) * propagator;

  LOG("HEDISPXSec", pINFO) << "d2xsec/dxdy[FreeN] (x= " << x  << ", y= " << y << ", Q2= " << Q2 << ") = " << xsec;

  // The algorithm computes d^2xsec/dxdy. Check whether variable tranformation is needed
  if( kps!=kPSxQ2fE ) xsec *= utils::kinematics::Jacobian(interaction,kPSxQ2fE,kps);

  // If requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must have been included in the structure functions)
  int NNucl = (pdg::IsProton(init_state.Tgt().HitNucPdg())) ? init_state.Tgt().Z() : init_state.Tgt().N(); 
  xsec *= NNucl; 

  return xsec;

}
//____________________________________________________________________________
double HEDISPXSec::ds_dxdy(SF_xQ2 sf, double x, double y ) const
{

    // We neglect F4 and F5 and higher order terms.
    double term1 = y * ( x*y );
    double term2 = ( 1 - y );
    double term3 = ( x*y*(1-y/2) );

    LOG("HEDISPXSec", pDEBUG) << sf.F1 << "  " << sf.F2 << "  " << sf.F3;
    LOG("HEDISPXSec", pDEBUG) << term1*sf.F1 + term2*sf.F2 + term3*sf.F3;

    return fmax( term1*sf.F1 + term2*sf.F2 + term3*sf.F3 , 0.);

}
//____________________________________________________________________________
double HEDISPXSec::ds_dxdy_mass(SF_xQ2 sf, double x, double y, double e, double mt, double ml2 ) const
{

    double term1 = y * ( x*y + ml2/2/e/mt );
    double term2 = ( 1 - y - mt*x*y/2/e - ml2/4/e/e );
    double term3 = (x*y*(1-y/2) - y*ml2/4/mt/e);
    double term4 = x*y*ml2/2/mt/e + ml2*ml2/4/mt/mt/e/e;
    double term5 = -1.*ml2/2/mt/e;

    double F4 = 0.;
    double F5 = sf.F2/x;

    LOG("HEDISPXSec", pDEBUG) << sf.F1 << "  " << sf.F2 << "  " << sf.F3;
    LOG("HEDISPXSec", pDEBUG) << term1*sf.F1 + term2*sf.F2 + term3*sf.F3;

    return fmax( term1*sf.F1 + term2*sf.F2 + term3*sf.F3 + term4*F4 + term5*F5 , 0.);

}
//____________________________________________________________________________
double HEDISPXSec::Integral(const Interaction * interaction) const
{

  return fXSecIntegrator->Integrate(this,interaction);

}
//____________________________________________________________________________
bool HEDISPXSec::ValidProcess(const Interaction * interaction) const
{

  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsDeepInelastic()) return false;

  const InitialState & init_state = interaction -> InitState();
  int probe_pdg = init_state.ProbePdg();
  if(!pdg::IsLepton(probe_pdg)) return false;

  if(! init_state.Tgt().HitNucIsSet()) return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if(!pdg::IsNeutronOrProton(hitnuc_pdg)) return false;

  return true;
}
//____________________________________________________________________________
void HEDISPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISPXSec::LoadConfig(void)
{

  //-- load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // Minimum value of W (typically driven by hadronization limitation)
  GetParam("Xsec-Wmin",       fWmin);
  GetParam("Mass-Terms", fMassTerms);

  // Information about Structure Functions
  GetParam("LHAPDF-set",      fSFinfo.LHAPDFset    );
  GetParam("LHAPDF-member",   fSFinfo.LHAPDFmember );
  GetParam("Is-NLO",          fSFinfo.IsNLO        );
  GetParam("Scheme",          fSFinfo.Scheme       );
  GetParam("Quark-Threshold", fSFinfo.QrkThrs      );
  GetParam("NGridX",          fSFinfo.NGridX       );
  GetParam("NGridQ2",         fSFinfo.NGridQ2      );
  GetParam("XGrid-Min",       fSFinfo.XGridMin     );
  GetParam("Q2Grid-Min",      fSFinfo.Q2GridMin    );
  GetParam("Q2Grid-Max",      fSFinfo.Q2GridMax    );
  GetParam("MassW",           fSFinfo.MassW        );
  GetParam("MassZ",           fSFinfo.MassZ        );
  GetParam("Rho",             fSFinfo.Rho          );
  GetParam("Sin2ThW",         fSFinfo.Sin2ThW      );
  GetParam("Mud",             fSFinfo.Vud          );
  GetParam("Mus",             fSFinfo.Vus          );
  GetParam("Mub",             fSFinfo.Vub          );
  GetParam("Mcd",             fSFinfo.Vcd          );
  GetParam("Mcs",             fSFinfo.Vcs          );
  GetParam("Mcb",             fSFinfo.Vcb          );
  GetParam("Mtd",             fSFinfo.Vtd          );
  GetParam("Mts",             fSFinfo.Vts          );
  GetParam("Mtb",             fSFinfo.Vtb          );

}