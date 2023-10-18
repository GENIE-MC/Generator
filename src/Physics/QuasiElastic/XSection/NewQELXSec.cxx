//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory 
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/QuasiElastic/XSection/NewQELXSec.h"

#include "Physics/NuclearState/NuclearModelI.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"
#include "Physics/Common/VertexGenerator.h"
#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelMap.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::gsl;

//____________________________________________________________________________
NewQELXSec::NewQELXSec() : XSecIntegratorI("genie::NewQELXSec")
{

}
//____________________________________________________________________________
NewQELXSec::NewQELXSec(std::string config) : XSecIntegratorI("genie::NewQELXSec", config)
{

}
//____________________________________________________________________________
double NewQELXSec::Integrate(const XSecAlgorithmI* model, const Interaction* in) const
{
  LOG("NewQELXSec",pDEBUG) << "Beginning integrate";
  if ( !model->ValidProcess(in) ) return 0.;

  Interaction* interaction = new Interaction( *in );
  interaction->SetBit( kISkipProcessChk );
  //interaction->SetBit( kISkipKinematicChk );

  const NuclearModelI* nucl_model = dynamic_cast<const NuclearModelI*>(
    model->SubAlg("IntegralNuclearModel") );
  assert( nucl_model );

  AlgFactory* algf = AlgFactory::Instance();
  const VertexGenerator* vtx_gen = dynamic_cast<const VertexGenerator*>(
    algf->GetAlgorithm(fVertexGenID) );
  assert( vtx_gen );

  // Determine the appropriate binding energy mode to use.
  // The default given here is for the case of a free nucleon.
  QELEvGen_BindingMode_t bind_mode = kOnShell;
  Target* tgt = interaction->InitState().TgtPtr();
  if ( tgt->IsNucleus() ) {
    std::string bind_mode_str = model->GetConfig()
      .GetString( "IntegralNucleonBindingMode" );
    bind_mode = genie::utils::StringToQELBindingMode( bind_mode_str );
  }

  utils::gsl::FullQELdXSec* func = new utils::gsl::FullQELdXSec(model,
    interaction, bind_mode, fMinAngleEM);
  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString( fGSLIntgType );

  // Switch to using the copy of the interaction in the integrator rather than
  // the copy that we made in this function
  delete interaction;
  interaction = func->GetInteractionPtr();

  // Also update the pointer to the Target
  tgt = interaction->InitState().TgtPtr();

  double abstol = 1e-16; // We mostly care about relative tolerance
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);

  // Integration ranges for the lepton COM frame scattering angles (in the
  // kPSQELEvGen phase space, these are measured with respect to the COM
  // velocity as observed in the lab frame)
  Range1D_t cos_theta_0_lim( -1., 1. );
  Range1D_t phi_0_lim( 0., 2.*kPi );

  double kine_min[2] = { cos_theta_0_lim.min, phi_0_lim.min };
  double kine_max[2] = { cos_theta_0_lim.max, phi_0_lim.max };

  // If averaging over the initial nucleon distribution has been
  // disabled, just integrate over angles and return the result.
  if ( !fAverageOverNucleons ) {
    double xsec_total = ig.Integral(kine_min, kine_max);
    delete func;
    return xsec_total;
  }

  // For a free nucleon target (hit nucleon is at rest in the lab frame), we
  // don't need to do an MC integration over the initial state variables. In
  // this case, just set up the nucleon at the origin, on-shell, and at rest,
  // then integrate over the angles and return the result.

  // Also use this approach if we're over the "nuclear influence" cutoff
  // energy for the probe. Beyond the cutoff, the effects of Fermi motion
  // and the removal energy are assumed to be small enough to be neglected
  double E_lab_cutoff = model->GetConfig()
    .GetDouble("IntegralNuclearInfluenceCutoffEnergy");

  double probeE = interaction->InitState().ProbeE( kRfLab );
  if ( !tgt->IsNucleus() || probeE > E_lab_cutoff ) {
    tgt->SetHitNucPosition(0.);

    if ( tgt->IsNucleus() ) nucl_model->GenerateNucleon(*tgt, 0.);
    else {
      nucl_model->SetRemovalEnergy(0.);
      interaction->SetBit( kIAssumeFreeNucleon );
    }

    nucl_model->SetMomentum3( TVector3(0., 0., 0.) );
    double xsec_total = ig.Integral(kine_min, kine_max);
    delete func;
    return xsec_total;
  }

  // For a nuclear target, we need to loop over a bunch of nucleons sampled
  // from the nuclear model (with positions sampled from the vertex generator
  // to allow for using the local Fermi gas model). The MC estimator for the
  // total cross section is simply the mean of ig.Integral() for all of the
  // sampled nucleons.
  double xsec_sum = 0.;
  for (int n = 0; n < fNumNucleonThrows; ++n) {

    // Select a new position for the initial hit nucleon (needed for the local
    // Fermi gas model, but other than slowing things down a bit, it doesn't
    // hurt to do this for other models)
    TVector3 vertex_pos = vtx_gen->GenerateVertex( interaction, tgt->A() );
    double radius = vertex_pos.Mag();
    tgt->SetHitNucPosition( radius );

    // Sample a new nucleon 3-momentum and removal energy (this will be applied
    // to the nucleon via a call to genie::utils::ComputeFullQELPXSec(), so
    // there's no need to mess with its 4-momentum here)
    nucl_model->GenerateNucleon(*tgt, radius);

    // The initial state variables have all been defined, so integrate over
    // the final lepton angles.
    double xsec = ig.Integral(kine_min, kine_max);

    xsec_sum += xsec;
  }

  delete func;

  // MC estimator of the total cross section is the mean of the xsec values
  double xsec_mean = xsec_sum / fNumNucleonThrows;

  return xsec_mean;
}
//____________________________________________________________________________
void NewQELXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NewQELXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NewQELXSec::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, std::string("adaptive") ) ;
  GetParamDef( "gsl-relative-tolerance", fGSLRelTol, 1e-2 ) ;
  int max;
  GetParamDef( "gsl-max-eval", max, 500000 ) ;
  fGSLMaxEval  = static_cast<unsigned int>( max );

  RgAlg vertexGenID;
  GetParamDef( "VertexGenAlg", vertexGenID, RgAlg("genie::VertexGenerator", "Default") );
  fVertexGenID = AlgId( vertexGenID );

  GetParamDef( "NumNucleonThrows", fNumNucleonThrows, 5000 );

  GetParamDef( "SF_MinAngleEM", fMinAngleEM, 0.);

  // If true, then the integration of the total cross section will include an
  // MC integration over the initial state nuclear model
  GetParamDef( "AverageOverNucleons", fAverageOverNucleons, true );
}

genie::utils::gsl::FullQELdXSec::FullQELdXSec(const XSecAlgorithmI* xsec_model,
  const Interaction* interaction, QELEvGen_BindingMode_t binding_mode, double min_angle_EM)
  : fXSecModel( xsec_model ), fInteraction( new Interaction(*interaction) ),
  fHitNucleonBindingMode( binding_mode ), fMinAngleEM( min_angle_EM )
{
  fNuclModel = dynamic_cast<const NuclearModelI*>( fXSecModel->SubAlg("IntegralNuclearModel") );
  assert( fNuclModel );
}

genie::utils::gsl::FullQELdXSec::~FullQELdXSec()
{
  delete fInteraction;
}

Interaction* genie::utils::gsl::FullQELdXSec::GetInteractionPtr()
{
  return fInteraction;
}

const Interaction& genie::utils::gsl::FullQELdXSec::GetInteraction() const
{
  return *fInteraction;
}

ROOT::Math::IBaseFunctionMultiDim* genie::utils::gsl::FullQELdXSec::Clone(void) const
{
  return new FullQELdXSec(fXSecModel, fInteraction, fHitNucleonBindingMode, fMinAngleEM);
}

unsigned int genie::utils::gsl::FullQELdXSec::NDim(void) const
{
  return 2;
}

double genie::utils::gsl::FullQELdXSec::DoEval(const double* xin) const
{
  // Elements of "xin"
  //
  // element 0: "cos_theta0" = Cosine of theta0, the angle between the COM frame
  //                           3-momentum of the outgoing lepton and the COM frame velocity
  //                           as measured in the laboratory frame
  // element 1: "phi_theta0" = Azimuthal angle of the COM frame 3-momentum of the
  //                           outgoing lepton measured with respect to the COM frame
  //                           velocity as measured in the laboratory frame

  double cos_theta0 = xin[0];
  double phi0 = xin[1];

  // Dummy storage for the binding energy of the hit nucleon
  double dummy_Eb = 0.;

  // Compute the full differential cross section
  double xsec = genie::utils::ComputeFullQELPXSec(fInteraction, fNuclModel,
    fXSecModel, cos_theta0, phi0, dummy_Eb, fHitNucleonBindingMode, fMinAngleEM, true);

  return xsec;
}
