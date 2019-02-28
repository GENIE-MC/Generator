//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Steven Gardiner <gardiner \at fnal.gov>
         Fermi National Accelerator Laboratory - Feb 26, 2019

 For the class documentation see the corresponding header file.

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
#include "Physics/NuclearState/PauliBlocker.h"
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

  // We're doing an MC integration over the nucleon momentum distribution
  // ourselves (including Pauli blocking), so we don't want to apply the
  // nuclear suppression factor. To turn it off, we'll set the "assume free
  // nucleon" flag.
  interaction->SetBit( kIAssumeFreeNucleon );

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

  ROOT::Math::IBaseFunctionMultiDim* func = new utils::gsl::FullQELdXSec(model,
    interaction, fPauliBlock, fPauliBlockID, bind_mode);
  ROOT::Math::IntegrationMultiDim::Type ig_type =
    utils::gsl::IntegrationNDimTypeFromString( fGSLIntgType );

  double abstol = 1e-16; // We mostly care about relative tolerance
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);

  // TODO: perhaps use genie::utils::CosTheta0Max() to speed things up here?
  Range1D_t cos_theta_0_lim( -1., 1. );
  Range1D_t phi_0_lim( 0., 2.*kPi );

  double kine_min[2] = { cos_theta_0_lim.min, phi_0_lim.min };
  double kine_max[2] = { cos_theta_0_lim.max, phi_0_lim.max };

  // For a free nucleon target (hit nucleon is at rest in the lab frame), we
  // don't need to do an MC integration over the initial state variables. In
  // this case, just set up the nucleon at the origin, on-shell, and at rest,
  // then integrate over the angles and return the result.
  if ( !tgt->IsNucleus() ) {
    tgt->SetHitNucPosition(0.);
    nucl_model->SetMomentum3( TVector3(0., 0., 0.) );
    nucl_model->SetRemovalEnergy(0.);
    double xsec_total = ig.Integral(kine_min, kine_max);
    return xsec_total;
  }

  // For a nuclear target, we need to loop over a bunch of nucleons sampled
  // from the nuclear model (with positions sampled from the vertex generator
  // to allow for using the local Fermi gas model). The MC estimator for the
  // total cross section is simply the mean of ig.Integral() for all of the
  // sampled nucleons.
  double xsec_sum = 0.;
  for (int n = 0; n < fNumNucleonThrows; ++n) {

    // Select a new position for the initial hit nucleon (needed for
    // the local Fermi gas model)
    TVector3 vertex_pos = vtx_gen->GenerateVertex( interaction, tgt->A() );
    double radius = vertex_pos.Mag();
    tgt->SetHitNucPosition( radius );

    // Sample a new nucleon 3-momentum and removal energy (this will be applied
    // to the nucleon via a call to genie::utils::ComputeFullQELPXSec(), so
    // there's no need to mess with its 4-momentum here)
    nucl_model->GenerateNucleon(*tgt, radius);

    // The initial state variables have all been defined, so integrate over
    // the final lepton angles. Note that we handle Pauli blocking and
    // multiplying by the number of active nucleons within the call
    // to ig.Integral()
    double xsec = ig.Integral(kine_min, kine_max);

    xsec_sum += xsec;
  }

  delete func;
  delete interaction;

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

  GetParamDef( "DoPauliBlocking", fPauliBlock, true );

  RgAlg pauliBlockID;
  GetParamDef( "PauliBlockerAlg", pauliBlockID, RgAlg("genie::PauliBlocker", "Default") );
  fPauliBlockID = AlgId( pauliBlockID );

  GetParamDef( "NumNucleonThrows", fNumNucleonThrows, 5000 );
}

genie::utils::gsl::FullQELdXSec::FullQELdXSec(const XSecAlgorithmI* xsec_model,
  const Interaction* interaction, bool do_Pauli_blocking,
  const AlgId& pauli_blocker_ID, QELEvGen_BindingMode_t binding_mode)
  : fXSecModel( xsec_model ), fInteraction( new Interaction(*interaction) ),
  fDoPauliBlocking( do_Pauli_blocking ), fHitNucleonBindingMode( binding_mode )
{
  fNuclModel = dynamic_cast<const NuclearModelI*>( fXSecModel->SubAlg("IntegralNuclearModel") );
  assert( fNuclModel );

  AlgFactory* algf = AlgFactory::Instance();
  fPauliBlocker = dynamic_cast<const PauliBlocker*>( algf->GetAlgorithm(pauli_blocker_ID) );
  assert( fPauliBlocker );
}

genie::utils::gsl::FullQELdXSec::~FullQELdXSec()
{
  delete fInteraction;
}

ROOT::Math::IBaseFunctionMultiDim* genie::utils::gsl::FullQELdXSec::Clone(void) const
{
  return new FullQELdXSec(fXSecModel, fInteraction, fDoPauliBlocking,
    fPauliBlocker->Id(), fHitNucleonBindingMode);
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

  // TODO: min_angle_EM!!!!
  double min_angle_EM = 0.;

  // Compute the full differential cross section
  double xsec = genie::utils::ComputeFullQELPXSec(fInteraction, fNuclModel,
    fXSecModel, cos_theta0, phi0, dummy_Eb, fHitNucleonBindingMode, min_angle_EM, true);

  // ComputeFullQELPXSec() sets the final nucleon's lab frame 4-momentum via
  // interaction->KinePtr()->SetHadSystP4(), so we can now check for Pauli
  // blocking

  // If the final nucleon would be Pauli blocked, then return zero immediately
  const Target& tgt = fInteraction->InitState().Tgt();
  if ( fDoPauliBlocking && tgt.IsNucleus() ) {
    double kF = fPauliBlocker->GetFermiMomentum(tgt, fInteraction->RecoilNucleonPdg(),
      tgt.HitNucPosition());
    double pNf = fInteraction->Kine().HadSystP4().P();
    if ( pNf < kF ) return 0.;
  }

  // PDG code for the initial hit nucleon
  int pdg_Ni = tgt.HitNucPdg();
  // Number of active nucleons in the target
  int num_active_nucleons = genie::pdg::IsProton( pdg_Ni ) ? tgt.Z() : tgt.N();

  xsec *= num_active_nucleons;

  return xsec;
}
