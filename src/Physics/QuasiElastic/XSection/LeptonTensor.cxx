#include "Framework/Conventions/RefFrame.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/QuasiElastic/XSection/LeptonTensor.h"

genie::LeptonTensor::LeptonTensor(const genie::Interaction& interaction, bool SF)
{
  // Make copies of the initial and final lepton 4-momenta in the lab frame
  TLorentzVector* temp_probeP4 = interaction.InitState()
    .GetProbeP4( genie::kRfLab );
  fProbeP4 = *temp_probeP4;
  delete temp_probeP4;

  fFSLepP4 = interaction.Kine().FSLeptonP4();

  fInitialLeptonPDG = interaction.InitState().ProbePdg();
  fInteractionType = interaction.ProcInfo().InteractionTypeId();

  // Initial lepton mass (needed for EM channel)
  fMLep2 = std::pow(interaction.InitState().Probe()->Mass(), 2);
 
  // Need to use a different tensor convention for SF cross sections
  SFtensor = SF;
}

genie::LeptonTensor::LeptonTensor(const TLorentzVector& p4Probe,
  const TLorentzVector& p4Lep, int probe_pdg, InteractionType_t type, bool SF)
{
  fProbeP4 = p4Probe;
  fFSLepP4 = p4Lep;

  fInitialLeptonPDG = probe_pdg;
  fInteractionType = type;

  genie::PDGLibrary* pdg_lib = genie::PDGLibrary::Instance();
  double probe_mass = pdg_lib->Find( probe_pdg )->Mass();
  fMLep2 = std::pow( probe_mass, 2 );

  SFtensor = SF;
}

std::complex<double> genie::LeptonTensor::operator()(genie::TensorIndex_t mu,
  genie::TensorIndex_t nu) const
{
  double g_mu_nu = Metric(mu, nu);
  double common_terms = fProbeP4[mu]*fFSLepP4[nu] + fProbeP4[nu]*fFSLepP4[mu]
    - g_mu_nu*(fProbeP4 * fFSLepP4);

  std::complex<double> result = 0.;
  if ( fInteractionType == genie::kIntEM ) {
    result = 2. * (fMLep2*g_mu_nu + common_terms);
  }
  else if ( fInteractionType == genie::kIntWeakCC
    || fInteractionType == genie::kIntWeakNC )
  {
    // NS: Tricky but Noemi's SF conventions
    // Require sign(nu) = +1, etc. and a different LeviCivita
    // I'll choose Steven's sign as the default 
    double sign = -1.;
    if(SFtensor) sign =  +1.;

    if ( genie::pdg::IsAntiNeutrino(fInitialLeptonPDG) ) sign = (!SFtensor) ? 1. : -1.;
    else if ( ! genie::pdg::IsNeutrino(fInitialLeptonPDG) ) {
      LOG("LeptonTensor", pFATAL) << "Invalid probe PDG code "
        << fInitialLeptonPDG << " encountered for a weak interaction";
      std::exit(1);
    }
    if (SFtensor) result = 8. * ( common_terms + sign*LeviCivitaProductSF(mu, nu) );
    else result = 8. * ( common_terms + sign*LeviCivitaProduct(mu, nu) );
    
  }
  else {
    LOG("LeptonTensor", pFATAL) << "Invalid interaction type "
      << fInteractionType << " encountered by genie::LeptonTensor";
    std::exit(1);
  }

  return result;
}


std::complex<double> genie::LeptonTensor::LeviCivitaProduct(
  genie::TensorIndex_t mu, genie::TensorIndex_t nu) const
{
  return std::complex<double>( 0.,
    Rank2LorentzTensorI::LeviCivitaProduct(mu, nu, fProbeP4, fFSLepP4) );
}

std::complex<double> genie::LeptonTensor::LeviCivitaProductSF(
  genie::TensorIndex_t mu, genie::TensorIndex_t nu) const
{
  return std::complex<double> (0., 
    Rank2LorentzTensorI::LeviCivitaProductSF(mu, nu, fProbeP4, fFSLepP4) );
}
 

