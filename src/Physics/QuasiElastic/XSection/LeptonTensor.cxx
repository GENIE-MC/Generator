#include "Framework/Conventions/RefFrame.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/QuasiElastic/XSection/FreeNucleonTensor.h"
#include "Physics/QuasiElastic/XSection/LeptonTensor.h"

genie::LeptonTensor::LeptonTensor(const genie::Interaction& interaction)
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
    double sign = 1.;
    if ( genie::pdg::IsAntiNeutrino(fInitialLeptonPDG) ) sign = -1.;
    else if ( ! genie::pdg::IsNeutrino(fInitialLeptonPDG) ) {
      LOG("LeptonTensor", pFATAL) << "Invalid probe PDG code "
        << fInitialLeptonPDG << " encountered for a weak interaction";
      std::exit(1);
    }
    result = 8. * ( common_terms + sign*LeviCivitaProduct(mu, nu) );
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

std::complex<double> genie::LeptonTensor::Contract(
  const Rank2LorentzTensorI& other) const
{
  // If we're dealing with the contraction of the EM nucleon current with the
  // leptonic current, then we need to add a couple of terms to correct for
  // current conservation. Code to do that lives in
  // FreeNucleonTensor::Contract(), so call that method if needed.
  const genie::FreeNucleonTensor* fnt = dynamic_cast<const
    genie::FreeNucleonTensor*>( &other );
  if ( fnt && fInteractionType == kIntEM ) return fnt->Contract( *this );
  else return Rank2LorentzTensorI::Contract( other );
}
