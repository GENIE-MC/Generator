//_________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 For the class documentation see the corresponding header file.
*/
//_________________________________________________________________________

#include "Framework/EventGen/HybridXSecAlgorithm.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"

using namespace genie;

//_________________________________________________________________________
HybridXSecAlgorithm::HybridXSecAlgorithm() : XSecAlgorithmI("genie::HybridXSecAlgorithm")
{
}
//_________________________________________________________________________
HybridXSecAlgorithm::HybridXSecAlgorithm(string config)
  : XSecAlgorithmI("genie::HybridXSecAlgorithm", config)
{
}
//_________________________________________________________________________
HybridXSecAlgorithm::~HybridXSecAlgorithm()
{
}
//_________________________________________________________________________
const XSecAlgorithmI* HybridXSecAlgorithm::ChooseXSecAlg(
  const Interaction& interaction) const
{
  std::string inter_str = interaction.AsString();
  RgKey key = "XSecAlg@Interaction=" + inter_str;

  std::map<std::string, const XSecAlgorithmI*>::const_iterator
    cend = fXSecAlgMap.cend();

  std::map<std::string, const XSecAlgorithmI*>::const_iterator
    citer = fXSecAlgMap.find( key );

  // If the algorithm doesn't appear in the map, try to load it.
  if ( citer == cend ) {

    // If a key exists for the algorithm in the registry, load it
    // and store it for rapid retrieval later
    const Registry& temp_reg = this->GetConfig();
    if ( temp_reg.Exists(key) ) {
      const XSecAlgorithmI* temp_alg = dynamic_cast< const XSecAlgorithmI* >(
        this->SubAlg(key) );
      assert( temp_alg );

      fXSecAlgMap[ inter_str ] = temp_alg;
      return temp_alg;
    }

    // Otherwise, if the user has specified a default algorithm, then store a
    // new entry for that one in the map
    else if ( fDefaultXSecAlg ) {
      fXSecAlgMap[ inter_str ] = fDefaultXSecAlg;
      return fDefaultXSecAlg;
    }

    // Otherwise, store and return a null pointer. No suitable algorithm could
    // be found for the requested interaction.
    else {
      fXSecAlgMap[ inter_str ] = NULL;
      return NULL;
    }

  }

  // If an entry was found in the map, then just use that
  else return citer->second;
}
//_________________________________________________________________________
double HybridXSecAlgorithm::XSec(const Interaction* interaction,
  KinePhaseSpace_t kps) const
{
  const XSecAlgorithmI* alg_to_use = this->ChooseXSecAlg( *interaction );

  if ( !alg_to_use ) return 0.;
  // Special case, when the process is EM scattering on free nucleon.
  // In this case genie::RosenbluthPXSec is used to calculate cross-section.
  // However, input kps is equal to either kPSTlctl for SuSAv2 
  // or kPSQELEvGen for LwlynSmithQELCCPXSec (2d phase spaces), 
  // while RosenbluthPXSec uses kPSQ2fE (1d phase space)
  // and there is no way to transform 1d to 2d phase space.
  if ( interaction->ProcInfo().IsEM() && 
      (interaction->InitState().Tgt().Pdg() == kPdgTgtFreeP || 
       interaction->InitState().Tgt().Pdg() == kPdgTgtFreeN) 
     )
        return alg_to_use->XSec( interaction, kPSQ2fE );
  return alg_to_use->XSec( interaction, kps );
  
}
//_________________________________________________________________________
double HybridXSecAlgorithm::Integral(const Interaction* interaction) const
{
  const XSecAlgorithmI* alg_to_use = this->ChooseXSecAlg( *interaction );

  if ( !alg_to_use ) return 0.;
  else return alg_to_use->Integral( interaction );
}
//_________________________________________________________________________
bool HybridXSecAlgorithm::ValidProcess(const Interaction* interaction) const
{
  if ( interaction->TestBit(kISkipProcessChk) ) return true;

  const XSecAlgorithmI* alg_to_use = this->ChooseXSecAlg( *interaction );
  if ( !alg_to_use ) return false;
  else return alg_to_use->ValidProcess( interaction );
}
//_________________________________________________________________________
void HybridXSecAlgorithm::Configure(const Registry& config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HybridXSecAlgorithm::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//_________________________________________________________________________
void HybridXSecAlgorithm::LoadConfig(void)
{
  fDefaultXSecAlg = NULL;

  // The user can optionally configure a cross section algorithm
  // to use by default (i.e., whenever one wasn't explicitly specified
  // for an input interaction). If one was given, then configure it.
  // Handling of the interaction-specific algorithms is done via
  // lazy initialization in ChooseXSecAlg().
  const Registry& temp_reg = this->GetConfig();
  if ( temp_reg.Exists("DefaultXSecAlg") ) {
    fDefaultXSecAlg = dynamic_cast< const XSecAlgorithmI* >(
      this->SubAlg("DefaultXSecAlg") );
    assert( fDefaultXSecAlg );
  }
}
//____________________________________________________________________________
const TVector3 & HybridXSecAlgorithm::FinalLeptonPolarization (const Interaction* interaction) const
{
  const XSecAlgorithmI* alg_to_use = this->ChooseXSecAlg( *interaction );

  if ( !alg_to_use ) return XSecAlgorithmI::FinalLeptonPolarization(interaction);
  
  return alg_to_use->FinalLeptonPolarization( interaction );
}
//____________________________________________________________________________
