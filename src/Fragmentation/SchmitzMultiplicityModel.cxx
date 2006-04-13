//____________________________________________________________________________
/*!

\class    genie::SchmitzMultiplicityModel

\brief    The 'Schmitz' multiplicity probability model as used in NeuGEN.
          Is a concerete implementation of the MultiplicityProbModelI interface.
          
\ref      N. Schmitz, Proc. Intl. Symp. on Lepton & Photon Interactions at
          High Energies, Bonn 1981 p.527

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Fragmentation/SchmitzMultiplicityModel.h"
#include "Fragmentation/KNODistribution.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
SchmitzMultiplicityModel::SchmitzMultiplicityModel() :
MultiplicityProbModelI("genie::SchmitzMultiplicityModel")
{                

}
//____________________________________________________________________________
SchmitzMultiplicityModel::SchmitzMultiplicityModel(string config) :
MultiplicityProbModelI("genie::SchmitzMultiplicityModel", config)
{                 

}
//____________________________________________________________________________
SchmitzMultiplicityModel::~SchmitzMultiplicityModel()
{

}
//____________________________________________________________________________
TH1D * SchmitzMultiplicityModel::ProbabilityDistribution(
                                        const Interaction * interaction) const
{
  // Compute the average multiplicity for the given interaction:
  // <n> = a + b * ln(W^2)
  
  double alpha = this->SelectOffset(interaction);
  double W     = interaction->GetKinematics().W();

  assert(W>kNeutronMass+kPionMass);
  
  double avn = alpha + fB * 2*TMath::Log(W);
  if(avn < 1) {
      LOG("Schmitz", pWARN) << "Average multiplicity too small: " << avn;
      return 0;
  }

  // Create a multiplicity probability distribution

  TH1D * prob = new TH1D("", "", kMaxMultiplicity, 0, kMaxMultiplicity);

  for(int n = 0; n < kMaxMultiplicity; n++) {
     // KNO distribution is <n>*P(n) vs n/<n>
     double n_avn = (double)n / avn;      // n/<n>
     double avnP  = fKNO->Value(n_avn);   // <n>*P(n)
     double P     = avnP / avn;           // P(n)

     LOG("Schmitz", pDEBUG)
          << "W = " << W << ", <n> = " << avn << ", n/<n> = " << n_avn
          << ", <n>*P = " << avnP << ", P = " << P;

     prob->Fill( (double)n, P);
  }
  //----- Normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  return prob; // Note: The calling function adopts the object
}
//____________________________________________________________________________
double SchmitzMultiplicityModel::SelectOffset(
                                        const Interaction * interaction) const
{
  const InitialState & init_state = interaction->GetInitialState();
  int nu_pdg  = init_state.GetProbePDGCode();
  int nuc_pdg = init_state.GetTarget().StruckNucleonPDGCode();

  if( pdg::IsNeutrino( nu_pdg ) ) {
      if ( pdg::IsProton(nuc_pdg)  )  return fAvp;
      if ( pdg::IsNeutron(nuc_pdg) )  return fAvn;
      else {
         LOG("Schmitz", pERROR)
                          << "PDG-Code = " << nuc_pdg << " is not a nucleon!";
      }
  } else  if (  pdg::IsAntiNeutrino(nu_pdg) ) {
      if ( pdg::IsProton(nuc_pdg)  )  return fAvbp;
      if ( pdg::IsNeutron(nuc_pdg) )  return fAvbn;
      else {
         LOG("Schmitz", pERROR)
                          << "PDG-Code = " << nuc_pdg << " is not a nucleon!";
      }
  } else {
    LOG("Schmitz", pERROR)<< "PDG-Code = " << nu_pdg << " is not a neutrino!";
  }
  return 0;        
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SchmitzMultiplicityModel::LoadConfig(void)
{
  // Load config parameters
  fAvp  = fConfig->GetDouble("alpha-vp");
  fAvn  = fConfig->GetDouble("alpha-vn");
  fAvbp = fConfig->GetDouble("alpha-vbp");
  fAvbn = fConfig->GetDouble("alpha-vbn");
  fB    = fConfig->GetDouble("beta");
  fKNOParamSet = fConfig->GetString("kno-param-set");

  // Get the KNO Distribution
  AlgFactory * algf = AlgFactory::Instance();
  fKNO = dynamic_cast<const KNODistribution *> (
                 algf->GetAlgorithm("genie::KNODistribution", fKNOParamSet));
}
//____________________________________________________________________________

