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
  //----- Make sure it has the configuration parameters for computing an
  //      average multiplicity

  assert(
         fConfig->Exists("beta")      &&
         fConfig->Exists("alpha-vp")  && fConfig->Exists("alpha-vn")  &&
         fConfig->Exists("alpha-vbp") && fConfig->Exists("alpha-vbn") 
        );
  
  //----- Make sure it knows which parameter set to use for the KNO
  //      distribution and get its name
  
  assert( fConfig->Exists("kno-param-set") );

  string kno_param_set = fConfig->GetString("kno-param-set");

  //----- Get the KNO Distribution

  //-- Get an instance of the AlgFactory

  AlgFactory * algf = AlgFactory::Instance();

  //-- Request the specified KNO distribution algorithm

  const Algorithm * algbase =
                 algf->GetAlgorithm("genie::KNODistribution", kno_param_set);

  const KNODistribution * kno =
                             dynamic_cast<const KNODistribution *> (algbase);

  //----- Compute the average multiplicity for the given interaction:
  //      <n> = a + b * ln(W^2)
  
  double beta  = fConfig->GetDouble("beta");
  double alpha = SelectOffset(interaction);

  double W = interaction->GetKinematics().W();

  assert( W > 0 );
  
  double avn   = alpha + beta * 2 * log(W);

  if(avn < 1) {
      LOG("Schmitz", pWARN) << "Average multiplicity too small: " << avn;
      return 0;
  }

  //----- Create a multiplicity probability distribution

  TH1D * prob = new TH1D("", "", kMaxMultiplicity, 0, kMaxMultiplicity);

  for(int n = 0; n < kMaxMultiplicity; n++) {

     // KNO distribution is <n>*P(n) vs n/<n>

     double n_avn = n / avn;              // n/<n>
     double avnP  = kno->Value(n_avn);    // <n>*P(n)
     double P     = avnP / avn;           // P(n)

     LOG("Schmitz", pDEBUG)
             << "W = " << W << ", <n> = " << avn << ", n/<n> = " << n_avn
             << ", <n>*P = " << avnP << ", P = " << P;

     prob->Fill(n, P);
  }

  //----- Normalize the probability distribution

  prob->Scale( 1.0 / prob->Integral("width") );

  return prob; // Note: The calling function adopts the object
}
//____________________________________________________________________________
double SchmitzMultiplicityModel::SelectOffset(
                                        const Interaction * interaction) const
{
  //----- Make sure we know which nucleon type was hit

  const InitialState & init_state = interaction->GetInitialState();
  
  int nu_pdg  = init_state.GetProbePDGCode();
  int nuc_pdg = init_state.GetTarget().StruckNucleonPDGCode();

  if( pdg::IsNeutrino( nu_pdg ) ) {

      if ( pdg::IsProton(nuc_pdg)  )    return fConfig->GetDouble("alpha-vp");
      if ( pdg::IsNeutron(nuc_pdg) )    return fConfig->GetDouble("alpha-vn");
      else {
         LOG("Schmitz", pERROR)
                          << "PDG-Code = " << nuc_pdg << " is not a nucleon!";
      }

  } else  if (  pdg::IsAntiNeutrino(nu_pdg) ) {

      if ( pdg::IsProton(nuc_pdg)  )   return fConfig->GetDouble("alpha-vbp");
      if ( pdg::IsNeutron(nuc_pdg) )   return fConfig->GetDouble("alpha-vbn");
      else {
         LOG("Schmitz", pERROR)
                          << "PDG-Code = " << nuc_pdg << " is not a nucleon!";
      }
  } else {
    LOG("Schmitz", pERROR) << "PDG-Code = " << nu_pdg << " is not a neutrino!";
  }

  return 0;        
}
//____________________________________________________________________________

