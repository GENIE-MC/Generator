//____________________________________________________________________________
/*!

\class    genie::KovalenkoQELCharmPXSec

\brief    Computes the QEL Charm Production Differential Cross Section
          using \b Kovalenko's duality model approach. 

          The computed differential cross section is the Dxsec = dxsec/dQ^2
          where \n
            \li \c Q2 is the momentum transfer.

          It models the differential cross sections for: \n
             \li v + n \rightarrow mu- + Lambda_{c}^{+} (2285)
             \li v + n \rightarrow mu- + Sigma_{c}^{+}  (2455)
             \li v + p \rightarrow mu- + Sigma_{c}^{++} (2455)

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      S.G.Kovalenko, Sov.J.Nucl.Phys.52:934 (1990)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#include <iostream>

#include "AlgFactory/AlgFactory.h"
#include "Charm/KovalenkoQELCharmPXSec.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "PDF/PDF.h"
#include "PDF/PDFModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
KovalenkoQELCharmPXSec::KovalenkoQELCharmPXSec() :
XSecAlgorithmI()
{
  fName     = "genie::KovalenkoQELCharmPXSec";
}
//____________________________________________________________________________
KovalenkoQELCharmPXSec::KovalenkoQELCharmPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::KovalenkoQELCharmPXSec";

  FindConfig();
}
//____________________________________________________________________________
KovalenkoQELCharmPXSec::~KovalenkoQELCharmPXSec()
{

}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::XSec(const Interaction * interaction) const
{
  LOG("CharmXSec", pDEBUG) << *fConfig;

  this->AssertProcessValidity(interaction);

  //----- get scattering params & init state - compute auxiliary vars

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();
  
  //final state primary lepton & nucleon mass
  double ml    = interaction -> GetFSPrimaryLepton() -> Mass();
  double Mnuc  = kNucleonMass; 
  double Mnuc2 = Mnuc * Mnuc;
  
  //neutrino energy & momentum transfer  

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);
  
  double E     = p4->Energy();
  double E2    = E * E;
  double Q2    = sc_params.Q2();

  delete p4;
  
  //resonance mass
  double MR    = this -> MRes  (interaction);
  double MR2   = this -> MRes2 (interaction);
  
  //resonance threshold
  double ER = ( pow(MR+ml,2) - Mnuc2 ) / (2*Mnuc);

  //----- check for user cuts on Q2;

  double Q2min, Q2max;

  this->Q2Cuts(Q2min, Q2max);

  if(Q2 > Q2max || Q2 < Q2min) return 0;

  //----- Calculate the differential cross section dxsec/dQ^2

  double xsec = 0;
  
  if(E > ER) {

    double Gf  = kGF_2 / (2*kPi);

    double vR = (MR2 - Mnuc2 + Q2) / (2*Mnuc);

    double xiR       = this->xiBar(interaction, vR);
    double vR2       = vR*vR;
    double vR_E      = vR/E;
    double Q2_4E2    = Q2/(4*E2);
    double Q2_2MExiR = Q2/(2*Mnuc*E*xiR);

    double Z = this->ZR(interaction);
    double D = this->DR(interaction);
    
    LOG("CharmXSec", pDEBUG) << "Z = " << Z << ", D = " << D;

    xsec = Gf * Z * D * ( 1 - vR_E + Q2_4E2 + Q2_2MExiR ) * sqrt(vR2 + Q2) / (vR*xiR);    
  }

  return xsec;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::ZR(const Interaction * interaction) const
{
  double Mo    = this->MoScale();
  double Mo2   = Mo*Mo;  
  double Mnuc2 = kNucleonMass_2;
  double MR2   = this->MRes2(interaction);

  double D0    = this->DR(interaction, true); // D^R(Q^2=0)
  double sumF2 = this->SumF2(interaction);    // FA^2+F1^2
  
  double Z  = 2*Mo2*kSin8c_2 * sumF2 / (D0 * (MR2-Mnuc2));

  return Z;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::DR(
                             const Interaction * interaction, bool norm) const
{
  //----- get the requested PDF model & attach it to a PDF object

  const Algorithm * algbase = this->SubAlg("pdf-alg-name", "pdf-param-set");
  
  const PDFModelI * pdf_model = dynamic_cast<const PDFModelI *> (algbase);

  PDF pdfs;
  pdfs.SetModel(pdf_model);   // <-- attach algorithm

  //----- compute integration area = [xi_bar_plus, xi_bar_minus]

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();

  double Q2     = sc_params.Q2();
  double Mnuc   = kNucleonMass;
  double Mnuc2  = kNucleonMass_2;
  double MR     = this->MRes(interaction);
  double DeltaR = this->ResDM(interaction);

  double vR_minus  = ( pow(MR-DeltaR,2) - Mnuc2 + Q2 ) / (2*Mnuc);
  double vR_plus   = ( pow(MR+DeltaR,2) - Mnuc2 + Q2 ) / (2*Mnuc);

  LOG("CharmXSec", pDEBUG)
            << "vR = [plus: " << vR_plus << ", minus: " << vR_minus << "]";
  
  double xi_bar_minus = this->xiBar(interaction, vR_minus);
  double xi_bar_plus  = this->xiBar(interaction, vR_plus);

  LOG("CharmXSec", pDEBUG) << "Integration limits = ["
                             << xi_bar_plus << ", " << xi_bar_minus << "]";  

  //----- define the integration grid & instantiate a FunctionMap

  UnifGrid grid;

  int nbins = (fConfig->Exists("nbins")) ? fConfig->GetInt("nbins") : 201;
  
  grid.AddDimension(nbins, xi_bar_plus, xi_bar_minus);

  FunctionMap fmap(grid);


  //----- auxiliary variables

  const InitialState & init_state = interaction -> GetInitialState();

  double delta_xi_bar = (xi_bar_plus - xi_bar_minus) / (nbins - 1);

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  assert(isP || isN);

  //----- loop over x range (at fixed Q^2) & compute the function map

  for(int i = 0; i < nbins; i++) {

     double t = xi_bar_plus + i * delta_xi_bar;
     
     if( t<0 || t>1) fmap.AddPoint( 0., i );
     else {

       if(norm) pdfs.Calculate(t, 0.);
       else     pdfs.Calculate(t, Q2);

       double f = (isP) ? ( pdfs.DownValence() + pdfs.DownSea() ):
                          ( pdfs.UpValence()   + pdfs.UpSea()   );

       fmap.AddPoint( f, i );

       LOG("CharmXSec", pDEBUG)
          << "point...." << i+1 << "/" << nbins << " : "
          << "x*pdf(Q^2 = " << Q2 << ", x = " << t << ") = " << f;
    }
  }

  //----- Numerical integration

  const IntegratorI * integrator = this->Integrator();

  double D = integrator->Integrate(fmap);

  return D;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::xiBar(
                              const Interaction * interaction, double v) const
{
  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();

  double Q2     = sc_params.Q2();
  double Mnuc   = kNucleonMass;
  double Mo     = this->MoScale();
  double Mo2    = Mo*Mo;
  double v2     = v *v;

  LOG("CharmXSec", pDEBUG)
                     << "Q2 = " << Q2 << ", Mo = " << Mo << ", v = " << v;

  double xi  = ( Q2/Mnuc ) / ( v + sqrt(v2+Q2) );

  double xi_bar = xi * ( 1 + (1 + Mo2/(Q2+Mo2))*Mo2/Q2 );

  return xi_bar;
}
//____________________________________________________________________________
void KovalenkoQELCharmPXSec::Q2Cuts(double & Q2min, double & Q2max) const
{
  // init:
  Q2min = -999999;
  Q2max =  999999;

  // read from config (if they exist):
  if( fConfig->Exists("Q2min") ) fConfig->Get("Q2min", Q2min);
  if( fConfig->Exists("Q2max") ) fConfig->Get("Q2max", Q2max);

  assert(Q2min < Q2max);
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::MoScale(void) const
{
// Return the 'proper scale of internal nucleon dynamics'.
// In the original paper Mo = 0.08 +/- 0.02 GeV.
// Try to get this value from the algorithm config registry, and if it is not
// exists set it to the value from Eq.(20) in Sov.J.Nucl.Phys.52:934 (1990)

  double Mo = (fConfig->Exists("Mo")) ? fConfig->GetDouble("Mo") : 0.1;

  return Mo;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::ResDM(const Interaction * interaction) const
{
// Resonance Delta M obeys the constraint DM(R+/-) <= |M(R+/-) - M(R)|
// where M(R-) <= M(R) <= M(R+) are the masses of the neighboring
// resonances R+, R-.
// Get the values from the algorithm conf. registry, and if they do not exist
// set them to default values (Eq.(20) in Sov.J.Nucl.Phys.52:934 (1990)

  double ResDMLambda = (fConfig->Exists("Res-DeltaM-Lambda")) ?
                     fConfig->GetDouble("Res-DeltaM-Lambda") : 0.56; /*GeV*/
  double ResDMSigma  = (fConfig->Exists("Res-DeltaM-Sigma")) ?
                     fConfig->GetDouble("Res-DeltaM-Sigma")  : 0.20; /*GeV*/

  //----- get final state charm particle & initial state (struck) nucleon

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();

  int pdgc = sc_params.GetInt("fs-charm-hadron-pdgc");

  bool isLambda = (pdgc == kPdgLambdacP);
  bool isSigma  = (pdgc == kPdgSigmacP || pdgc == kPdgSigmacPP);

  if      ( isLambda ) return ResDMLambda;
  else if ( isSigma  ) return ResDMSigma;
  else
    assert(1);

  return 0;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::MRes(const Interaction * interaction) const
{
  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();

  int pdgc = sc_params.GetInt("fs-charm-hadron-pdgc");

  double MR = PDGLibrary::Instance()->Find(pdgc)->Mass();

  return MR;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::MRes2(const Interaction * interaction) const
{
  double MR  = this->MRes(interaction);

  double MR2 = MR*MR;

  return MR2;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::vR_minus(const Interaction * interaction) const
{
  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();

  double Q2     = sc_params.GetDouble("Q2");
  double DeltaR = this->ResDM(interaction);
  double MR     = MRes(interaction);
  double MN     = kNucleonMass;
  double MN2    = kNucleonMass_2;

  double vR = ( pow(MR-DeltaR,2) - MN2 + Q2 ) / (2*MN);

  return vR;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::vR_plus(const Interaction * interaction) const
{
  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();

  double Q2     = sc_params.GetDouble("Q2");
  double DeltaR = this->ResDM(interaction);
  double MR     = MRes(interaction);
  double MN     = kNucleonMass;
  double MN2    = kNucleonMass_2;

  double vR = ( pow(MR+DeltaR,2) - MN2 + Q2 ) / (2*MN);

  return vR;
}
//____________________________________________________________________________
const IntegratorI * KovalenkoQELCharmPXSec::Integrator(void) const
{
  //----- get the specified integrator

  string integrator_name;

  if( fConfig->Exists("integrator") )
                          fConfig->Get("integrator", integrator_name );
  else integrator_name = "genie::Simpson1D";

  //-- ask the AlgFactory for the integrator algorithm

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * intg_alg_base = algf->GetAlgorithm(integrator_name);

  const IntegratorI * integrator =
                           dynamic_cast<const IntegratorI *> (intg_alg_base);

  assert( integrator != 0 );

  return integrator;
}
//____________________________________________________________________________
double KovalenkoQELCharmPXSec::SumF2(const Interaction * interaction) const
{
// Returns F1^2 (Q^2=0) + FA^2 (Q^2 = 0) for the normalization factor.
// Get the values from the algorithm conf. registry, and if they do not exist
// set them to default values I computed using Sov.J.Nucl.Phys.52:934 (1990)

  double F2LambdaP = (fConfig->Exists("F1^2+FA^2-LambdaP")) ?
                            fConfig->GetDouble("F1^2+FA^2-LambdaP") : 2.07;
  double F2SigmaP  = (fConfig->Exists("F1^2+FA^2-SigmaP")) ?
                             fConfig->GetDouble("F1^2+FA^2-SigmaP") : 0.71;
  double F2SigmaPP = (fConfig->Exists("F1^2+FA^2-SigmaPP")) ?
                            fConfig->GetDouble("F1^2+FA^2-SigmaPP") : 1.42;

  //----- get final state charm particle & initial state (struck) nucleon

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();

  int pdgc = sc_params.GetInt("fs-charm-hadron-pdgc");

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  if      ( pdgc == kPdgLambdacP && isN ) return F2LambdaP;
  else if ( pdgc == kPdgSigmacP  && isN ) return F2SigmaP;
  else if ( pdgc == kPdgSigmacPP && isP ) return F2SigmaPP;
  else
    assert(1);

  return 0;
}
//____________________________________________________________________________
void KovalenkoQELCharmPXSec::AssertProcessValidity(
                                        const Interaction * interaction) const
{
  //----- get final state charm particle & initial state (struck) nucleon

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();

  assert( sc_params.Exists("fs-charm-hadron-pdgc") );

  int pdgc = sc_params.GetInt("fs-charm-hadron-pdgc");

  bool isP = pdg::IsProton ( init_state.GetTarget().StruckNucleonPDGCode() );
  bool isN = pdg::IsNeutron( init_state.GetTarget().StruckNucleonPDGCode() );

  assert(isP || isN);

  //----- make sure we are dealing with one of the following channels:
  //
  //   v + n --> mu- + Lambda_{c}^{+} (2285)
  //   v + n --> mu- + Sigma_{c}^{+} (2455)
  //   v + p --> mu- + Sigma_{c}^{++} (2455)

  assert(
         (pdgc == kPdgLambdacP && isN) || /* v + n -> l + #Lambda_{c}^{+} */
         (pdgc == kPdgSigmacP  && isN) || /* v + n -> l + #Sigma_{c}^{+}  */
         (pdgc == kPdgSigmacPP && isP)    /* v + p -> l + #Sigma_{c}^{++} */
       );
}
//____________________________________________________________________________
