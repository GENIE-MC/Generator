//____________________________________________________________________________
/*!

\class    genie::QELXSec

\brief    Computes the Quasi Elastic (QEL) cross section.

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include <iostream>

#include "AlgFactory/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/QELXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QELXSec::QELXSec() :
XSecAlgorithmI()
{
  fName = "genie::QELXSec";
}
//____________________________________________________________________________
QELXSec::QELXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::QELXSec";

  FindConfig();
}
//____________________________________________________________________________
QELXSec::~QELXSec()
{

}
//____________________________________________________________________________
double QELXSec::XSec(const Interaction * interaction) const
{
  //----- Get an algorithm to calculate differential cross sections dxsec/dQ2
  
  const XSecAlgorithmI * partial_xsec_alg = this->PartialXSecModel();

  //----- get the input number of dxsec/dlogQ2 points for num. integration
  //      or use a default if no number is specified

  int nbins; // <-- must be odd number (for Simpson rule)

  if ( fConfig->Exists("N-logQ2-bins") )
                                     fConfig->Get("N-logQ2-bins", nbins );
  else nbins = 101;
  
  LOG("QELXSec", pDEBUG)
                     << "Number of integration (logQ^2) bins = " << nbins;
  
  //----- get initial & final state information

  const InitialState & init_state = interaction -> GetInitialState();

  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double E    = nu_p4->Energy();
  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double Mnuc = init_state.GetTarget().StruckNucleonMass();

  double ml2   = ml   * ml; 
  double Mnuc2 = Mnuc * Mnuc;

  delete nu_p4;
        
  //----- "phase-space" cut

  double Ethr = (ml2 + 2*ml*Mnuc) / (2*Mnuc);

  LOG("QELXSec", pDEBUG) << "Neutrino Energy = " << E
                          << " and Neutrino Energy Threshold = " << Ethr;
  
  if(E < Ethr) return 0;
  
  //----- estimate the integration limits & step

  double s     = 2*Mnuc*E + Mnuc2;
  double tmp   = s + ml2 -Mnuc2;
  double tmp2  = tmp*tmp - 4*s*ml2;
  
  //if(tmp2<0) tmp2=0;
  tmp2 = sqrt(tmp2);
  
  double Q2min = -(ml2 - 0.5 * (s-Mnuc2) * (tmp-tmp2) /s);
  double Q2max = -(ml2 - 0.5 * (s-Mnuc2) * (tmp+tmp2) /s);

  LOG("QELXSec", pDEBUG) << "Q2 integration range = ("
                                       << Q2min << ", " << Q2max << ")";
  
  //if( TMath::Abs(q2max) < 1E-12 ) Q2max = -1E-12;

  double logQ2max = TMath::Log(Q2max);
  double logQ2min = TMath::Log(Q2min);
    
  double dlogQ2   = (logQ2max - logQ2min) / (nbins - 1);

  //----- define the integration grid & instantiate a FunctionMap
  
  UnifGrid grid;

  grid.AddDimension(nbins, logQ2min, logQ2max);

  FunctionMap Q2dxsec_dQ2(grid);

  //----- loop over logQ2 range, estimate/store Q2*dxsec/dQ2
  
  for(int iQ2 = 0; iQ2 < nbins; iQ2++) {
        
     double Q2 = TMath::Exp(logQ2min + iQ2 * dlogQ2);

     interaction->GetScatParamsPtr()->Set("Q2", Q2); // <-- update Q2 
         
     double partial_xsec = partial_xsec_alg->XSec(interaction);
     
     Q2dxsec_dQ2.AddPoint( Q2 * partial_xsec, iQ2 );

     LOG("QELXSec", pDEBUG) 
          << "point...." << iQ2+1 << "/" << nbins << " : "
          << "dxsec/dQ^2 (Q^2 = " << Q2 << " ) = " << partial_xsec;
  }
  
  //----- Numerical integration

  const IntegratorI * integrator = this->Integrator();
  
  double sQEtot = integrator->Integrate(Q2dxsec_dQ2);
            
  return sQEtot;
}
//____________________________________________________________________________
const XSecAlgorithmI * QELXSec::PartialXSecModel(void) const
{
// Get an algorithm to calculate differential cross sections (dxsec/dQ2)

  assert( fConfig->Exists("partial-xsec-alg-name")  &&
                                  fConfig->Exists("partial-xsec-param-set") );
 
  //-- get the QEL alg-name & param-set from the configuration registry

  string pxsec_alg_name, pxsec_param_set;

  fConfig->Get("partial-xsec-alg-name",  pxsec_alg_name  );
  fConfig->Get("partial-xsec-param-set", pxsec_param_set );

  //-- get the specified algorithm from the AlgFactory

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase =
                       algf->GetAlgorithm(pxsec_alg_name, pxsec_param_set);

  const XSecAlgorithmI * partial_xsec_alg =
                            dynamic_cast<const XSecAlgorithmI *> (algbase);

  assert(partial_xsec_alg);

  return partial_xsec_alg;
}
//____________________________________________________________________________
const IntegratorI * QELXSec::Integrator(void) const
{
// Get the integration algorithm specified in the configuration registry.
// Use Simpson1D if no one else is defined

  string integrator_name;

  if( fConfig->Exists("integrator") )
                          fConfig->Get("integrator", integrator_name );
  else integrator_name = "genie::Simpson1D";

  //-- ask the AlgFactory for the integrator algorithm

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * intg_alg_base = algf->GetAlgorithm(integrator_name);

  const IntegratorI * integrator =
                           dynamic_cast<const IntegratorI *> (intg_alg_base);

  assert(integrator);

  return integrator;
}
//____________________________________________________________________________
