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

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/QELXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
#include "Numerical/IntegratorI.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QELXSec::QELXSec() :
XSecAlgorithmI("genie::QELXSec")
{

}
//____________________________________________________________________________
QELXSec::QELXSec(string config) :
XSecAlgorithmI("genie::QELXSec", config)
{

}
//____________________________________________________________________________
QELXSec::~QELXSec()
{

}
//____________________________________________________________________________
double QELXSec::XSec(const Interaction * interaction) const
{
  //----- Get an algorithm to calculate differential cross sections dxsec/dQ2
  const XSecAlgorithmI * partial_xsec_alg =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));

  //----- get the input number of dxsec/dlogQ2 points for num. integration
  //      or use a default if no number is specified
  //      (must be odd number for the Simpson rule)
  int nbins = fConfig->GetIntDef("N-logQ2-bins", 101);
  LOG("QELXSec", pDEBUG)
                 << "Number of integration (logQ^2) bins = " << nbins;

  //----- get initial & final state information

  const InitialState & init_state = interaction -> GetInitialState();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double Mnuc = init_state.GetTarget().StruckNucleonMass();
  double ml2  = ml * ml;

  //----- "phase-space" cut
  double Ethr = (ml2 + 2*ml*Mnuc) / (2*Mnuc);

  LOG("QELXSec", pDEBUG)
       << "Computing QEL XSec for Ev = " << E
                            << " / Neutrino Energy Threshold = " << Ethr;
  if(E <= Ethr) {
     LOG("QELXSec", pINFO) << "Ev = " << E << " <= Ethreshold = "<< Ethr;
     return 0;
  }

  //----- estimate the integration limits & step
  Range1D_t  rQ2 = utils::kinematics::Q2Range_M(interaction);

  LOG("QELXSec", pDEBUG) << "Q2 integration range = ("
                                    << rQ2.min << ", " << rQ2.max << ")";

  double logQ2max = TMath::Log(rQ2.max);
  double logQ2min = TMath::Log(rQ2.min);
  double dlogQ2   = (logQ2max - logQ2min) / (nbins - 1);

  //----- define the integration grid & instantiate a FunctionMap
  UnifGrid grid;
  grid.AddDimension(nbins, logQ2min, logQ2max);

  FunctionMap Q2dxsec_dQ2(grid);

  //----- loop over logQ2 range, estimate/store Q2*dxsec/dQ2

  for(int iQ2 = 0; iQ2 < nbins; iQ2++) {

     double Q2 = TMath::Exp(logQ2min + iQ2 * dlogQ2);
     interaction->GetKinematicsPtr()->SetQ2(Q2);

     double partial_xsec = partial_xsec_alg->XSec(interaction);
     Q2dxsec_dQ2.AddPoint( Q2 * partial_xsec, iQ2 );

     LOG("QELXSec", pDEBUG)
          << "point...." << iQ2+1 << "/" << nbins << " : "
          << "dxsec/dQ^2 (Q^2 = " << Q2 << " ) = " << partial_xsec;
  }

  //----- Numerical integration
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");

  AlgFactory * algf = AlgFactory::Instance();
  const IntegratorI * integrator =
          dynamic_cast<const IntegratorI *> (algf->GetAlgorithm(intgr));

  double sQEtot = integrator->Integrate(Q2dxsec_dQ2);
  return sQEtot;
}
//____________________________________________________________________________
