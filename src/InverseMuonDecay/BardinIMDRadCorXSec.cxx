//____________________________________________________________________________
/*!

\class    genie::BardinIMDRadCorXSec

\brief    Computes the Inverse Muon Decay cross section using the Bardin -
          Dokuchaeva model which includes all 1-loop radiative corrections. \n

          This algorithm merely integrates the Bardin differential IMD cross
          section. The specific differential cross section algorithm is
          specified in this algorithm's XML config file.

          The exact 'type' of the cross section depends on the specified
          differential IMD cross section algorithm. It can be a 'trully'
          inclusive IMD cross section or a cross section where part of the
          brem cross section contribution is subtracted
          (for futher details, see the documentation of the Bardin-Dokuchaeva
          model diffential IMD cross section algorithms and the Bardin paper,
          cited below).
          
          BardinIMDRadCorXSec is a concrete implementation of the
          XSecAlgorithmI interface. \n

\ref      D.Yu.Bardin and V.A.Dokuchaeva, Nucl.Phys.B287:839 (1987)
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  Fabruary 14, 2005

*/
//____________________________________________________________________________

#include <iostream>

#include "AlgFactory/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "InverseMuonDecay/BardinIMDRadCorXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BardinIMDRadCorXSec::BardinIMDRadCorXSec() :
XSecAlgorithmI()
{
  fName = "genie::BardinIMDRadCorXSec";
}
//____________________________________________________________________________
BardinIMDRadCorXSec::BardinIMDRadCorXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::BardinIMDRadCorXSec";

  FindConfig();
}
//____________________________________________________________________________
BardinIMDRadCorXSec::~BardinIMDRadCorXSec()
{

}
//____________________________________________________________________________
double BardinIMDRadCorXSec::XSec(const Interaction * interaction) const
{
  const XSecAlgorithmI * pxsec      = this->DiffXSec();
  const IntegratorI *    integrator = this->Integrator();

  const int    nsteps  = 201;
  const double min     = 0;
  const double max     = 0.999;
  const double step    = (max-min)/(nsteps-1);

  UnifGrid grid;

  grid.AddDimension(nsteps, min, max);

  FunctionMap fmap(grid);

  //-- all kinematical cuts (energy threshold, physical y range) are applied
  //   within the differential cross section algorithm - returns 0 if kinematic
  //   params are not valid.
  
  for(int i = 0; i < nsteps; i++) {

    double y  = min + i * step;

    interaction->GetScatParamsPtr()->Set("y",y);
    
    double dsig_dy  = pxsec->XSec(interaction);

    fmap.AddPoint(dsig_dy, i);
  }

  double sig = integrator->Integrate(fmap);

  LOG("InverseMuDecay", pDEBUG) << "*** xsec[IMD] = " << sig;

  return sig;
}
//____________________________________________________________________________
const XSecAlgorithmI * BardinIMDRadCorXSec::DiffXSec(void) const
{
 LOG("InverseMuDecay", pINFO)
                           << "Getting requested differential xsec algorithm";

 assert(fConfig->Exists("partial-xsec-alg-name") &&
                                  fConfig->Exists("partial-xsec-param-set"));

 string alg_name  = fConfig->GetString("partial-xsec-alg-name");
 string param_set = fConfig->GetString("partial-xsec-param-set");

 //-- Get the requested cross section calculator from the AlgFactory

 AlgFactory * algf = AlgFactory::Instance();

 const Algorithm * alg_base = algf->GetAlgorithm(alg_name, param_set);

 const XSecAlgorithmI * xs = dynamic_cast<const XSecAlgorithmI *> (alg_base);

 assert(xs);

 return xs;
}
//____________________________________________________________________________
const IntegratorI * BardinIMDRadCorXSec::Integrator(void) const
{
 LOG("InverseMuDecay", pINFO) << "Getting requested integrator algorithm";

 // read integrator name from config or set default
 string integrator_name = ( fConfig->Exists("integrator-name") ) ?
                  fConfig->GetString("integrator-name") : "genie::Simpson1D";
 
 // get integrator algorithm

 AlgFactory * algf = AlgFactory::Instance();

 const Algorithm * alg_base = algf->GetAlgorithm(integrator_name);

 const IntegratorI * intgr = dynamic_cast<const IntegratorI *> (alg_base);

 assert(intgr);

 return intgr;
}
//____________________________________________________________________________
