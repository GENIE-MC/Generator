//____________________________________________________________________________
/*!

\class    genie::KNODistribution

\brief    Describes a KNO Distribution.
          Its data are loaded from its XML configuration file.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

*/
//____________________________________________________________________________

#include <sstream>

#include "Fragmentation/KNODistribution.h"
#include "Messenger/Messenger.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________
KNODistribution::KNODistribution() :
Algorithm()
{
  fName = "genie::KNODistribution";

  fKNOSpline = 0;
  fMaxScaledMultiplicity = 0;
}
//____________________________________________________________________________
KNODistribution::KNODistribution(const char * param_set) :
Algorithm(param_set)
{
  fName = "genie::KNODistribution";

  fKNOSpline = 0;
  fMaxScaledMultiplicity = 0;
  
  FindConfig();
}
//____________________________________________________________________________
KNODistribution::~KNODistribution()
{
  if (fKNOSpline) delete fKNOSpline;
}
//____________________________________________________________________________
double KNODistribution::Value(double n_avn) const
{
  if (fKNOSpline) {

      if(n_avn < fMaxScaledMultiplicity) return fKNOSpline->Eval(n_avn);

      else {        
          LOG("KNO", pDEBUG) 
                  << "n/<n> = " << n_avn << " > maximum scaled multiplicity";
          return 0;
      }

  } else {

      LOG("KNO", pERROR) << "KNO spline is not built!";
      return 0;
  }

  return 0;
}
//____________________________________________________________________________
void KNODistribution::KNOFromXmlConfig2Spline(void)
{
  LOG("KNO", pDEBUG) << "Loading KNO data";

  if (fKNOSpline) delete fKNOSpline;

  assert( fConfig->Exists("n-bins") );
  
  int nbins = fConfig->GetInt("n-bins");

  double x[nbins], y[nbins], xmax = -1;

  for(int ibin = 0; ibin < nbins; ibin++) {

     ostringstream x_key, y_key;

     x_key << "n/avn--bin="    << ibin;
     y_key << "P(n)*avn--bin=" << ibin;
  
     assert( fConfig->Exists(x_key.str()) && fConfig->Exists(y_key.str()) );

     x[ibin] = fConfig->GetDouble(x_key.str());
     y[ibin] = fConfig->GetDouble(y_key.str());

     xmax = TMath::Max(xmax, x[ibin]);

     LOG("KNO", pINFO)
                  << "n/<n> = " << x[ibin] << " --> <n>*P(n) = " << y[ibin];
  }

  fKNOSpline = new TSpline3("fKNOSpline", x, y, nbins, "0", 0, xmax);

  fMaxScaledMultiplicity = xmax;

  LOG("KNO", pINFO)
                << "Maximum scaled multiplicity = " << fMaxScaledMultiplicity;
}
//____________________________________________________________________________
void KNODistribution::Configure(const Registry & config)
{
  Algorithm::Configure(config);

  KNOFromXmlConfig2Spline();
}
//____________________________________________________________________________
void KNODistribution::Configure(string param_set)
{
  Algorithm::Configure(param_set);

  KNOFromXmlConfig2Spline();
}
//____________________________________________________________________________
