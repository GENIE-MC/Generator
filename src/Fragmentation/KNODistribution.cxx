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

#include <TSystem.h>

#include "Fragmentation/KNODistribution.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"

using std::ostringstream;
using namespace genie;

//____________________________________________________________________________
KNODistribution::KNODistribution() :
Algorithm("genie::KNODistribution")
{
  fKNOSpline = 0;
}
//____________________________________________________________________________
KNODistribution::KNODistribution(string config) :
Algorithm("genie::KNODistribution", config)
{
  fKNOSpline = 0;
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
      bool inrange = n_avn<fKNOSpline->XMax() && n_avn>fKNOSpline->XMin();
      if(inrange) 
          return fKNOSpline->Evaluate(n_avn);
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
void KNODistribution::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadKNO();
}
//____________________________________________________________________________
void KNODistribution::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadKNO();
}
//____________________________________________________________________________
void KNODistribution::LoadKNO(void)
{
  assert(gSystem->Getenv("GENIE"));

  string basedir    = gSystem->Getenv("GENIE");
  string defknodata = basedir + "/data/kno/KNO.dat";
  string knodata    = fConfig->GetStringDef("kno-data",defknodata);

  fKNOSpline = new Spline(knodata);
}
//____________________________________________________________________________

