//_____________________________________________________________________________
/*!

\class    genie::nuvld::BeamFluxBin

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "BeamFluxBin.h"

namespace genie {
namespace nuvld {
  
//______________________________________________________________________________________
ostream & operator << (ostream & stream, const BeamFluxBin & bin)
{
  stream << "[Emin = " << bin._Emin << ", Emax = " << bin._Emax << "] " 
         << " - <E> = " << bin._E << " ---> flux = " << bin._Flux 
         << " (+" << bin._Flux_perr << ", -" << bin._Flux_nerr << ")" << endl;
	 
  return stream;
}
//______________________________________________________________________________________
BeamFluxBin::BeamFluxBin() :
_E(""),
_Emin(""),
_Emax(""),
_Flux(""),
_Flux_perr(""),
_Flux_nerr("")
{

}
//______________________________________________________________________________________
BeamFluxBin::BeamFluxBin(
                 string E, string Emin, string Emax, string F, string dFp, string dFn) :
_E(E),
_Emin(Emin),
_Emax(Emax),
_Flux(F),
_Flux_perr(dFp),
_Flux_nerr(dFn)
{

}
//______________________________________________________________________________________
BeamFluxBin::BeamFluxBin(const BeamFluxBin & /*bin*/)
{

}
//______________________________________________________________________________________

} // nuvld namespace
} // genie namespace
