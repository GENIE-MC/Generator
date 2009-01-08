//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlBeamFluxBin

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _BEAM_FLUX_BIN_H_
#define _BEAM_FLUX_BIN_H_

#include <string>
#include <iostream>

using std::string;
using std::ostream;
using std::cout;
using std::endl;

namespace genie {
namespace nuvld {
  
class XmlBeamFluxBin
{
public:

     XmlBeamFluxBin();
     XmlBeamFluxBin(string E, string Emin, string Emax, string F, string dFp, string dFn);
     XmlBeamFluxBin(const XmlBeamFluxBin & info);
     
     virtual ~XmlBeamFluxBin() { };

     friend ostream & operator <<(ostream & stream, const XmlBeamFluxBin & info);

     virtual const string MeanEnergy (void) const { return _E;         }
     virtual const string MinEnergy  (void) const { return _Emin;      }
     virtual const string MaxEnergy  (void) const { return _Emax;      }
     virtual const string Flux       (void) const { return _Flux;      }
     virtual const string FluxPErr   (void) const { return _Flux_perr; }
     virtual const string FluxNErr   (void) const { return _Flux_nerr; }

protected:

     string  _E;
     string  _Emin;
     string  _Emax;
     string  _Flux;
     string  _Flux_perr;
     string  _Flux_nerr;
};

} // nuvld namespace
} // genie namespace

#endif // _BEAM_FLUX_BIN_H_
