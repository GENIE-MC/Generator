//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlBeamFluxSpectrum

\brief    Encapsulates an XML data file beam flux

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _BEAM_FLUX_SPECTRUM_H_
#define _BEAM_FLUX_SPECTRUM_H_

#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>

#include "ValidationTools/NuVld/XmlBeamFluxBin.h"

using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::binary_function;

namespace genie {
namespace nuvld {

struct less_energy :
              public binary_function <XmlBeamFluxBin *, XmlBeamFluxBin *, bool>  {

   bool operator()(XmlBeamFluxBin * a, XmlBeamFluxBin * b) {
     
      return atof(a->MeanEnergy().c_str()) < atof(b->MeanEnergy().c_str()); 
   }
};

class XmlBeamFluxSpectrum
{
public:

  XmlBeamFluxSpectrum();
  XmlBeamFluxSpectrum(const XmlBeamFluxSpectrum & spectrum);
     
  virtual ~XmlBeamFluxSpectrum() { delete _spectrum; };

  void   Add(XmlBeamFluxBin * bin);

  const vector<XmlBeamFluxBin *> & GetSpectrum(void) const { return *_spectrum; }

  virtual string FluxUnits(void)   const { return _flux_units; }
  virtual string EnergyUnits(void) const { return _E_units;    }
  virtual string EnergyFrame(void) const { return _E_frame;    }

  virtual void SetFluxUnits(string units)   { _flux_units = units; }
  virtual void SetEnergyUnits(string units) { _E_units    = units; }
  virtual void SetEnergyFrame(string frame) { _E_frame    = frame; }
          
  friend ostream & operator <<(ostream & stream, const XmlBeamFluxSpectrum & spectrum);

protected:

  string _flux_units;
  string _E_units;
  string _E_frame;          
     
  vector<XmlBeamFluxBin *> *  _spectrum;
};

} // nuvld namespace
} // genie namespace

#endif // _BEAM_FLUX_SPECTRUM_H_
