//_____________________________________________________________________________
/*!

\class    genie::nuvld::BeamFluxSpectrum

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "BeamFluxSpectrum.h"

namespace genie {
namespace nuvld {
  
//______________________________________________________________________________________
ostream & operator << (ostream & stream, const BeamFluxSpectrum & spectrum)
{
  stream << endl;
  stream << "----------------- Printing beam flux spectrum -------------------" << endl;
  
  stream << " --> flux units   = " << spectrum._flux_units << endl;
  stream << " --> energy units = " << spectrum._E_units    << endl;
  stream << " --> energy frame = " << spectrum._E_frame    << endl;

  vector<BeamFluxBin *>::iterator bin_iter;
  for(bin_iter = spectrum._spectrum->begin(); 
                         bin_iter != spectrum._spectrum->end(); ++bin_iter)
                                                                 stream << *(*bin_iter);
  return stream;
}
//______________________________________________________________________________________
BeamFluxSpectrum::BeamFluxSpectrum()
{
  _flux_units = "cm-2 sec-1 sr-1 GeV-1";
  _E_units    = "GeV";
  _E_frame    = "lab";
  
  _spectrum = new vector<BeamFluxBin *>;
}
//______________________________________________________________________________________
BeamFluxSpectrum::BeamFluxSpectrum(const BeamFluxSpectrum & spectrum)
{

}
//______________________________________________________________________________________
void BeamFluxSpectrum::Add(BeamFluxBin * bin)
{
  _spectrum->push_back(bin);

  sort( _spectrum->begin(), _spectrum->end(), less_energy() );
}
//______________________________________________________________________________________

} // nuvld namespace
} // genie namespace
