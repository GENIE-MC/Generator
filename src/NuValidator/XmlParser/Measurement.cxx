//_____________________________________________________________________________
/*!

\class    genie::nuvld::Measurement

\brief

\author  Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created August 2003         
*/
//_____________________________________________________________________________

#include "Measurement.h"

namespace genie {
namespace nuvld {
  
//_______________________________________________________________________________
Measurement::Measurement()
{
  _data = new vector<RecordBase *>;
} 
//_______________________________________________________________________________
Measurement::Measurement(const Measurement & meas)
{

}
//_______________________________________________________________________________
void Measurement::Add(MeasurementHeader * header)
{
  _header = header;
}
//_______________________________________________________________________________
void Measurement::Add(RecordBase * rec)
{
  _data->push_back(rec);
}
//_______________________________________________________________________________

} // nuvld namespace
} // genie namespace
