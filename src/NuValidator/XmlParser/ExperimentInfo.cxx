//_____________________________________________________________________________
/*!

\class    genie::nuvld::ExperimentInfo

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "ExperimentInfo.h"

namespace genie {
namespace nuvld {
  
//______________________________________________________________________________________
ostream & operator << (ostream & stream, const ExperimentInfo & info)
{
  stream << endl;
  stream << "printing experiment info: " << endl;

  stream << "experiment name..: " << info._name        << endl;
  stream << "facility.........: " << info._facility    << endl;
  stream << "detector.........: " << info._detector    << endl;
  stream << "beam.............: " << info._beam        << endl;
  stream << "target...........: " << info._target      << endl;
  stream << "year [begin].....: " << info._year_begin  << endl;
  stream << "year [end].......: " << info._year_end    << endl;
  stream << "exposure.........: " << info._exposure    << endl;
  stream << "energy [min].....: " << info._energy_min  << endl;
  stream << "energy [max].....: " << info._energy_max  << endl;
  stream << "comment..........: " << info._comment     << endl;

  return stream;
}
//______________________________________________________________________________________
ExperimentInfo::ExperimentInfo()
{
  _name           = "";
  _comment        = "";
  _facility       = "";
  _detector       = "";
  _beam           = "";
  _target         = "";
  _exposure       = "";
  _exposure_units = "";
  _year_begin     = "";
  _year_end       = "";
  _energy_min     = "";
  _energy_max     = "";
  _energy_units   = "";
  _energy_frame   = "";
}
//______________________________________________________________________________________
ExperimentInfo::ExperimentInfo(const ExperimentInfo & /*info*/)
{

}
//______________________________________________________________________________________
void ExperimentInfo::Add(string  key, string value)
{
   if (key.compare("name") == 0)            _name     = value;
   else if (key.compare("comment") == 0)    _comment  = value;
   else if (key.compare("facility") == 0)   _facility = value;
   else if (key.compare("detector") == 0)   _detector = value;
   else if (key.compare("target") == 0)     _target   = value;
   else if (key.compare("beam") == 0)       _beam     = value;
   else if (key.compare("exposure") == 0)   _exposure = value;
   else if (key.compare("year") == 0)       ProcessYearField(value);
   else if (key.compare("E") == 0)          ProcessEnergyField(value);
   else 
       cerr << "ExperimentInfo::add(string, string) : unknown key" << endl;
}
//______________________________________________________________________________________
void ExperimentInfo::ProcessYearField(string value)
{
  unsigned int ipos = value.find('-');

  if(ipos < value.length()) {
     _year_begin.assign( value, 0,      ipos           );
     _year_end.assign(   value, ipos+1, value.length() );
  } else {
     _year_begin = value;
     _year_end   = value;
  }
}
//______________________________________________________________________________________
void ExperimentInfo::ProcessEnergyField(string value)
{
  unsigned int ipos = value.find('-');

  if(ipos < value.length()) {
     _energy_min.assign( value, 0,      ipos           );
     _energy_max.assign( value, ipos+1, value.length() );
     
  } else {
     _energy_min = value;
     _energy_max   = value;
  }
}
//______________________________________________________________________________________

} // nuvld namespace
} // genie namespace
