//_____________________________________________________________________________
/*!

\class    genie::nuvld::Citation

\brief    Encapsulates an XML data file citation

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "Citation.h"

namespace genie {
namespace nuvld {
  
//___________________________________________________________________________
ostream & operator << (ostream & stream, const Citation & ref)
{
  stream << "[-] Reference " << endl;
  stream << " |-> author...........: " << ref._author   << endl;
  stream << " |-> journal..........: " << ref._journal  << endl;
  stream << " |-> year.............: " << ref._year     << endl;

  return stream;
}
//___________________________________________________________________________
Citation::Citation() :
_author(string("")),
_journal(string("")),
_year(string(""))
{

}
//___________________________________________________________________________
Citation::Citation(const Citation & /*ref*/)
{

}
//___________________________________________________________________________
void Citation::Add(string key, string value)
{
  if      (key.compare("author")   == 0) _author   = value;
  else if (key.compare("journal")  == 0) _journal  = value;
  else if (key.compare("year")     == 0) _year     = value;
  else cerr << "Citation::add(string, string): unknown key" << endl;
}
//___________________________________________________________________________

} // nuvld namespace
} // genie namespace

