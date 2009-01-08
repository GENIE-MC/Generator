//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 01, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/XmlCitation.h"

namespace genie {
namespace nuvld {
  
//___________________________________________________________________________
ostream & operator << (ostream & stream, const XmlCitation & ref)
{
  stream << "[-] Reference " << endl;
  stream << " |-> author...........: " << ref._author   << endl;
  stream << " |-> journal..........: " << ref._journal  << endl;
  stream << " |-> year.............: " << ref._year     << endl;

  return stream;
}
//___________________________________________________________________________
XmlCitation::XmlCitation() :
_author(string("")),
_journal(string("")),
_year(string(""))
{

}
//___________________________________________________________________________
XmlCitation::XmlCitation(const XmlCitation & /*ref*/)
{

}
//___________________________________________________________________________
void XmlCitation::Add(string key, string value)
{
  if      (key.compare("author")   == 0) _author   = value;
  else if (key.compare("journal")  == 0) _journal  = value;
  else if (key.compare("year")     == 0) _year     = value;
  else cerr << "XmlCitation::add(string, string): unknown key" << endl;
}
//___________________________________________________________________________

} // nuvld namespace
} // genie namespace

