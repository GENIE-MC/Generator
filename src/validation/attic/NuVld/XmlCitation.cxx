//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
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

#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/XmlCitation.h"

namespace genie {
namespace nuvld {
  
//___________________________________________________________________________
ostream & operator << (ostream & stream, const XmlCitation & ref)
{
  stream << "[-] Reference " << endl;
  stream << " |-> author...........: " << ref.fAuthor   << endl;
  stream << " |-> journal..........: " << ref.fJournal  << endl;
  stream << " |-> year.............: " << ref.fYear     << endl;

  return stream;
}
//___________________________________________________________________________
XmlCitation::XmlCitation() :
fAuthor (string("")),
fJournal(string("")),
fYear   (string(""))
{

}
//___________________________________________________________________________
XmlCitation::XmlCitation(const XmlCitation & /*ref*/)
{

}
//___________________________________________________________________________
void XmlCitation::Add(string key, string value)
{
  if      (key.compare("author")   == 0) fAuthor   = value;
  else if (key.compare("journal")  == 0) fJournal  = value;
  else if (key.compare("year")     == 0) fYear     = value;
  else {
    LOG("NuVld", pERROR) << "Unknown key: " << key;
  }
}
//___________________________________________________________________________

} // nuvld namespace
} // genie namespace

