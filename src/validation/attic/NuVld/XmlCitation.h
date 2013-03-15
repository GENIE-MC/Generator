//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlCitation

\brief    A citation for a measurement stored in an NuValidator XML file.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _CITATION_H_
#define _CITATION_H_

#include <string>
#include <iostream>

using std::string;
using std::ostream;
using std::endl;
using std::cerr;

namespace genie {
namespace nuvld {

const int c_ref_ntags = 3;

const string c_ref_tag[c_ref_ntags] = {
   string("author"),
   string("journal"),
   string("year")
};

class XmlCitation
{
public:
  XmlCitation();
  XmlCitation(const XmlCitation & ref);
  virtual ~XmlCitation() { }

  void Add(string key, string value);

  virtual const string Author  (void) const { return fAuthor;  }     
  virtual const string Journal (void) const { return fJournal; }
  virtual const string Year    (void) const { return fYear;    }

  friend ostream & operator << (ostream & stream, const XmlCitation & ref);

protected:

  string  fAuthor;
  string  fJournal;
  string  fYear;
};

} // nuvld namespace
} // genie namespace

#endif // _CITATION_H_
