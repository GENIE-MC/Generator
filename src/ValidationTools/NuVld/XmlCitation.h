//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlCitation

\brief    Encapsulates an XML data file citation

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

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

  virtual const string Author  (void) const { return _author;  }     
  virtual const string Journal (void) const { return _journal; }
  virtual const string Year    (void) const { return _year;    }

  friend ostream & operator << (ostream & stream, const XmlCitation & ref);

protected:

  string  _author;
  string  _journal;
  string  _year;
};

} // nuvld namespace
} // genie namespace

#endif // _CITATION_H_
