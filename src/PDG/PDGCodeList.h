//____________________________________________________________________________
/*!

\class   genie::PDGCodeList

\brief   A list of PDG codes

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 13, 2005

*/
//____________________________________________________________________________

#ifndef _PDG_CODE_LIST_H_
#define _PDG_CODE_LIST_H_

#include <list>
#include <ostream>

using std::list;
using std::ostream;

namespace genie {

class PDGCodeList : public list<int> {

public :

  PDGCodeList();
  PDGCodeList(size_type n);
  ~PDGCodeList();

  //-- override all list<int> insertion methods to explicitly check for
  //   PDG code validity

  void push_front (int pdg_code);
  void push_back  (int pdg_code);

  //-- print list
  
  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const PDGCodeList & list);

private:

  bool ExistsInPDGLibrary(int pdg_code);
};

}      // genie namespace

#endif // _PDG_CODE_LIST_H_
