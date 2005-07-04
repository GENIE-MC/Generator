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

#include <vector>
#include <ostream>

using std::vector;
using std::ostream;

namespace genie {

class PDGCodeList : public vector<int> {

public :

  PDGCodeList();
  PDGCodeList(size_type n);
  ~PDGCodeList();

  //-- override the vector<int> insertion methods to explicitly check for
  //   PDG code validity and that no PDG code is listed more than once

  void push_back  (int pdg_code);
  void insert     (iterator pos, size_type n, const int& x);

  //-- print 
  
  void Print(ostream & stream) const;
  friend ostream & operator << (ostream & stream, const PDGCodeList & list);

private:

  bool CheckPDGCode        (int pdg_code);
  bool ExistsInPDGLibrary  (int pdg_code);
  bool ExistsInPDGCodeList (int pdg_code);
};

}      // genie namespace

#endif // _PDG_CODE_LIST_H_
