//____________________________________________________________________________
/*!

\class   genie::PDGCodeList

\brief   A list of PDG codes

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 13, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PDG_CODE_LIST_H_
#define _PDG_CODE_LIST_H_

#include <vector>
#include <ostream>

using std::vector;
using std::ostream;

namespace genie {

class PDGCodeList;
ostream & operator << (ostream & stream, const PDGCodeList & list);

class PDGCodeList : public vector<int> {

public :

  PDGCodeList(bool allowdup=false);
  PDGCodeList(size_type n, bool allowdup=false);
  PDGCodeList(const PDGCodeList & list);
  ~PDGCodeList();

  //! override the vector<int> insertion methods to explicitly check for
  //! PDG code validity and that no PDG code is listed more than once
  void push_back  (int pdg_code);
  void insert     (iterator pos, size_type n, const int& x);

  //! PDG code checks used by PDGCodeList
  bool CheckPDGCode        (int pdg_code) const;
  bool ExistsInPDGLibrary  (int pdg_code) const;
  bool ExistsInPDGCodeList (int pdg_code) const;

  //! copy / print
  void Copy  (const PDGCodeList & list);
  void Print (ostream & stream) const;

  //! check state
  bool DuplEntriesAllowed(void) const { return fAllowDuplicateEntries; }

  //! overloaded operators
  PDGCodeList &    operator =  (const PDGCodeList & list);
  friend ostream & operator << (ostream & stream, const PDGCodeList & list);

private:

  bool fAllowDuplicateEntries; ///< allow duplicate entries in the list?
};

}      // genie namespace

#endif // _PDG_CODE_LIST_H_
