//____________________________________________________________________________
/*!

\class    genie::CharmFractionTable

\brief    A Charm Fractions table.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _CHARM_FRACTION_TABLE_H_
#define _CHARM_FRACTION_TABLE_H_

namespace genie {

class CharmFractionTable {

public:

  CharmFractionTable();
  ~CharmFractionTable();

  int    GenerateCharmHadron (double E) const;    
  double Fraction            (double E, int HadronPdgCode) const;

private:
  
};

}         // genie namespace

#endif    // _CHARM_FRACTION_TABLE_H_

