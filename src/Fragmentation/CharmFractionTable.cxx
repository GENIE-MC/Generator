//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - November 22, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Fragmentation/CharmFractionTable.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"

using namespace genie;

//____________________________________________________________________________
CharmFractionTable::CharmFractionTable() 
{

}
//____________________________________________________________________________
CharmFractionTable::~CharmFractionTable()
{

}
//____________________________________________________________________________
int CharmFractionTable::GenerateCharmHadron(double E) const
{
// Generates a charmed hadron pdg code using the charm fraction table

  return 0;
}
//____________________________________________________________________________
double CharmFractionTable::Fraction(double E, int HadronPdgCode) const
{
  return 0;
}
//____________________________________________________________________________

