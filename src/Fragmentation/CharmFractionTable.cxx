//____________________________________________________________________________
/*!

\class    genie::CharmFractionTable

\brief    A Charm Fractions table.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  November 22, 2004

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

