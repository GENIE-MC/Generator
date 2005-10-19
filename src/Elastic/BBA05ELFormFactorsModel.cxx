//____________________________________________________________________________
/*!

\class    genie::BBA05ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2005 parameterization.

\ref

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 19, 2005

*/
//____________________________________________________________________________

#include "Elastic/BBA05ELFormFactorsModel.h"
#include "Interaction/Interaction.h"

using namespace genie;

//____________________________________________________________________________
BBA05ELFormFactorsModel::BBA05ELFormFactorsModel() :
ELFormFactorsModelI()
{

}
//____________________________________________________________________________
BBA05ELFormFactorsModel::BBA05ELFormFactorsModel(const char * param_set) :
ELFormFactorsModelI(param_set)
{
  fName = "genie::BBA05ELFormFactorsModel";

  this->FindConfig();
}
//____________________________________________________________________________
BBA05ELFormFactorsModel::~BBA05ELFormFactorsModel()
{

}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Ge(const Interaction * /*interaction*/) const
{
  return 0;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gm(const Interaction * /*interaction*/) const
{
  return 0;
}
//____________________________________________________________________________

