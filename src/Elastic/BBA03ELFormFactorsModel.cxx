//____________________________________________________________________________
/*!

\class    genie::BBA03ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2003 parameterization.

\ref

\author

\created  October 19, 2005

*/
//____________________________________________________________________________

#include "Elastic/BBA03ELFormFactorsModel.h"
#include "Interaction/Interaction.h"

using namespace genie;

//____________________________________________________________________________
BBA03ELFormFactorsModel::BBA03ELFormFactorsModel() :
ELFormFactorsModelI()
{
  fName = "genie::BBA03ELFormFactorsModel";
}
//____________________________________________________________________________
BBA03ELFormFactorsModel::BBA03ELFormFactorsModel(const char * param_set) :
ELFormFactorsModelI(param_set)
{
  fName = "genie::BBA03ELFormFactorsModel";

  this->FindConfig();
}
//____________________________________________________________________________
BBA03ELFormFactorsModel::~BBA03ELFormFactorsModel()
{

}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gep(double) const
{
  return 0;
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gmp(double) const
{
  return 0;
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gen(double) const
{
  return 0;
}
//____________________________________________________________________________
double BBA03ELFormFactorsModel::Gmn(double) const
{
  return 0;
}
//____________________________________________________________________________

