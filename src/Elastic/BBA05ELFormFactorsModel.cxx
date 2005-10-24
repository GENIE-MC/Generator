//____________________________________________________________________________
/*!

\class    genie::BBA05ELFormFactorsModel

\brief    Concrete implementation of the ELFormFactorsModelI interface.
          Computes elastic form factors using the BBA2005 parameterization.

\ref

\author

\created  October 19, 2005

*/
//____________________________________________________________________________

#include "Elastic/BBA05ELFormFactorsModel.h"
#include "Interaction/Interaction.h"

using namespace genie;

//____________________________________________________________________________
BBA05ELFormFactorsModel::BBA05ELFormFactorsModel() :
ELFormFactorsModelI("genie::BBA05ELFormFactorsModel")
{

}
//____________________________________________________________________________
BBA05ELFormFactorsModel::BBA05ELFormFactorsModel(string config) :
ELFormFactorsModelI("genie::BBA05ELFormFactorsModel", config)
{

}
//____________________________________________________________________________
BBA05ELFormFactorsModel::~BBA05ELFormFactorsModel()
{

}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gep(double) const
{
  return 0;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gmp(double) const
{
  return 0;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gen(double) const
{
  return 0;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gmn(double) const
{
  return 0;
}
//____________________________________________________________________________
