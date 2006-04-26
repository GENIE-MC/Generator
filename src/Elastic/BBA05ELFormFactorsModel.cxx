//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 19, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
double BBA05ELFormFactorsModel::Gep(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gmp(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gen(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
double BBA05ELFormFactorsModel::Gmn(const Interaction * ) const
{
  return 0;
}
//____________________________________________________________________________
