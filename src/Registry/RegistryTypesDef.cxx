//____________________________________________________________________________
/*!

\class    genie::RegistryItem

\brief    A templated concrete implementation of the RegistryItemI interface.
          Provides an arbitrary basic type (bool, int, double, string) value
          for RegistryI-type containers.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include "Registry/RegistryTypesDef.h"

using std::endl;

//____________________________________________________________________________
RgAlg::RgAlg() 
{

}
//____________________________________________________________________________
RgAlg::RgAlg(string n, string c) : 
name(n), 
config(c) 
{

}
//____________________________________________________________________________
RgAlg::~RgAlg()
{

}
//____________________________________________________________________________
ostream & operator << (ostream & stream, const RgAlg & alg)
{
  stream << alg.name << "/" << alg.config;
  return stream;
}
//____________________________________________________________________________
RgAlg & RgAlg::operator = (const RgAlg & alg)
{
  this->name   = alg.name;
  this->config = alg.config;
  return (*this);
}
//____________________________________________________________________________

