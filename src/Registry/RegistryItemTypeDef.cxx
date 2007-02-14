//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Registry/RegistryItemTypeDef.h"

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

