//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Framework/Registry/RegistryItemTypeDef.h"

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

