//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Jan 12, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/DBSFTableFields.h"

using namespace genie::nuvld;

ClassImp(DBSFTableFields)

//____________________________________________________________________________
DBSFTableFields::DBSFTableFields() :
DBTableFields()
{
  AddField("name");
  AddField("measurement_tag");
  AddField("sf");
  AddField("stat_err_p");
  AddField("stat_err_m");
  AddField("syst_err_p");
  AddField("syst_err_m");
  AddField("R");
  AddField("p");
  AddField("x");
  AddField("Q2");
}
//____________________________________________________________________________
DBSFTableFields::DBSFTableFields(const DBTableFields * fields):
DBTableFields(fields)
{

}
//____________________________________________________________________________
DBSFTableFields::~DBSFTableFields()
{

}
//____________________________________________________________________________


