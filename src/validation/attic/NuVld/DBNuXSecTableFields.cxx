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

#include "ValidationTools/NuVld/DBNuXSecTableFields.h"

using namespace genie::nuvld;

ClassImp(DBNuXSecTableFields)

//____________________________________________________________________________
DBNuXSecTableFields::DBNuXSecTableFields() :
DBTableFields()
{
  AddField("name");
  AddField("measurement_tag");
  AddField("xsec");
  AddField("stat_err_p");
  AddField("stat_err_m");
  AddField("syst_err_p");
  AddField("syst_err_m");
  AddField("xsec_units");
  AddField("xsec_norm");
  AddField("stat_err_type");
  AddField("syst_err_type");
  AddField("E");
  AddField("E_min");
  AddField("E_max");
  AddField("E_units");
  AddField("E_frame");
}
//____________________________________________________________________________
DBNuXSecTableFields::DBNuXSecTableFields(const DBTableFields * fields):
DBTableFields(fields)
{

}
//____________________________________________________________________________
DBNuXSecTableFields::~DBNuXSecTableFields()
{

}
//____________________________________________________________________________


