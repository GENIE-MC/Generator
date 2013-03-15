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

#include "ValidationTools/NuVld/DBElDiffXSecTableFields.h"

using namespace genie::nuvld;

ClassImp(DBElDiffXSecTableFields)

//____________________________________________________________________________
DBElDiffXSecTableFields::DBElDiffXSecTableFields () :
DBTableFields()
{
  AddField("name");
  AddField("measurement_tag");
  AddField("Sigma");
  AddField("Sigma_units");
  AddField("dSigma");
  AddField("E");
  AddField("E_units");
  AddField("EP");
  AddField("EP_units");
  AddField("Theta");
  AddField("Theta_units");
  AddField("Q2");
  AddField("Q2_units");
  AddField("W2");
  AddField("W2_units");
  AddField("Nu");
  AddField("Nu_units");
  AddField("Epsilon");
  AddField("Gamma");
  AddField("x");
}
//____________________________________________________________________________
DBElDiffXSecTableFields::DBElDiffXSecTableFields(const DBTableFields * fields):
DBTableFields(fields)
{

}
//____________________________________________________________________________
DBElDiffXSecTableFields::~DBElDiffXSecTableFields()
{

}
//____________________________________________________________________________


