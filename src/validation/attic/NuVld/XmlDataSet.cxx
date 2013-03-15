//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 01, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency.
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include "ValidationTools/NuVld/XmlDataSet.h" 

namespace genie {
namespace nuvld {
  
//__________________________________________________________________________
XmlDataSet::XmlDataSet()
{
   _data = new map<string, XmlExperimentMeasurements *>;
}
//__________________________________________________________________________
XmlDataSet::~XmlDataSet()
{

}
//__________________________________________________________________________
void XmlDataSet::Add(string str_unique_id, XmlExperimentMeasurements * mlist)
{
   _data->insert(
         map<string, XmlExperimentMeasurements *>::value_type(
                                                     str_unique_id, mlist)
   );
}
//__________________________________________________________________________
const map<string, XmlExperimentMeasurements *> & XmlDataSet::Get(void) const
{
   return *_data;
}
//__________________________________________________________________________

} // nuvld namespace
} // genie namespace
