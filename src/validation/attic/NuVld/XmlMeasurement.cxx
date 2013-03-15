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

#include "ValidationTools/NuVld/XmlMeasurement.h"

namespace genie {
namespace nuvld {
  
//_______________________________________________________________________________
XmlMeasurement::XmlMeasurement()
{
  _data = new vector<XmlRecordBase *>;
} 
//_______________________________________________________________________________
XmlMeasurement::XmlMeasurement(const XmlMeasurement & /*meas*/)
{

}
//_______________________________________________________________________________
void XmlMeasurement::Add(XmlMeasurementHeader * header)
{
  _header = header;
}
//_______________________________________________________________________________
void XmlMeasurement::Add(XmlRecordBase * rec)
{
  _data->push_back(rec);
}
//_______________________________________________________________________________

} // nuvld namespace
} // genie namespace
