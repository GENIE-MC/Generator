//_____________________________________________________________________________
/*!

\class    genie::nuvld::XmlDataSet

\brief    Holds the contents of a parsed NuValidator XML data file

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include "XmlDataSet.h" 

namespace genie {
namespace nuvld {
  
//__________________________________________________________________________
XmlDataSet::XmlDataSet()
{
   _data = new map<string, ExperimentMeasurements *>;
}
//__________________________________________________________________________
XmlDataSet::~XmlDataSet()
{

}
//__________________________________________________________________________
void XmlDataSet::Add(string str_unique_id, ExperimentMeasurements * mlist)
{
   _data->insert(
         map<string, ExperimentMeasurements *>::value_type(
                                                     str_unique_id, mlist)
   );
}
//__________________________________________________________________________
const map<string, ExperimentMeasurements *> & XmlDataSet::Get(void) const
{
   return *_data;
}
//__________________________________________________________________________

} // nuvld namespace
} // genie namespace
