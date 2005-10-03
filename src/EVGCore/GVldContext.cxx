//____________________________________________________________________________
/*!

\class   genie::GVldContext

\brief   Validity Context for an Event Generator

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 20, 2004

*/
//____________________________________________________________________________

#include <algorithm>

#include "EVGCore/GVldContext.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"

using std::count;

using namespace genie;

//___________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const GVldContext & vldc)
 {
   vldc.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
GVldContext::GVldContext()
{
  Init();
}
//___________________________________________________________________________
GVldContext::GVldContext(const GVldContext & validity_context)
{
  Init();
}
//___________________________________________________________________________
GVldContext::~GVldContext()
{
  if(fCurr)    delete fCurr;
  if(fProbes)  delete fProbes;
  if(fTargets) delete fTargets;
}
//___________________________________________________________________________
bool GVldContext::IsValid(const Interaction * interaction) const
{
  SLOG("VldContext", pINFO)
    << "Examining EvGen's VldContext against an interaction to be generated";

  const ProcessInfo & proc_info = interaction->GetProcessInfo();

  ScatteringType_t  sc_type_id = proc_info.ScatteringTypeId();
  InteractionType_t in_type_id = proc_info.InteractionTypeId();

  int matches = count(fCurr->begin(), fCurr->end(), in_type_id);

  if(fProc == sc_type_id && matches > 0) return true;
  else return false;
}
//___________________________________________________________________________
bool GVldContext::Clashes(const GVldContext & vld_context) const
{
  return true;
}
//___________________________________________________________________________
void GVldContext::Decode(string encoded_vld_context)
{
//example:
//  PROC:RES;CURR:CC,NC;PROBE:nue,nuebar,numu,numubar;TARGET:all;ENERGY:0-100

  vector<string> fields = utils::str::Split(encoded_vld_context, ";");

  vector<string>::const_iterator field_iter;

  for(field_iter = fields.begin(); field_iter != fields.end(); ++field_iter){

     SLOG("VldContext", pINFO) << " ************ " << *field_iter;

     vector<string> curr_field = utils::str::Split(*field_iter, ":");

     assert(curr_field.size() == 2);

     string name   = curr_field[0];
     string values = curr_field[1];

     //-- Make lowercase/uppercase irrelevant

     name = utils::str::ToUpper(name);

     //-- send the string to an appropriate decoder

     if      (name.find("PROC")   != string::npos) DecodePROC   (values);
     else if (name.find("CURR")   != string::npos) DecodeCURR   (values);
     else if (name.find("PROBE")  != string::npos) DecodePROBE  (values);
     else if (name.find("TARGET") != string::npos) DecodeTARGET (values);
     else if (name.find("ENERGY") != string::npos) DecodeENERGY (values);
     else {
         SLOG("VldContext", pWARN)
                      << "**** Unknown field named: "
                                        << name << " in encoded vld-context";
     }
  }
}
//___________________________________________________________________________
void GVldContext::DecodePROC(string encoded_proc)
{
// Decode process: QEL,DIS,RES,COH...
//
  SLOG("VldContext", pDEBUG) << "Decoding PROC: " << encoded_proc;

  ScatteringType_t type = ScatteringType::FromString(encoded_proc);

  assert( type != kScNull );

  fProc = type;
}
//___________________________________________________________________________
void GVldContext::DecodeCURR(string encoded_curr)
{
// Decode current: NC, CC, E/M

  SLOG("VldContext", pDEBUG) << "Decoding CURR: " << encoded_curr;

  vector<string> curr = utils::str::Split(encoded_curr, ",");

  vector<string>::const_iterator curr_iter;

  if(fCurr) delete fCurr;

  fCurr = new vector<InteractionType_t>;

  for(curr_iter = curr.begin(); curr_iter != curr.end(); ++curr_iter) {

     InteractionType_t type = InteractionType::FromString(*curr_iter);

     assert( type != kIntNull );

     fCurr->push_back(type);
  }
}
//___________________________________________________________________________
void GVldContext::DecodePROBE(string encoded_probe)
{
  SLOG("VldContext", pDEBUG) << "Decoding PROBE: " << encoded_probe;
}
//___________________________________________________________________________
void GVldContext::DecodeTARGET(string encoded_target)
{
  SLOG("VldContext", pDEBUG) << "Decoding TARGET: " << encoded_target;
}
//___________________________________________________________________________
void GVldContext::DecodeENERGY(string encoded_energy)
{
  SLOG("VldContext", pDEBUG) << "Decoding ENERGY: " << encoded_energy;

  vector<string> energy = utils::str::Split(encoded_energy, "-");

  assert ( energy.size() == 2 );

  fEmin = atof( energy[0].c_str() );
  fEmax = atof( energy[1].c_str() );
}
//___________________________________________________________________________
void GVldContext::Init(void)
{
  fProc    =  kScNull;
  fCurr    =  0;
  fEmin    = -1.0;
  fEmax    = -1.0;
  fProbes  =  0;
  fTargets =  0;
}
//___________________________________________________________________________
void GVldContext::Print(ostream & stream) const
{
  stream << "\n Scattering-type:...." << fProc;

  vector<int>::const_iterator int_iter;

  stream << "\n Interaction-types:..";
  for(int_iter = fProbes->begin();
           int_iter != fProbes->end(); ++int_iter)
                         stream << InteractionType::AsString(
                                      (InteractionType_t) *int_iter) << "  ";


  stream << "\n Energy range:......." << "[" << fEmin << ", " << fEmax << "]";

  stream << "\n";
}
//___________________________________________________________________________
