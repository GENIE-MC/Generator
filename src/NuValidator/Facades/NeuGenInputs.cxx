//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenInputs

\brief    Encapsulation of NeuGEN's Input Card

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#include <fstream>
#include <sstream>
#include <iomanip>

#include "Facades/NeuGenInputs.h"

using std::endl;
using std::setw;
using std::ios;
using std::setiosflags;
using std::setfill;
using std::ofstream;
using std::ostringstream;

using namespace genie::nuvld::facades;

ClassImp(NeuGenInputs)

//____________________________________________________________________________
namespace genie {
 namespace nuvld {
  namespace facades {
    ostream & operator << (ostream & stream, const NeuGenInputs & conf)
    {
      conf.Print(stream);
      return stream;
    }
  }
 }
}
//____________________________________________________________________________
NeuGenInputs::NeuGenInputs()
{
  this->Init();
}
//____________________________________________________________________________
NeuGenInputs::NeuGenInputs(const NeuGenInputs * inputs)
{
  this->Init();

   // neugen card params

  _nbins             = inputs->_nbins;       
  _xsec_type_code    = inputs->_xsec_type_code;
  _emin              = inputs->_emin;
  _emax              = inputs->_emax;
  _e                 = inputs->_e;
  _plot_var_code     = inputs->_plot_var_code;
  _flux_id_code      = inputs->_flux_id_code;
  _plot_range_code   = inputs->_plot_range_code;
  _plot_var_min      = inputs->_plot_var_min;
  _plot_var_max      = inputs->_plot_var_max;
  _nu_type_code      = inputs->_nu_type_code;
  _wcurrent_code     = inputs->_wcurrent_code;
  _target_code       = inputs->_target_code;
  _final_state_code  = inputs->_final_state_code;
  _init_state_code   = inputs->_init_state_code;
  _process_mask_code = inputs->_process_mask_code;
  _cut_var_code      = inputs->_cut_var_code;
  _cut_var_min       = inputs->_cut_var_min;
  _cut_var_max       = inputs->_cut_var_max;
  _qel_sum           = inputs->_qel_sum;
  _res_sum           = inputs->_res_sum;
  _dis_sum           = inputs->_dis_sum;

  _qel_bit_in_mask   = inputs->_qel_bit_in_mask;
  _res_bit_in_mask   = inputs->_res_bit_in_mask;
  _dis_bit_in_mask   = inputs->_dis_bit_in_mask;

   // aux

  _fin_p             = inputs->_fin_p;
  _fin_n             = inputs->_fin_n;       
  _fin_pi_plus       = inputs->_fin_pi_plus;   
  _fin_pi_0          = inputs->_fin_pi_0;   
  _fin_pi_minus      = inputs->_fin_pi_minus;  
  _xsec_type_str     = inputs->_xsec_type_str;
  _plot_var_str      = inputs->_plot_var_str;
  _flux_id_str       = inputs->_flux_id_str;
  _plot_range_str    = inputs->_plot_range_str;
  _nu_type_str       = inputs->_nu_type_str;
  _wcurrent_str      = inputs->_wcurrent_str;
  _target_str        = inputs->_target_str;
  _final_state_str   = inputs->_final_state_str;
  _init_state_str    = inputs->_init_state_str;
  _cut_var_str       = inputs->_cut_var_str;
}
//____________________________________________________________________________
NeuGenInputs::~NeuGenInputs()
{

}
//____________________________________________________________________________
void NeuGenInputs::WriteNeuGenInputCard(const char * filename) const
{
// This methods writes out the object state in the form of 'data cards' that
// NeuGEN understands.
//  
  ofstream data_card(filename);

  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _nbins
            << " \\\\ nbins:  number of points in plot "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _xsec_type_code
            << " \\\\ xsec type: 1=total, 2=differential "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _emin
            << " \\\\ xsec type=1 - Emin "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _emax
            << " \\\\             - Emax "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _e
            << " \\\\ xsec type=2 - E "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _plot_var_code
            << " \\\\ plot variable "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _flux_id_code
            << " \\\\ flux id (1=ANL, 2=GGM, 3=BNL, 4=BEBC) "
            << endl;;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _plot_range_code
            << " \\\\ plot  range (1=auto, 2=custom) "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _plot_var_min
            << " \\\\ plot range - min "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _plot_var_max
            << " \\\\            - max "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _nu_type_code
            << " \\\\ neutrino (nue/bar=5/6, numu/bar=7/8, nutau/bar=9/10 "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _wcurrent_code
            << " \\\\ weak current (1=CC, 2=NC, 3=BOTH) "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _target_code
            << " \\\\ target: nucleus / particle code, -1 for isoscalar "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _final_state_code
            << " \\\\ final state - in the form pn+-0 "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _init_state_code
            << " \\\\ initial state - (1=v-p, 2=v-n, 3=vbar-p, 4=vbar-n "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _process_mask_code
            << " \\\\ process mask: bits for qel, res, dis, coh "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _cut_var_code
            << " \\\\ cuts variable (0=none, 1=|q^2|, 2=W, 3=x, 4=y "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _cut_var_min
            << " \\\\ cut variable - min "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _cut_var_max
            << " \\\\             - max "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _qel_sum 
            << " \\\\ qelsum: >0 means add qel channel "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _res_sum 
            << " \\\\ ressum: >0 means add all res channels "
            << endl;
  data_card << setiosflags(ios::left) << setfill(' ') << setw(9)
            << _dis_sum 
            << " \\\\ dissum: >0 means add all dis "
            << endl;
}
//____________________________________________________________________________
NGInteraction NeuGenInputs::GetInteraction(void) const
{
  NGFlavor_t     f = NGFlavor::GetFromCode(_nu_type_code);
  NGNucleus_t    n = e_free;
  NGCcNc_t       c = NGCcNc::GetFromCode(_wcurrent_code);
  NGInitState_t  i = NGInitState::GetInitStateFromCode(_init_state_code);

  NGInteraction inter(f, n, c, i);

  return inter;
}
//____________________________________________________________________________
NGFinalState NeuGenInputs::GetFinalState(void) const
{
  NGFinalState state;

  state.SetFinalState(_fin_p, _fin_n, _fin_pi_plus, _fin_pi_minus, _fin_pi_0);

  return state;
}
//____________________________________________________________________________
NeuGenCuts NeuGenInputs::GetCuts(void) const
{
  NGKineVar_t kvid = NGKineVar::GetKineVarFromCode(_cut_var_code);

  bool sumQel = (_qel_sum == 1);
  bool sumRes = (_res_sum == 1);
  bool sumDis = (_dis_sum == 1);

  NeuGenCuts cuts(kvid, _cut_var_min, _cut_var_max, 
                                  _process_mask_code, sumQel, sumRes, sumDis);

  return cuts;
}
//____________________________________________________________________________
void NeuGenInputs::SetNBins(int nbins)
{
  _nbins = nbins;
}
//____________________________________________________________________________
void NeuGenInputs::SetQelBitInMask(bool on)
{
  ( on ) ? _qel_bit_in_mask = 1 : _qel_bit_in_mask = 0;

  this->ComputeProcessMask();
}
//____________________________________________________________________________
void NeuGenInputs::SetResBitInMask(bool on)
{
  ( on ) ? _res_bit_in_mask = 1 : _res_bit_in_mask = 0;

  this->ComputeProcessMask();
}
//____________________________________________________________________________
void NeuGenInputs::SetDisBitInMask(bool on)
{
  ( on ) ? _dis_bit_in_mask = 1 : _dis_bit_in_mask = 0;

  this->ComputeProcessMask();
}
//____________________________________________________________________________
void NeuGenInputs::SetXSecType(string xsec_type)
{
  _xsec_type_str  = xsec_type;
  _xsec_type_code = this->NeuGenXSecTypeCode(xsec_type);
}
//____________________________________________________________________________
void NeuGenInputs::SetEmin(float e_min)
{
  _emin = e_min;
}
//____________________________________________________________________________
void NeuGenInputs::SetEmax(float e_max)
{
  _emax = e_max;
}
//____________________________________________________________________________
void NeuGenInputs::SetE( float e)
{
  _e = e;
}
//____________________________________________________________________________
void NeuGenInputs::SetPlotVar(string plot_variable)
{
  _plot_var_str  = plot_variable;
  _plot_var_code = this->NeuGenVariableCode(plot_variable);
}
//____________________________________________________________________________
void NeuGenInputs::SetFlux(string flux)
{
  _flux_id_str  = flux;
  _flux_id_code = this->NeuGenFluxCode(flux);
}
//____________________________________________________________________________
void NeuGenInputs::SetRange(string range)
{
  _plot_range_str  = range;
  _plot_range_code = this->NeuGenPlotRangeCode(range);
}
//____________________________________________________________________________
void NeuGenInputs::SetPlotVarMin(float var_min)
{
  _plot_var_min = var_min;
}
//____________________________________________________________________________
void NeuGenInputs::SetPlotVarMax(float var_max)
{
  _plot_var_max = var_max;
}
//____________________________________________________________________________
void NeuGenInputs::SetNeutrino(string neutrino)
{
  _nu_type_str  = neutrino;
  _nu_type_code = this->NeuGenNeutrinoCode(neutrino);
}
//____________________________________________________________________________
void NeuGenInputs::SetWkCurrent(string wcurrent)
{
  _wcurrent_str  = wcurrent;
  _wcurrent_code = this->NeuGenWkCurrentCode(wcurrent);
}
//____________________________________________________________________________
void NeuGenInputs::SetTarget(string /*target*/)
{
  _target_str  = ""; // unused
  _target_code = 0;  // unused
}
//____________________________________________________________________________
void NeuGenInputs::SetCutVar(string cut_variable)
{
  _cut_var_str  = cut_variable;
  _cut_var_code = this->NeuGenVariableCode(cut_variable);
}
//____________________________________________________________________________
void NeuGenInputs::SetCutVarMin(float var_min)
{
  _cut_var_min = var_min;
}
//____________________________________________________________________________
void NeuGenInputs::SetCutVarMax(float var_max)
{
  _cut_var_max = var_max;
}
//____________________________________________________________________________
void NeuGenInputs::SetInclusive(bool on)
{
  _inclusive = on;
  
  this->SetQelSum(on);
  this->SetResSum(on);
  this->SetDisSum(on);
}
//____________________________________________________________________________
void NeuGenInputs::SetQelSum(bool on)
{
  _qel_sum  = this->Bool2Int(on);
}
//____________________________________________________________________________
void NeuGenInputs::SetResSum(bool on)
{
  _res_sum  = this->Bool2Int(on);
}
//____________________________________________________________________________
void NeuGenInputs::SetDisSum(bool on)
{
  _dis_sum  = this->Bool2Int(on);
}
//____________________________________________________________________________
void NeuGenInputs::SetFinalState(string fin_state)
{
  _final_state_str  = fin_state;
  _final_state_code = this->NeuGenFinalStateCode(fin_state);
}
//____________________________________________________________________________
void NeuGenInputs::SetInitialState(string init_state)
{
  _init_state_str  = init_state;
  _init_state_code = this->NeuGenInitialStateCode(init_state);
}  
//____________________________________________________________________________
int NeuGenInputs::NeuGenXSecTypeCode(string xsec_type)
{
  if (xsec_type.find("differential") != string::npos) return 2;
  else return 1;
}
//____________________________________________________________________________
int NeuGenInputs::NeuGenFluxCode(string flux)
{
  if      (flux.find("ANL")  != string::npos) return 1;
  else if (flux.find("GGM")  != string::npos) return 2;
  else if (flux.find("BNL")  != string::npos) return 3;
  else if (flux.find("BEBC") != string::npos) return 4;
  else return 0;  
}
//____________________________________________________________________________
int NeuGenInputs::NeuGenPlotRangeCode(string range)
{
  if (range.find("custom") != string::npos) return 2;
  else return 1;
}
//____________________________________________________________________________
int NeuGenInputs::NeuGenNeutrinoCode(string neutrino)
{
  if      (neutrino.find("nu_e")       != string::npos) return  5;
  else if (neutrino.find("nu_e_bar")   != string::npos) return  6;
  else if (neutrino.find("nu_mu")      != string::npos) return  7;
  else if (neutrino.find("nu_mu_bar")  != string::npos) return  8;
  else if (neutrino.find("nu_tau")     != string::npos) return  9;
  else if (neutrino.find("nu_tau_bar") != string::npos) return 10;
  else                                                  return  0;
}
//____________________________________________________________________________
int NeuGenInputs::NeuGenWkCurrentCode(string wcurrent)
{
  if (wcurrent.find("+") != string::npos)       return 3;
  else if (wcurrent.find("NC") != string::npos) return 2;
  else                                          return 1;
}
//____________________________________________________________________________
string NeuGenInputs::NeuGenFinalStateCode(string fin_state)
{
// build neugen's final state - in the form (pn+-0)

  _fin_p        = 0;   // init
  _fin_n        = 0;
  _fin_pi_plus  = 0;
  _fin_pi_0     = 0;
  _fin_pi_minus = 0;

  _fin_p        = this->NMatches(fin_state, "p ");
  _fin_n        = this->NMatches(fin_state, "n ");
  _fin_pi_plus  = this->NMatches(fin_state, "pi(+)");
  _fin_pi_0     = this->NMatches(fin_state, "pi(0)");
  _fin_pi_minus = this->NMatches(fin_state, "pi(-)");

  ostringstream code;

  code << _fin_p << _fin_n << _fin_pi_plus << _fin_pi_minus << _fin_pi_0;

  return code.str();
}
//____________________________________________________________________________
int NeuGenInputs::NeuGenInitialStateCode(string init_state)
{
  if      (init_state.find("nu + p")     != string::npos) return 1;
  else if (init_state.find("nu + n")     != string::npos) return 2;
  else if (init_state.find("nu_bar + p") != string::npos) return 3;
  else if (init_state.find("nu_bar + n") != string::npos) return 4;
  else if (init_state.find("nu + N")     != string::npos) return 5; // ????
  else if (init_state.find("nu_bar + N") != string::npos) return 6; // ????
  else                                              return 0;
}
//____________________________________________________________________________
int NeuGenInputs::NeuGenVariableCode(string cut_variable)
{
  if      (cut_variable.find("none")  != string::npos) return 0;
  else if (cut_variable.find("|q^2|") != string::npos) return 1;
  else if (cut_variable.find("W")     != string::npos) return 2;
  else if (cut_variable.find("x")     != string::npos) return 3;
  else if (cut_variable.find("y")     != string::npos) return 4;
  else                                                 return 0;
}      
//____________________________________________________________________________
void NeuGenInputs::ComputeProcessMask(void)
{
  int qel, dis, res;
  
  (_qel_bit_in_mask == 1) ? qel = 0 : qel = 1;
  (_res_bit_in_mask == 1) ? res = 0 : res = 1;
  (_dis_bit_in_mask == 1) ? dis = 0 : dis = 1;
  
  _process_mask_code = qel + 2 * res + 4 * dis;
}
//____________________________________________________________________________
int NeuGenInputs::NMatches(string input, string pattern)
{
  // if max = 1 then this is ok - do something more generic later
  
  if (input.find(pattern) != string::npos) return 1;
  else return 0;

  /*
  string::size_type pos = 0;
  int n=0;
  
  while( (pos = input.find_first_of(pattern, pos)) != string::npos ) {

      n++;
      input.erase(pos, pattern.length());
  }
  return n;
  */
}
//____________________________________________________________________________
int NeuGenInputs::Bool2Int(bool on)
{
  if (on) return 1;
  else    return 0;
}
//____________________________________________________________________________
void NeuGenInputs::Init(void)
{
  //-- init neugen cards variables
  
  _nbins             = 0;
  _xsec_type_code    = 0;
  _emin              = 0;
  _emax              = 0;
  _e                 = 0;
  _plot_var_code     = 0;
  _flux_id_code      = 0;
  _plot_range_code   = 0;
  _plot_var_min      = 0;
  _plot_var_max      = 0;
  _nu_type_code      = 0;
  _wcurrent_code     = 0;
  _target_code       = 0;
  _final_state_code  = "00000";
  _init_state_code   = 0;
  _process_mask_code = 0;
  _cut_var_code      = 0;
  _cut_var_min       = 0;
  _cut_var_max       = 0;
  _qel_sum           = 0;
  _res_sum           = 0;
  _dis_sum           = 0;

  _qel_bit_in_mask   = 0;
  _res_bit_in_mask   = 0;
  _dis_bit_in_mask   = 0;

  //-- init auxiliary variables

  _fin_p             = 0;
  _fin_n             = 0;
  _fin_pi_plus       = 0;
  _fin_pi_0          = 0;
  _fin_pi_minus      = 0;
  
  _xsec_type_str     = "";
  _plot_var_str      = "";
  _flux_id_str       = "";
  _plot_range_str    = "";
  _nu_type_str       = "";
  _wcurrent_str      = "";
  _target_str        = "";
  _final_state_str   = "";
  _init_state_str    = "";
  _cut_var_str       = "";    
}
//____________________________________________________________________________
void NeuGenInputs::Print(ostream & stream) const
{
  stream << "number of bins =  " << _nbins             << endl;
  stream << "xsec type =       " << _xsec_type_code    << endl;
  stream << "E min =           " << _emin              << endl;
  stream << "E max =           " << _emax              << endl;
  stream << "E =               " << _e                 << endl;
  stream << "plot var =        " << _plot_var_code     << endl;
  stream << "flux id =         " << _flux_id_code      << endl;
  stream << "plot range =      " << _plot_range_code   << endl;
  stream << "plot var - min =  " << _plot_var_min      << endl;
  stream << "plot var - max =  " << _plot_var_max      << endl;
  stream << "neutrino type =   " << _nu_type_code      << endl;
  stream << "weak current =    " << _wcurrent_code     << endl;
  stream << "target =          " << _target_code       << endl;
  stream << "final state =     " << _final_state_code  << endl;
  stream << "initial state =   " << _init_state_code   << endl;
  stream << "process mask =    " << _process_mask_code << endl;
  stream << "qel bit in mask = " << _qel_bit_in_mask   << endl;
  stream << "res bit in mask = " << _res_bit_in_mask   << endl;
  stream << "dis bit in mask = " << _dis_bit_in_mask   << endl;
  stream << "cut variable =    " << _cut_var_code      << endl;
  stream << "cut var - min =   " << _cut_var_min       << endl;
  stream << "cut var - max =   " << _cut_var_max       << endl;
  stream << "qel sum =         " << _qel_sum           << endl;
  stream << "res sum =         " << _res_sum           << endl;
  stream << "dis sum =         " << _dis_sum           << endl;
}
//____________________________________________________________________________

