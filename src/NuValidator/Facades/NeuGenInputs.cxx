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

#include "Messenger/Messenger.h"
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

  fNBins          = inputs->fNBins;       
  fPlotType       = inputs->fPlotType;
  fEmin           = inputs->fEmin;
  fEmax           = inputs->fEmax;
  fE              = inputs->fE;
  fPlotVarCode    = inputs->fPlotVarCode;
  fFluxCode       = inputs->fFluxCode;
  fPlotRangeCode  = inputs->fPlotRangeCode;
  fPlotVarMin     = inputs->fPlotVarMin;
  fPlotVarMax     = inputs->fPlotVarMax;
  fProbeCode      = inputs->fProbeCode;
  fWkCurr         = inputs->fWkCurr;
  fTgtCode        = inputs->fTgtCode;
  fA              = inputs->fA;
  fFinalStateCode = inputs->fFinalStateCode;
  fInitStateCode  = inputs->fInitStateCode;
  fCutVarCode     = inputs->fCutVarCode;
  fCutVarMin      = inputs->fCutVarMin;
  fCutVarMax      = inputs->fCutVarMax;
  fIncludeQel     = inputs->fIncludeQel;
  fIncludeRes     = inputs->fIncludeRes;
  fIncludeDis     = inputs->fIncludeDis;
  fInclusive      = inputs->fInclusive;
  fSFRawDis       = inputs->fSFRawDis;
  fSFCode         = inputs->fSFCode;
  fSFFixedVar     = inputs->fSFFixedVar;
  
   // aux

  fFsP            = inputs->fFsP;
  fFsN            = inputs->fFsN;       
  fFsPiPlus       = inputs->fFsPiPlus;   
  fFsPi0          = inputs->fFsPi0;   
  fFsPiMinus      = inputs->fFsPiMinus;  
  fPlotTypeStr    = inputs->fPlotTypeStr;
  fPlotVarStr     = inputs->fPlotVarStr;
  fFluxStr        = inputs->fFluxStr;
  fPlotRangeStr   = inputs->fPlotRangeStr;
  fProbeStr       = inputs->fProbeStr;
  fWkCurrStr      = inputs->fWkCurrStr;
  fTgtStr         = inputs->fTgtStr;
  fSFStr          = inputs->fSFStr;
  fFinalStateStr  = inputs->fFinalStateStr;
  fInitStateStr   = inputs->fInitStateStr;
  fCutVarStr      = inputs->fCutVarStr;
}
//____________________________________________________________________________
NeuGenInputs::~NeuGenInputs()
{

}
//____________________________________________________________________________
NGInteraction NeuGenInputs::GetInteraction(void) const
{
  NGFlavor_t     f = NGFlavor::GetFromCode(fProbeCode);
  NGNucleus_t    n = e_free;
  NGCcNc_t       c = NGCcNc::GetFromCode(fWkCurr);
  NGInitState_t  i = NGInitState::GetInitStateFromCode(fInitStateCode);

  NGInteraction inter(f, n, c, i);

  return inter;
}
//____________________________________________________________________________
NGFinalState NeuGenInputs::GetFinalState(void) const
{
  NGFinalState state;

  state.SetFinalState(fFsP, fFsN, fFsPiPlus, fFsPiMinus, fFsPi0);

  return state;
}
//____________________________________________________________________________
NeuGenCuts NeuGenInputs::GetCuts(void) const
{
  NGKineVar_t kvid = NGKineVar::GetKineVarFromCode(fCutVarCode);

  NeuGenCuts cuts(kvid, fCutVarMin, fCutVarMax,
                           fInclusive, fIncludeQel, fIncludeRes, fIncludeDis);

  LOG("NeuGen", pINFO) << cuts;
  
  return cuts;
}
//____________________________________________________________________________
void NeuGenInputs::SetNBins(int nbins)
{
  fNBins = nbins;
}
//____________________________________________________________________________
void NeuGenInputs::SetPlotType(string plot_type)
{
  fPlotTypeStr = plot_type;
  fPlotType    = this->NeuGenPlotType(plot_type);
}
//____________________________________________________________________________
void NeuGenInputs::SetEmin(float e_min)
{
  fEmin = e_min;
}
//____________________________________________________________________________
void NeuGenInputs::SetEmax(float e_max)
{
  fEmax = e_max;
}
//____________________________________________________________________________
void NeuGenInputs::SetE( float e)
{
  fE = e;
}
//____________________________________________________________________________
void NeuGenInputs::SetPlotVar(string plot_variable)
{
  fPlotVarStr  = plot_variable;
  fPlotVarCode = this->NeuGenVariableCode(plot_variable);
}
//____________________________________________________________________________
void NeuGenInputs::SetFlux(string flux)
{
  fFluxStr  = flux;
  fFluxCode = this->NeuGenFluxCode(flux);
}
//____________________________________________________________________________
void NeuGenInputs::SetRange(string range)
{
  fPlotRangeStr  = range;
  fPlotRangeCode = this->NeuGenPlotRangeCode(range);
}
//____________________________________________________________________________
void NeuGenInputs::SetPlotVarMin(float var_min)
{
  fPlotVarMin = var_min;
}
//____________________________________________________________________________
void NeuGenInputs::SetPlotVarMax(float var_max)
{
  fPlotVarMax = var_max;
}
//____________________________________________________________________________
void NeuGenInputs::SetNeutrino(string neutrino)
{
  fProbeStr  = neutrino;
  fProbeCode = this->NeuGenNeutrinoCode(neutrino);
}
//____________________________________________________________________________
void NeuGenInputs::SetWkCurrent(string wcurrent)
{
  fWkCurrStr = wcurrent;
  fWkCurr    = this->NeuGenWkCurrentCode(wcurrent);
}
//____________________________________________________________________________
void NeuGenInputs::SetTarget(string /*target*/)
{
  fTgtStr  = ""; // unused
  fTgtCode = 0;  // unused
}
//____________________________________________________________________________
void NeuGenInputs::SetA(int A)
{
  fA = A;
}
//____________________________________________________________________________
void NeuGenInputs::SetCutVar(string cut_variable)
{
  fCutVarStr  = cut_variable;
  fCutVarCode = this->NeuGenVariableCode(cut_variable);
}
//____________________________________________________________________________
void NeuGenInputs::SetCutVarMin(float var_min)
{
  fCutVarMin = var_min;
}
//____________________________________________________________________________
void NeuGenInputs::SetCutVarMax(float var_max)
{
  fCutVarMax = var_max;
}
//____________________________________________________________________________
void NeuGenInputs::SetSFFixedVar(float var)
{
  fSFFixedVar = var;
}
//____________________________________________________________________________
void NeuGenInputs::SetInclusive(bool on)
{
  fInclusive = on;
}
//____________________________________________________________________________
void NeuGenInputs::SetIncludeQel(bool on)
{
  fIncludeQel = on;
}
//____________________________________________________________________________
void NeuGenInputs::SetIncludeRes(bool on)
{
  fIncludeRes = on;
}
//____________________________________________________________________________
void NeuGenInputs::SetIncludeDis(bool on)
{
  fIncludeDis = on;
}
//____________________________________________________________________________
void NeuGenInputs::SetSFRawDis(bool on)
{
  fSFRawDis = on;
}
//____________________________________________________________________________
void NeuGenInputs::SetFinalState(string fin_state)
{
  fFinalStateStr  = fin_state;
  fFinalStateCode = this->NeuGenFinalStateCode(fin_state);
}
//____________________________________________________________________________
void NeuGenInputs::SetInitialState(string init_state)
{
  fInitStateStr  = init_state;
  fInitStateCode = this->NeuGenInitialStateCode(init_state);
}  
//____________________________________________________________________________
void NeuGenInputs::SetSF(string sf)
{
  fSFStr  = sf;
  fSFCode = this->NeuGenSFCode(sf);
}
//____________________________________________________________________________
NGPlotType_t NeuGenInputs::NeuGenPlotType(string plot_type)
{
  if      (plot_type.find("diff")  != string::npos) return e_DiffXSec;
  else if (plot_type.find("struc") != string::npos) return e_SF; 
  else                                              return e_XSec;
}
//____________________________________________________________________________
int NeuGenInputs::SFRawDisCode(void) const
{
  if(fSFRawDis) return 2;
  else          return 1;
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

  fFsP       = 0;   // init
  fFsN       = 0;
  fFsPiPlus  = 0;
  fFsPi0     = 0;
  fFsPiMinus = 0;

  fFsP       = this->NMatches(fin_state, "p ");
  fFsN       = this->NMatches(fin_state, "n ");
  fFsPiPlus  = this->NMatches(fin_state, "pi(+)");
  fFsPi0     = this->NMatches(fin_state, "pi(0)");
  fFsPiMinus = this->NMatches(fin_state, "pi(-)");

  ostringstream code;

  code << fFsP << fFsN << fFsPiPlus << fFsPiMinus << fFsPi0;

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
int NeuGenInputs::NeuGenSFCode(string sf)
{
  if      (sf.find("xF3")  != string::npos) return 6; 
  else if (sf.find("F1")   != string::npos) return 0;
  else if (sf.find("F2")   != string::npos) return 1; 
  else if (sf.find("F3")   != string::npos) return 2;
  else if (sf.find("F4")   != string::npos) return 3;
  else if (sf.find("F5")   != string::npos) return 4;
  else if (sf.find("F6")   != string::npos) return 5;
  else return 0;
}
//____________________________________________________________________________
NGSF_t NeuGenInputs::SF(void) const
{
  return NGSF::GetSFFromCode(fSFCode);
}
//____________________________________________________________________________
NGKineVar_t NeuGenInputs::PlotVar(void) const
{
  return NGKineVar::GetKineVarFromCode(fPlotVarCode);
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
  
  fNBins            = 0;
  fPlotType         = e_UndefinedPlotType;
  fEmin             = 0;
  fEmax             = 0;
  fE                = 0;
  fPlotVarCode      = 0;
  fFluxCode         = 0;
  fPlotRangeCode    = 0;
  fPlotVarMin       = 0;
  fPlotVarMax       = 0;
  fProbeCode        = 0;
  fWkCurr           = 0;
  fTgtCode          = 0;
  fA                = 1;
  fFinalStateCode   = "00000";
  fInitStateCode    = 0;
  fCutVarCode       = 0;
  fCutVarMin        = 0;
  fCutVarMax        = 0;
  fIncludeQel       = 0;
  fIncludeRes       = 0;
  fIncludeDis       = 0;
  fInclusive        = true;
  fSFRawDis         = false;
  fSFCode           = -1;
  fSFFixedVar       = 0;
  
  //-- init auxiliary variables

  fFsP              = 0;
  fFsN              = 0;
  fFsPiPlus         = 0;
  fFsPi0            = 0;
  fFsPiMinus        = 0;
  
  fPlotTypeStr      = "";
  fPlotVarStr       = "";
  fFluxStr          = "";
  fPlotRangeStr     = "";
  fProbeStr         = "";
  fWkCurrStr        = "";
  fTgtStr           = "";
  fSFStr            = "";
  fFinalStateStr    = "";
  fInitStateStr     = "";
  fCutVarStr        = "";    
}
//____________________________________________________________________________
void NeuGenInputs::Print(ostream & stream) const
{
  stream << "number of bins =  " << fNBins                          << endl;
  stream << "plot type =       " << NGPlotType::AsString(fPlotType) << endl;
  stream << "E min =           " << fEmin             << endl;
  stream << "E max =           " << fEmax             << endl;
  stream << "E =               " << fE                << endl;
  stream << "plot var =        " << fPlotVarCode      << endl;
  stream << "flux id =         " << fFluxCode         << endl;
  stream << "plot range =      " << fPlotRangeCode    << endl;
  stream << "plot var - min =  " << fPlotVarMin       << endl;
  stream << "plot var - max =  " << fPlotVarMax       << endl;
  stream << "probe type =      " << fProbeCode        << endl;
  stream << "weak current =    " << fWkCurr           << endl;
  stream << "target =          " << fTgtCode          << endl;
  stream << "A =               " << fA                << endl;
  stream << "final state =     " << fFinalStateCode   << endl;
  stream << "initial state =   " << fInitStateCode    << endl;
  stream << "cut variable =    " << fCutVarCode       << endl;
  stream << "cut var - min =   " << fCutVarMin        << endl;
  stream << "cut var - max =   " << fCutVarMax        << endl;
  stream << "include qel =     " << fIncludeQel       << endl;
  stream << "include res =     " << fIncludeRes       << endl;
  stream << "include dis =     " << fIncludeDis       << endl;
  stream << "inclusive =       " << fInclusive        << endl;
  stream << "SF raw dis =      " << fSFRawDis         << endl;
  stream << "SF code =         " << fSFCode           << endl;
  stream << "SF fixed var =    " << fSFFixedVar       << endl;
}
//____________________________________________________________________________

