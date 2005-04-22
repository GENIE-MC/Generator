//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenInputs

\brief    Encapsulation of NeuGEN's Input Card

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_INPUTS_H_
#define _NEUGEN_INPUTS_H_

#include <string>
#include <ostream>

#include <TObject.h>

#include "Facades/NeuGenCuts.h"       
#include "Facades/NGInitState.h"       
#include "Facades/NGFinalState.h" 
#include "Facades/NGInteraction.h"

using std::string;
using std::ostream;

namespace genie   {
namespace nuvld   {
namespace facades {

class NeuGenInputs : public TObject
{
public:

  NeuGenInputs();
  NeuGenInputs(const NeuGenInputs * inputs);
  ~NeuGenInputs();

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NeuGenInputs & inputs);

  void WriteNeuGenInputCard(const char * filename) const;

  NGInteraction GetInteraction (void) const;
  NGFinalState  GetFinalState  (void) const;
  NeuGenCuts    GetCuts        (void) const;

  bool   QelSum           (void) const { return _qel_sum == 1;      }
  bool   ResSum           (void) const { return _res_sum == 1;      }
  bool   DisSum           (void) const { return _dis_sum == 1;      }
  bool   Inclusive        (void) const { return _inclusive;         }
  int    NBins            (void) const { return _nbins;             }
  int    XSecTypeCode     (void) const { return _xsec_type_code;    }
  int    PlotVarCode      (void) const { return _plot_var_code;     }
  int    FluxIdCode       (void) const { return _flux_id_code;      }
  int    PlotRangeCode    (void) const { return _plot_range_code;   }
  int    NuTypeCode       (void) const { return _nu_type_code;      }
  int    WkCurrentCode    (void) const { return _wcurrent_code;     }
  int    TargetCode       (void) const { return _target_code;       }
  int    InitStateCode    (void) const { return _init_state_code;   }
  int    CutVarCode       (void) const { return _cut_var_code;      }
  float  EnergyMin        (void) const { return _emin;              }
  float  EnergyMax        (void) const { return _emax;              }
  float  Energy           (void) const { return _e;                 }
  float  PlotVarMin       (void) const { return _plot_var_min;      }
  float  PlotVarMax       (void) const { return _plot_var_max;      }
  float  CutVarMin        (void) const { return _cut_var_min;       }
  float  CutVarMax        (void) const { return _cut_var_max;       }
  string FinalStateCode   (void) const { return _final_state_code;  }
  string XSecTypeString   (void) const { return _xsec_type_str;     }
  string PlotVarString    (void) const { return _plot_var_str;      }
  string FluxIdString     (void) const { return _flux_id_str;       }
  string PlotRangeString  (void) const { return _plot_range_str;    }
  string NuTypeString     (void) const { return _nu_type_str;       }
  string WkCurrentString  (void) const { return _wcurrent_str;      }
  string TargetString     (void) const { return _target_str;        }
  string FinalStateString (void) const { return _final_state_str;   }
  string InitStateString  (void) const { return _init_state_str;    }
  string CutVarString     (void) const { return _cut_var_str;       }

  void SetNBins         ( int    nbins      );
  void SetXSecType      ( string xsec_type  );
  void SetEmin          ( float e_min       );
  void SetEmax          ( float e_max       );
  void SetE             ( float e           );
  void SetPlotVar       ( string variable   );
  void SetFlux          ( string flux       );
  void SetRange         ( string range      );
  void SetPlotVarMin    ( float var_min     );
  void SetPlotVarMax    ( float var_max     );
  void SetNeutrino      ( string neutrino   );
  void SetWkCurrent     ( string wcurrent   );
  void SetTarget        ( string target     );
  void SetFinalState    ( string fin_state  );
  void SetInitialState  ( string init_state );
  void SetCutVar        ( string variable   );
  void SetCutVarMin     ( float var_min     );
  void SetCutVarMax     ( float var_max     );
  void SetQelSum        ( bool on           );
  void SetResSum        ( bool on           );
  void SetDisSum        ( bool on           );
  void SetInclusive     ( bool on           );

private:

  void   Init (void);

  int    NeuGenXSecTypeCode     (string xsec_type);
  int    NeuGenFluxCode         (string flux);
  int    NeuGenPlotRangeCode    (string range);
  int    NeuGenNeutrinoCode     (string neutrino);
  int    NeuGenWkCurrentCode    (string wcurrent);
  string NeuGenFinalStateCode   (string fin_state);
  int    NeuGenInitialStateCode (string init_state);
  int    NeuGenVariableCode     (string plot_variable);
  int    Bool2Int               (bool on);
  int    NMatches               (string input, string pattern);

  // NeuGEN input variables - in the order they appear in the NeuGEN input card:

  int    _nbins;
  int    _xsec_type_code;
  float  _emin;
  float  _emax;
  float  _e;
  int    _plot_var_code;
  int    _flux_id_code;
  int    _plot_range_code;
  float  _plot_var_min;
  float  _plot_var_max;
  int    _nu_type_code;
  int    _wcurrent_code;
  int    _target_code;
  string _final_state_code;
  int    _init_state_code;
  int    _cut_var_code;
  float  _cut_var_min;
  float  _cut_var_max;
  int    _qel_sum;
  int    _res_sum;
  int    _dis_sum;
  bool   _inclusive;
  
  // auxiliary variables:
  
  int    _fin_p;          // final state: number of p's
  int    _fin_n;          //              number of n's
  int    _fin_pi_plus;    //              number of pi+'s
  int    _fin_pi_0;       //              number of pi0's
  int    _fin_pi_minus;   //              number of pi-'s

  string _xsec_type_str;
  string _plot_var_str;
  string _flux_id_str;
  string _plot_range_str;
  string _nu_type_str;
  string _wcurrent_str; 
  string _target_str;   
  string _final_state_str;
  string _init_state_str;
  string _cut_var_str;

  ClassDef(NeuGenInputs, 1)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif
