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
#include "Facades/NGSF.h"
#include "Facades/NGPlotType.h"

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

  NGInteraction GetInteraction (void) const;
  NGFinalState  GetFinalState  (void) const;
  NeuGenCuts    GetCuts        (void) const;

  bool   IncludeQel       (void) const { return fIncludeQel;       }
  bool   IncludeRes       (void) const { return fIncludeRes;       }
  bool   IncludeDis       (void) const { return fIncludeDis;       }
  bool   Inclusive        (void) const { return fInclusive;        }
  bool   SFRawDis         (void) const { return fSFRawDis;         }
  int    NBins            (void) const { return fNBins;            }
  int    PlotVarCode      (void) const { return fPlotVarCode;      }
  int    FluxIdCode       (void) const { return fFluxCode;         }
  int    PlotRangeCode    (void) const { return fPlotRangeCode;    }
  int    ProbeCode        (void) const { return fProbeCode;        }
  int    WkCurrentCode    (void) const { return fWkCurr;           }
  int    TargetCode       (void) const { return fTgtCode;          }
  int    InitStateCode    (void) const { return fInitStateCode;    }
  int    CutVarCode       (void) const { return fCutVarCode;       }
  int    A                (void) const { return fA;                }
  int    SFCode           (void) const { return fSFCode;           }
  float  EnergyMin        (void) const { return fEmin;             }
  float  EnergyMax        (void) const { return fEmax;             }
  float  Energy           (void) const { return fE;                }
  float  PlotVarMin       (void) const { return fPlotVarMin;       }
  float  PlotVarMax       (void) const { return fPlotVarMax;       }
  float  CutVarMin        (void) const { return fCutVarMin;        }
  float  CutVarMax        (void) const { return fCutVarMax;        }
  float  SFFixedVar       (void) const { return fSFFixedVar;       }
  string FinalStateCode   (void) const { return fFinalStateCode;   }
  string PlotTypeString   (void) const { return fPlotTypeStr;      }
  string PlotVarString    (void) const { return fPlotVarStr;       }
  string FluxIdString     (void) const { return fFluxStr;          }
  string PlotRangeString  (void) const { return fPlotRangeStr;     }
  string ProbeString      (void) const { return fProbeStr;         }
  string WkCurrentString  (void) const { return fWkCurrStr;        }
  string TargetString     (void) const { return fTgtStr;           }
  string SFString         (void) const { return fSFStr;            }
  string FinalStateString (void) const { return fFinalStateStr;    }
  string InitStateString  (void) const { return fInitStateStr;     }
  string CutVarString     (void) const { return fCutVarStr;        }

  int          SFRawDisCode (void) const;
  NGSF_t       SF           (void) const;
  NGKineVar_t  PlotVar      (void) const;
  NGPlotType_t PlotType     (void) const { return fPlotType; }
  
  void SetNBins         ( int    nbins      );
  void SetPlotType      ( string plot_type  );
  void SetEmin          ( float  e_min      );
  void SetEmax          ( float  e_max      );
  void SetE             ( float  e          );
  void SetPlotVar       ( string variable   );
  void SetFlux          ( string flux       );
  void SetRange         ( string range      );
  void SetPlotVarMin    ( float  var_min    );
  void SetPlotVarMax    ( float  var_max    );
  void SetNeutrino      ( string neutrino   );
  void SetWkCurrent     ( string wcurrent   );
  void SetTarget        ( string target     );
  void SetSF            ( string sf         );
  void SetA             ( int A             );
  void SetFinalState    ( string fin_state  );
  void SetInitialState  ( string init_state );
  void SetCutVar        ( string variable   );
  void SetCutVarMin     ( float  var_min    );
  void SetCutVarMax     ( float  var_max    );
  void SetSFFixedVar    ( float  var        );
  void SetIncludeQel    ( bool   on         );
  void SetIncludeRes    ( bool   on         );
  void SetIncludeDis    ( bool   on         );
  void SetInclusive     ( bool   on         );
  void SetSFRawDis      ( bool   on         );

private:

  void         Init                   (void);
  NGPlotType_t NeuGenPlotType         (string plot_type);
  int          NeuGenFluxCode         (string flux);
  int          NeuGenPlotRangeCode    (string range);
  int          NeuGenNeutrinoCode     (string neutrino);
  int          NeuGenWkCurrentCode    (string wcurrent);
  string       NeuGenFinalStateCode   (string fin_state);
  int          NeuGenInitialStateCode (string init_state);
  int          NeuGenVariableCode     (string plot_variable);
  int          NeuGenSFCode           (string sf);
  int          Bool2Int               (bool on);
  int          NMatches               (string input, string pattern);

  // NeuGEN input variables - in the order they appear in the NeuGEN input card:
  int          fNBins;
  float        fEmin;
  float        fEmax;
  float        fE;
  int          fPlotVarCode;
  int          fFluxCode;
  int          fPlotRangeCode;
  float        fPlotVarMin;
  float        fPlotVarMax;
  int          fProbeCode;
  int          fWkCurr;
  int          fTgtCode;
  int          fA;
  string       fFinalStateCode;
  int          fInitStateCode;
  int          fCutVarCode;
  float        fCutVarMin;
  float        fCutVarMax;
  float        fSFFixedVar;
  int          fSFCode;
  bool         fIncludeQel;
  bool         fIncludeRes;
  bool         fIncludeDis;
  bool         fInclusive;
  bool         fSFRawDis;
  NGPlotType_t fPlotType;
  
  // auxiliary variables:  
  int          fFsP;        ///< final state: number of p's
  int          fFsN;        ///<              number of n's
  int          fFsPiPlus;   ///<              number of pi+'s
  int          fFsPi0;      ///<              number of pi0's
  int          fFsPiMinus;  ///<              number of pi-'s
  string       fPlotTypeStr;
  string       fPlotVarStr;
  string       fFluxStr;
  string       fPlotRangeStr;
  string       fProbeStr;
  string       fWkCurrStr; 
  string       fTgtStr;
  string       fSFStr;
  string       fFinalStateStr;
  string       fInitStateStr;
  string       fCutVarStr;

  ClassDef(NeuGenInputs, 1)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif
