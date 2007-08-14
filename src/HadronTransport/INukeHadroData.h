//____________________________________________________________________________
/*!

\class    genie::INukeHadroData

\brief    Singleton class to load & serve hadron x-section splines used by
          GENIE's version of the INTRANUKE cascade MC.

          See $GENIE/src/HadronTransport/Intranuke.h for more details on the
          INTRANUKE cascade MC developed primarity by S.Dytman and H.Gallagher
          continuing older work from R.Edgecock, G.F.Pearce, W.A.Mann, 
          R.Merenyi and others.

          The hadron x-section data used to build the x-section splines stored
          at this singleton are provided & maintained by Steve Dytman.
          See the data files in $GENIE/data/hadron_xsec/ for more details on
          Steve's data sources and applied corrections. 
          In a nutshell:
          The h+N x-sections come mostly from the SAID (Arndt et al.) PWA fit
          while the h+A x-sections come from a combination of Ashery, Carroll 
          data and extrapolations, and INC model results from Mashnik et al.
          for h+Fe56.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.
          Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.

\created  February 01, 2007

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_HADRON_CROSS_SECTIONS_H_
#define _INTRANUKE_HADRON_CROSS_SECTIONS_H_

#include "HadronTransport/INukeHadroFates.h"

namespace genie {

class Spline;

class INukeHadroData
{
public:
  static INukeHadroData * Instance (void);

// Note that, unlike most the rest of GENIE where everything is expressed
// in natural units, all x-section splines included here are evaluated in
// kinetic energies given in MeV and return the x-section value in mbarns

  double XSec (int hpdgc, INukeFateHA_t fate, double ke) const;
  double XSec (int hpdgc, INukeFateHN_t fate, double ke) const;
  double Frac (int hpdgc, INukeFateHA_t fate, double ke) const;
  double Frac (int hpdgc, INukeFateHN_t fate, double ke) const;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // hN mode hadron x-section splines
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const Spline * XSecPN_Tot       (void) const { return fXSecPN_Tot;       }
  const Spline * XSecPN_Elas      (void) const { return fXSecPN_Elas;      }
  const Spline * XSecPN_Reac      (void) const { return fXSecPN_Reac;      }
  const Spline * XSecNN_Tot       (void) const { return fXSecNN_Tot;       }    
  const Spline * XSecNN_Elas      (void) const { return fXSecNN_Elas;      } 
  const Spline * XSecNN_Reac      (void) const { return fXSecNN_Reac;      }
  const Spline * XSecPipN_Tot     (void) const { return fXSecPipN_Tot;     }
  const Spline * XSecPipN_CEx     (void) const { return fXSecPipN_CEx;     }
  const Spline * XSecPipN_Elas    (void) const { return fXSecPipN_Elas;    }
  const Spline * XSecPipN_Reac    (void) const { return fXSecPipN_Reac;    }
  const Spline * XSecPipN_Abs     (void) const { return fXSecPipN_Abs;     }
  const Spline * XSecPimN_Tot     (void) const { return fXSecPimN_Tot;     }
  const Spline * XSecPimN_CEx     (void) const { return fXSecPimN_CEx;     }
  const Spline * XSecPimN_Elas    (void) const { return fXSecPimN_Elas;    }
  const Spline * XSecPimN_Reac    (void) const { return fXSecPimN_Reac;    }
  const Spline * XSecPimN_Abs     (void) const { return fXSecPimN_Abs;     }
  const Spline * XSecPi0N_Tot     (void) const { return fXSecPi0N_Tot;     }
  const Spline * XSecPi0N_CEx     (void) const { return fXSecPi0N_CEx;     }
  const Spline * XSecPi0N_Elas    (void) const { return fXSecPi0N_Elas;    }
  const Spline * XSecPi0N_Reac    (void) const { return fXSecPi0N_Reac;    }
  const Spline * XSecPi0N_Abs     (void) const { return fXSecPi0N_Abs;     }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // hA mode hadron x-section splines
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const Spline * XSecPA_Tot       (void) const { return fXSecPA_Tot;       }
  const Spline * XSecPA_Elas      (void) const { return fXSecPA_Elas;      }
  const Spline * XSecPA_Inel      (void) const { return fXSecPA_Inel;      }
  const Spline * XSecPA_CEx       (void) const { return fXSecPA_CEx;       }
  const Spline * XSecPA_NP        (void) const { return fXSecPA_NP;        }
  const Spline * XSecPA_PP        (void) const { return fXSecPA_PP;        }
  const Spline * XSecPA_NPP       (void) const { return fXSecPA_NPP;       }
  const Spline * XSecPA_NNP       (void) const { return fXSecPA_NNP;       }
  const Spline * XSecPA_NNPPP     (void) const { return fXSecPA_NNPPP;     }
  const Spline * XSecPA_NPip      (void) const { return fXSecPA_NPip;      }
  const Spline * XSecPA_NPipPi0   (void) const { return fXSecPA_NPipPi0;   }
  const Spline * XSecNA_Tot       (void) const { return fXSecNA_Tot;       }
  const Spline * XSecNA_Elas      (void) const { return fXSecNA_Elas;      }
  const Spline * XSecNA_Inel      (void) const { return fXSecNA_Inel;      }
  const Spline * XSecNA_CEx       (void) const { return fXSecNA_CEx;       }
  const Spline * XSecNA_NP        (void) const { return fXSecNA_NP;        }
  const Spline * XSecNA_PP        (void) const { return fXSecNA_PP;        }
  const Spline * XSecNA_NPP       (void) const { return fXSecNA_NPP;       }
  const Spline * XSecNA_NNP       (void) const { return fXSecNA_NNP;       }
  const Spline * XSecNA_NNPPP     (void) const { return fXSecNA_NNPPP;     }
  const Spline * XSecNA_NPip      (void) const { return fXSecNA_NPip;      }
  const Spline * XSecNA_NPipPi0   (void) const { return fXSecNA_NPipPi0;   }
  const Spline * XSecPipA_Tot     (void) const { return fXSecPipA_Tot;     }
  const Spline * XSecPipA_Elas    (void) const { return fXSecPipA_Elas;    }
  const Spline * XSecPipA_Inel    (void) const { return fXSecPipA_Inel;    }
  const Spline * XSecPipA_CEx     (void) const { return fXSecPipA_CEx;     }
  const Spline * XSecPipA_NP      (void) const { return fXSecPipA_NP;      }
  const Spline * XSecPipA_PP      (void) const { return fXSecPipA_PP;      }
  const Spline * XSecPipA_NPP     (void) const { return fXSecPipA_NPP;     }
  const Spline * XSecPipA_NNP     (void) const { return fXSecPipA_NNP;     }
  const Spline * XSecPipA_NNPP    (void) const { return fXSecPipA_NNPP;    }
  const Spline * XSecPipA_NPipPi0 (void) const { return fXSecPipA_NPipPi0; }
  const Spline * XSecPimA_Tot     (void) const { return fXSecPimA_Tot;     }
  const Spline * XSecPimA_Elas    (void) const { return fXSecPimA_Elas;    }
  const Spline * XSecPimA_Inel    (void) const { return fXSecPimA_Inel;    }
  const Spline * XSecPimA_CEx     (void) const { return fXSecPimA_CEx;     }
  const Spline * XSecPimA_NP      (void) const { return fXSecPimA_NP;      }
  const Spline * XSecPimA_PP      (void) const { return fXSecPimA_PP;      }
  const Spline * XSecPimA_NPP     (void) const { return fXSecPimA_NPP;     }
  const Spline * XSecPimA_NNP     (void) const { return fXSecPimA_NNP;     }
  const Spline * XSecPimA_NNPP    (void) const { return fXSecPimA_NNPP;    }
  const Spline * XSecPimA_NPipPi0 (void) const { return fXSecPimA_NPipPi0; }
  const Spline * XSecPi0A_Tot     (void) const { return fXSecPi0A_Tot;     }
  const Spline * XSecPi0A_Elas    (void) const { return fXSecPi0A_Elas;    }
  const Spline * XSecPi0A_Inel    (void) const { return fXSecPi0A_Inel;    }
  const Spline * XSecPi0A_CEx     (void) const { return fXSecPi0A_CEx;     }
  const Spline * XSecPi0A_NP      (void) const { return fXSecPi0A_NP;      }
  const Spline * XSecPi0A_PP      (void) const { return fXSecPi0A_PP;      }
  const Spline * XSecPi0A_NPP     (void) const { return fXSecPi0A_NPP;     }
  const Spline * XSecPi0A_NNP     (void) const { return fXSecPi0A_NNP;     }
  const Spline * XSecPi0A_NNPP    (void) const { return fXSecPi0A_NNPP;    }
  const Spline * XSecPi0A_NPipPi0 (void) const { return fXSecPi0A_NPipPi0; }

  static double fMinKinEnergy;   ///<
  static double fMaxKinEnergyHA; ///<
  static double fMaxKinEnergyHN; ///<

private:
  INukeHadroData();
  INukeHadroData(const INukeHadroData & shx);
 ~INukeHadroData();

  void LoadCrossSections(void); 

  static INukeHadroData * fInstance;

  Spline * fXSecPN_Tot;        ///< N+N x-section splines
  Spline * fXSecPN_Elas;       ///<
  Spline * fXSecPN_Reac;       ///<
  Spline * fXSecNN_Tot;        ///<
  Spline * fXSecNN_Elas;       ///<
  Spline * fXSecNN_Reac;       ///<
  Spline * fXSecPipN_Tot;      ///< pi+N x-section splines
  Spline * fXSecPipN_CEx;      ///<
  Spline * fXSecPipN_Elas;     ///<
  Spline * fXSecPipN_Reac;     ///<
  Spline * fXSecPipN_Abs;      ///<
  Spline * fXSecPimN_Tot;      ///<
  Spline * fXSecPimN_CEx;      ///<	
  Spline * fXSecPimN_Elas;     ///<
  Spline * fXSecPimN_Reac;     ///<
  Spline * fXSecPimN_Abs;      ///<
  Spline * fXSecPi0N_Tot;      ///<
  Spline * fXSecPi0N_CEx;      ///<
  Spline * fXSecPi0N_Elas;     ///<
  Spline * fXSecPi0N_Reac;     ///<
  Spline * fXSecPi0N_Abs;      ///<
  Spline * fXSecPA_Tot;        ///< N+A x-section splines
  Spline * fXSecPA_Elas;       ///<
  Spline * fXSecPA_Inel;       ///<
  Spline * fXSecPA_CEx;        ///<
  Spline * fXSecPA_NP;         ///<
  Spline * fXSecPA_PP;         ///<
  Spline * fXSecPA_NPP;        ///<
  Spline * fXSecPA_NNP;        ///<
  Spline * fXSecPA_NNPPP;      ///<
  Spline * fXSecPA_NPip;       ///<
  Spline * fXSecPA_NPipPi0;    ///<
  Spline * fXSecNA_Tot;        ///<
  Spline * fXSecNA_Elas;       ///<
  Spline * fXSecNA_Inel;       ///<
  Spline * fXSecNA_CEx;        ///<
  Spline * fXSecNA_NP;         ///<
  Spline * fXSecNA_PP;         ///<
  Spline * fXSecNA_NPP;        ///<
  Spline * fXSecNA_NNP;        ///<
  Spline * fXSecNA_NNPPP;      ///<
  Spline * fXSecNA_NPip;       ///<
  Spline * fXSecNA_NPipPi0;    ///<
  Spline * fXSecPipA_Tot;      ///< pi+A x-section splines
  Spline * fXSecPipA_Elas;     ///<
  Spline * fXSecPipA_Inel;     ///<
  Spline * fXSecPipA_CEx;      ///<
  Spline * fXSecPipA_NP;       ///<
  Spline * fXSecPipA_PP;       ///<
  Spline * fXSecPipA_NPP;      ///<
  Spline * fXSecPipA_NNP;      ///<
  Spline * fXSecPipA_NNPP;     ///<
  Spline * fXSecPipA_NPipPi0;  ///<
  Spline * fXSecPimA_Tot;      ///<
  Spline * fXSecPimA_Elas;     ///<
  Spline * fXSecPimA_Inel;     ///<
  Spline * fXSecPimA_CEx;      ///<
  Spline * fXSecPimA_NP;       ///<
  Spline * fXSecPimA_PP;       ///<
  Spline * fXSecPimA_NPP;      ///<
  Spline * fXSecPimA_NNP;      ///<
  Spline * fXSecPimA_NNPP;     ///<
  Spline * fXSecPimA_NPipPi0;  ///<
  Spline * fXSecPi0A_Tot;      ///<
  Spline * fXSecPi0A_Elas;     ///<
  Spline * fXSecPi0A_Inel;     ///<
  Spline * fXSecPi0A_CEx;      ///<
  Spline * fXSecPi0A_NP;       ///<
  Spline * fXSecPi0A_PP;       ///<
  Spline * fXSecPi0A_NPP;      ///<
  Spline * fXSecPi0A_NNP;      ///<
  Spline * fXSecPi0A_NNPP;     ///<
  Spline * fXSecPi0A_NPipPi0;  ///<

  //-- Sinleton cleaner
  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (INukeHadroData::fInstance !=0) {
            delete INukeHadroData::fInstance;
            INukeHadroData::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace
#endif //_INTRANUKE_HADRON_CROSS_SECTIONS_H_

