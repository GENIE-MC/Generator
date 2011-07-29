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

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>, Rutherford Lab.
          Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
	  Aaron Meyer <asm58@pitt.edu>, Pittsburgh Univ.
	  Alex Bell, Pittsburgh Univ.

\created  February 01, 2007

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_HADRON_CROSS_SECTIONS_H_
#define _INTRANUKE_HADRON_CROSS_SECTIONS_H_

#include "HadronTransport/INukeHadroFates.h"
#include "GHEP/GHepParticle.h"
#include "Numerical/BLI2D.h"

class TGraph2D;

namespace genie {

class Spline;

class INukeHadroData
{
public:
  static INukeHadroData * Instance (void);

// Note that, unlike most the rest of GENIE where everything is expressed
// in natural units, all x-section splines included here are evaluated in
// kinetic energies given in MeV and return the x-section value in mbarns

  double XSec (int hpdgc, int tgt, int nprod, INukeFateHN_t rxnType, double ke, double costh) const;

  double XSec (int hpdgc, INukeFateHA_t fate, double ke) const;
  double XSec (int hpdgc, INukeFateHN_t fate, double ke, int targA, int targZ) const;
  double Frac (int hpdgc, INukeFateHA_t fate, double ke) const;
  double Frac (int hpdgc, INukeFateHN_t fate, double ke, int targA=0, int targZ=0) const;
  double IntBounce       (const GHepParticle* p, int target, int s1, INukeFateHN_t fate);
  //int    AngleAndProduct (const GHepParticle* p, int target, double &angle, INukeFateHN_t fate);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // hN mode hadron x-section splines
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const Spline * XSecPipn_Tot     (void) const { return fXSecPipn_Tot;     }
  const Spline * XSecPipn_CEx     (void) const { return fXSecPipn_CEx;     }
  const Spline * XSecPipn_Elas    (void) const { return fXSecPipn_Elas;    }
  const Spline * XSecPipn_Reac    (void) const { return fXSecPipn_Reac;    }
  const Spline * XSecPipp_Tot     (void) const { return fXSecPipp_Tot;     }
  const Spline * XSecPipp_CEx     (void) const { return fXSecPipp_CEx;     }
  const Spline * XSecPipp_Elas    (void) const { return fXSecPipp_Elas;    }
  const Spline * XSecPipp_Reac    (void) const { return fXSecPipp_Reac;    }
  const Spline * XSecPipp_Abs     (void) const { return fXSecPipd_Abs;     }
  const Spline * XSecPi0n_Tot     (void) const { return fXSecPi0n_Tot;     }
  const Spline * XSecPi0n_CEx     (void) const { return fXSecPi0n_CEx;     }
  const Spline * XSecPi0n_Elas    (void) const { return fXSecPi0n_Elas;    }
  const Spline * XSecPi0n_Reac    (void) const { return fXSecPi0n_Reac;    }
  const Spline * XSecPi0p_Tot     (void) const { return fXSecPi0p_Tot;     }
  const Spline * XSecPi0p_CEx     (void) const { return fXSecPi0p_CEx;     }
  const Spline * XSecPi0p_Elas    (void) const { return fXSecPi0p_Elas;    }
  const Spline * XSecPi0p_Reac    (void) const { return fXSecPi0p_Reac;    }
  const Spline * XSecPi0p_Abs     (void) const { return fXSecPi0d_Abs;     }
  const Spline * XSecPp_Tot       (void) const { return fXSecPp_Tot;       }
  const Spline * XSecPp_Elas      (void) const { return fXSecPp_Elas;      }
  const Spline * XSecPp_Reac      (void) const { return fXSecPp_Reac;      }
  const Spline * XSecPn_Tot       (void) const { return fXSecPn_Tot;       }
  const Spline * XSecPn_Elas      (void) const { return fXSecPn_Elas;      } 
  const Spline * XSecPn_Reac      (void) const { return fXSecPn_Reac;      }
  const Spline * XSecNn_Tot       (void) const { return fXSecNn_Tot;       }
  const Spline * XSecNn_Elas      (void) const { return fXSecNn_Elas;      } 
  const Spline * XSecNn_Reac      (void) const { return fXSecNn_Reac;      }
  const Spline * XSecKpn_Elas     (void) const { return fXSecKpn_Elas;     }
  const Spline * XSecKpp_Elas     (void) const { return fXSecKpp_Elas;     }
  const Spline * XSecKpN_Tot      (void) const { return fXSecKpN_Tot;      }
  const Spline * XSecGamp_fs      (void) const { return fXSecGamp_fs;      }  
  const Spline * XSecGamn_fs      (void) const { return fXSecGamn_fs;      }
  const Spline * XSecGamN_Tot     (void) const { return fXSecGamN_Tot;     }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // hA mode hadron x-section splines
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const Spline * XSecPA_Tot       (void) const { return fXSecPA_Tot;       }
  const Spline * XSecPA_Elas      (void) const { return fXSecPA_Elas;      }
  const Spline * XSecPA_Inel      (void) const { return fXSecPA_Inel;      }
  const Spline * XSecPA_CEx       (void) const { return fXSecPA_CEx;       }
  const Spline * XSecPA_Abs       (void) const { return fXSecPA_Abs;       }
  const Spline * XSecPA_Pipro     (void) const { return fXSecPA_Pipro;     }
  const Spline * XSecNA_Tot       (void) const { return fXSecNA_Tot;       }
  const Spline * XSecNA_Elas      (void) const { return fXSecNA_Elas;      }
  const Spline * XSecNA_Inel      (void) const { return fXSecNA_Inel;      }
  const Spline * XSecNA_CEx       (void) const { return fXSecNA_CEx;       }
  const Spline * XSecNA_Abs       (void) const { return fXSecNA_Abs;       }
  const Spline * XSecNA_Pipro     (void) const { return fXSecNA_Pipro;     }
  const Spline * XSecPipA_Tot     (void) const { return fXSecPipA_Tot;     }
  const Spline * XSecPipA_Elas    (void) const { return fXSecPipA_Elas;    }
  const Spline * XSecPipA_Inel    (void) const { return fXSecPipA_Inel;    }
  const Spline * XSecPipA_CEx     (void) const { return fXSecPipA_CEx;     }
  const Spline * XSecPipA_Abs     (void) const { return fXSecPipA_Abs;     }
  const Spline * XSecPipA_NPipPi0 (void) const { return fXSecPipA_NPipPi0; }
  const Spline * XSecPimA_Tot     (void) const { return fXSecPimA_Tot;     }
  const Spline * XSecPimA_Elas    (void) const { return fXSecPimA_Elas;    }
  const Spline * XSecPimA_Inel    (void) const { return fXSecPimA_Inel;    }
  const Spline * XSecPimA_CEx     (void) const { return fXSecPimA_CEx;     }
  const Spline * XSecPimA_Abs     (void) const { return fXSecPimA_Abs;     }
  const Spline * XSecPimA_NPipPi0 (void) const { return fXSecPimA_NPipPi0; }
  const Spline * XSecPi0A_Tot     (void) const { return fXSecPi0A_Tot;     }
  const Spline * XSecPi0A_Elas    (void) const { return fXSecPi0A_Elas;    }
  const Spline * XSecPi0A_Inel    (void) const { return fXSecPi0A_Inel;    }
  const Spline * XSecPi0A_CEx     (void) const { return fXSecPi0A_CEx;     }
  const Spline * XSecPi0A_Abs     (void) const { return fXSecPi0A_Abs;     }
  const Spline * XSecPi0A_NPipPi0 (void) const { return fXSecPi0A_NPipPi0; }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // hN mode TGraph2D XSec objects
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /*const TGraph2D * const hN2dXSecPP_Elas          (void) const { return fhN2dXSecPP_Elas;          }
  const TGraph2D * const hN2dXSecNP_Elas          (void) const { return fhN2dXSecNP_Elas;          }
  const TGraph2D * const hN2dXSecPipN_Elas        (void) const { return fhN2dXSecPipN_Elas;        }
  const TGraph2D * const hN2dXSecPi0N_Elas        (void) const { return fhN2dXSecPi0N_Elas;        }
  const TGraph2D * const hN2dXSecPimN_Elas        (void) const { return fhN2dXSecPimN_Elas;        }
  const TGraph2D * const hN2dXSecKpN_Elas         (void) const { return fhN2dXSecKpN_Elas;         }
  const TGraph2D * const hN2dXSecKpP_Elas         (void) const { return fhN2dXSecKpP_Elas;         }
  const TGraph2D * const hN2dXSecPiN_CEx          (void) const { return fhN2dXSecPiN_CEx;          }
  const TGraph2D * const hN2dXSecPiN_Abs          (void) const { return fhN2dXSecPiN_Abs;          }
  const TGraph2D * const hN2dXSecGamPi0P_Inelas   (void) const { return fhN2dXSecGamPi0P_Inelas;   }
  const TGraph2D * const hN2dXSecGamPi0N_Inelas   (void) const { return fhN2dXSecGamPi0N_Inelas;   }
  const TGraph2D * const hN2dXSecGamPipN_Inelas   (void) const { return fhN2dXSecGamPipN_Inelas;   }
  const TGraph2D * const hN2dXSecGamPimP_Inelas   (void) const { return fhN2dXSecGamPimP_Inelas;   }*/

  const BLI2DNonUnifGrid * const hN2dXSecPP_Elas          (void) const { return fhN2dXSecPP_Elas;          }
  const BLI2DNonUnifGrid * const hN2dXSecNP_Elas          (void) const { return fhN2dXSecNP_Elas;          }
  const BLI2DNonUnifGrid * const hN2dXSecPipN_Elas        (void) const { return fhN2dXSecPipN_Elas;        }
  const BLI2DNonUnifGrid * const hN2dXSecPi0N_Elas        (void) const { return fhN2dXSecPi0N_Elas;        }
  const BLI2DNonUnifGrid * const hN2dXSecPimN_Elas        (void) const { return fhN2dXSecPimN_Elas;        }
  const BLI2DNonUnifGrid * const hN2dXSecKpN_Elas         (void) const { return fhN2dXSecKpN_Elas;         }
  const BLI2DNonUnifGrid * const hN2dXSecKpP_Elas         (void) const { return fhN2dXSecKpP_Elas;         }
  const BLI2DNonUnifGrid * const hN2dXSecPiN_CEx          (void) const { return fhN2dXSecPiN_CEx;          }
  const BLI2DNonUnifGrid * const hN2dXSecPiN_Abs          (void) const { return fhN2dXSecPiN_Abs;          }
  const BLI2DNonUnifGrid * const hN2dXSecGamPi0P_Inelas   (void) const { return fhN2dXSecGamPi0P_Inelas;   }
  const BLI2DNonUnifGrid * const hN2dXSecGamPi0N_Inelas   (void) const { return fhN2dXSecGamPi0N_Inelas;   }
  const BLI2DNonUnifGrid * const hN2dXSecGamPipN_Inelas   (void) const { return fhN2dXSecGamPipN_Inelas;   }
  const BLI2DNonUnifGrid * const hN2dXSecGamPimP_Inelas   (void) const { return fhN2dXSecGamPimP_Inelas;   }

  static double fMinKinEnergy;   ///<
  static double fMaxKinEnergyHA; ///<
  static double fMaxKinEnergyHN; ///<

private:
  INukeHadroData();
  INukeHadroData(const INukeHadroData & shx);
 ~INukeHadroData();

  void LoadCrossSections(void); 

  void ReadhNFile(
         string filename, double ke, int npoints, int & curr_point,
         /*double * ke_array,*/ double * costh_array, double * xsec_array, int cols);

  static INukeHadroData * fInstance;

  Spline * fXSecPipn_Tot;      ///< pi+n hA x-section splines
  Spline * fXSecPipn_CEx;      ///<
  Spline * fXSecPipn_Elas;     ///<
  Spline * fXSecPipn_Reac;     ///<
  Spline * fXSecPipp_Tot;      ///< pi+p hA x-section splines
  Spline * fXSecPipp_CEx;      ///<
  Spline * fXSecPipp_Elas;     ///<
  Spline * fXSecPipp_Reac;     ///<
  Spline * fXSecPipd_Abs;      ///<
  Spline * fXSecPi0n_Tot;      ///< pi0n hA x-section splines
  Spline * fXSecPi0n_CEx;      ///<
  Spline * fXSecPi0n_Elas;     ///<
  Spline * fXSecPi0n_Reac;     ///<
  Spline * fXSecPi0p_Tot;      ///< pi0p hA x-section splines
  Spline * fXSecPi0p_CEx;      ///<
  Spline * fXSecPi0p_Elas;     ///<
  Spline * fXSecPi0p_Reac;     ///<
  Spline * fXSecPi0d_Abs;      ///<
  Spline * fXSecPp_Tot;        ///< p/nN x-section splines
  Spline * fXSecPp_Elas;       ///<
  Spline * fXSecPp_Reac;       ///<
  Spline * fXSecPn_Tot;        ///<
  Spline * fXSecPn_Elas;       ///<
  Spline * fXSecPn_Reac;       ///<
  Spline * fXSecNn_Tot;        ///<
  Spline * fXSecNn_Elas;       ///<
  Spline * fXSecNn_Reac;       ///<
  Spline * fXSecKpn_Elas;      ///< K+N x-section splines
  Spline * fXSecKpp_Elas;      ///<
  Spline * fXSecKpN_Tot;       ///<
  Spline * fXSecPA_Tot;        ///< N+A x-section splines
  Spline * fXSecPA_Elas;       ///<
  Spline * fXSecPA_Inel;       ///<
  Spline * fXSecPA_CEx;        ///<
  Spline * fXSecPA_Abs;        ///<
  Spline * fXSecPA_Pipro;       ///<
  //  Spline * fXSecPA_NPipPi0;    ///<
  Spline * fXSecNA_Tot;        ///<
  Spline * fXSecNA_Elas;       ///<
  Spline * fXSecNA_Inel;       ///<
  Spline * fXSecNA_CEx;        ///<
  Spline * fXSecNA_Abs;        ///<
  Spline * fXSecNA_Pipro;       ///<
  //  Spline * fXSecNA_NPipPi0;    ///<
  Spline * fXSecPipA_Tot;      ///< pi+A x-section splines
  Spline * fXSecPipA_Elas;     ///<
  Spline * fXSecPipA_Inel;     ///<
  Spline * fXSecPipA_CEx;      ///<
  Spline * fXSecPipA_Abs;      ///<
  Spline * fXSecPipA_NPipPi0;  ///<
  Spline * fXSecPimA_Tot;      ///<
  Spline * fXSecPimA_Elas;     ///<
  Spline * fXSecPimA_Inel;     ///<
  Spline * fXSecPimA_CEx;      ///<
  Spline * fXSecPimA_Abs;      ///<
  Spline * fXSecPimA_NPipPi0;  ///<
  Spline * fXSecPi0A_Tot;      ///<
  Spline * fXSecPi0A_Elas;     ///<
  Spline * fXSecPi0A_Inel;     ///<
  Spline * fXSecPi0A_CEx;      ///<
  Spline * fXSecPi0A_Abs;      ///<
  Spline * fXSecPi0A_NPipPi0;  ///<
  Spline * fXSecKpA_Tot;       ///< K+A x-section splines
  Spline * fXSecKpA_Elas;      ///<
  Spline * fXSecGamp_fs;       ///< gamma A x-section splines
  Spline * fXSecGamn_fs;       ///<
  Spline * fXSecGamN_Tot;      ///<

  /*TGraph2D * fhN2dXSecPP_Elas;
  TGraph2D * fhN2dXSecNP_Elas;
  TGraph2D * fhN2dXSecPipN_Elas;
  TGraph2D * fhN2dXSecPi0N_Elas;
  TGraph2D * fhN2dXSecPimN_Elas;
  TGraph2D * fhN2dXSecKpN_Elas;
  TGraph2D * fhN2dXSecKpP_Elas;
  TGraph2D * fhN2dXSecPiN_CEx;
  TGraph2D * fhN2dXSecPiN_Abs;
  TGraph2D * fhN2dXSecGamPi0P_Inelas;
  TGraph2D * fhN2dXSecGamPi0N_Inelas;
  TGraph2D * fhN2dXSecGamPipN_Inelas;
  TGraph2D * fhN2dXSecGamPimP_Inelas;*/
  BLI2DNonUnifGrid * fhN2dXSecPP_Elas;
  BLI2DNonUnifGrid * fhN2dXSecNP_Elas;
  BLI2DNonUnifGrid * fhN2dXSecPipN_Elas;
  BLI2DNonUnifGrid * fhN2dXSecPi0N_Elas;
  BLI2DNonUnifGrid * fhN2dXSecPimN_Elas;
  BLI2DNonUnifGrid * fhN2dXSecKpN_Elas;
  BLI2DNonUnifGrid * fhN2dXSecKpP_Elas;
  BLI2DNonUnifGrid * fhN2dXSecPiN_CEx;
  BLI2DNonUnifGrid * fhN2dXSecPiN_Abs;
  BLI2DNonUnifGrid * fhN2dXSecGamPi0P_Inelas;
  BLI2DNonUnifGrid * fhN2dXSecGamPi0N_Inelas;
  BLI2DNonUnifGrid * fhN2dXSecGamPipN_Inelas;
  BLI2DNonUnifGrid * fhN2dXSecGamPimP_Inelas;

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

