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

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_HADRON_CROSS_SECTIONS_2018_H_
#define _INTRANUKE_HADRON_CROSS_SECTIONS_2018_H_

#include "Physics/HadronTransport/INukeHadroFates2018.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Numerical/BLI2D.h"

class TGraph2D;

namespace genie {

class Spline;

class INukeHadroData2018
{
public:
  static INukeHadroData2018 * Instance (void);

// Note that, unlike most the rest of GENIE where everything is expressed
// in natural units, all x-section splines included here are evaluated in
// kinetic energies given in MeV and return the x-section value in mbarns

  double XSec (int hpdgc, int tgt, int nprod, INukeFateHN_t rxnType, double ke, double costh) const;
  double XSec (int hpdgc, INukeFateHN_t fate, double ke, int targA, int targZ) const;
  double FracADep (int hpdgc, INukeFateHA_t fate, double ke, int targA) const;
  double FracAIndep (int hpdgc, INukeFateHA_t fate, double ke) const;
  double Frac (int hpdgc, INukeFateHN_t fate, double ke, int targA=0, int targZ=0) const;
  double IntBounce       (const GHepParticle* p, int target, int s1, INukeFateHN_t fate);


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
  const Spline * XSecKpn_CEx      (void) const { return fXSecKpn_CEx;      }
  const Spline * XSecKpp_Elas     (void) const { return fXSecKpp_Elas;     }
  const Spline * XSecKpN_Abs      (void) const { return fXSecKpN_Abs;      } //not implemented
  const Spline * XSecKpN_Tot      (void) const { return fXSecKpN_Tot;      }
  const Spline * XSecGamp_fs      (void) const { return fXSecGamp_fs;      }
  const Spline * XSecGamn_fs      (void) const { return fXSecGamn_fs;      }
  const Spline * XSecGamN_Tot     (void) const { return fXSecGamN_Tot;     }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // hA mode hadron x-section splines
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const Spline * FracPA_Tot       (void) const { return fFracPA_Tot;       }
  const Spline * FracPA_Elas      (void) const { return fFracPA_Elas;      }
  const Spline * FracPA_Inel      (void) const { return fFracPA_Inel;      }
  const Spline * FracPA_CEx       (void) const { return fFracPA_CEx;       }
  const Spline * FracPA_Abs       (void) const { return fFracPA_Abs;       }
  const Spline * FracPA_PiPro     (void) const { return fFracPA_PiPro;     }
  const Spline * FracNA_Tot       (void) const { return fFracNA_Tot;       }
  const Spline * FracNA_Elas      (void) const { return fFracNA_Elas;      }
  const Spline * FracNA_Inel      (void) const { return fFracNA_Inel;      }
  const Spline * FracNA_CEx       (void) const { return fFracNA_CEx;       }
  const Spline * FracNA_Abs       (void) const { return fFracNA_Abs;       }
  const Spline * FracNA_PiPro     (void) const { return fFracNA_PiPro;     }
  const Spline * FracKA_Tot       (void) const { return fFracKA_Tot;       }
  const Spline * FracKA_Elas      (void) const { return fFracKA_Elas;      }
  const Spline * FracKA_Inel      (void) const { return fFracKA_Inel;      }
  const Spline * FracKA_CEx       (void) const { return fFracKA_CEx;       }
  const Spline * FracKA_Abs       (void) const { return fFracKA_Abs;       }

  const Spline * FracPA_Cmp       (void) const { return fFracPA_Cmp;      }
  const Spline * FracNA_Cmp       (void) const { return fFracNA_Cmp;      } //suarez


  const BLI2DNonUnifGrid * hN2dXSecPP_Elas          (void) const { return fhN2dXSecPP_Elas;          }
  const BLI2DNonUnifGrid * hN2dXSecNP_Elas          (void) const { return fhN2dXSecNP_Elas;          }
  const BLI2DNonUnifGrid * hN2dXSecPipN_Elas        (void) const { return fhN2dXSecPipN_Elas;        }
  const BLI2DNonUnifGrid * hN2dXSecPi0N_Elas        (void) const { return fhN2dXSecPi0N_Elas;        }
  const BLI2DNonUnifGrid * hN2dXSecPimN_Elas        (void) const { return fhN2dXSecPimN_Elas;        }
  const BLI2DNonUnifGrid * hN2dXSecKpN_Elas         (void) const { return fhN2dXSecKpN_Elas;         }
  const BLI2DNonUnifGrid * hN2dXSecKpP_Elas         (void) const { return fhN2dXSecKpP_Elas;         }
  const BLI2DNonUnifGrid * hN2dXSecPiN_CEx          (void) const { return fhN2dXSecPiN_CEx;          }
  const BLI2DNonUnifGrid * hN2dXSecPiN_Abs          (void) const { return fhN2dXSecPiN_Abs;          }
  const BLI2DNonUnifGrid * hN2dXSecGamPi0P_Inelas   (void) const { return fhN2dXSecGamPi0P_Inelas;   }
  const BLI2DNonUnifGrid * hN2dXSecGamPi0N_Inelas   (void) const { return fhN2dXSecGamPi0N_Inelas;   }
  const BLI2DNonUnifGrid * hN2dXSecGamPipN_Inelas   (void) const { return fhN2dXSecGamPipN_Inelas;   }
  const BLI2DNonUnifGrid * hN2dXSecGamPimP_Inelas   (void) const { return fhN2dXSecGamPimP_Inelas;   }

  static double fMinKinEnergy;   ///<
  static double fMaxKinEnergyHA; ///<
  static double fMaxKinEnergyHN; ///<

private:
  INukeHadroData2018();
  INukeHadroData2018(const INukeHadroData2018 & shx);
 ~INukeHadroData2018();

  void LoadCrossSections(void);

  void ReadhNFile(
         string filename, double ke, int npoints, int & curr_point,
         /*double * ke_array,*/ double * costh_array, double * xsec_array, int cols);

  static INukeHadroData2018 * fInstance;

  Spline * fXSecPipn_Tot;      ///< pi+n hN x-section splines
  Spline * fXSecPipn_CEx;      ///<
  Spline * fXSecPipn_Elas;     ///<
  Spline * fXSecPipn_Reac;     ///<
  Spline * fXSecPipp_Tot;      ///< pi+p hN x-section splines
  Spline * fXSecPipp_CEx;      ///<
  Spline * fXSecPipp_Elas;     ///<
  Spline * fXSecPipp_Reac;     ///<
  Spline * fXSecPipd_Abs;      ///<
  Spline * fXSecPi0n_Tot;      ///< pi0n hN x-section splines
  Spline * fXSecPi0n_CEx;      ///<
  Spline * fXSecPi0n_Elas;     ///<
  Spline * fXSecPi0n_Reac;     ///<
  Spline * fXSecPi0p_Tot;      ///< pi0p hN x-section splines
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
  Spline * fXSecKpn_CEx;       ///<
  Spline * fXSecKpN_Abs;        ///<
  Spline * fXSecKpN_Tot;       ///<
  Spline * fFracPA_Tot;        ///< N+A x-section splines
  Spline * fFracPA_Elas;       ///<
  Spline * fFracPA_Inel;       ///<
  Spline * fFracPA_CEx;        ///<
  Spline * fFracPA_Abs;        ///<
  Spline * fFracPA_PiPro;       ///<
  Spline * fFracNA_Tot;        ///<
  Spline * fFracNA_Elas;       ///<
  Spline * fFracNA_Inel;       ///<
  Spline * fFracNA_CEx;        ///<
  Spline * fFracNA_Abs;        ///<
  Spline * fFracNA_PiPro;       ///<
  Spline * fFracKA_Tot;        ///< K+A x-section splines
  Spline * fFracKA_Elas;       ///<
  Spline * fFracKA_CEx;        ///<
  Spline * fFracKA_Inel;       ///<
  Spline * fFracKA_Abs;        ///<
  Spline * fXSecGamp_fs;       ///< gamma A x-section splines
  Spline * fXSecGamn_fs;       ///<
  Spline * fXSecGamN_Tot;      ///<
  Spline * fXSecPp_Cmp;        ///< NN cmp (compound nucleus) fate
  Spline * fXSecPn_Cmp;
  Spline * fXSecNn_Cmp;
  Spline * fFracPA_Cmp;
  Spline * fFracNA_Cmp;

//  TGraph2D * TPipA_Tot;
  TGraph2D * TfracPipA_CEx;
//  TGraph2D * TfracPipA_Elas;
  TGraph2D * TfracPipA_Inelas;
  TGraph2D * TfracPipA_Abs;
  TGraph2D * TfracPipA_PiPro;

  BLI2DNonUnifGrid * fhN2dXSecPP_Elas;
  BLI2DNonUnifGrid * fhN2dXSecNP_Elas;
  BLI2DNonUnifGrid * fhN2dXSecPipN_Elas;
  BLI2DNonUnifGrid * fhN2dXSecPi0N_Elas;
  BLI2DNonUnifGrid * fhN2dXSecPimN_Elas;
  BLI2DNonUnifGrid * fhN2dXSecKpN_Elas;
  BLI2DNonUnifGrid * fhN2dXSecKpP_Elas;
  BLI2DNonUnifGrid * fhN2dXSecKpN_CEx;
  BLI2DNonUnifGrid * fhN2dXSecKpN_Abs;
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
         if (INukeHadroData2018::fInstance !=0) {
            delete INukeHadroData2018::fInstance;
            INukeHadroData2018::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace
#endif //_INTRANUKE_HADRON_CROSS_SECTIONS_2015_H_

