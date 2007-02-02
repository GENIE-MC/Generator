//____________________________________________________________________________
/*!

\class    genie::INukeHadroData

\brief    Singleton class to load & serve INTRANUKE hadron x-section splines.
          h+N x-sections come from SAID (Arndt, Workman, Strakovsky) PWA fit.
          p+Fe, pi+Fe x-sections from the Mashnik calculation.

	  Adapted from the inc/hadron_data.inc, gen/read_hadron_data.F and
          gen/calc_hadron_data.F from the fortran vrs of INTRANUKE.

\author   Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
          Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
          Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.

\created  October 03, 2006

*/
//____________________________________________________________________________

#ifndef _INTRANUKE_HADRON_CROSS_SECTIONS_H_
#define _INTRANUKE_HADRON_CROSS_SECTIONS_H_

namespace genie {

class Spline;

class INukeHadroData
{
public:
  static INukeHadroData * Instance (void);

  // ----
  // Note that, unlike most the rest of GENIE where everything is expressed
  // in natural units, all x-section splines included here are evaluated in
  // kinetic energies given in MeV and return the x-section value in mbarn
  // ----

  // Cross sections from SAID (Arndt, Workman, Strakovsky) PWA fit.
  //
  const Spline * XSecPipP_Elas (void) const { return fXSecPipP_Elas; }
  const Spline * XSecPipP_Reac (void) const { return fXSecPipP_Reac; }
  const Spline * XSecPipD_Abs  (void) const { return fXSecPipD_Abs;  }
  const Spline * XSecPimP_Elas (void) const { return fXSecPimP_Elas; }
  const Spline * XSecPimP_Reac (void) const { return fXSecPimP_Reac; }
  const Spline * XSecPimP_CEx  (void) const { return fXSecPimP_CEx;  }
  const Spline * XSecPP_Elas   (void) const { return fXSecPP_Elas;   }
  const Spline * XSecPP_Reac   (void) const { return fXSecPP_Reac;   }
  const Spline * XSecNP_Elas   (void) const { return fXSecNP_Elas;   }
  const Spline * XSecNP_Reac   (void) const { return fXSecNP_Reac;   }

  // Cross sections from Mashnik's calculations (p+Fe)
  //
  const Spline * XSecPFe_Elas  (void) const { return fXSecPFe_Elas; }
  const Spline * XSecPFe_Reac  (void) const { return fXSecPFe_Reac; }
  const Spline * XSecPFe_N     (void) const { return fXSecPFe_N;    }
  const Spline * XSecPFe_P     (void) const { return fXSecPFe_P;    }
  const Spline * XSecPFe_NP    (void) const { return fXSecPFe_NP;   }
  const Spline * XSecPFe_PP    (void) const { return fXSecPFe_PP;   }
  const Spline * XSecPFe_NPP   (void) const { return fXSecPFe_NPP;  }
  const Spline * XSecPFe_NNP   (void) const { return fXSecPFe_NNP;  }
  const Spline * XSecPFe_NNPP  (void) const { return fXSecPFe_NNPP; }
  const Spline * XSecPFe_Pim   (void) const { return fXSecPFe_Pim;  }
  const Spline * XSecPFe_Pi0   (void) const { return fXSecPFe_Pi0;  }
  const Spline * XSecPFe_Pip   (void) const { return fXSecPFe_Pip;  }

  // Cross sections from Mashnik's calculations (pi+Fe)
  //
  const Spline * XSecPiFe_P    (void) const { return fXSecPiFe_P;    }
  const Spline * XSecPiFe_PP   (void) const { return fXSecPiFe_PP;   }
  const Spline * XSecPiFe_PPP  (void) const { return fXSecPiFe_PPP;  }
  const Spline * XSecPiFe_N    (void) const { return fXSecPiFe_N;    }
  const Spline * XSecPiFe_NN   (void) const { return fXSecPiFe_NN;   }
  const Spline * XSecPiFe_NNN  (void) const { return fXSecPiFe_NNN;  }
  const Spline * XSecPiFe_NP   (void) const { return fXSecPiFe_NP;   }
  const Spline * XSecPiFe_NPP  (void) const { return fXSecPiFe_NPP;  }
  const Spline * XSecPiFe_NPPP (void) const { return fXSecPiFe_NPPP; }
  const Spline * XSecPiFe_NNP  (void) const { return fXSecPiFe_NNP;  }
  const Spline * XSecPiFe_NNPP (void) const { return fXSecPiFe_NNPP; }
  const Spline * XSecPiFe_Pi0  (void) const { return fXSecPiFe_Pi0;  }
  
  // Cross sections from data (Ashery, Carrol) and extrapolations
  //
  const Spline * XSecAshPiFe_Abs  (void) const { return fXSecAshPiFe_Abs;  }
  const Spline * XSecAshPiFe_Reac (void) const { return fXSecAshPiFe_Reac; }
  const Spline * XSecCarPiFe_Tot  (void) const { return fXSecCarPiFe_Tot;  }

  // Total x-sections needed from mean free path
  //
  const Spline * XSecPip_Tot   (void) const { return fXSecPip_Tot; }
  const Spline * XSecPim_Tot   (void) const { return fXSecPim_Tot; }
  const Spline * XSecPi0_Tot   (void) const { return fXSecPi0_Tot; }
  const Spline * XSecP_Tot     (void) const { return fXSecP_Tot;   }
  const Spline * XSecN_Tot     (void) const { return fXSecN_Tot;   }

  // Fractions of total cross sections needed for deciding particle fates
  //
  const Spline * FracPip_CEx   (void) const { return fFracPip_CEx;    }
  const Spline * FracPip_Elas  (void) const { return fFracPip_Elas;   }
  const Spline * FracPip_Reac  (void) const { return fFracPip_Reac;   }
  const Spline * FracPip_Abs   (void) const { return fFracPip_Abs;    }
  const Spline * FracPim_CEx   (void) const { return fFracPim_CEx;    }
  const Spline * FracPim_Elas  (void) const { return fFracPim_Elas;   }
  const Spline * FracPim_Reac  (void) const { return fFracPim_Reac;   }
  const Spline * FracPim_Abs   (void) const { return fFracPim_Abs;    }
  const Spline * FracPi0_CEx   (void) const { return fFracPi0_CEx;    }
  const Spline * FracPi0_Elas  (void) const { return fFracPi0_Elas;   }
  const Spline * FracPi0_Reac  (void) const { return fFracPi0_Reac;   }
  const Spline * FracPi0_Abs   (void) const { return fFracPi0_Abs;    }
  const Spline * FracP_Reac    (void) const { return fFracP_Reac;     }
  const Spline * FracN_Reac    (void) const { return fFracN_Reac;     }
  const Spline * FracPiA_Elas  (void) const { return fFracPiA_Elas;   }
  const Spline * FracPiA_Inel  (void) const { return fFracPiA_Inel;   }
  const Spline * FracPiA_CEx   (void) const { return fFracPiA_CEx;    }
  const Spline * FracPiA_Abs   (void) const { return fFracPiA_Abs;    }
  const Spline * FracPiA_PP    (void) const { return fFracPiA_PP;     }
  const Spline * FracPiA_NPP   (void) const { return fFracPiA_NPP;    }
  const Spline * FracPiA_NNP   (void) const { return fFracPiA_NNP;    }
  const Spline * FracPiA_4N4P  (void) const { return fFracPiA_4N4P;   }
  const Spline * FracPiA_PiProd(void) const { return fFracPiA_PiProd; }
  const Spline * FracPA_Elas   (void) const { return fFracPA_Elas;    }
  const Spline * FracPA_Inel   (void) const { return fFracPA_Inel;    }
  const Spline * FracPA_Abs    (void) const { return fFracPA_Abs;     }
  const Spline * FracPA_PP     (void) const { return fFracPA_PP;      }
  const Spline * FracPA_NPP    (void) const { return fFracPA_NPP;     }
  const Spline * FracPA_NNP    (void) const { return fFracPA_NNP;     }
  const Spline * FracPA_4N4P   (void) const { return fFracPA_4N4P;    }
  const Spline * FracPA_PiProd (void) const { return fFracPA_PiProd;  }

private:
  INukeHadroData();
  INukeHadroData(const INukeHadroData & shx);
 ~INukeHadroData();

  void LoadData (void); ///< adapted from neugen3's read_hadron_data.F 
  void CalcData (void); ///< adapted from neugen3's calc_hadron_data.F 

  static INukeHadroData * fInstance;

  //-- loaded cross section splines

  // SAID h+N x-section splines
  Spline * fXSecPipP_Elas;  ///< pi+p ELAS x-section
  Spline * fXSecPipP_Reac;  ///< pi+p REAC x-section
  Spline * fXSecPipD_Abs;   ///< pi+d ABS  x-section
  Spline * fXSecPimP_Elas;  ///< pi-p ELAS x-section
  Spline * fXSecPimP_Reac;  ///< pi-p REAC x-section
  Spline * fXSecPimP_CEx;   ///< pi-p CEX  x-section
  Spline * fXSecPP_Elas;    ///< pp   ELAS x-section
  Spline * fXSecPP_Reac;    ///< pp   REAC x-section
  Spline * fXSecNP_Elas;    ///< np   ELAS x-section
  Spline * fXSecNP_Reac;    ///< np   REAC x-section //?????????

  // Mashnik p+Fe x-section splines
  Spline * fXSecPFe_Elas;   ///< pFe ELAS   x-section
  Spline * fXSecPFe_Reac;   ///< pFe REAC   x-section
  Spline * fXSecPFe_P;      ///< pFe->pX    x-section
  Spline * fXSecPFe_N;      ///< pFe->nX    x-section
  Spline * fXSecPFe_NP;     ///< pFe->npX   x-section
  Spline * fXSecPFe_PP;     ///< pFe->ppX   x-section
  Spline * fXSecPFe_NPP;    ///< pFe->nppX  x-section
  Spline * fXSecPFe_NNP;    ///< pFe->nnpX  x-section
  Spline * fXSecPFe_NNPP;   ///< pFe->nnppX x-section
  Spline * fXSecPFe_Pim;    ///< pFe->pi-X  x-section
  Spline * fXSecPFe_Pi0;    ///< pFe->pi0X  x-section
  Spline * fXSecPFe_Pip;    ///< pFe->pi+X  x-section

  // Mashnik pi+Fe x-section splines
  Spline * fXSecPiFe_P;     ///< piFe->pX    x-section
  Spline * fXSecPiFe_PP;    ///< piFe->ppX   x-section
  Spline * fXSecPiFe_PPP;   ///< piFe->pppX  x-section
  Spline * fXSecPiFe_N;     ///< piFe->nX    x-section
  Spline * fXSecPiFe_NN;    ///< piFe->nnX   x-section
  Spline * fXSecPiFe_NNN;   ///< piFe->nnnX  x-section
  Spline * fXSecPiFe_NP;    ///< piFe->npX   x-section
  Spline * fXSecPiFe_NPP;   ///< piFe->nppX  x-section
  Spline * fXSecPiFe_NPPP;  ///< piFe->npppX x-section
  Spline * fXSecPiFe_NNP;   ///< piFe->nnpX  x-section
  Spline * fXSecPiFe_NNPP;  ///< piFe->nnppX x-section
  Spline * fXSecPiFe_Pi0;   ///< piFe->pi0X  x-section

  // Hadronic x-sections from data (Ash & Carrol)
  Spline * fXSecAshPiFe_Abs;  ///<
  Spline * fXSecAshPiFe_Reac; ///<
  Spline * fXSecCarPiFe_Tot;  ///<

  // Total x-sections needed from mean free path
  Spline * fXSecPip_Tot;     ///<
  Spline * fXSecPim_Tot;     ///<
  Spline * fXSecPi0_Tot;     ///<
  Spline * fXSecP_Tot;       ///<
  Spline * fXSecN_Tot;       ///<

  // Fractions of total x-sections needed for particle fates
  Spline * fFracPip_CEx;    ///< fraction of CEx for pi+
  Spline * fFracPip_Elas;   ///< fraction of CEx+Elas for pi+
  Spline * fFracPip_Reac;   ///< fraction of CEx+Elas+PiProd for pi+
  Spline * fFracPip_Abs;    ///< fraction of CEx+Elas+PiProd+Abs for pi+
  Spline * fFracPim_CEx;    ///<
  Spline * fFracPim_Elas;   ///<
  Spline * fFracPim_Reac;   ///<
  Spline * fFracPim_Abs;    ///<
  Spline * fFracPi0_CEx;    ///<
  Spline * fFracPi0_Elas;   ///<
  Spline * fFracPi0_Reac;   ///<
  Spline * fFracPi0_Abs;    ///<
  Spline * fFracP_Reac;     ///<
  Spline * fFracN_Reac;     ///<
  Spline * fFracPiA_Elas;   ///<
  Spline * fFracPiA_Inel;   ///<
  Spline * fFracPiA_CEx;    ///<
  Spline * fFracPiA_Abs;    ///<
  Spline * fFracPiA_PP;     ///<
  Spline * fFracPiA_NPP;    ///<
  Spline * fFracPiA_NNP;    ///<
  Spline * fFracPiA_4N4P;   ///<
  Spline * fFracPiA_PiProd; ///<
  Spline * fFracPA_Elas;    ///<
  Spline * fFracPA_Inel;    ///<
  Spline * fFracPA_Abs;     ///<
  Spline * fFracPA_PP;      ///<
  Spline * fFracPA_NPP;     ///<
  Spline * fFracPA_NNP;     ///<
  Spline * fFracPA_4N4P;    ///<
  Spline * fFracPA_PiProd;  ///<

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

