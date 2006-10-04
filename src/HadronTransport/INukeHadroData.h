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

  //-- cross section splines

  // Cross sections from SAID (Arndt, Workman, Strakovsky) PWA fit.
  //
  const Spline * XsPipPElas (void) const { return fXsPipPElas; }
  const Spline * XsPipPReac (void) const { return fXsPipPReac; }
  const Spline * XsPipDAbs  (void) const { return fXsPipDAbs;  }
  const Spline * XsPimPElas (void) const { return fXsPimPElas; }
  const Spline * XsPimPReac (void) const { return fXsPimPReac; }
  const Spline * XsPimPCEx  (void) const { return fXsPimPCEx;  }
  const Spline * XsPPElas   (void) const { return fXsPPElas;   }
  const Spline * XsPPReac   (void) const { return fXsPPReac;   }
  const Spline * XsNPElas   (void) const { return fXsNPElas;   }
  const Spline * XsNPReac   (void) const { return fXsNPReac;   }

  // Cross sections from Mashnik's calculations (p+Fe)
  //
  const Spline * XsPFeElas  (void) const { return fXsPFeElas; }
  const Spline * XsPFeReac  (void) const { return fXsPFeReac; }
  const Spline * XsPFeP     (void) const { return fXsPFeP;    }
  const Spline * XsPFePP    (void) const { return fXsPFePP;   }
  const Spline * XsPFeNPP   (void) const { return fXsPFeNPP;  }
  const Spline * XsPFeNNP   (void) const { return fXsPFeNNP;  }
  const Spline * XsPFeNNPP  (void) const { return fXsPFeNNPP; }
  const Spline * XsPFePim   (void) const { return fXsPFePim;  }
  const Spline * XsPFePi0   (void) const { return fXsPFePi0;  }
  const Spline * XsPFePip   (void) const { return fXsPFePip;  }

  // Cross sections from Mashnik's calculations (pi+Fe)
  //
  const Spline * XsPiFeP    (void) const { return fXsPiFeP;    }
  const Spline * XsPiFePP   (void) const { return fXsPiFePP;   }
  const Spline * XsPiFePPP  (void) const { return fXsPiFePPP;  }
  const Spline * XsPiFeN    (void) const { return fXsPiFeN;    }
  const Spline * XsPiFeNN   (void) const { return fXsPiFeNN;   }
  const Spline * XsPiFeNNN  (void) const { return fXsPiFeNNN;  }
  const Spline * XsPiFeNP   (void) const { return fXsPiFeNP;   }
  const Spline * XsPiFeNPP  (void) const { return fXsPiFeNPP;  }
  const Spline * XsPiFeNPPP (void) const { return fXsPiFeNPPP; }
  const Spline * XsPiFeNNP  (void) const { return fXsPiFeNNP;  }
  const Spline * XsPiFeNNPP (void) const { return fXsPiFeNNPP; }
  const Spline * XsPiFePi0  (void) const { return fXsPiFePi0;  }
  
  // Cross sections from data (Ashery, Carrol) and extrapolations
  //
  const Spline * XsAshPiFeAbs  (void) const { return fXsAshPiFeAbs;  }
  const Spline * XsAshPiFeReac (void) const { return fXsAshPiFeReac; }
  const Spline * XsCarPiFeTot  (void) const { return fXsCarPiFeTot;  }

  // Total x-sections needed from mean free path
  //
  const Spline * XsPipTot   (void) const { return fXsPipTot; }
  const Spline * XsPimTot   (void) const { return fXsPimTot; }
  const Spline * XsPi0Tot   (void) const { return fXsPi0Tot; }
  const Spline * XsPTot     (void) const { return fXsPTot;   }
  const Spline * XsNTot     (void) const { return fXsNTot;   }

  // Fractions of total cross sections needed for deciding particle fates
  //
  const Spline * FrPipCEx   (void) const { return fFrPipCEx;    }
  const Spline * FrPipElas  (void) const { return fFrPipElas;   }
  const Spline * FrPipReac  (void) const { return fFrPipReac;   }
  const Spline * FrPipAbs   (void) const { return fFrPipAbs;    }
  const Spline * FrPimCEx   (void) const { return fFrPimCEx;    }
  const Spline * FrPimElas  (void) const { return fFrPimElas;   }
  const Spline * FrPimReac  (void) const { return fFrPimReac;   }
  const Spline * FrPimAbs   (void) const { return fFrPimAbs;    }
  const Spline * FrPi0CEx   (void) const { return fFrPi0CEx;    }
  const Spline * FrPi0Elas  (void) const { return fFrPi0Elas;   }
  const Spline * FrPi0Reac  (void) const { return fFrPi0Reac;   }
  const Spline * FrPi0Abs   (void) const { return fFrPi0Abs;    }
  const Spline * FrPReac    (void) const { return fFrPReac;     }
  const Spline * FrNReac    (void) const { return fFrNReac;     }
  const Spline * FrPiAElas  (void) const { return fFrPiAElas;   }
  const Spline * FrPiAInel  (void) const { return fFrPiAInel;   }
  const Spline * FrPiACEx   (void) const { return fFrPiACEx;    }
  const Spline * FrPiAAbs   (void) const { return fFrPiAAbs;    }
  const Spline * FrPiAPP    (void) const { return fFrPiAPP;     }
  const Spline * FrPiANPP   (void) const { return fFrPiANPP;    }
  const Spline * FrPiANNP   (void) const { return fFrPiANNP;    }
  const Spline * FrPiA4N4P  (void) const { return fFrPiA4N4P;   }
  const Spline * FrPiAPiProd(void) const { return fFrPiAPiProd; }
  const Spline * FrPAElas   (void) const { return fFrPAElas;    }
  const Spline * FrPAInel   (void) const { return fFrPAInel;    }
  const Spline * FrPAAbs    (void) const { return fFrPAAbs;     }
  const Spline * FrPAPP     (void) const { return fFrPAPP;      }
  const Spline * FrPANPP    (void) const { return fFrPANPP;     }
  const Spline * FrPANNP    (void) const { return fFrPANNP;     }
  const Spline * FrPA4N4P   (void) const { return fFrPA4N4P;    }
  const Spline * FrPAPiProd (void) const { return fFrPAPiProd;  }

private:
  INukeHadroData();
  INukeHadroData(const INukeHadroData & shx);
  ~INukeHadroData();

  void LoadXsData    (void);
  void CalcXsData    (void);
  void CalcFractions (void);

  static INukeHadroData * fInstance;

  //-- loaded cross section splines

  // SAID h+N x-section splines
  Spline * fXsPipPElas;  ///< pi+p ELAS x-section
  Spline * fXsPipPReac;  ///< pi+p REAC x-section
  Spline * fXsPipDAbs;   ///< pi+d ABS  x-section
  Spline * fXsPimPElas;  ///< pi-p ELAS x-section
  Spline * fXsPimPReac;  ///< pi-p REAC x-section
  Spline * fXsPimPCEx;   ///< pi-p CEX  x-section
  Spline * fXsPPElas;    ///< pp   ELAS x-section
  Spline * fXsPPReac;    ///< pp   REAC x-section
  Spline * fXsNPElas;    ///< np   ELAS x-section
  Spline * fXsNPReac;    ///< np   REAC x-section //?????????

  // Mashnik p+Fe x-section splines
  Spline * fXsPFeElas;   ///< pFe ELAS   x-section
  Spline * fXsPFeReac;   ///< pFe REAC   x-section
  Spline * fXsPFeP;      ///< pFe->pX    x-section
  Spline * fXsPFePP;     ///< pFe->ppX   x-section
  Spline * fXsPFeNPP;    ///< pFe->nppX  x-section
  Spline * fXsPFeNNP;    ///< pFe->nnpX  x-section
  Spline * fXsPFeNNPP;   ///< pFe->nnppX x-section
  Spline * fXsPFePim;    ///< pFe->pi-X  x-section
  Spline * fXsPFePi0;    ///< pFe->pi0X  x-section
  Spline * fXsPFePip;    ///< pFe->pi+X  x-section

  // Mashnik pi+Fe x-section splines
  Spline * fXsPiFeP;     ///< piFe->pX    x-section
  Spline * fXsPiFePP;    ///< piFe->ppX   x-section
  Spline * fXsPiFePPP;   ///< piFe->pppX  x-section
  Spline * fXsPiFeN;     ///< piFe->nX    x-section
  Spline * fXsPiFeNN;    ///< piFe->nnX   x-section
  Spline * fXsPiFeNNN;   ///< piFe->nnnX  x-section
  Spline * fXsPiFeNP;    ///< piFe->npX   x-section
  Spline * fXsPiFeNPP;   ///< piFe->nppX  x-section
  Spline * fXsPiFeNPPP;  ///< piFe->npppX x-section
  Spline * fXsPiFeNNP;   ///< piFe->nnpX  x-section
  Spline * fXsPiFeNNPP;  ///< piFe->nnppX x-section
  Spline * fXsPiFePi0;   ///< piFe->pi0X  x-section

  // Hadronic x-sections from data (Ash & Carrol)
  Spline * fXsAshPiFeAbs;  ///<
  Spline * fXsAshPiFeReac; ///<
  Spline * fXsCarPiFeTot;  ///<

  // Total x-sections needed from mean free path
  Spline * fXsPipTot;     ///<
  Spline * fXsPimTot;     ///<
  Spline * fXsPi0Tot;     ///<
  Spline * fXsPTot;       ///<
  Spline * fXsNTot;       ///<

  // Fractions of total x-sections needed for particle fates
  Spline * fFrPipCEx;    ///< fraction of CEx for pi+
  Spline * fFrPipElas;   ///< fraction of CEx+Elas for pi+
  Spline * fFrPipReac;   ///< fraction of CEx+Elas+PiProd for pi+
  Spline * fFrPipAbs;    ///< fraction of CEx+Elas+PiProd+Abs for pi+
  Spline * fFrPimCEx;    ///<
  Spline * fFrPimElas;   ///<
  Spline * fFrPimReac;   ///<
  Spline * fFrPimAbs;    ///<
  Spline * fFrPi0CEx;    ///<
  Spline * fFrPi0Elas;   ///<
  Spline * fFrPi0Reac;   ///<
  Spline * fFrPi0Abs;    ///<
  Spline * fFrPReac;     ///<
  Spline * fFrNReac;     ///<
  Spline * fFrPiAElas;   ///<
  Spline * fFrPiAInel;   ///<
  Spline * fFrPiACEx;    ///<
  Spline * fFrPiAAbs;    ///<
  Spline * fFrPiAPP;     ///<
  Spline * fFrPiANPP;    ///<
  Spline * fFrPiANNP;    ///<
  Spline * fFrPiA4N4P;   ///<
  Spline * fFrPiAPiProd; ///<
  Spline * fFrPAElas;    ///<
  Spline * fFrPAInel;    ///<
  Spline * fFrPAAbs;     ///<
  Spline * fFrPAPP;      ///<
  Spline * fFrPANPP;     ///<
  Spline * fFrPANNP;     ///<
  Spline * fFrPA4N4P;    ///<
  Spline * fFrPAPiProd;  ///<

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

