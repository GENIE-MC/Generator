//____________________________________________________________________________
/*!

\class    genie::FKR

\brief    A class holding the Feynmann-Kislinger-Ravndall (FKR) baryon
          excitation model parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _FKR_H_
#define _FKR_H_

#include <iostream>

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResParams.h"
#include "Interaction/Interaction.h"

using std::ostream;

namespace genie {

class FKR {

public:

  FKR();
  ~FKR();

  void Calculate  (const Interaction * interaction);
  void Initialize (void);
  void Print      (ostream & stream) const;

  //-- options that need to be set by the Rein-Seghal xsec algorithm

  void SetZeta      (double zeta );
  void SetOmega     (double omega);
  void SetResParams (const BaryonResParams & rp);
  void SetMa2       (double ma2);
  void SetMv2       (double mv2);

  //-- Feynmann-Kislinger-Ravndall (FKR) parameters

  const double Lamda    (void) const { return fLamda;        }
  const double Tv       (void) const { return fTv;           }
  const double Rv       (void) const { return fRv;           }
  const double S        (void) const { return fS;            }
  const double Ta       (void) const { return fTa;           }
  const double Ra       (void) const { return fRa;           }
  const double B        (void) const { return fB;            }
  const double C        (void) const { return fC;            }
  const double R        (void) const { return fR;            }
  const double T        (void) const { return fT;            }
  const double Tplus    (void) const { return fTplus;        }
  const double Tminus   (void) const { return fTminus;       }
  const double Rplus    (void) const { return fRplus;        }
  const double Rminus   (void) const { return fRminus;       }

  //-- frequently used combinations of FKR parameters
  //   ( w = - sin^2(theta-weinberg) )

  const double LRminus     (void) const { return fLamdaRminus;  }
  const double LRplus      (void) const { return fLamdaRplus;   }
  const double LTminus     (void) const { return fLamdaTminus;  }
  const double LTplus      (void) const { return fLamdaTplus;   }
  const double LC          (void) const { return fLC;           }
  const double LS          (void) const { return fLS;           }
  const double LR          (void) const { return fLR;           }
  const double LT          (void) const { return fLT;           }
  const double LC_2B       (void) const { return fLC_2B;        }
  const double LC_3B       (void) const { return fLC_3B;        }
  const double LC_5B       (void) const { return fLC_5B;        }
  const double L2          (void) const { return fLamda2;       }
  const double L2Rminus    (void) const { return fLamda2Rminus; }
  const double L2Rplus     (void) const { return fLamda2Rplus;  }
  const double L2C         (void) const { return fLamda2C;      }
  const double L2S         (void) const { return fLamda2S;      }
  const double L2R         (void) const { return fLamda2R;      }
  const double Rminus_wR   (void) const { return fRminus_wR;    }
  const double Rminus_2wR  (void) const { return fRminus_2wR;   }
  const double Rminus_3wR  (void) const { return fRminus_3wR;   }
  const double Rminus_4wR  (void) const { return fRminus_4wR;   }
  const double Rplus_wR    (void) const { return fRplus_wR;     }
  const double Rplus_2wR   (void) const { return fRplus_2wR;    }
  const double Rplus_3wR   (void) const { return fRplus_3wR;    }
  const double Rplus_4wR   (void) const { return fRplus_4wR;    }
  const double Tminus_wTv  (void) const { return fTminus_wTv;   }
  const double Tminus_2wTv (void) const { return fTminus_2wTv;  }
  const double Tminus_3wTv (void) const { return fTminus_3wTv;  }
  const double Tminus_4wTv (void) const { return fTminus_4wTv;  }
  const double Tplus_wTv   (void) const { return fTplus_wTv;    }
  const double Tplus_2wTv  (void) const { return fTplus_2wTv;   }
  const double Tplus_3wTv  (void) const { return fTplus_3wTv;   }
  const double Tplus_4wTv  (void) const { return fTplus_4wTv;   }

  friend ostream & operator<< (ostream & stream, const FKR & parameters);

private:

  double fLamda;
  double fTv;
  double fRv;
  double fS;
  double fTa;
  double fRa;
  double fB;
  double fC;
  double fR;
  double fT;
  double fTplus;
  double fTminus;
  double fRplus;
  double fRminus;

  //-- recurring FKR param cobinations

  double fLamdaRminus;
  double fLamdaRplus;
  double fLamdaTminus;
  double fLamdaTplus;
  double fLC;
  double fLS;
  double fLR;
  double fLT;
  double fLC_2B;
  double fLC_3B;
  double fLC_5B;
  double fLamda2;
  double fLamda2Rminus;
  double fLamda2Rplus;
  double fLamda2C;
  double fLamda2S;
  double fLamda2R;
  double fRminus_wR;
  double fRminus_2wR;
  double fRminus_3wR;
  double fRminus_4wR;
  double fRplus_wR;
  double fRplus_2wR;
  double fRplus_3wR;
  double fRplus_4wR;
  double fTminus_wTv;
  double fTminus_2wTv;
  double fTminus_3wTv;
  double fTminus_4wTv;
  double fTplus_wTv;
  double fTplus_2wTv;
  double fTplus_3wTv;
  double fTplus_4wTv;

  //-- options
  double fZeta;  ///< parameter Z used in computing the FKRs
  double fOmega; ///< parameter Omega used in computing the FKRs
  double fMa2;   ///< resonance Ma^2 used
  double fMv2;   ///< resonance Mv^2 used
  const BaryonResParams * fResParams;
};

}        // genie namespace

#endif   // _FKR_PARAMS_H_
