//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenConfig

\brief    Encapsulation of NeuGEN's Configuration Card

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_CONFIG_H_
#define _NEUGEN_CONFIG_H_

#include <string>
#include <ostream>

#include <TObject.h>

#include "Facades/NGInitState.h"

using std::string;
using std::ostream;

using namespace genie::nuvld::facades;

const unsigned int kNMlt  = 2; ///< n. hadronic multiplicities
const unsigned int kNInSt = 4; ///< number of initial states (vp,vn,vbp,vbn)

//! Default (best) NeuGEN parameters

const int    kNGDefPdfGrp   =  5;
const int    kNGDefPdfSet   = 12;
const float  kNGDefMaQel    =  1.032;
const float  kNGDefMaRes    =  1.032;
const float  kNGDefMaCoh    =  1.000;    
const float  kNGDefFa0Qel   = -1.260;   
const float  kNGDefEtaQel   =  0.120;   
const float  kNGDefOmegaRes =  1.050; 
const float  kNGDefZRes     =  0.750;     
const float  kNGDefR0Coh    =  1.000;    
const float  kNGDefREICoh   =  0.300;   
const float  kNGDefKnoB     =  1.300;

const float kNGDefKnoA[kNInSt] = {
                 /* v+p   v+n   vb+p  vb+n */
                    0.50, 0.00, 0.20, 0.20
};
const float kNGDefDisRes[kNMlt][kNInSt] = {
                 /* v+p   v+n   vb+p  vb+n */
                   {0.20, 0.20, 0.20, 0.20},    /* multiplicity = 2 */
                   {1.00, 1.00, 1.00, 1.00}     /* multiplicity = 3 */
};

namespace genie   {
namespace nuvld   {
namespace facades {

class NeuGenConfig : public TObject
{

public:

  NeuGenConfig();
  NeuGenConfig(const char * name);
  NeuGenConfig(const NeuGenConfig * config);
  ~NeuGenConfig();

  string Name     (void) const { return fName;     }
  int    PdfGroup (void) const { return fPdfGrp;   }
  int    PdfSet   (void) const { return fPdfSet;   }
  float  MaQel    (void) const { return fMaQel;    }
  float  MaRes    (void) const { return fMaRes;    }
  float  MaCoh    (void) const { return fMaCoh;    }
  float  Fa0Qel   (void) const { return fFa0Qel;   }
  float  EtaQel   (void) const { return fEtaQel;   }
  float  OmegaRes (void) const { return fOmegaRes; }
  float  ZRes     (void) const { return fZRes;     }
  float  R0Coh    (void) const { return fR0Coh;    }
  float  REICoh   (void) const { return fREICoh;   }
  float  KnoB     (void) const { return fKnoB;     }  
  float  KnoA     (NGInitState_t init) const;
  float  DisRes   (unsigned int multiplicity, NGInitState_t init) const;

  void SetBestParameters(void);

  void SetPdfGroup (int   pdf_grp  ) { fPdfGrp   = pdf_grp;   }
  void SetPdfSet   (int   pdf_set  ) { fPdfSet   = pdf_set;   }
  void SetMaQel    (float ma_qel   ) { fMaQel    = ma_qel;    }
  void SetMaRes    (float ma_res   ) { fMaRes    = ma_res;    }
  void SetMaCoh    (float ma_coh   ) { fMaCoh    = ma_coh;    }
  void SetFa0Qel   (float fa0_qel  ) { fFa0Qel   = fa0_qel;   }
  void SetEtaQel   (float eta_qel  ) { fEtaQel   = eta_qel;   }
  void SetOmegaRes (float omega_res) { fOmegaRes = omega_res; }
  void SetZRes     (float z_res    ) { fZRes     = z_res;     }
  void SetR0Coh    (float r0_coh   ) { fR0Coh    = r0_coh;    }
  void SetREICoh   (float rei_coh  ) { fREICoh   = rei_coh;   }
  void SetKnoB     (float kno_b    ) { fKnoB     = kno_b;     }
  void SetKnoA     (NGInitState_t init, float kno);    
  void SetDisRes   (unsigned int multiplicity, NGInitState_t init, float dis_res);

  void Copy  (const NeuGenConfig * config);
  void Print (ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NeuGenConfig & conf);

private:

  int   Multiplicity2IPos   (unsigned int multiplicity) const;
  int   InitState2IPos      (NGInitState_t init)        const;
  bool  IsValidMultiplicity (unsigned int multiplicity) const;
  bool  IsValidInitState    (NGInitState_t init)        const;

  string fName;                   ///< name of NeuGEN configuration set
  int    fPdfGrp;                 ///< PDFLIB Group
  int    fPdfSet;                 ///< PDFLIB Set
  float  fMaQel;                  ///< Axial mass for QEL scattering
  float  fMaRes;                  ///< Axial mass for RES scattering
  float  fMaCoh;                  ///< Axial mass for COH scattering
  float  fFa0Qel;                 ///< QEL axial form factor at Q^2 = 0
  float  fEtaQel;                 ///< Elastic scattering parameter
  float  fOmegaRes;               ///< R-S Model parameter Omega
  float  fZRes;                   ///< R-S Model parameter Z
  float  fR0Coh;                  ///< Nuclear size scale param in COH scattering
  float  fREICoh;                 ///< Re/Im for pion scattering amplitude  
  float  fKnoB;                   ///< KNO Hadronization parameter B,  <n> = A + B*ln(W^2)
  float  fKnoA[kNInSt];           ///< KNO Hadronization parameters A, B,  <n> = A + B*ln(W^2)
  float  fDisRes[kNMlt][kNInSt];  ///< DIS/RES tuning factors
  
ClassDef(NeuGenConfig, 1)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif
