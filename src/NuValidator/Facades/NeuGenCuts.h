//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenCuts

\brief    NeuGEN's cuts

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_CUTS_H_
#define _NEUGEN_CUTS_H_

#include <string>
#include <ostream>

#include <TObject.h>

#include "Facades/NGKineVar.h"

using std::string;
using std::ostream;

namespace genie   {
namespace nuvld   {
namespace facades {

class NeuGenCuts : public TObject {

public:

  NeuGenCuts();
  NeuGenCuts(const char * name);
  NeuGenCuts(NGKineVar_t kvid, float kvmin, float kvmax,
                      bool inclusive, bool incQel, bool incRes, bool incDis);
  ~NeuGenCuts();

  void SetCut(NGKineVar_t kvid, float kvmin, float kvmax);

  void SetInclusive  (bool on);
  void SetIncludeQel (bool on);
  void SetIncludeRes (bool on);
  void SetIncludeDis (bool on);

  bool         IsCutSet   (void) const { return _cutset;    }
  NGKineVar_t  CutKVId    (void) const { return _kvid;      } 
  float        KVMin      (void) const { return _kvmin;     }
  float        KVMax      (void) const { return _kvmax;     }
  int          ProcMask   (void) const { return _procmask;  }
  bool         Inclusive  (void) const { return _inclusive; }
  bool         IncludeQel (void) const { return _incqel;    }
  bool         IncludeRes (void) const { return _incres;    }
  bool         IncludeDis (void) const { return _incdis;    }
  bool         SumQel     (void) const;
  bool         SumRes     (void) const;
  bool         SumDis     (void) const;

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NeuGenCuts & final);
  
private:

  void UpdateProcessMask (void);

  string       _name;
  bool         _cutset;
  NGKineVar_t  _kvid;
  float        _kvmin;
  float        _kvmax;
  int          _procmask;
  bool         _inclusive;
  bool         _incqel;
  bool         _incres;
  bool         _incdis;

ClassDef(NeuGenCuts, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

