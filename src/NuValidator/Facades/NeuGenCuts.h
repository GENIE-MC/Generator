//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NeuGenCuts

\brief    NeuGEN's cuts

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
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
                   int procmask, bool sumQel, bool sumRes, bool sumDis);
  ~NeuGenCuts();

  void SetCut(NGKineVar_t kvid, float kvmin, float kvmax);

  void SetProcmask (int in)   { _procmask = in;  }
  void SetSumQel   (bool sum) { _sumqel   = sum; }
  void SetSumRes   (bool sum) { _sumres   = sum; }
  void SetSumDis   (bool sum) { _sumdis   = sum; }

  bool         IsCutSet (void) const { return _cutset;   }
  NGKineVar_t  CutKVId  (void) const { return _kvid;     } 
  float        KVMin    (void) const { return _kvmin;    }
  float        KVMax    (void) const { return _kvmax;    }
  int          ProcMask (void) const { return _procmask; }
  bool         SumQel   (void) const { return _sumqel;   }
  bool         SumRes   (void) const { return _sumres;   }
  bool         SumDis   (void) const { return _sumdis;   }

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NeuGenCuts & final);
  
private:

  string       _name;
  bool         _cutset;
  NGKineVar_t  _kvid;
  float        _kvmin;
  float        _kvmax;
  int          _procmask;
  bool         _sumqel;
  bool         _sumres;
  bool         _sumdis;

ClassDef(NeuGenCuts, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

