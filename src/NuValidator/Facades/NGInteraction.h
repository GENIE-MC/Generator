//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGInteraction

\brief    NeuGEN's Interaction

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _INTERACTION_H_
#define _INTERACTION_H_

#include <string>
#include <ostream>

#include <TObject.h>

#include "Facades/NGFlavor.h"
#include "Facades/NGCcNc.h"
#include "Facades/NGInitState.h"
#include "Facades/NGFinalState.h"
#include "Facades/NGNucleus.h"

using std::string;
using std::ostream;

extern "C" {
  void makestate_(int*,int*,int*,int*,int*,int*,int*,int*,int*);
  void writestate_(int*);
}

namespace genie   {
namespace nuvld   {
namespace facades {

class NGInteraction : public TObject {

public:

  friend ostream & operator << (ostream & stream, const NGInteraction & inter);

  NGInteraction();
  NGInteraction(const char *name);
  NGInteraction(const NGInteraction * inter);
  NGInteraction(NGFlavor_t f, NGNucleus_t n, NGCcNc_t c, NGInitState_t in);
  ~NGInteraction();
  
  const char * Name (void) const { return _name.c_str(); }

  void SetFlavor    (NGFlavor_t f   )  { _flavor     = f; }
  void SetNucleus   (NGNucleus_t n  )  { _nucleus    = n; }
  void SetCCNC      (NGCcNc_t c     )  { _ccnc       = c; }
  void SetInitState (NGInitState_t i)  { _init_state = i; }

  NGFlavor_t     GetFlavor    (void) const { return _flavor;     }
  NGNucleus_t    GetNucleus   (void) const { return _nucleus;    }
  NGCcNc_t       GetCCNC      (void) const { return _ccnc;       }
  NGInitState_t  GetInitState (void) const { return _init_state; }
  int            GetProcess   (void) const;
  int            GetProcess   (NGFinalState * final) const;

  void Print(ostream & stream) const;

private:

  string          _name;
  NGFlavor_t      _flavor;
  NGNucleus_t     _nucleus;
  NGCcNc_t        _ccnc;
  NGInitState_t  _init_state;

ClassDef(NGInteraction, 0)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

