//_____________________________________________________________________________
/*!

\class    genie::nuvld::facades::NGFinalState

\brief    NeuGEN's Final State information

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>
          Hugh Gallagher      (Tufts University) <gallag@minos.phy.tufts.edu>

\created  August 2004
*/
//_____________________________________________________________________________

#ifndef _NG_FINAL_STATE_H_
#define _NG_FINAL_STATE_H_

#include <string>
#include <ostream>

#include <TObject.h>

using std::string;
using std::ostream;

namespace genie   {
namespace nuvld   {
namespace facades {


class NGFinalState : public TObject {

public:

  NGFinalState();
  NGFinalState(const char * name);
  NGFinalState(int proton, int neutron, int piplus, int piminus, int pizero);
  ~NGFinalState();

  void SetFinalState(int proton, int neutron, int piplus, int piminus, int pizero);

  int GetProton  (void)  const { return _proton;  }
  int GetNeutron (void)  const { return _neutron; }
  int GetPiplus  (void)  const { return _piplus;  }
  int GetPiminus (void)  const { return _piminus; }
  int GetPizero  (void)  const { return _pizero;  }

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const NGFinalState & final);

private:

  string _name;
  int    _proton;
  int    _neutron;
  int    _piplus;
  int    _piminus;
  int    _pizero;

ClassDef(NGFinalState, 1)
};

} // facades namespace
} // nuvld   namespace
} // genie   namespace

#endif

