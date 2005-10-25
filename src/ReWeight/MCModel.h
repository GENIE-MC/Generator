//____________________________________________________________________________
/*!

\class   genie::MCModel

\brief   A collection of cross section models

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 22, 2005

*/
//____________________________________________________________________________

#ifndef _MC_MODEL_H_
#define _MC_MODEL_H_

#include <string>
#include <map>

#include "Algorithm/AlgId.h"
#include "Interaction/InitialState.h"
#include "Interaction/ProcessInfo.h"

using std::string;
using std::map;

namespace genie {

class XSecAlgorithmI;
class Interaction;

class MCModel {

public :

  MCModel();
  MCModel(string name);
  MCModel(const MCModel & model);
  ~MCModel();

  void UseXSecAlg(const ProcessInfo & proc, const AlgId & algid);
  void UseXSecAlg(const ProcessInfo & proc, const InitialState & init, const AlgId & algid);

  const XSecAlgorithmI * XSecAlg(const Interaction * interaction) const;

  void Copy(const MCModel & model);

private :

  void   Initialize();
  string BuildKey(const ProcessInfo & proc) const;
  string BuildKey(const ProcessInfo & proc, const InitialState & init) const;

  string fName;

  map<string, const XSecAlgorithmI *> fXSecModelList;
};

}      // genie namespace

#endif // _MC_MODEL_H_
