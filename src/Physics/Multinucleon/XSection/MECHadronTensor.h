//____________________________________________________________________________
/*!

\class    genie::MECHadronTensor

\brief    Singleton class to load and store MEC hadron tensor tables,
          to aid in the implementation (and improve the CPU efficiency of)
          MEC cross-section models.

\author   Code contributed by Jackie Schwehr
          Substantial refactorization by the core GENIE group.

\ref      Hadron tensors used here are those computed by the following models:
          
          J. Nieves, I. Ruiz Simo, M.J. Vicente Vacas,
          Inclusive quasi-elastic neutrino reactions, PRC 83 (2011) 045501

\created  September 12, 2014

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_HADRON_TENSOR_H_
#define _MEC_HADRON_TENSOR_H_

#include <map>
#include <vector>
#include <string>

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#include "Framework/Numerical/BLI2D.h"

using std::map;
using std::vector;
using std::string;

namespace genie {

class MECHadronTensor
{
public:

  // ................................................................
  // MEC hadron tensor type
  //

  typedef enum EMECHadronTensorType {  
    kMHTUndefined = -1,
    kMHTValenciaFullAll, 
    kMHTValenciaFullpn, 
    kMHTValenciaDeltaAll, 
    kMHTValenciaDeltapn
  } 
  MECHadronTensorType_t;

  // ................................................................
  // MEC hadron tensor table
  //

  class MECHadronTensorTable 
  {
  public:
     MECHadronTensorTable() { }
    ~MECHadronTensorTable() { /* note: should delete the grids! */ }
     map<MECHadronTensor::MECHadronTensorType_t, vector<genie::BLI2DNonUnifGrid *> > Table;
  };

  // ................................................................

  static MECHadronTensor * Instance();

  // method to return whether the targetpdg is in fKnownTargets
  bool KnownTarget(int targetpdg);

  // method to return whether the targetpdg is in fKnownTensors
  bool KnownTensor(int targetpdg);

  // method to access a specific set of tables
  const vector<genie::BLI2DNonUnifGrid *> &
     TensorTable(int targetpdg, MECHadronTensor::MECHadronTensorType_t type);

private:

  // Ctors & dtor
  MECHadronTensor();
  MECHadronTensor(const MECHadronTensor &);
 ~MECHadronTensor();

  // Self
  static MECHadronTensor * fgInstance;

  // Load available hadron tensor tables.
  // NOTES: All tables are loaded by the ctor - Consider loading them on demand.
  //        This will also need to be extended to load tensors for requested model.
  void LoadTensorTables(int targetpdg);

  // This map holds all known tensor tables (target PDG code is the key)
  std::map<int, MECHadronTensorTable> fTargetTensorTables;

  // List of targets for which we can provide a calculation
  // some known targets use scale from the tensor table from another target.
  std::vector<int> fKnownTensors;

  // that array size 14400 is 120x120.  Could go dynamic in the future, or not.
  void ReadHadTensorqzq0File(string filename, int nwpoints, int nqzpoints, int nq0points, double hadtensor_w_array[][14400]);

  // singleton cleaner
  struct Cleaner {
    void DummyMethodAndSilentCompiler(){}
    ~Cleaner(){
      if (MECHadronTensor::fgInstance !=0){
	delete MECHadronTensor::fgInstance;
	MECHadronTensor::fgInstance = 0;
      }
    }
  };
  friend struct Cleaner;
};

} // genie namespace

#endif // _MEC_HARDRON_TENSOR_H_
