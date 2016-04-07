//____________________________________________________________________________
/*!

\class    genie::INukeHadroData

\brief    Singleton class to load & save mec hadron tensor tables - heavily
          based on MECLoadXSecFiles (which is based on INukeHadroData.cxx)

\author   Jackie Schwehr

\created  September 12, 2014

\notes    January 14, 2015 - Generalized for all neutrino flavors, Oxygen and Carbon


*/
//_____________________________________________________________________

#ifndef _MEC_LOAD_HAD_TENSOR_H_
#define _MEC_LOAD_HAD_TENSOR_H_

#include "Numerical/BLI2D.h"
#include "GHEP/GHepParticle.h"

#include <map>

namespace genie {

class MECLoadHadTensor
{
 public:
  static MECLoadHadTensor * Instance();
  
  // public functions

  // the main call to get a cross section.
  // need target nucleus, neutrino, neutrino energy, desired muon KE and angle.
  // returns a four element array with cross sections for all and pn fraction
  // default is to also get information to tag the delta (its XS and pn fraction)
  // but the spline call to this function turns that off.
  void XSecFour(int targetpdg, int nupdg, double Enu, double Tmu, double CosTheta, double* FourXSec, bool delta = true);

  // integrate to get the total cross section, needed for splines
  double TotalXsecAtE(int targetpdg, int nupdg, double Enu);

  // kinematic calculation.  does this exist elsewhere in GENIE in this form?
  // give q0q3 (and Enu and lepton mass) return muon KE and costheta, and jacobian area.
  double GetTmuCostFromq0q3(double dq0, double dq3, double Enu, double lmass, double &tmu, double &cost, double &area);

  // Public Variables

  // The tensor tables for one target
  // four hadron tensors in a map, with an enum as the key.
  struct HadTensorTables
  {
    enum EHadTensorType { kFullAll=0, kFullpn, kDeltaAll, kDeltapn, kNHadTensorTypes };
    std::map<EHadTensorType, std::vector<genie::BLI2DNonUnifGrid *> > fTables;
  };

  // This map holds all the different target nuclei
  // the pdg code is the key, and the above struct is the data.
  std::map<int, HadTensorTables> fTargetTensorTables;

  // list of targets for which we can provide a calculation
  // some known targets use scale from the tensor table from another target.
  std::vector<int> fKnownTensors;

  // method to return whether the targetpdg is in fKnownTargets
  bool KnownTarget(int targetpdg);

  // method to return whether the targetpdg is in fKnownTensors
  bool KnownTensor(int targetpdg);

  // utility that encodes the Qvalues for the kinematic calculation.
  // this is used in the code that contracts the hadron tensor with the lepton tensor
  double Qvalue(int targetpdg, int nupdg);

 private:

  // These methods return one cross section directly from one table
  // They don't know how to scale from table to other nuclei, so keep them private.
  // These wrap a call the XSec function for each of the four tensor types.
  double XSecFullAll(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta);
  double XSecFullpn(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta);
  double XSecDeltaAll(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta);
  double XSecDeltapn(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double CosTheta);

  // pull a single cross section out of the table, contract it with the lepton tensor.
  // this implements Juan's and Manolo's code, but is probably mostly generic actually.
  double XSec(int targetpdg, int tensorpdg, int nupdg, double Enu, double Tmu, double Costheta,  vector <genie::BLI2DNonUnifGrid *> HadTensor);

  
  // private offical stuff
  MECLoadHadTensor();
  MECLoadHadTensor(const MECLoadHadTensor & shx);
  ~MECLoadHadTensor();
  static MECLoadHadTensor * fInstance;

  // private functions

  // the constructor does this for each table.
  // (could load on demand...)
  void LoadTensorTables(int targetpdg);

  // that array size 14400 is 120x120.  Could go dynamic in the future, or not.
  void ReadHadTensorqzq0File(string filename, int nwpoints, int nqzpoints, int nq0points, double hadtensor_w_array[][14400]);



  // private varables


  // singleton cleaner
  struct Cleaner {
    void DummyMethodAndSilentCompiler(){}
    ~Cleaner(){
      if (MECLoadHadTensor::fInstance !=0){
	delete MECLoadHadTensor::fInstance;
	MECLoadHadTensor::fInstance = 0;
      }
    }
  };
  friend struct Cleaner;
};

} // genie namespace

#endif // _MEC_LOAD_HAD_TENSOR_H_
