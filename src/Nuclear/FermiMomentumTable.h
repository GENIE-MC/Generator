//____________________________________________________________________________
/*!

\class    genie::FermiMomentumTable

\brief    A table of Fermi momentum constants

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 18, 2005

*/
//____________________________________________________________________________

#ifndef _FERMI_MOMENTUM_TABLE_H_
#define _FERMI_MOMENTUM_TABLE_H_

#include <map>

using std::map;

namespace genie {

typedef struct EKF_t {
  double p; // Fermi momentum for proton
  double n; // Fermi momentum for neutron
} KF_t;

class FermiMomentumTable
{
public:
  FermiMomentumTable();
  FermiMomentumTable(const FermiMomentumTable & fmt);
  virtual ~FermiMomentumTable();

  double FindClosestKF (int target_pdgc, int nucleon_pdgc) const;
  void   AddTableEntry (int target_pdgc, KF_t kf);

private:
  map<int, KF_t> fKFSets; // the actual Fermi momenta table
};

}      // genie namespace

#endif // _FERMI_MOMENTUM_TABLE_H_
