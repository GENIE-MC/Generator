//____________________________________________________________________________
/*!

\class    genie::FermiMomentumTable

\brief    A table of Fermi momentum constants

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 18, 2005

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Messenger/Messenger.h"
#include "Nuclear/FermiMomentumTable.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
FermiMomentumTable::FermiMomentumTable()
{
}
//____________________________________________________________________________
FermiMomentumTable::FermiMomentumTable(const FermiMomentumTable & fmt)
{

}
//____________________________________________________________________________
FermiMomentumTable::~FermiMomentumTable()
{

}
//____________________________________________________________________________
void FermiMomentumTable::AddTableEntry(int tgt_pdgc, KF_t kf)
{
  fKFSets.insert(map<int, KF_t>::value_type(tgt_pdgc, kf));
}
//____________________________________________________________________________
double FermiMomentumTable::FindClosestKF(int tgt_pdgc, int nucleon_pdgc) const
{
  LOG("FermiP", pDEBUG)
       << "Finding Fermi momenta table entry for nucleus closest to pdgc = "
       << tgt_pdgc;

  if(fKFSets.size()==0) {
      LOG("FermiP", pWARN)
         << "The Fermi momenta table is empty! Returning kf(tgt = "
                  << tgt_pdgc << ", nucl = " << nucleon_pdgc << ") = 0";
      return 0;
  }

  if(fKFSets.count(tgt_pdgc) == 1) {
     LOG("FermiP", pDEBUG) << "Got exact match in Fermi momenta table";
     map<int, KF_t>::const_iterator table_iter = fKFSets.find(tgt_pdgc);
     return table_iter->second.p;
  }
  LOG("FermiP", pDEBUG) << "Couldn't find exact match in Fermi momenta table";

  int  Z   = pdg::IonPdgCodeToZ(tgt_pdgc);
  bool isp = pdg::IsProton(nucleon_pdgc);

  int    Ac=9999, Zc=9999, dZmin=9999;
  double kf=0;
  map<int, KF_t>::const_iterator kfiter;
  for(kfiter=fKFSets.begin(); kfiter!=fKFSets.end(); ++kfiter) {
    int pdgc = kfiter->first;
    int Zt = pdg::IonPdgCodeToZ(pdgc);
    int dZ = TMath::Abs(Zt-Z);
    if(dZ<dZmin) {
      dZmin = dZ;
      Zc = Zt;
      Ac = pdg::IonPdgCodeToA(pdgc);
      KF_t kft  = kfiter->second;
      if(isp) kf=kft.p;
      else    kf=kft.n;
    }
  }
  LOG("FermiP", pDEBUG)
       << "The closest nucleus in table is pdgc = " << pdg::IonPdgCode(Ac,Zc);
  return kf;
}
//____________________________________________________________________________
