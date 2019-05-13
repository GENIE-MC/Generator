//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//____________________________________________________________________________
FermiMomentumTable::FermiMomentumTable()
{
}
//____________________________________________________________________________
FermiMomentumTable::FermiMomentumTable(const FermiMomentumTable & )
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
  LOG("FermiP", pINFO)
       << "Finding Fermi momenta table entry for (tgt = "
                      << tgt_pdgc << ", nucl = " << nucleon_pdgc << ")";

  if(fKFSets.size()==0) {
      LOG("FermiP", pWARN)
         << "The Fermi momenta table is empty! Returning kf(tgt = "
                  << tgt_pdgc << ", nucl = " << nucleon_pdgc << ") = 0";
      return 0;
  }

  double kf=0;
  bool isp = pdg::IsProton(nucleon_pdgc);

  if(fKFSets.count(tgt_pdgc) == 1) {
     LOG("FermiP", pDEBUG) << "Got exact match in Fermi momenta table";
     map<int, KF_t>::const_iterator table_iter = fKFSets.find(tgt_pdgc);
     if(isp) kf = table_iter->second.p;
     else    kf = table_iter->second.n;
     LOG("FermiP", pINFO) << "kF = " << kf;
     return kf;
  }
  LOG("FermiP", pINFO) << "Couldn't find exact match in Fermi momenta table";

  int  Z   = pdg::IonPdgCodeToZ(tgt_pdgc);
  int    Ac=9999, Zc=9999, dZmin=9999;
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
  LOG("FermiP", pINFO)
       << "The closest nucleus in table is pdgc = " << pdg::IonPdgCode(Ac,Zc);
  LOG("FermiP", pINFO) << "kF = " << kf;
  return kf;
}
//____________________________________________________________________________
