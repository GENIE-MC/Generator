//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <iomanip>

#include "EVGCore/EventRecord.h"
#include "Interaction/Interaction.h"
#include "Ntuple/NtpMCSummary.h"
#include "GHEP/GHepSummaryBuilder.h"

using namespace genie;

ClassImp(NtpMCSummary)

using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCSummary & ntpp)
  {
     ntpp.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCSummary::NtpMCSummary() :
TObject()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCSummary::NtpMCSummary(const NtpMCSummary & mcs) :
TObject()
{
  this->Copy(mcs);
}
//____________________________________________________________________________
NtpMCSummary::~NtpMCSummary()
{

}
//____________________________________________________________________________
void NtpMCSummary::PrintToStream(ostream & stream) const
{
  stream << endl;
  stream << "Probe PDG code:............... " << this->probe  << endl;
  stream << "F/S primary lepton PDG code:.. " << this->fsl    << endl;
  stream << "Nuclear target PDG code:......." << this->tgt    << endl;
  stream << "Hit nucleon PDG code:.........." << this->nucl   << endl;
  stream << "Hit quark PDG code:............" << this->iqrk   << endl;
  stream << "Outgoing quark PDG code:......." << this->fqrk   << endl;
  stream << "Resonance PDG code:............" << this->res    << endl;
  stream << "Charm hadron PDG code:........." << this->ch     << endl;
  stream << "Target Z:......................" << this->Z      << endl;
  stream << "Target A:......................" << this->A      << endl;
  stream << "Target N:......................" << this->N      << endl;
  stream << "Scattering Type:..............." << this->scat   << endl;
  stream << "Interaction Type:.............." << this->proc   << endl;
  stream << "Bjorken x:....................." << this->x      << endl;
  stream << "Inelasticity y:................" << this->y      << endl;
  stream << "Fragmentation parameter z:....." << this->z      << endl;
  stream << "Momentum transfer Q^2:........." << this->Q2     << endl;
  stream << "Invariant mass W:.............." << this->W      << endl;
  stream << "E/M fraction:.................:" << this->emfrac << endl;

  stream << "Vertex (x,y,z,t)..............: ("
         << this->v[0] << ", " << this->v[1] << ", "
         << this->v[2] << ", " << this->v[3] << ")"
         << endl;
  stream << "Probe 4P (px,py,pz,E).........: ("
         << this->p4p[0] << ", " << this->p4p[1] << ", "
         << this->p4p[2] << ", " << this->p4p[3] << ")"
         << endl;
  stream << "Hit nucleon 4P (px,py,pz,E)...: ("
         << this->p4nucl[0] << ", " << this->p4nucl[1] << ", "
         << this->p4nucl[2] << ", " << this->p4nucl[3] << ")"
         << endl;
  stream << "F/S lepton 4P (px,py,pz,E)...: ("
         << this->p4fsl[0] << ", " << this->p4fsl[1] << ", "
         << this->p4fsl[2] << ", " << this->p4fsl[3] << ")"
         << endl;
  stream << "Num (p,n)....................: ("
         << this->np << ", " << this->nn << ")"
         << endl;
  stream << "Num (pi^0,pi^+,pi^-).........: ("
         << this->npi0 << ", " << this->npip << ", " << this->npim << ")"
         << endl;
  stream << "Num (K^0,K^+,K^-)............: ("
         << this->nK0 << ", " << this->nKp << ", " << this->nKm << ")"
         << endl;

  stream << "Cross Section (E):............." << this->xsec  << endl;
  stream << "Cross Section (E, Kinematics).:" << this->dxsec << endl;
  stream << "Event weight..................:" << this->wght  << endl;
}
//____________________________________________________________________________
void NtpMCSummary::Copy(const EventRecord & evrec)
{
  GHepSummaryBuilder sbld;

  sbld.AnalyzeEventRecord(evrec);

  this->probe     = sbld.ProbePdgC();
  this->fsl       = sbld.FslPdgC();
  this->tgt       = sbld.TgtPdgC();
  this->nucl      = sbld.HitNuclPdgC();
  this->iqrk      = sbld.HitQuarkPdgC();
  this->fqrk      = sbld.OutQuarkPdgC();
  this->res       = sbld.ResPdgC();
  this->ch        = sbld.CharmHadPdgC();
  this->Z         = sbld.NuclTgtZ();
  this->A         = sbld.NuclTgtA();
  this->N         = sbld.NuclTgtN();
  this->scat      = int(sbld.ScatType());
  this->proc      = int(sbld.ProcType());
  this->x         = sbld.KineX();
  this->y         = sbld.KineY();
  this->z         = sbld.FragmZ();
  this->Q2        = sbld.KineQ2();
  this->W         = sbld.KineW();
  this->emfrac    = sbld.EmFrac();

  this->v[0]      = sbld.Vtx().X();
  this->v[1]      = sbld.Vtx().Y();
  this->v[2]      = sbld.Vtx().Z();
  this->v[3]      = sbld.Vtx().T();

  this->p4p[0]    = sbld.Probe4P().Px();
  this->p4p[1]    = sbld.Probe4P().Py();
  this->p4p[2]    = sbld.Probe4P().Pz();
  this->p4p[3]    = sbld.Probe4P().E();

  this->p4nucl[0] = sbld.HitNucl4P().Px();
  this->p4nucl[1] = sbld.HitNucl4P().Py();
  this->p4nucl[2] = sbld.HitNucl4P().Pz();
  this->p4nucl[3] = sbld.HitNucl4P().E();

  this->p4fsl[0]  = sbld.Fsl4P().Px();
  this->p4fsl[1]  = sbld.Fsl4P().Py();
  this->p4fsl[2]  = sbld.Fsl4P().Pz();
  this->p4fsl[3]  = sbld.Fsl4P().E();

  this->p4fsh[0]  = sbld.HadShw4P().Px();
  this->p4fsh[1]  = sbld.HadShw4P().Py();
  this->p4fsh[2]  = sbld.HadShw4P().Pz();
  this->p4fsh[3]  = sbld.HadShw4P().E();

  this->p4fsl2[0] = 0;
  this->p4fsl2[1] = 0;
  this->p4fsl2[2] = 0;
  this->p4fsl2[3] = 0;

  this->np   = sbld.NProton();
  this->nn   = sbld.NNeutron();
  this->npi0 = sbld.NPi0();
  this->npip = sbld.NPiPlus();
  this->npim = sbld.NPiMinus();
  this->nK0  = sbld.NK0();
  this->nKp  = sbld.NKPlus();
  this->nKm  = sbld.NKMinus();

  this->xsec  = evrec.XSec();
  this->dxsec = evrec.DiffXSec();
  this->wght  = evrec.Weight();
}
//____________________________________________________________________________
void NtpMCSummary::Copy(const NtpMCSummary & mcs)
{
  this->probe     =  mcs.probe;
  this->fsl       =  mcs.fsl;
  this->tgt       =  mcs.tgt;
  this->nucl      =  mcs.nucl;
  this->iqrk      =  mcs.iqrk;
  this->fqrk      =  mcs.fqrk;
  this->res       =  mcs.res;
  this->ch        =  mcs.ch;
  this->Z         =  mcs.Z;
  this->A         =  mcs.A;
  this->N         =  mcs.N;
  this->scat      =  mcs.scat;
  this->proc      =  mcs.proc;
  this->x         =  mcs.x;
  this->y         =  mcs.y;
  this->z         =  mcs.z;
  this->Q2        =  mcs.Q2;
  this->W         =  mcs.W;
  this->emfrac    =  mcs.emfrac;
  this->p4p[0]    =  mcs.p4p[0];
  this->p4p[1]    =  mcs.p4p[1];
  this->p4p[2]    =  mcs.p4p[2];
  this->p4p[3]    =  mcs.p4p[3];
  this->p4nucl[0] =  mcs.p4nucl[0];
  this->p4nucl[1] =  mcs.p4nucl[1];
  this->p4nucl[2] =  mcs.p4nucl[2];
  this->p4nucl[3] =  mcs.p4nucl[3];
  this->p4fsl[0]  =  mcs.p4fsl[0];
  this->p4fsl[1]  =  mcs.p4fsl[1];
  this->p4fsl[2]  =  mcs.p4fsl[2];
  this->p4fsl[3]  =  mcs.p4fsl[3];
  this->p4fsl2[0] =  mcs.p4fsl2[0];
  this->p4fsl2[1] =  mcs.p4fsl2[1];
  this->p4fsl2[2] =  mcs.p4fsl2[2];
  this->p4fsl2[3] =  mcs.p4fsl2[3];
  this->p4fsh[0]  =  mcs.p4fsh[0];
  this->p4fsh[1]  =  mcs.p4fsh[1];
  this->p4fsh[2]  =  mcs.p4fsh[2];
  this->p4fsh[3]  =  mcs.p4fsh[3];
  this->np        =  mcs.np;
  this->nn        =  mcs.nn;
  this->npi0      =  mcs.npi0;
  this->npip      =  mcs.npip;
  this->npim      =  mcs.npim;
  this->nK0       =  mcs.nK0;
  this->nKp       =  mcs.nKp;
  this->nKm       =  mcs.nKm;
  this->xsec      =  mcs.xsec;
  this->dxsec     =  mcs.dxsec;
  this->wght      =  mcs.wght;
}
//____________________________________________________________________________
void NtpMCSummary::Init(void)
{
  this->probe     =  0;
  this->fsl       =  0;
  this->tgt       =  0;
  this->nucl      =  0;
  this->iqrk      =  0;
  this->fqrk      =  0;
  this->res       =  0;
  this->ch        =  0;
  this->Z         =  0;
  this->A         =  0;
  this->N         =  0;
  this->scat      =  0;
  this->proc      =  0;
  this->x         =  0.;
  this->y         =  0.;
  this->z         =  0.;
  this->Q2        =  0.;
  this->W         =  0.;
  this->emfrac    =  0.;
  this->p4p[0]    =  0.;
  this->p4p[1]    =  0.;
  this->p4p[2]    =  0.;
  this->p4p[3]    =  0.;
  this->p4nucl[0] =  0.;
  this->p4nucl[1] =  0.;
  this->p4nucl[2] =  0.;
  this->p4nucl[3] =  0.;
  this->p4fsl[0]  =  0.;
  this->p4fsl[1]  =  0.;
  this->p4fsl[2]  =  0.;
  this->p4fsl[3]  =  0.;
  this->p4fsl2[0] =  0.;
  this->p4fsl2[1] =  0.;
  this->p4fsl2[2] =  0.;
  this->p4fsl2[3] =  0.;
  this->p4fsh[0]  =  0.;
  this->p4fsh[1]  =  0.;
  this->p4fsh[2]  =  0.;
  this->p4fsh[3]  =  0.;
  this->np        =  0;
  this->nn        =  0;
  this->npi0      =  0;
  this->npip      =  0;
  this->npim      =  0;
  this->nK0       =  0;
  this->nKp       =  0;
  this->nKm       =  0;
  this->xsec      =  0.;
  this->dxsec     =  0.;
  this->wght      =  0.;
}
//____________________________________________________________________________
