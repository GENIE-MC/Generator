//____________________________________________________________________________
/*!

\class   genie::NtpMCSummary

\brief   MINOS-style Ntuple Class to hold a MC Summary Information

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#include <iomanip>

#include "EventGeneration/EventRecord.h"
#include "Interaction/Interaction.h"
#include "Ntuple/NtpMCSummary.h"
#include "Ntuple/NtpGHepAnalyzer.h"

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
NtpMCSummary::NtpMCSummary()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCSummary::NtpMCSummary(const NtpMCSummary & mcs)
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

}
//____________________________________________________________________________
void NtpMCSummary::Copy(const EventRecord & evrec)
{
  NtpGHepAnalyzer ghep_analyzer;

  ghep_analyzer.AnalyzeEventRecord(evrec);

  this->probe     = ghep_analyzer.ProbePdgCode();
  this->fsl       = ghep_analyzer.FslPdgCode();
  this->tgt       = ghep_analyzer.TgtPdgCode();
  this->nucl      = ghep_analyzer.HitNuclPdgCode();
  this->iqrk      = ghep_analyzer.HitQuarkPdgCode();
  this->fqrk      = ghep_analyzer.OutgQuarkPdgCode();
  this->res       = ghep_analyzer.ResPdgCode();
  this->ch        = ghep_analyzer.CharmHadPdgCode();
  this->Z         = ghep_analyzer.NuclTgtZ();
  this->A         = ghep_analyzer.NuclTgtA();
  this->N         = ghep_analyzer.NuclTgtN();
  this->scat      = ghep_analyzer.ScatType();
  this->proc      = ghep_analyzer.ProcType();
  this->x         = ghep_analyzer.KineX();
  this->y         = ghep_analyzer.KineY();
  this->z         = ghep_analyzer.FragmZ();
  this->Q2        = ghep_analyzer.KineQ2();
  this->W         = ghep_analyzer.KineW();
  this->xsec      = ghep_analyzer.XSec();
  this->dxsec     = ghep_analyzer.dXSec();
  this->emfrac    = ghep_analyzer.EmFrac();

  this->v[0]      = ghep_analyzer.Vtx().X();
  this->v[1]      = ghep_analyzer.Vtx().Y();
  this->v[2]      = ghep_analyzer.Vtx().Z();
  this->v[3]      = ghep_analyzer.Vtx().T();

  this->p4p[0]    = ghep_analyzer.Probe4P().Px();
  this->p4p[1]    = ghep_analyzer.Probe4P().Py();
  this->p4p[2]    = ghep_analyzer.Probe4P().Pz();
  this->p4p[3]    = ghep_analyzer.Probe4P().E();

  this->p4nucl[0] = ghep_analyzer.HitNucl4P().Px();
  this->p4nucl[1] = ghep_analyzer.HitNucl4P().Py();
  this->p4nucl[2] = ghep_analyzer.HitNucl4P().Pz();
  this->p4nucl[3] = ghep_analyzer.HitNucl4P().E();

  this->p4fsl[0]  = ghep_analyzer.Fsl4P().Px();
  this->p4fsl[1]  = ghep_analyzer.Fsl4P().Py();
  this->p4fsl[2]  = ghep_analyzer.Fsl4P().Pz();
  this->p4fsl[3]  = ghep_analyzer.Fsl4P().E();

  this->p4fsh[0]  = ghep_analyzer.HadShw4P().Px();
  this->p4fsh[1]  = ghep_analyzer.HadShw4P().Py();
  this->p4fsh[2]  = ghep_analyzer.HadShw4P().Pz();
  this->p4fsh[3]  = ghep_analyzer.HadShw4P().E();

  this->p4fsl2[0] = 0;
  this->p4fsl2[1] = 0;
  this->p4fsl2[2] = 0;
  this->p4fsl2[3] = 0;

  this->np   = ghep_analyzer.NumProton();
  this->nn   = ghep_analyzer.NumNeutron();
  this->npi0 = ghep_analyzer.NumPi0();
  this->npip = ghep_analyzer.NumPiPlus();
  this->npim = ghep_analyzer.NumPiMinus();
  this->nK0  = ghep_analyzer.NumK0();
  this->nKp  = ghep_analyzer.NumKPlus();
  this->nKm  = ghep_analyzer.NumKMinus();
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
  this->xsec      =  mcs.xsec;
  this->dxsec     =  mcs.dxsec;
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
  this->xsec      =  0.;
  this->dxsec     =  0.;
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
}
//____________________________________________________________________________
