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
  const Interaction * interaction = evrec.GetInteraction();

  const InitialState &     init_state = interaction->GetInitialState();
  const ProcessInfo &      proc_info  = interaction->GetProcessInfo();
  const ScatteringParams & scp        = interaction->GetScatteringParams();
  const Target &           tgt        = init_state.GetTarget();

  TLorentzVector * p4lprobe = init_state.GetProbeP4(kRfLab);
  TLorentzVector * p4probe  = init_state.GetProbeP4(kRfStruckNucAtRest);

  this->probe     = init_state.GetProbePDGCode();
  this->fsl       = 0;
  this->tgt       = tgt.PDGCode();
  this->nucl      = tgt.StruckNucleonPDGCode();
  this->iqrk      = 0;
  this->fqrk      = 0;
  this->Z         = tgt.Z();
  this->A         = tgt.A();
  this->N         = tgt.N();  
  this->scat      = proc_info.ScatteringTypeId();
  this->proc      = proc_info.InteractionTypeId();
  this->x         = (scp.Exists("x"))  ? scp.GetDouble("x")  : 0.;
  this->y         = (scp.Exists("y"))  ? scp.GetDouble("y")  : 0.;
  this->z         = 0;
  this->Q2        = (scp.Exists("Q2")) ? scp.GetDouble("Q2") : 0.;
  this->W         = (scp.Exists("W"))  ? scp.GetDouble("W")  : 0.;
  this->xsec      = 0;
  this->p4p[0]    = p4probe->Px();
  this->p4p[1]    = p4probe->Py();
  this->p4p[2]    = p4probe->Pz();
  this->p4p[3]    = p4probe->Energy();
  this->p4nucl[0] = tgt.StruckNucleonP4()->Px();
  this->p4nucl[1] = tgt.StruckNucleonP4()->Py();
  this->p4nucl[2] = tgt.StruckNucleonP4()->Pz();
  this->p4nucl[3] = tgt.StruckNucleonP4()->Energy();
  this->p4fsl[0]  =  0;
  this->p4fsl[1]  =  0;
  this->p4fsl[2]  =  0;
  this->p4fsl[3]  =  0;
  this->p4fsl2[0] =  0;
  this->p4fsl2[1] =  0;
  this->p4fsl2[2] =  0;
  this->p4fsl2[3] =  0;
  this->p4fsh[0]  =  0;
  this->p4fsh[1]  =  0;
  this->p4fsh[2]  =  0;
  this->p4fsh[3]  =  0;

  delete p4lprobe;
  delete p4probe;
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
}
//____________________________________________________________________________
