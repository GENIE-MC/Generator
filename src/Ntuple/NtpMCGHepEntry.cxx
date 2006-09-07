//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 01, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <iomanip>

#include "Ntuple/NtpMCGHepEntry.h"
#include "GHEP/GHepParticle.h"

using namespace genie;

ClassImp(NtpMCGHepEntry)

using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream& stream, const NtpMCGHepEntry & ntpp)
  {
     ntpp.PrintToStream(stream);
     return stream;
  }
}
//____________________________________________________________________________
NtpMCGHepEntry::NtpMCGHepEntry()
{
  this->Init();
}
//____________________________________________________________________________
NtpMCGHepEntry::NtpMCGHepEntry(const NtpMCGHepEntry & particle)
{
  this->Copy(particle);
}
//____________________________________________________________________________
NtpMCGHepEntry::~NtpMCGHepEntry()
{

}
//____________________________________________________________________________
void NtpMCGHepEntry::Copy(const GHepParticle & particle)
{
  this->idx    = -1;
  this->ist    = particle.Status();
  this->pdgc   = particle.Pdg();
  this->mom[0] = particle.FirstMother();
  this->mom[1] = particle.LastMother();
  this->dtr[0] = particle.FirstDaughter();
  this->dtr[1] = particle.LastDaughter();
  this->mass   = particle.Mass();
  this->p4[0]  = particle.Px();
  this->p4[1]  = particle.Py();
  this->p4[2]  = particle.Pz();
  this->p4[3]  = particle.E();
  this->v4[0]  = particle.Vx();
  this->v4[1]  = particle.Vy();
  this->v4[2]  = particle.Vz();
  this->v4[3]  = particle.Vt();
}
//____________________________________________________________________________
void NtpMCGHepEntry::Copy(const NtpMCGHepEntry & particle)
{
  this->idx    = particle.idx;
  this->ist    = particle.ist;
  this->pdgc   = particle.pdgc;
  this->mom[0] = particle.mom[0];
  this->mom[1] = particle.mom[1];
  this->dtr[0] = particle.dtr[0];
  this->dtr[1] = particle.dtr[1];
  this->mass   = particle.mass;
  this->p4[0]  = particle.p4[0];
  this->p4[1]  = particle.p4[1];
  this->p4[2]  = particle.p4[2];
  this->p4[3]  = particle.p4[3];
  this->v4[0]  = particle.v4[0];
  this->v4[1]  = particle.v4[1];
  this->v4[2]  = particle.v4[2];
  this->v4[3]  = particle.v4[3];
}
//____________________________________________________________________________
void NtpMCGHepEntry::PrintToStream(ostream & stream) const
{
  stream << "|";
  stream << setfill(' ') << setw(3)  << this->idx    << " | ";
  stream << setfill(' ') << setw(3)  << this->ist    << " | ";
  stream << setfill(' ') << setw(11) << this->pdgc   << " | ";
  stream << setfill(' ') << setw(3)  << this->mom[0] << " | ";
  stream << setfill(' ') << setw(3)  << this->dtr[0] << " | ";
  stream << setfill(' ') << setw(3)  << this->dtr[1] << " | ";

  stream << setiosflags(ios::fixed)  << setprecision(3);

  stream << setfill(' ') << setw(7)  << this->p4[0]  << " | ";
  stream << setfill(' ') << setw(7)  << this->p4[1]  << " | ";
  stream << setfill(' ') << setw(7)  << this->p4[2]  << " | ";
  stream << setfill(' ') << setw(7)  << this->p4[3]  << " | ";
}
//____________________________________________________________________________
void NtpMCGHepEntry::Init(void)
{
  this->idx    = -1;
  this->ist    = -1;
  this->pdgc   = -1;
  this->mom[0] = -1;
  this->mom[1] = -1;
  this->dtr[0] = -1;
  this->dtr[0] = -1;
  this->mass   = -1.;
  this->p4[0]  =  0.;
  this->p4[1]  =  0.;
  this->p4[2]  =  0.;
  this->p4[3]  =  0.;
  this->v4[0]  =  0.;
  this->v4[1]  =  0.;
  this->v4[2]  =  0.;
  this->v4[3]  =  0.;
}
//____________________________________________________________________________
