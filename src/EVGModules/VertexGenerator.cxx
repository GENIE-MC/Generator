//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "EVGModules/VertexGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::utils;

//___________________________________________________________________________
VertexGenerator::VertexGenerator() :
EventRecordVisitorI("genie::VertexGenerator")
{

}
//___________________________________________________________________________
VertexGenerator::VertexGenerator(string config) :
EventRecordVisitorI("genie::VertexGenerator", config)
{

}
//___________________________________________________________________________
VertexGenerator::~VertexGenerator()
{

}
//___________________________________________________________________________
void VertexGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// generate a vtx and set it to all GHEP physical particles

  RandomGen * rnd = RandomGen::Instance();
  TVector3 vtx(999999.,999999.,999999.);

  GHepParticle * nucltgt = evrec->TargetNucleus();

  if(!nucltgt) {
    vtx.SetXYZ(0.,0.,0.);
  } 
  else {
    double A = nucltgt->A();
    double R = fR0 * TMath::Power(A, 1./3.);

    while(vtx.Mag() > R) {
      vtx.SetX(-R + 2*R * rnd->RndFsi().Rndm());
      vtx.SetY(-R + 2*R * rnd->RndFsi().Rndm());
      vtx.SetZ(-R + 2*R * rnd->RndFsi().Rndm());
    }
  }

  LOG("Vtx", pNOTICE) 
     << "Generated vtx @ r = " << vtx.Mag() << " fm / " 
                                          << print::Vec3AsString(&vtx);
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(p->IsFake()   ) continue;
    if(p->IsNucleus()) continue;
    LOG("Vtx", pNOTICE) << "Setting vertex position for: " << p->Name();
    p->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);
  }
}
//___________________________________________________________________________
void VertexGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void VertexGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void VertexGenerator::LoadConfig(void)
{
 AlgConfigPool * confp = AlgConfigPool::Instance();
 const Registry * gc = confp->GlobalParameterList();

 fR0 = fConfig->GetDoubleDef ("R0", gc->GetDouble("NUCL-R0")); // fm
}
//____________________________________________________________________________

