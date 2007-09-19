//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Sep 19, 2007 - CA
   Added code to generate the interaction according to a realistic nuclear
   density and made that the default setting.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/VertexGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::controls;
using namespace genie::constants;

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

  bool uniform   = fVtxGenMethod==0;
  bool realistic = !uniform;;

  if(!nucltgt) {
    vtx.SetXYZ(0.,0.,0.);
  } 
  else {
    double A = nucltgt->A();
    double R = fR0 * TMath::Power(A, 1./3.);

    // Generate the vertex uniformly within the nucleus
    //
    if(uniform) {
       LOG("Vtx", pINFO) 
             << "Generating intranuclear vertex uniformly in volume";
       while(vtx.Mag() > R) {
           vtx.SetX(-R + 2*R * rnd->RndFsi().Rndm());
           vtx.SetY(-R + 2*R * rnd->RndFsi().Rndm());
           vtx.SetZ(-R + 2*R * rnd->RndFsi().Rndm());
       }
    } 

    // Generate the vertex using a realistic nuclear density
    //
    if(realistic) {
       LOG("Vtx", pINFO) 
           << "Generating vertex according to a realistic nuclear density profile";
        // get inputs to the rejection method
        double ymax = -1;
        double rmax = 3*R;
        double dr   = R/40.;
        for(double r = 0; r < rmax; r+=dr) { 
           ymax = TMath::Max(ymax, r*r * utils::nuclear::Density(r,A)); 
        }
        ymax *= 1.2;
        
        // select a vertex using the rejection method
        unsigned int iter = 0;
        while(1) {
          iter++;
          
          // throw an exception if it hasn't find a solution after many attempts
          if(iter > kRjMaxIterations) {
               LOG("Vtx", pWARN)
                 << "Couldn't generate a vertex position after " << iter << " iterations";
               evrec->EventFlags()->SetBitNumber(kGenericErr, true);
               genie::exceptions::EVGThreadException exception;
               exception.SetReason("Couldn't generate vertex");
               exception.SwitchOnFastForward();
               throw exception;
           }

           double r = rmax * rnd->RndFsi().Rndm();
           double t = ymax * rnd->RndFsi().Rndm();
           double y = r*r * utils::nuclear::Density(r,A);
           if(y > ymax) {
               LOG("Vtx", pERROR)
                 << "y = " << y << " > ymax = " << ymax 
                              << " for r = " << r << ", A = " << A;
           }
           bool accept = (t < y);
           if(accept) {
             double phi      = 2*kPi * rnd->RndFsi().Rndm();
             double cosphi   = TMath::Cos(phi);
             double sinphi   = TMath::Sin(phi);
             double costheta = -1 + 2 * rnd->RndFsi().Rndm();
             double sintheta = TMath::Sqrt(1-costheta*costheta);
             vtx.SetX(r*sintheta*cosphi);
             vtx.SetY(r*sintheta*sinphi);
             vtx.SetZ(r*costheta);
             break;
           }
        }//w(1)
    } //use density?

  } // nuclear target ?

  LOG("Vtx", pNOTICE) 
     << "Generated vtx @ r = " << vtx.Mag() << " fm / " 
                                          << print::Vec3AsString(&vtx);

  // Copy the vertex info to the particles already in the event  record
  //
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

 fVtxGenMethod = fConfig->GetIntDef (
           "VtxGenerationMethod", gc->GetInt("NUCL-VtxGenerationMethod")); 

 fR0 = fConfig->GetDoubleDef (
           "R0", gc->GetDouble("NUCL-R0")); // fm
}
//____________________________________________________________________________

