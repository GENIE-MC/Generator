//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2007 - CA
   Added code to generate the interaction according to a realistic nuclear
   density and made that the default setting.
 @ Dec 01, 2007 - CA
   For COH and ve- interactions setting the vertex on the nuclear boundary
   rather than inside the nucleus.
 @ Sep 15, 2009 - CA
   IsFake() and IsNucleus() are no longer available in GHepParticle. 
   Use pdg::IsPseudoParticle() and pdg::IsIon().
 @ Feb 12, 2013 - CA (code from Rosen Matev)
   Handle the IMD annihilation channel.
 @ Jun 03, 2016 - JJ (SD)
   Move code to generate the position to a public method, so this code
   can easily be reused by other classes. (Specifically LwlynSmithQELCCPXSec
   and NievesQELCCPXSec to generate a position before calculating the xsec
   when making splines).
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Common/VertexGenerator.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

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
  Interaction * interaction = evrec->Summary();
  GHepParticle * nucltgt = evrec->TargetNucleus();
  TVector3 vtx(9999999.,999999.,999999.);
  if(!nucltgt){
    vtx.SetXYZ(0.,0.,0.);
  }else{
    double A = nucltgt->A();
    vtx = GenerateVertex(interaction,A);
  }

  // Copy the vertex info to the particles already in the event  record
  //
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(pdg::IsPseudoParticle(p->Pdg())) continue;
    if(pdg::IsIon           (p->Pdg())) continue;

    LOG("Vtx", pINFO) << "Setting vertex position for: " << p->Name();
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

  GetParam( "VtxGenerationMethod", fVtxGenMethod ) ;
  GetParam( "NUCL-R0",             fR0 ) ;  //fm

}
//____________________________________________________________________________
TVector3 VertexGenerator::GenerateVertex(const Interaction * interaction,
					 double A) const{
  RandomGen * rnd = RandomGen::Instance();
  TVector3 vtx(999999.,999999.,999999.);

  //GHepParticle * nucltgt = evrec->TargetNucleus();

  bool uniform   = fVtxGenMethod==0;
  bool realistic = !uniform;

  //if(!nucltgt) {
  //vtx.SetXYZ(0.,0.,0.);
  //} 
  //else {
    //double A = nucltgt->A();
  double R = fR0 * TMath::Power(A, 1./3.);
  
  //Interaction * interaction = evrec->Summary();
  const ProcessInfo & proc_info = interaction->ProcInfo();
  bool is_coh = proc_info.IsCoherent() || proc_info.IsCoherentElas();
  bool is_ve  = proc_info.IsInverseMuDecay() ||
    proc_info.IsIMDAnnihilation() ||
    proc_info.IsNuElectronElastic();
  
  if(is_coh||is_ve) {
    // ** For COH or ve- set a vertex positon on the nuclear boundary
    //
    LOG("Vtx", pINFO)  << "Setting vertex on the nuclear boundary";
    double phi      = 2*kPi * rnd->RndFsi().Rndm();
    double cosphi   = TMath::Cos(phi);
    double sinphi   = TMath::Sin(phi);
    double costheta = -1 + 2 * rnd->RndFsi().Rndm();
    double sintheta = TMath::Sqrt(1-costheta*costheta);
    vtx.SetX(R*sintheta*cosphi);
    vtx.SetY(R*sintheta*sinphi);
    vtx.SetZ(R*costheta);
  }
  else {  
    // ** For other events on nuclear targets set the interaction vertex
    // ** using the requested method: either using a realistic nuclear 
    // ** density or by setting it uniformly within the nucleus
    //
    if(realistic) {
      // Generate the vertex using a realistic nuclear density
      //
      LOG("Vtx", pINFO) 
	<< "Generating vertex according to a realistic nuclear density profile";
      // get inputs to the rejection method
      double ymax = -1;
      double rmax = 3*R;
      double dr   = R/40.;
      for(double r = 0; r < rmax; r+=dr) { 
	ymax = TMath::Max(ymax, r*r * utils::nuclear::Density(r,(int)A)); 
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
	  //evrec->EventFlags()->SetBitNumber(kGenericErr, true);
	  genie::exceptions::EVGThreadException exception;
	  exception.SetReason("Couldn't generate vertex");
	  exception.SwitchOnFastForward();
	  throw exception;
	}
	
	double r = rmax * rnd->RndFsi().Rndm();
	double t = ymax * rnd->RndFsi().Rndm();
	double y = r*r * utils::nuclear::Density(r,(int)A);
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
    
    if(uniform) {
      // Generate the vertex uniformly within the nucleus
      //
      LOG("Vtx", pINFO) 
	<< "Generating intranuclear vertex uniformly in volume";
      while(vtx.Mag() > R) {
	vtx.SetX(-R + 2*R * rnd->RndFsi().Rndm());
	vtx.SetY(-R + 2*R * rnd->RndFsi().Rndm());
	vtx.SetZ(-R + 2*R * rnd->RndFsi().Rndm());
      }
    }// uniform? 
    
  } // coh or ve-?
    //} // nuclear target ?
  
  LOG("Vtx", pINFO) 
    << "Generated vtx @ r = " << vtx.Mag() << " fm / " 
    << print::Vec3AsString(&vtx);
  return vtx;
}
