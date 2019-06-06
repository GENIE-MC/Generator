//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________

#include "ARSampledNucleus.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <istream>
#include <complex>

#include <TSystem.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/GaussLegendreIntegrator.h>

#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Numerical/IntegrationTools.h"
#include "Physics/Coherent/XSection/ARConstants.h"

//
// Equation/Table numbers refer to:
// J. Nieves and E. Oset
// A Theoretical approach to pionic atoms and the problem of anomalies
// Nuclear Physics A 554 (1993) 509-553
//

using namespace genie::constants;

namespace genie {
namespace alvarezruso {

double ARSampledNucleus::mean_radius_squared = 0.69; // fm (CA)

ARSampledNucleus::ARSampledNucleus(unsigned int ZNumber, unsigned int ANumber, unsigned int fSampling_):
  fZ(ZNumber),
  fA(ANumber),
  fSampling(fSampling_)
{
  fNDensities = 2*fSampling;
  
  fDensities = NULL;
  fDensitiesOfCentres = NULL;
  fRadii = NULL;
  fSample_points_1  = NULL;
  fSample_points_2  = NULL;
  fSample_weights_1 = NULL;
  fSample_weights_2 = NULL;
  
  if(fA>20) { // fermi gas
    if      (fA ==  27) { fNucRadius = 3.07; fDiffuseness = 0.52; }  // aluminum
    else if (fA ==  28) { fNucRadius = 3.07; fDiffuseness = 0.54; }  // silicon
    else if (fA ==  40) { fNucRadius = 3.53; fDiffuseness = 0.54; }  // argon
    else if (fA ==  56) { fNucRadius = 4.10; fDiffuseness = 0.56; }  // iron
    else if (fA == 208) { fNucRadius = 6.62; fDiffuseness = 0.55; }  // lead
    else { 
       fNucRadius = pow(fA,0.35); fDiffuseness = 0.54; 
    } //others
    fUseHarmonicOscillator = false;
  }
  else {
    if      (fA ==  7) { fNucRadius = 1.77   ; fDiffuseness = 0.327;} // lithium
    else if (fA == 12) { fNucRadius = 1.692  ; fDiffuseness = 1.083;} // carbon
    else if (fA == 14) { fNucRadius = 1.76   ; fDiffuseness = 1.23; } // nitrogen
    else if (fA == 16) { fNucRadius = 1.83   ; fDiffuseness = 1.54; } // oxygen
    else if (fA <= 4 ) { fNucRadius = 1.344  ; fDiffuseness = 0.00; } // helium
    else {
      fNucRadius=1.75; fDiffuseness=-0.4+.12*fA; 
    }  //others- fDiffuseness=0.08 if A=4
    
    fUseHarmonicOscillator = true;
  }
  
  fNucRadiusSq = fNucRadius*fNucRadius;
  
  // Calculate number densities (density of 'centres'):
  double diff2 = fDiffuseness*fDiffuseness;
  if (fUseHarmonicOscillator) {
    fRadiusCentres =  TMath::Sqrt( (fNucRadiusSq) - (mean_radius_squared / 1.5) );
    double x = (fDiffuseness * fNucRadiusSq) / ((1 + (1.5*fDiffuseness)) * (fRadiusCentres*fRadiusCentres));
    fDiffusenessCentres =   2.*x / (2. - 3.*x ) ;
  }
  else { //fermi liquid
    double numerator = 5.0 * mean_radius_squared * fNucRadius;
    double denominator = (15.0 * (fNucRadiusSq)) + (7.0 * kPi2 * fDiffuseness*fDiffuseness);
    fRadiusCentres = fNucRadius + (numerator / denominator);
    
    numerator = (fNucRadiusSq*fNucRadius) + (kPi2 * diff2 * fRadiusCentres) - (fRadiusCentres*fRadiusCentres*fRadiusCentres);
    denominator = kPi2 * fRadiusCentres;
    
    fDiffusenessCentres = sqrt( numerator / denominator );
  }
  
  this->Fill();
}

ARSampledNucleus::~ARSampledNucleus()
{
  for(unsigned int i = 0; i != fNDensities; ++i)
  {
    if (fDensities          && fDensities         [i]) delete[] fDensities[i];
    if (fDensitiesOfCentres && fDensitiesOfCentres[i]) delete[] fDensitiesOfCentres[i];
    if (fRadii              && fRadii             [i]) delete[] fRadii    [i];
  }
  
  if (fDensities         ) delete[] fDensities;
  if (fDensitiesOfCentres) delete[] fDensitiesOfCentres;
  if (fRadii             ) delete[] fRadii    ;
  if (fSample_points_1   ) delete[] fSample_points_1;
  if (fSample_points_2   ) delete[] fSample_points_2;
  if (fSample_weights_1  ) delete[] fSample_weights_1;
  if (fSample_weights_2  ) delete[] fSample_weights_2;
}
void ARSampledNucleus::Fill()
{
  
  fR_max = 3.0 * TMath::Power(this->A(), (1.0/3.0));
  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("AR_PiWFunction_Table", pDEBUG)<< "N:: fR_max = " << fR_max
  << "N:: z = " << fZ
  << "N:: a = " << fA
#endif

  this->FillSamplePoints();
  this->FillDensities();
}

double ARSampledNucleus::Density(const int i, const int j) const
{
  return fDensities[i][j];
}

double ARSampledNucleus::DensityOfCentres(const int i, const int j) const
{
  return fDensitiesOfCentres[i][j];
}

double ARSampledNucleus::Radius(const int i, const int j) const
{
  return fRadii[i][j];
}

void ARSampledNucleus::FillSamplePoints()
{
  if (fSample_points_1)  delete[] fSample_points_1;
  if (fSample_points_2)  delete[] fSample_points_2;
  if (fSample_weights_1) delete[] fSample_weights_1;
  if (fSample_weights_2) delete[] fSample_weights_2;
  
  fSample_points_1 = new double[fNDensities];
  fSample_points_2 = new double[fNDensities];
  
  fSample_weights_1 = new double[fNDensities];
  fSample_weights_2 = new double[fNDensities];
  unsigned int decoy;
  
  integrationtools::SGNR(    0.0, fR_max, 2, fSampling, fSample_points_1, decoy, fSample_weights_1);
  integrationtools::SGNR(-fR_max, fR_max, 2, fSampling, fSample_points_2, decoy, fSample_weights_2);
  
}

void ARSampledNucleus::FillDensities()
{
  double r;
  
  for(unsigned int i = 0; i != fNDensities; ++i)
  {
    if (fDensities          && fDensities         [i]) delete[] fDensities         [i];
    if (fDensitiesOfCentres && fDensitiesOfCentres[i]) delete[] fDensitiesOfCentres[i];
    if (fRadii              && fRadii             [i]) delete[] fRadii             [i];
  }
  if (fDensities         ) delete[] fDensities;
  if (fDensitiesOfCentres) delete[] fDensitiesOfCentres;
  if (fRadii             ) delete[] fRadii    ;
  
  fDensities          = new double*[fNDensities];
  fDensitiesOfCentres = new double*[fNDensities];
  fRadii              = new double*[fNDensities];
  
  for(unsigned int i = 0; i != fNDensities; ++i)
  {
    fDensities         [i] = new double[fNDensities];
    fDensitiesOfCentres[i] = new double[fNDensities];
    fRadii             [i] = new double[fNDensities];
    
    for(unsigned int j = 0; j != fNDensities; ++j)
    {
      r = TMath::Sqrt( fSample_points_1[i]*fSample_points_1[i] + fSample_points_2[j]*fSample_points_2[j] );
      fRadii[i][j] = r;
      fDensities         [i][j]=this->CalcMatterDensity(r);
      fDensitiesOfCentres[i][j]=this->CalcNumberDensity(r);
    }
  }
}

unsigned int ARSampledNucleus::GetSampling  (void) const
{
  return fSampling;
}
unsigned int ARSampledNucleus::GetNDensities(void) const
{
  return fNDensities;
}

double ARSampledNucleus::CalcDensity(double r, double nuc_rad, double nuc_diff) const
{
  double dens_rel;
  if (fUseHarmonicOscillator) {
    dens_rel = (1.0 + nuc_diff*r*r/(nuc_rad*nuc_rad)) * exp(-r*r/(nuc_rad*nuc_rad));
  }
  else { //fermi liquid
    dens_rel = 1.0 / (1.0 + exp((r - nuc_rad)/nuc_diff));
  }
  
  //~ double dens_0  = utils::nuclear::Density(nuc_rad, fA);
  double dens_0  = Density0(fA,nuc_rad,nuc_diff);
  
  return dens_0 * dens_rel;
}

double ARSampledNucleus::Density0(
  unsigned int number, // A
  double radius,  // R
  double diffuseness  // a
  ) const
{
  double result = 0.0;
  if( fUseHarmonicOscillator )
  {
    double u = fR_max/radius;
    double dterm = (3*diffuseness+2);
    double term1 = TMath::Sqrt(TMath::Pi())*dterm*TMath::Power(radius,3)*TMath::Exp(u*u)*TMath::Erf(u);
    double term2 = 2*fR_max*(dterm*radius*radius+2*diffuseness*radius*radius);
    result = (0.5)*TMath::Pi()*TMath::Exp(-u*u)*(term1 - term2);
  }
  else
  {
    //Probably faster to do this via the integral for the moment as ROOT doesn't have a builtin
    //PolyLog function
    TF1* f = this->Density0Function();
    f->SetParameter(0, diffuseness);
    f->SetParameter(1, radius);
    result = f->Integral(0.0,fR_max);
    delete f;
  }
  
  return ( number / result );
}

TF1* ARSampledNucleus::Density0Function() const
{
  return (new TF1("density0function", Density0FunctionFermiLiquid, 0.0, fR_max, 2));
}

Double_t ARSampledNucleus::Density0FunctionFermiLiquid(Double_t* r_, Double_t* parameters)
{
  double r = (*r_);
  double diffuseness = parameters[0];
  double radius = parameters[1];
  
  double dens = 1.0 / (1.0 + exp((r - radius)/diffuseness));
  
  return 4.0*TMath::Pi()*r*r*dens;
}

double ARSampledNucleus::CalcMatterDensity(double r) const
{
  return this->CalcDensity(r,fNucRadius,fDiffuseness);
}

double ARSampledNucleus::CalcNumberDensity(double r) const
{
  return this->CalcDensity(r,fRadiusCentres,fDiffusenessCentres);
}
  
} //namespace alvarezruso
} //namespace genie
