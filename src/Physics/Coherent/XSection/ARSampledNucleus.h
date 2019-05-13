//____________________________________________________________________________
/*!

\class    genie::alvarezruso::ARSampledNucleus

\brief    Nucleus class for Alvarez-Ruso Coherent Pion Production xsec

\ref      

\author   Steve Dennis
          University of Warwick, Rutherford Appleton Laboratory

\created  05/12/2013

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _AR_NUCLEUS_H_
#define _AR_NUCLEUS_H_

#include <TF1.h>

namespace genie
{
namespace alvarezruso
{

class ARSampledNucleus
{
  public:
    
    ARSampledNucleus(unsigned int ZNumber, unsigned int ANumber, unsigned int sampling = 20);
    
    ~ARSampledNucleus();
    
    unsigned int A() const  {  return fA;  }
    
    unsigned int Z() const  {  return fZ;  }
    
    unsigned int N() const  {  return (fA-fZ);  }
    
    double Density         (const int i, const int j) const;
    double DensityOfCentres(const int i, const int j) const;
    double Radius          (const int i, const int j) const;
    
    double RadiusMax() const
    {
      return fR_max;
    }
    double SamplePoint1(const unsigned int i) const // absib in original fortran
    {
      return fSample_points_1[i];
    }
    double SamplePoint2(const unsigned int i) const // absiz in original fortran
    {
      return fSample_points_2[i];
    }
    //~ double SampleWeight1(const unsigned int i) const
    //~ {
      //~ return fSample_weights_1[i];
    //~ }
    //~ double SampleWeight2(const unsigned int i) const
    //~ {
      //~ return fSample_weights_2[i];
    //~ }
    
    unsigned int GetSampling  (void) const;
    unsigned int GetNDensities(void) const;
    
    double CalcMatterDensity(double r) const;
    double CalcNumberDensity(double r) const;
    
  private:
  
    void Fill();
    void FillSamplePoints();
    void FillDensities();
    // Members
    
    const unsigned int fZ;
    const unsigned int fA;
    unsigned int fSampling;
    
    unsigned int fNDensities;
    
    double fR_max;
    double** fRadii;
    double** fDensities;
    double** fDensitiesOfCentres;
    double* fSample_points_1;  // absib: 0 < r < r_max
    double* fSample_points_2;  // absiz: -r_max < r < r_max
    double* fSample_weights_1;
    double* fSample_weights_2;
    
    double fDiffuseness;
    double fNucRadius;
    double fNucRadiusSq;
    double fDiffusenessCentres;
    double fRadiusCentres;
  //double fRadiusCentresSq;
    double fUseHarmonicOscillator;
    
    double CalcDensity(double radius, double nuc_rad, double nuc_diff) const;
    
    // warning: in-class initializer for static data member of type 'const double' is a GNU extension [-Wgnu-static-float-init]
    // static const double mean_radius_squared = 0.69; // in fermi
    static double mean_radius_squared;
    
    double Density0(unsigned int number, double diffuseness, double radius) const;
    TF1* Density0Function() const;
    static Double_t Density0FunctionFermiLiquid(Double_t* r, Double_t* parameters);

};

} //namespace alvarezruso
} //namespace genie

#endif
