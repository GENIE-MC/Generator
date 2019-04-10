//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Daniel Scully ( d.i.scully \at warwick.ac.uk)
   University of Warwick

*/
//____________________________________________________________________________


#include <TMath.h>

#include <cstdlib>

#include "Framework/Numerical/IntegrationTools.h"
#include "Physics/Coherent/XSection/AREikonalSolution.h"
#include "Physics/NuclearState/NuclearUtils.h"

typedef std::complex<double> cdouble;

namespace genie {
namespace alvarezruso {

cdouble AREikonalSolution::Element(const double radius, const double cosine_rz, 
                                   const double e_pion)
{
  
  const double mpik = this->Parent()->GetPiMass();
  const double mpi = this->Con()->PiPMass();
  const double hb = this->Con()->HBar() * 1000.0;
  
  const double r = (radius);
  
  const double ekin = (e_pion - mpik) * hb;
  
  const double cosa = cosine_rz;

  const double rmax = this->Nucleus()->RadiusMax();
  
  const unsigned int nz = 1;
  
  const double za = r * cosa;
  const double be = r * TMath::Sqrt(1.0 - cosa*cosa);
  const double omepi = ekin / hb + mpi;
  const double ppim = TMath::Sqrt(omepi*omepi - mpi*mpi);
  unsigned int sampling = (this->Nucleus())->GetSampling();
  
  double absiz[sampling];
  double decoy[sampling];
  
  unsigned int junk;
  
  integrationtools::SGNR(za, rmax, nz, sampling, absiz, junk, decoy);
  
  //do i=1,nzs
//  cdouble ordez[sampling];
  cdouble * ordez = new cdouble[sampling]; // CA
  double zp, rp;
  cdouble piself;
  
  unsigned int A = fNucleus->A();
  unsigned int Z = fNucleus->Z();

  for(unsigned int i = 0; i != sampling; ++i)
  {
    // Sample point in nucleus
    zp = absiz[i];
    
    // Radius in nucleus
    rp = TMath::Sqrt( be*be + zp*zp );
    
    // Get nuclear densities
    double dens_cent = fNucleus->CalcNumberDensity(rp);
    double dens_p_cent = dens_cent * Z / A ;
    double dens_n_cent = dens_cent * (A-Z)/A;
    
    // Calculate pion self energy
    piself = this->PionSelfEnergy(dens_p_cent, dens_n_cent, omepi, ppim);
    
    // Optical potential at each point in the nucleus
    ordez[i] = piself / 2.0 / ppim;
    
  }

  //Integrate the optical potential through the nucleus
  cdouble resu = integrationtools::RGN1D(za, rmax, nz, sampling, ordez);
    
  // Eikonal approximation to the wave function  
  cdouble uwaveik = exp( - cdouble(0,1) * ( ppim*za + resu ) );
  
  delete [] ordez; // CA

  return uwaveik;
}


cdouble AREikonalSolution::PionSelfEnergy(const double rhop_cent, const double rhon_cent, const double omepi, const double ppim)
{
  const double rho0 = this->Con()->Rho0();
  const double mn = this->Con()->NucleonMass();
  const double hb = this->Con()->HBar() * 1000.0;
  const double pi = constants::kPi;
  const double fs = this->Con()->DeltaNCoupling();
  const double mdel = this->Con()->DeltaPMass();
  const double mpi = this->Con()->PiPMass();
  const double fs_mpi2 = fs*fs/(mpi*mpi);
  
  
  const cdouble ui(0,1);
  
  const double resig = -53.0/hb;
  const double rho = rhop_cent + rhon_cent;
  const double rat = rho / rho0;
  const double pf = TMath::Power( (3.*pi*pi/2.*rho), (1.0/3.0) );
            
  const double sdel = mn*mn + mpi*mpi + 2.*omepi*(mn+3./5.*pf*pf/2./mn);
  const double sqsdel = TMath::Sqrt(sdel);
  
  double gamdpb, imsig;
  this->Deltamed(sdel, pf, rat, gamdpb, imsig, ppim, omepi);
  
  const cdouble pe = -1./6./pi*fs_mpi2*
    ( rhop_cent/(sqsdel-mdel-resig*(2.*rhon_cent/rho0)+ui*(gamdpb/2.-imsig)) +
    1./3.*rhon_cent/(sqsdel-mdel-resig*2./3.*(2.*rhon_cent+rhop_cent)/rho0+ui*(gamdpb/2.-imsig)) +
    rhon_cent/(-sqsdel-mdel+2.*mn-resig*(2.*rhop_cent/rho0))+
    1./3.*rhop_cent/(-sqsdel-mdel+2.*mn-resig*2./3.*(2.*rhop_cent+rhon_cent)/rho0) );

  const cdouble efe = 4.*pi*mn*mn/sdel*pe/(1.+4.*pi*0.63*pe);

  const cdouble piself = -1.0*efe*(1.-1./2.*omepi/mn)*ppim*ppim/(1.+efe*(1.-1./2.*omepi/mn));
  
  return piself;
}


void AREikonalSolution::Deltamed(const double sdel, const double pf, const double rat, double& gamdpb, double& imsig, const double ppim, const double omepi)
{
  unsigned int iapr = 1; // approximation chosen to calculate gamdpb
           
  // Calculation of the pauli-blocked width
  double gamdfree = this->Gamd(sdel);
                      
  if( gamdfree == 0.0)
  {
    gamdpb = 0.0;
  }
  else
  {
    double f;
    if( iapr == 1 )
    {
      // Approximation from Nieves et al. NPA 554(93)554           
      const double r = this->Qcm(sdel) / pf;
      if( r > 1.0 )
      {  
                                // f=1.+(-2./5./r**2+9./35./r**4-2./21./r**6)
        f = 1.0 + ( -2.0 / 5.0 / (r*r) + 9.0 / 35.0 / (r*r*r*r) - 2.0 / 21.0 / (r*r*r*r*r*r) );
      }
      else
      {
        //f=34./35.*r-22./105*r**3
        f = 34.0 / 35.0 * r - 22.0 / 105 * (r*r*r);
      }
    }
    else
    {
      //Approximation from Garcia-Recio, NPA 526(91)685
      
      const double mn = this->Con()->NucleonMass();
      const double wd = TMath::Sqrt(sdel); // Delta inv. mass
      const double ef = TMath::Sqrt(mn*mn + pf*pf); // Fermi energy
      const double kd = ppim; // modulus of the Delta 3-momentum in Lab.
      const double pn = this->Qcm(sdel); // nucleon(and pion) momentum in C.M.
      const double en = TMath::Sqrt(mn*mn + pn*pn);

      f = ( kd * pn + TMath::Sqrt(sdel+kd*kd) * en - ef * wd ) / (2.0 * kd * pn);
      if (f < 0.0) f = 0.0;
      if (f > 1.0) f = 1.0;
    }
    gamdpb = gamdfree * f;
  }
  
  //Calculation of the delta selfenergy

  // Imaginary part: using Oset, Salcedo, NPA 468(87)631
  // Using eq. (3.5) to relate the energy of the delta with the pion energy used 
  // in the parametrization

  // Prescriptions for the effective pion energy
  // nucleon at rest
  // ! ome=p(0)-mn
  // ! ome=(sdel-mn**2-mpi**2)/2./mn
  // nucleon with an average momentum
  // ! ekp=3./5.*pf**2/2./mn
  // ! ome=p(0)-mn-ekp    
  // ! ome=(sdel-mn**2-mpi**2-ekp**2)/2./(mn+ekp)
  double ome = omepi;
  double mpi = this->Parent()->GetPiMass();

  // The parameterization is valid for 85 MeV < tpi < 315. outside  we take a contant values
  const double hb = this->Con()->HBar() * 1000.0;
  if( ome <= (mpi + 85.0 / hb) ) ome = mpi + 85.0 / hb;
  if( ome >= (mpi + 315.0 / hb) ) ome = mpi + 315.0 / hb;
           
  // The parameterization of Oset, Salcedo, with ca3 extrapolated to zero at low kin. energies
  const double cq   = this->Cc(-5.19,15.35,2.06,ome)/hb;
  const double ca2  = this->Cc(1.06,-6.64,22.66,ome)/hb;
  double ca3  = this->Cc(-13.46,46.17,-20.34,ome)/hb;
  const double alpha= this->Cc(0.382,-1.322,1.466,ome);
  const double beta = this->Cc(-0.038,0.204,0.613,ome);
  const double gamma=2.*beta;
  if( ome <= (mpi + 85.0/hb) ) ca3 = this->Cc(-13.46,46.17,-20.34,(mpi+85./hb))/85.*(ome-mpi);

  imsig = - ( cq*TMath::Power(rat,alpha) + ca2*TMath::Power(rat, beta) + ca3*TMath::Power(rat,gamma) );
}


double AREikonalSolution::Cc(const double a, const double b, const double c, const double ome)
{
  double mpi = this->Parent()->GetPiMass();
  const double x = ome / mpi - 1.0;
  return (a*x*x + b*x + c);
}


double AREikonalSolution::Gamd(const double s)
{
  // Delta -> N pi
  const double mpi = this->Con()->PiPMass();
  const double mn = this->Con()->NucleonMass();
  if( s <= (mn+mpi)*(mn+mpi) )
  {
    return 0.0;
  }
  else
  {
    double fs_mpi2 = this->Con()->DeltaNCoupling()/mpi;
    fs_mpi2 *= fs_mpi2;
    const double qcm = this->Qcm(s);
    return 1.0 / (6.0*constants::kPi)*fs_mpi2*mn*(qcm*qcm*qcm)/TMath::Sqrt(s);
  }
}


double AREikonalSolution::Qcm(const double s)
{
  // Returns the 3-momentum of the pion formed after the decay of a
  // resonance (R->N pi) of inv. mass s in the rest frame of the resonance
  const double mpi = this->Con()->PiPMass();;
  const double mn = this->Con()->NucleonMass();
  return TMath::Sqrt((s-mpi*mpi-mn*mn)*(s-mpi*mpi-mn*mn) - 4.*mpi*mpi*mn*mn)/2.0/TMath::Sqrt(s);
}

AREikonalSolution::AREikonalSolution(bool debug, AlvarezRusoCOHPiPDXSec* parent): ARWFSolution(debug), parent_(parent)
{
  if( debug_ ) std::cerr << "AREikonalSolution::AREikonalSolution" << std::endl;
  this->fNucleus   = &(this->Parent()->GetNucleus());
  this->constants_ = &(this->Parent()->GetConstants());
  owns_constants = false;
}

AREikonalSolution::AREikonalSolution(bool debug, ARSampledNucleus* nucl): ARWFSolution(debug), fNucleus(nucl)
{
  if( debug_ ) std::cerr << "AREikonalSolution::AREikonalSolution" << std::endl;
  this->constants_ = new ARConstants();
  owns_constants = true;
}

AREikonalSolution::~AREikonalSolution(){
  if (owns_constants) delete constants_;
}
             
void AREikonalSolution::Solve()
{
  if(false) std::cout << "Hi!" << std::endl;
}

} //namespace alvarezruso
} //namespace genie
