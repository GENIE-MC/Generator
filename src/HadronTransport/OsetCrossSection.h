#ifndef OSET_CROSS_SECTION_H
#define OSET_CROSS_SECTION_H

#include <cmath>
#include <limits>

class OsetCrossSection // abstract oset class
{
  public:

  OsetCrossSection ();

  // use to set up Oset class (assign pion Tk, nuclear density etc)
  virtual void setupOset (const double &density, const double &pionTk,
                          const int &pionPDG,
                          const double &protonFraction) = 0;

  // get total = (qel+cex+abs) cross section 
  inline double getTotalCrossSection  () const
  {
    return totalCrossSection;
  }
  // get cex cross section 
  inline double getCexCrossSection () const
  {
    return cexCrossSection;
  }
  // get absorption cross section
  inline double getAbsorptionCrossSection () const
  {
    return absorptionCrossSection;
  }
  // get fraction of cex events
  inline double getCexFraction () const
  {
    return cexCrossSection / totalCrossSection;
  }
  // get fraction of absorption events
  inline double getAbsorptionFraction () const
  {
    return absorptionCrossSection / totalCrossSection;
  }

  protected:

  double nuclearDensity;    // nuclear density in fm-3
  double pionKineticEnergy; // pion kinetic energy in MeV

  // cross section averaged over proton / neutron fraction
  double totalCrossSection;      // el+cex+abs cross section
  double cexCrossSection;        //        cex cross section
  double absorptionCrossSection; // absorption cross section

  // cross section per channel 
  // if (pi0) channel = 2;
  // else channel = [(10 * pip + pim) == (10 * p + n)]
  // 0 -> pi+n or pi-p, 1 -> pi+p or pi-n, 2 -> pi0  double
  static const unsigned int nChannels = 3; // pi+n, pi+p, pi0
  double qelCrossSections[nChannels]; // qel = el + cex
  double cexCrossSections[nChannels]; // only cex
  
  // set up cross sections for all channels
  virtual void setCrossSections () = 0;
  // calculate avg cross sections according to proton / neutron fraction
  void setCrossSections (const int &pionPDG,
                         const double &protonFraction);

  static const int kPdgPiP =  211;
  static const int kPdgPiM = -211;
  static const int kPdgPi0 =  111;  
};

namespace osetUtils
{
  // workaround to get access to last instance
  extern OsetCrossSection* currentInstance; 
  // check if double == double with defined accuracy
  inline bool isEqual (const double &x, const double &y,
                       const double &epsilon)
  {
    return std::abs(x - y) < epsilon;
  }
  // check if double == double with best accuracy
  inline bool isEqual (const double &x, const double &y)
  {
    static const double epsilon = std::numeric_limits<double>::epsilon();
    return isEqual (x, y, epsilon);
  }
  // x -> variale, a->coefficients
  inline double quadraticFunction (const double &x, const double *a)
  {
    return a[0] * x * x + a[1] * x + a[2];
  }
} // namespace osetUtils

#endif
