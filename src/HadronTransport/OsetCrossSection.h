/**
 * @brief Oset model handler (abstract class)
 * 
 * @author Tomasz Golan
 * @date 2015
 * @warning Applicable for pion with Tk < 350 MeV
 * @remarks Based on E. Oset et al., Nucl. Phys. A484 (1988) 557-592
 * 
*/

#ifndef OSET_CROSS_SECTION_H
#define OSET_CROSS_SECTION_H

#include <cmath>
#include <limits>

#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Conventions/Constants.h"

using namespace genie;
using namespace genie::constants;

class OsetCrossSection
{
  public:

  OsetCrossSection (); //!< contructor

  //! use to set up Oset class (assign pion Tk, nuclear density etc)
  virtual void setupOset (const double &density, const double &pionTk,
                          const int &pionPDG,
                          const double &protonFraction) = 0;

  //! return total = (qel+cex+abs) cross section 
  inline double getTotalCrossSection  () const
  {
    return totalCrossSection;
  }
  
  //! return cex cross section 
  inline double getCexCrossSection () const
  {
    return cexCrossSection;
  }
  
  //! return absorption cross section
  inline double getAbsorptionCrossSection () const
  {
    return absorptionCrossSection;
  }
  
  //! return fraction of cex events
  inline double getCexFraction () const
  {
    return cexCrossSection / totalCrossSection;
  }
  
  //! return fraction of absorption events
  inline double getAbsorptionFraction () const
  {
    return absorptionCrossSection / totalCrossSection;
  }

  protected:

  double nuclearDensity;    //!< nuclear density in fm-3
  double pionKineticEnergy; //!< pion kinetic energy in MeV

  //! el+cex+abs cross section (averaged over proton / neutron fraction)
  double totalCrossSection;      
  //! cex cross section (averaged over proton / neutron fraction)
  double cexCrossSection;        
  //! absorption cross section (averaged over proton / neutron fraction)
  double absorptionCrossSection; 

  //! number of possible channels: pi+n, pi+p, pi0
  /*! if (pi0) channel = 2 \n
   *  else channel = [(10 * pip + pim) == (10 * p + n)] \n \n
   *  0 -> pi+n or pi-p, 1 -> pi+p or pi-n, 2 -> pi0
   */
  static const unsigned int nChannels = 3; 
  
  //! total qel (el+cex) cross section for each channel
  double qelCrossSections[nChannels]; 
  double cexCrossSections[nChannels]; //!< cex cross section for each channel
  
  //! calculalte cross sections for each channel
  virtual void setCrossSections () = 0;
  
  //! calculate avg cross sections according to proton / neutron fraction
  void setCrossSections (const int &pionPDG,
                         const double &protonFraction);
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
