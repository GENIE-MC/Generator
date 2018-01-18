/**
 * @brief Oset model handler (abstract class)
 * 
 * @author Tomasz Golan
 * @date 2015
 * @warning Applicable for pion with Tk < 350 MeV
 * @remarks Based on E. Oset et al., Nucl. Phys. A484 (1988) 557-592
 * 
*/

#ifndef INUKE_OSET_H
#define INUKE_OSET_H

#include <cmath>
#include <limits>

#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Conventions/Constants.h"

using namespace genie;
using namespace genie::constants;

class INukeOset
{
  public:

  INukeOset (); //!< contructor

  //! use to set up Oset class (assign pion Tk, nuclear density etc)
  virtual void setupOset (const double &density, const double &pionTk, const int &pionPDG,
                          const double &protonFraction) = 0;

  //! return total = (qel+cex+abs) cross section 
  inline double getTotalCrossSection  () const
  {
    return fTotalCrossSection;
  }
  
  //! return cex cross section 
  inline double getCexCrossSection () const
  {
    return fCexCrossSection;
  }
  
  //! return absorption cross section
  inline double getAbsorptionCrossSection () const
  {
    return fAbsorptionCrossSection;
  }
  
  //! return fraction of cex events
  inline double getCexFraction () const
  {
    return fCexCrossSection / fTotalCrossSection;
  }
  
  //! return fraction of absorption events
  inline double getAbsorptionFraction () const
  {
    return fAbsorptionCrossSection / fTotalCrossSection;
  }

  protected:

  double fNuclearDensity;    //!< nuclear density in fm-3
  double fPionKineticEnergy; //!< pion kinetic energy in MeV

  //! el+cex+abs cross section (averaged over proton / neutron fraction)
  double fTotalCrossSection;      
  //! cex cross section (averaged over proton / neutron fraction)
  double fCexCrossSection;        
  //! absorption cross section (averaged over proton / neutron fraction)
  double fAbsorptionCrossSection; 

  //! number of possible channels: pi+n, pi+p, pi0
  /*! if (pi0) channel = 2 \n
   *  else channel = [(10 * pip + pim) == (10 * p + n)] \n \n
   *  0 -> pi+n or pi-p, 1 -> pi+p or pi-n, 2 -> pi0
   */
  static const unsigned int fNChannels = 3; 
  
  //! total qel (el+cex) cross section for each channel
  double fQelCrossSections[fNChannels]; 
  double fCexCrossSections[fNChannels]; //!< cex cross section for each channel
  
  //! calculalte cross sections for each channel
  virtual void setCrossSections () = 0;
  
  //! calculate avg cross sections according to proton / neutron fraction
  void setCrossSections (const int &pionPDG,
                         const double &protonFraction);
};

namespace osetUtils
{
  // workaround to get access to last instance
  extern INukeOset* currentInstance; 
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

#endif // INUKE_OSET_H
