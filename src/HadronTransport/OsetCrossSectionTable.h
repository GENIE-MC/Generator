/*
 * Author: T. Golan
 * Tables taken from NuWro's implementation of:
 * E. Oset et al., Nucl. Phys. A484 (1988) 557-592
 * 
 * This implementation is kept from historical reasons.
 * Due to different approach to cascade in GENIE and NuWro there are some normalization issues.
 * Default Oset model can be found in OsetCrossSectionFormula
 * 
*/

#ifndef OSET_CROSS_SECTION_TABLE_H
#define OSET_CROSS_SECTION_TABLE_H

#include "OsetCrossSection.h"
#include <vector>
#include <string>

// handle tables with Oset cross sections
class OsetCrossSectionTable : public OsetCrossSection
{
  public:

  OsetCrossSectionTable (const char* filename); // load table from file

  // use to set up Oset class (assign pion Tk, nuclear density etc)
  void setupOset (const double &density, const double &pionTk,
                  const int &pionPDG, const double &protonFraction);

  private:
  
  // contain values in the following order:
  // d0 e0, d0 e1, ... , d0 en, d1 e0 ...
  std::vector <double> qelCrossSectionTable [nChannels];
  std::vector <double> cexCrossSectionTable [nChannels];
  std::vector <double> absorptionCrossSectionTable;

  unsigned int nDensityBins; // number of denisty bins
  unsigned int nEnergyBins;  // number of energy bins
  double densityBinWidth;    // density step (must be fixed)
  double energyBinWidth;     // energy step (must be fixed)

  // make interpolation between four points around (density, energy)
  double interpolate (const std::vector<double> &data) const;

  // process single line from table file, push values to proper vector
  int processLine (const std::string &line);
  // check if data in file is consistent (bin width etc)
  int checkIntegrity (const double &densityValue, 
                      const double &energyValue);
  // stop program and through an error if input file is corrupted
  void badFile (const char* file, const int &errorCode,
                const int &line = 0) const;

  // set up cross sections for all channels
  void setCrossSections ();

  struct PointHandler  // handle index&weights for given density/energy 
  {
    double value;      // exact value
    double lowWeight;  // distance from high boundary
    double highWeight; // distance from low boundary
    int index;         // point index = index of low boundary
    bool isEdge;       // true if value is on edge of table
                       // should never happen (just in case)
    double binWidth;   // bin width used to calculate distances
    int nBins;         // nBins to check isEdge

    PointHandler () : value (-1.0) {};
    
    inline void setHandler (const double &width, const int &bins)
    {
      binWidth = width;
      nBins    = bins;
    }
    
    void update (const double &newValue); // update point if changed
  };

  PointHandler densityHandler; // nuclear density handler
  PointHandler energyHandler;  // pion kinetic energy handler
};

#endif
