#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "INukeOsetTable.h"

using namespace osetUtils;

//! load tables from file and check its integrity
INukeOsetTable :: INukeOsetTable (const char *filename) : fNDensityBins (0), fNEnergyBins (0), 
                                                          fDensityBinWidth (0.0), fEnergyBinWidth (0.0)
{
  // open file with Oset table
  std::ifstream fileWithTables (filename);
  // check if file is open
  if (not fileWithTables.is_open()) badFile (filename, 1);
  // temporary string for getline
  std::string tempLine;
  // line counter (used in the case of error)
  int lineCounter = 0;

  // skip comments lines at the beginning of the file
  while (fileWithTables.get() == '#' and ++lineCounter)
    getline (fileWithTables, tempLine);
  // process other lines
  while (getline (fileWithTables, tempLine) and ++lineCounter)
    if (int error = processLine (tempLine))   // get exit code
      badFile (filename, error, lineCounter); // stop if exit code != 0

  // reset counter (in the case other file will be loaded)
  checkIntegrity (-1.0, -1.0);

  fDensityHandler.setHandler (fDensityBinWidth, fNDensityBins);
  fEnergyHandler.setHandler  (fEnergyBinWidth,  fNEnergyBins);
}

// process single line from table file, push values to proper vector
/*! <ul>
 * <li> split line read from the file; the following structure is used:
 * \n \n 
 * d; Tk; pi+n: tot fabs fcex; pi+p: tot fabs fcex; pi0: tot fabs fcex
 * \n \n
 * <li> check its integrity
 * <li> save values into proper variables
 * </ul> 
 */
int INukeOsetTable :: processLine (const std::string &line)
{
  std::istringstream splitLine (line); // use istringstream to split string
  // line = d; KE; pi+n: tot fabs fcex; pi+p: tot fabs fcex; pi0: tot fabs fcex

  double currentDensityValue, currentEnergyValue; // will be extracted from line

  // get current density and energy value
  splitLine >> currentDensityValue >> currentEnergyValue;

  if (int exitCode = checkIntegrity (currentDensityValue, currentEnergyValue))
    return exitCode; // stop in the case of error (exit code != 0)

  for (unsigned int i = 0; i < fNChannels; i++) // channel loop
  {
    // get qel and cex cross sections for i-th channel
    double xsecQel, xsecCex;
    splitLine >> xsecQel >> xsecCex;
    // save them in proper vector
    fQelCrossSectionTable[i].push_back (xsecQel);
    fCexCrossSectionTable[i].push_back (xsecCex);
  } // channel loop

  // get absorption cross section
  double absorption;
  splitLine >> absorption;
  // save it in proper vector
  fAbsorptionCrossSectionTable.push_back (absorption);

  return 0; // no errors
}

//! return 0 if all bins widths are constistent or return error code
int INukeOsetTable :: checkIntegrity (const double &densityValue, const double &energyValue)
{
  static unsigned int energyBinCounter = 0; // #energy bins for current density
  
  if (densityValue < 0) // checkIntegrity(-1.0, -1.0) is called to reset counter
    energyBinCounter = -1;
  else if (not isEqual (fNuclearDensity, densityValue)) // on new density value
  {
    fNDensityBins++; // update nof density bins

    if (fNuclearDensity >= 0) // skip first entry (density = -1.0)
    {
      // calculate current bin width
      double currentDensityBinWidth = densityValue - fNuclearDensity;

      // assign first density bin width
      // or check if next bins are the same as first
      // stop program with corresponding error message (exit code 2) if not ==
      if (fDensityBinWidth == 0.0) fDensityBinWidth = currentDensityBinWidth;
      else if (not isEqual (fDensityBinWidth, currentDensityBinWidth)) return 2;
      
      // assign energy bins for first density value
      // or check if next bins are the same as first
      if (fNEnergyBins == 0) fNEnergyBins = energyBinCounter;
      else if (fNEnergyBins != energyBinCounter) return 3; // exit code 3

      // reset counter and last energy
      energyBinCounter = 0;
      fPionKineticEnergy  = -1;
    }
  }
  else if (fPionKineticEnergy >= 0) // check energy width consistency (skip first)
  {
    // calculate current energy bin width
    double currentEnergyBinWidth = energyValue - fPionKineticEnergy;

    // assign first energy bin width
    // or check if next bins are the same as first
    // stop program with corresponding error message (exit code 4) if not ==
    if (fEnergyBinWidth == 0.0) fEnergyBinWidth = currentEnergyBinWidth;
    else if (not isEqual(fEnergyBinWidth, currentEnergyBinWidth)) return 4;
  }

  // increase energy bin counter (for current density)
  energyBinCounter++;
  // update last values
  fPionKineticEnergy  = energyValue;
  fNuclearDensity = densityValue;

  return 0; // no errors
}

/*! possible errors:
 * 
 *  <ul>
 *  <li> the file with tables could not be loaded
 *  <li> density bin width is not constant
 *  <li> there is different number of energy point per density value
 *  <li> energy bin width is not constant
 *  </ul> 
 */ 
void INukeOsetTable :: badFile (const char* file, const int &errorCode, const int &line) const
{
  std::cerr << "\nERROR: " << file << " is corrupted (";
  
  switch (errorCode)
  {
    case 1: std::cerr << "could not be opened).\n"; exit(errorCode);
    case 2: std::cerr << "wrong density"; break;
    case 3: std::cerr << "different number of energy bins per density"; break;
    case 4: std::cerr << "wrong energy"; break;
    default: std::cerr << "undefined error";
  }

  std::cerr << " in line " << line << ").\n";
  exit(errorCode);
}

//! make bilinear interpolation between four points around (density, energy)
double INukeOsetTable :: interpolate (const std::vector<double> &data) const
{
  // take four points adjacent to (density, energy) = (d,E):
  // (d0, E0), (d1, E0), (d0, E1), (d1, E1)
  // where d0 < d < d1, E0 < E < E1
  // each point goes in with weight = proper distance

  // total bin width for normalization
  static const double totalBinWidth = fDensityBinWidth * fEnergyBinWidth;
  
  // extract low boundary value from table
  double lowValue = data[fEnergyHandler.index +
                         fDensityHandler.index * fNEnergyBins]; // (d0, E0)
  // high boundaries = low boundary if on edge
  double highDensityValue = lowValue; // (d1, E0)
  double highEnergyValue  = lowValue; // (d0, E1)
  double highValue        = lowValue; // (d1, E1)

  // extract high boundaries values from table (if not on edge)
  if (not fDensityHandler.isEdge)
  {
    highDensityValue = data[fEnergyHandler.index +
                            (fDensityHandler.index + 1) * fNEnergyBins];

    if (not fEnergyHandler.isEdge)
    {
      highEnergyValue = data[fEnergyHandler.index + 1 +
                             fDensityHandler.index * fNEnergyBins];
      highValue       = data[fEnergyHandler.index + 1 +
                             (fDensityHandler.index + 1) * fNEnergyBins];
    }
  }
  else if (not fEnergyHandler.isEdge)
    highEnergyValue = data[fEnergyHandler.index + 1 +
                           fDensityHandler.index * fNEnergyBins];

  lowValue         *= fDensityHandler.lowWeight  * fEnergyHandler.lowWeight;
  highDensityValue *= fDensityHandler.highWeight * fEnergyHandler.lowWeight;
  highEnergyValue  *= fDensityHandler.lowWeight  * fEnergyHandler.highWeight;
  highValue        *= fDensityHandler.highWeight * fEnergyHandler.highWeight;
  
  return (lowValue + highDensityValue + highEnergyValue + highValue) /
          totalBinWidth;
}

//! set up table index and weights for given point
void INukeOsetTable :: PointHandler :: update (const double &newValue)
{
  value = newValue;         // update value
  index = value / binWidth; // update index

  // in the case value > max value use max; check if it is on edge
  if (index >= nBins - 1)
  {
    index  = nBins - 1;
    isEdge = true;
  }
  else isEdge = false;
  // calcualte weights for boundary values
  lowWeight  = (index + 1) * binWidth - value;
  highWeight = value - index * binWidth;
}


/*! <ul>
 *  <li> set up density, pion Tk and their handlers (for interpolation)
 *  <li> get interpolated cross sections
 *  <li> set up proper cross section variables
 *  </ul> 
 */ 
void INukeOsetTable :: setupOset (const double &density, const double &pionTk, const int &pionPDG,
                                  const double &protonFraction)
{
  fNuclearDensity    = density;
  fPionKineticEnergy = pionTk;
  fDensityHandler.update (fNuclearDensity);    // update density index / weights
  fEnergyHandler.update  (fPionKineticEnergy); // update energy index / weights
  setCrossSections();
  INukeOset::setCrossSections (pionPDG, protonFraction);  
}

/*! assign cross sections values to proper variables
 * using bilinear interpolation
 */ 
void INukeOsetTable :: setCrossSections ()
{
    for (unsigned int i = 0; i < fNChannels; i++) // channel loop
    {
      fQelCrossSections[i] = interpolate (fQelCrossSectionTable[i]);
      fCexCrossSections[i] = interpolate (fCexCrossSectionTable[i]);
    }

    fAbsorptionCrossSection = interpolate (fAbsorptionCrossSectionTable); 
}
