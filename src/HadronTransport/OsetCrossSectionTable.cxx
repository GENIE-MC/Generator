#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "OsetCrossSectionTable.h"

using namespace osetUtils;

// constructor -> load table from file
OsetCrossSectionTable :: OsetCrossSectionTable (const char *filename) :
                          nDensityBins (0), nEnergyBins (0),
                          densityBinWidth (0.0), energyBinWidth (0.0)
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

  densityHandler.setHandler (densityBinWidth, nDensityBins);
  energyHandler.setHandler  (energyBinWidth,  nEnergyBins);
}
// process single line from table file, push values to proper vector
int OsetCrossSectionTable :: processLine (const std::string &line)
{
  std::istringstream splitLine (line); // use istringstream to split string
  // line = d; KE; pi+n: tot fabs fcex; pi+p: tot fabs fcex; pi0: tot fabs fcex

  double currentDensityValue, currentEnergyValue; // will be extracted from line

  // get current density and energy value
  splitLine >> currentDensityValue >> currentEnergyValue;

  if (int exitCode = checkIntegrity (currentDensityValue, currentEnergyValue))
    return exitCode; // stop in the case of error (exit code != 0)

  for (unsigned int i = 0; i < nChannels; i++) // channel loop
  {
    // get qel and cex cross sections for i-th channel
    double xsecQel, xsecCex;
    splitLine >> xsecQel >> xsecCex;
    // save them in proper vector
    qelCrossSectionTable[i].push_back (xsecQel);
    cexCrossSectionTable[i].push_back (xsecCex);
  } // channel loop

  // get absorption cross section
  double absorption;
  splitLine >> absorption;
  // save it in proper vector
  absorptionCrossSectionTable.push_back (absorption);

  return 0; // no errors
}

// check if data in file is consistent (bin width etc)
int OsetCrossSectionTable :: checkIntegrity (const double &densityValue,
                                              const double &energyValue)
{
  static unsigned int energyBinCounter = 0; // #energy bins for current density
  
  if (densityValue < 0) // checkIntegrity(-1.0, -1.0) is called to reset counter
    energyBinCounter = -1;
  else if (not isEqual (nuclearDensity, densityValue)) // on new density value
  {
    nDensityBins++; // update nof density bins

    if (nuclearDensity >= 0) // skip first entry (density = -1.0)
    {
      // calculate current bin width
      double currentDensityBinWidth = densityValue - nuclearDensity;

      // assign first density bin width
      // or check if next bins are the same as first
      // stop program with corresponding error message (exit code 2) if not ==
      if (densityBinWidth == 0.0) densityBinWidth = currentDensityBinWidth;
      else if (not isEqual (densityBinWidth, currentDensityBinWidth)) return 2;
      
      // assign energy bins for first density value
      // or check if next bins are the same as first
      if (nEnergyBins == 0) nEnergyBins = energyBinCounter;
      else if (nEnergyBins != energyBinCounter) return 3; // exit code 3

      // reset counter and last energy
      energyBinCounter = 0;
      pionKineticEnergy  = -1;
    }
  }
  else if (pionKineticEnergy >= 0) // check energy width consistency (skip first)
  {
    // calculate current energy bin width
    double currentEnergyBinWidth = energyValue - pionKineticEnergy;

    // assign first energy bin width
    // or check if next bins are the same as first
    // stop program with corresponding error message (exit code 4) if not ==
    if (energyBinWidth == 0.0) energyBinWidth = currentEnergyBinWidth;
    else if (not isEqual(energyBinWidth, currentEnergyBinWidth)) return 4;
  }

  // increase energy bin counter (for current density)
  energyBinCounter++;
  // update last values
  pionKineticEnergy  = energyValue;
  nuclearDensity = densityValue;

  return 0; // no errors
}

// stop program and through an error if input file is corrupted
void OsetCrossSectionTable :: badFile (const char* file, const int &errorCode,
                                        const int &line) const
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

// make interpolation between four points around (density, energy)  
double OsetCrossSectionTable :: interpolate (const std::vector<double> &data)
                                              const
{
  // take four points adjacent to (density, energy) = (d,E):
  // (d0, E0), (d1, E0), (d0, E1), (d1, E1)
  // where d0 < d < d1, E0 < E < E1
  // each point goes in with weight = proper distance

  // total bin width for normalization
  static const double totalBinWidth = densityBinWidth * energyBinWidth;
  
  // extract low boundary value from table
  double lowValue = data[energyHandler.index +
                         densityHandler.index * nEnergyBins]; // (d0, E0)
  // high boundaries = low boundary if on edge
  double highDensityValue = lowValue; // (d1, E0)
  double highEnergyValue  = lowValue; // (d0, E1)
  double highValue        = lowValue; // (d1, E1)

  // extract high boundaries values from table (if not on edge)
  if (not densityHandler.isEdge)
  {
    highDensityValue = data[energyHandler.index +
                            (densityHandler.index + 1) * nEnergyBins];

    if (not energyHandler.isEdge)
    {
      highEnergyValue = data[energyHandler.index + 1 +
                             densityHandler.index * nEnergyBins];
      highValue       = data[energyHandler.index + 1 +
                             (densityHandler.index + 1) * nEnergyBins];
    }
  }
  else if (not energyHandler.isEdge)
    highEnergyValue = data[energyHandler.index + 1 +
                           densityHandler.index * nEnergyBins];

  lowValue         *= densityHandler.lowWeight * energyHandler.lowWeight;
  highDensityValue *= densityHandler.highWeight * energyHandler.lowWeight;
  highEnergyValue  *= densityHandler.lowWeight * energyHandler.highWeight;
  highValue        *= densityHandler.highWeight * energyHandler.highWeight;
  
  return (lowValue + highDensityValue + highEnergyValue + highValue) /
          totalBinWidth;
}

// update point if value has changed
void OsetCrossSectionTable :: PointHandler :: update (const double &newValue)
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


// set up density, pion Tk and their handlers (use for interpolation)
void OsetCrossSectionTable :: setupOset (const double &density,
                                         const double &pionTk,
                                         const int &pionPDG,
                                         const double &protonFraction)
{
  nuclearDensity    = density;
  pionKineticEnergy = pionTk;
  densityHandler.update (nuclearDensity);    // update density index / weights
  energyHandler.update  (pionKineticEnergy); // update energy index / weights
  setCrossSections();
  OsetCrossSection::setCrossSections (pionPDG, protonFraction);  
}

// calculate cross sections
void OsetCrossSectionTable :: setCrossSections ()
{
    for (int i = 0; i < nChannels; i++) // channel loop
    {
      qelCrossSections[i] = interpolate (qelCrossSectionTable[i]);
      cexCrossSections[i] = interpolate (cexCrossSectionTable[i]);
    }

    absorptionCrossSection = interpolate (absorptionCrossSectionTable); 
}
