//____________________________________________________________________________
/*!

\class    genie::SampleStats

\brief    A simple class to calculate running values for the mean,
          sample variance, and standard error for a set of input values.
          Non-unit event weights are allowed.

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  January 3, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _SAMPLESTATS_H_
#define _SAMPLESTATS_H_

// Standard library includes
#include <cmath> // std::sqrt

namespace genie {

// Computes running values of sample statistics for weighted events.
// Based on an algorithm by West (https://doi.org/10.1145%2F359146.359153).
// See also http://tinyurl.com/mean-var-onl-alg.
class SampleStats {

  public:

    SampleStats() : fSumW(0.), fSumW2(0.), fN(0), fMean(0.), fS(0.), fT(0.) {}

    inline void AddValue(double value, double weight = 1.) {
      ++fN;

      // Skip the rest for events with zero weight (this also prevents
      // numerical issues caused by dividing by zero)
      if ( weight == 0. ) return;

      fSumW += weight;
      fSumW2 += weight * weight;

      double delta = value - fMean;

      fMean += weight * delta / fSumW;

      // note that fMean has been updated in the line above, so
      // delta != (value - fMean) on the following line
      double term = weight * delta * ( value - fMean );
      fS += term;
      fT += weight * term;
    }

    inline double Mean() const { return fMean; }

    // Subtracting one here applies Bessel's correction for weighted samples.
    // Note also that we assume that frequency weights are used as opposed
    // to reliability weights
    inline double SampleVariance() const { return fS / (fSumW - 1.); }

    inline double SampleStdDeviation() const
      { return std::sqrt( this->SampleVariance() ); }

    inline int SampleSize() const { return fN; }
    inline double SumOfWeights() const { return fSumW; }
    inline double SumOfSquaredWeights() const { return fSumW2; }

    // For a weighted arithmetic mean, there is no universally agreed upon
    // estimator for the standard error. Here we use a formula based on
    // the bootstrapping study reported in
    // https://doi.org/10.1016/1352-2310(94)00210-C.
    inline double SquareOfStdErrorOnMean() const {
      return fN * fT / ( fSumW * fSumW * (fN - 1) );
    }
    // This version uses another estimator given in
    // http://www.analyticalgroup.com/download/WEIGHTED_MEAN.pdf
    //inline double SquareOfStdErrorOnMean() const
    //  { return this->SampleVariance() * fSumW2 / ( fSumW * fSumW ); }

    inline double StdErrorOnMean() const
      { return std::sqrt( this->SquareOfStdErrorOnMean() ); }

  protected:

    double fSumW; // sum of event weights
    double fSumW2; // sum of squared event weights
    int fN; // sample size (unweighted)
    double fMean; // arithmetic mean
    double fS; // running value needed to calculate the sample variance
    double fT; // running value needed to calculate the standard error
};

} // genie namespace

#endif // _SAMPLESTATS_H_
