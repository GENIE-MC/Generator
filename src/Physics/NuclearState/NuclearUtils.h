//____________________________________________________________________________
/*!

\namespace  genie::utils::nuclear

\brief      Simple nuclear physics empirical formulas (densities, radii, ...)
            and empirical nuclear corrections

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    May 06, 2004

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUCLEAR_UTILS_H_
#define _NUCLEAR_UTILS_H_

#include <string>

#include "Framework/Conventions/Constants.h"

using std::string;

namespace genie {

class Target;
class Interaction;

namespace utils {

namespace nuclear
{
  double BindEnergy             (const Target & target);
  double BindEnergy             (int nucA, int nucZ);
  double BindEnergyPerNucleon   (const Target & target);
  double BindEnergyLastNucleon  (const Target & target);
  double Radius                 (int A, double Ro=constants::kNucRo);

  double NuclQELXSecSuppression (string kftable, double pmax, const Interaction * in);

  double RQEFG_generic (
            double q2, double Mn, double kFi, double kFf, double pmax);

  double FmI1   (double alpha, double beta,
                  double a, double b,  double kFi, double kFf, double q);
  double FmI2   (double alpha, double beta,
                  double a, double b,  double kFi, double kFf, double q);
  double FmArea (double alpha, double beta, double kf, double pmax);

  double DISNuclFactor (double x, int A);

  double Density           (double r, int A, double ring=0.);
  double DensityGaus       (double r, double ap, double alf, double ring=0.);
  double DensityWoodsSaxon (double r, double c, double z, double ring=0.);

  double BindEnergyPerNucleonParametrization(const Target & target);
  double FermiMomentumForIsoscalarNucleonParametrization(const Target & target);


} // nuclear namespace
} // utils   namespace
} // genie   namespace

#endif // _NUCL_UTILS_H_
