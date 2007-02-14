//____________________________________________________________________________
/*!

\namespace genie::constants::bba2003

\brief     Coefficients for the BBA2003 (Budd,Bodek,Arrington) elastic form
           factor parameterization.
           Note: these values will be used as 'defaults' if no other values
           are given at the BBA03 model XML config file.

\ref       H.Budd, NuINT02 proceedings

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   October 20, 2005

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _BBA2003_CONSTANTS_H_
#define _BBA2003_CONSTANTS_H_

namespace genie {
namespace constants {
namespace bba2003 {

// Parameters for the Krutov Gen parameterization used in BBA2003
static const double kGen_a   =  0.942;
static const double kGen_b   =  4.610;

// Coefficients for the BBA2003 inverse polynomial fit for Gep,Gmp,Gmn
static const double kGep_a2  =  3.253;
static const double kGep_a4  =  1.422;
static const double kGep_a6  =  0.08582;
static const double kGep_a8  =  0.3318;
static const double kGep_a10 = -0.09371;
static const double kGep_a12 =  0.01076;

static const double kGmp_a2  =  3.104;
static const double kGmp_a4  =  1.428;
static const double kGmp_a6  =  0.1112;
static const double kGmp_a8  = -0.006981;
static const double kGmp_a10 =  0.0003705;
static const double kGmp_a12 = -0.7063E-5;

static const double kGmn_a2  =  3.043;
static const double kGmn_a4  =  0.8548;
static const double kGmn_a6  =  0.6806;
static const double kGmn_a8  = -0.1287;
static const double kGmn_a10 =  0.008912;
static const double kGmn_a12 =  0.;

// Gep/Gmp assummed const for Q2 > Q2Max
static const double kQ2Max   =  6.0;

} // namespace bba2003
} // namespace constants
} // namespace genie

#endif // _BBA2003_CONSTANTS_H_



