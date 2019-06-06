//____________________________________________________________________________
/*!

\namespace  genie::units

\brief      Physical System of Units

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    May 03, 2004

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _UNITS_H_
#define _UNITS_H_

namespace genie {

namespace units {

//-- Basic unit

static const double gigaelectronvolt  = 1.;
static const double GeV               = gigaelectronvolt;

//-- Conversion of conventionl [L], [M], [T] units in physical units

static const double meter    = 5.07e+15 / GeV;
static const double kilogram = 5.61e+26 * GeV;
static const double second   = 1.52e+24 / GeV;

// GeV^-2  -> mbarns : x 0.389;
// mbarns  -> cm^2   : x 1.00E-27;
// m       -> GeV^-1 : x 5.07E+15;
// cm      -> GeV^-1 : x 5.07E+13;
// kgr     -> GeV    : x 5.61E+26;
// gr      -> GeV    : x 5.61E+23;
// sec     -> GeV^-1 : x 1.52E+24;
// gr/cm^3 -> GeV^4  : x 4.30466E-18;

//-- [L: length],[S: area],[V: volume]

static const double kilometer   = 1000.*meter;
static const double millimeter  = 0.001*meter;
static const double millimeter2 = millimeter*millimeter;
static const double millimeter3 = millimeter*millimeter2;
static const double centimeter  = 0.01*meter;
static const double centimeter2 = centimeter*centimeter;
static const double centimeter3 = centimeter*centimeter2;
static const double decimeter   = 0.1*meter;
static const double decimeter2  = decimeter*decimeter;
static const double decimeter3  = decimeter*decimeter2;
static const double meter2      = meter*meter;
static const double meter3      = meter*meter2;
static const double micrometer  = 1.e-6 *meter;
static const double nanometer   = 1.e-9 *meter;
static const double angstrom    = 1.e-10*meter;
static const double fermi       = 1.e-15*meter;
static const double fermi2      = fermi*fermi;
static const double fermi3      = fermi*fermi2;
static const double barn        = 1.e-28*meter2;
static const double millibarn   = 1.e-3 *barn;
static const double microbarn   = 1.e-6 *barn;
static const double nanobarn    = 1.e-9 *barn;
static const double picobarn    = 1.e-12*barn;

static const double km  = kilometer;
static const double mm  = millimeter;
static const double mm2 = millimeter2;
static const double mm3 = millimeter3;
static const double cm  = centimeter;
static const double cm2 = centimeter2;
static const double cm3 = centimeter3;
static const double m   = meter;
static const double m2  = meter2;
static const double m3  = meter3;
static const double A   = angstrom;
static const double fm  = fermi;
static const double fm2 = fermi2;
static const double fm3 = fermi3;
static const double b   = barn;
static const double mb  = millibarn;
static const double ub  = microbarn;
static const double nb  = nanobarn;
static const double pb  = picobarn;

//-- [T: time]

static const double millisecond = 1.e-03 *second;
static const double microsecond = 1.e-06 *second;
static const double nanosecond  = 1.e-09 *second;
static const double picosecond  = 1.e-12 *second;

static const double s  = second;
static const double ms = millisecond;
static const double us = microsecond;
static const double ns = nanosecond;
static const double ps = picosecond;

static const double hertz     = 1./second;
static const double kilohertz = 1.e+3*hertz;
static const double megahertz = 1.e+6*hertz;
static const double gigahertz = 1.e+9*hertz;

static const double  Hz  = hertz;
static const double  kHz = kilohertz;
static const double  MHz = megahertz;
static const double  GHz = gigahertz;

//-- [Q: Charge]

static const double qe          = 1.;
static const double qe_coulomb  = 1.60217733e-19;

//-- [E: Energy]

static const double     electronvolt = 1.e-09 *GeV;
static const double kiloelectronvolt = 1.e+03 *electronvolt;
static const double megaelectronvolt = 1.e+06 *electronvolt ;
static const double teraelectronvolt = 1.e+12 *electronvolt;
static const double petaelectronvolt = 1.e+15 *electronvolt;

static const double  eV = electronvolt;
static const double keV = kiloelectronvolt;
static const double MeV = megaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;

static const double GeV2 = GeV * GeV;
static const double GeV3 = GeV * GeV2;
static const double GeV4 = GeV * GeV3;
static const double GeV5 = GeV * GeV4;

//-- [M: Mass]

static const double      gram = 1.e-3 *kilogram;
static const double milligram = 1.e-3 *gram;

static const double  kg = kilogram;
static const double   g = gram;
static const double  mg = milligram;

//-- [Density]

static const double  kilogram_meter3  = kilogram / meter3;
static const double  gram_centimeter3 = gram     / centimeter3;

static const double kg_m3 = kilogram_meter3;
static const double g_cm3 = gram_centimeter3;

//-- [Dimensionless quantities]

// Angle

static const double radian      = 1.;
static const double milliradian = 1.e-3*radian;
static const double degree      = (3.14159265358979323846/180.0)*radian;
static const double steradian   = 1.;

static const double rad  = radian;
static const double mrad = milliradian;
static const double sr   = steradian;
static const double deg  = degree;

//-- [Etc]

static const double clhep_def_density_unit = g_cm3/(0.62415185185E+19);

} // namespace units
} // namespace genie

#endif // _UNITS_H_

