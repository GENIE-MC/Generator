//____________________________________________________________________________
/*!

\namespace  genie::units

\brief      Physical System of Units

\author     Costas Andreopoulos <c.andreopoulos \at cern.ch>
            University of Liverpool

\created    May 03, 2004

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _UNITS_H_
#define _UNITS_H_

namespace genie {

namespace units {

//-- Basic unit

static constexpr double gigaelectronvolt  = 1.;
static constexpr double GeV               = gigaelectronvolt;

//-- Conversion of conventionl [L], [M], [T] units in physical units

static constexpr double kSpeedOfLight = 2.99792458e+08;   // [m/s]  exact by definition
static constexpr double qe_coulomb    = 1.602176634e-19; // electron charge magnitude, exact by definition since 20 May 2019
static constexpr double hbarc         = 1.973269804e-16;  // [GeV*m] exact by definition although approxicamted here
static constexpr double meter         = 1. / (hbarc * GeV);
static constexpr double kilogram      = kSpeedOfLight * kSpeedOfLight * 10e-9 * GeV / qe_coulomb ;
static constexpr double second        = meter * kSpeedOfLight / GeV;


static constexpr double kilometer   = 1000.*meter;
static constexpr double millimeter  = 0.001*meter;
static constexpr double millimeter2 = millimeter*millimeter;
static constexpr double millimeter3 = millimeter*millimeter2;
static constexpr double centimeter  = 0.01*meter;
static constexpr double centimeter2 = centimeter*centimeter;
static constexpr double centimeter3 = centimeter*centimeter2;
static constexpr double decimeter   = 0.1*meter;
static constexpr double decimeter2  = decimeter*decimeter;
static constexpr double decimeter3  = decimeter*decimeter2;
static constexpr double meter2      = meter*meter;
static constexpr double meter3      = meter*meter2;
static constexpr double micrometer  = 1.e-6 *meter;
static constexpr double nanometer   = 1.e-9 *meter;
static constexpr double angstrom    = 1.e-10*meter;
static constexpr double fermi       = 1.e-15*meter;
static constexpr double fermi2      = fermi*fermi;
static constexpr double fermi3      = fermi*fermi2;
static constexpr double barn        = 1.e-28*meter2;
static constexpr double millibarn   = 1.e-3 *barn;
static constexpr double microbarn   = 1.e-6 *barn;
static constexpr double nanobarn    = 1.e-9 *barn;
static constexpr double picobarn    = 1.e-12*barn;

static constexpr double km  = kilometer;
static constexpr double mm  = millimeter;
static constexpr double mm2 = millimeter2;
static constexpr double mm3 = millimeter3;
static constexpr double cm  = centimeter;
static constexpr double cm2 = centimeter2;
static constexpr double cm3 = centimeter3;
static constexpr double m   = meter;
static constexpr double m2  = meter2;
static constexpr double m3  = meter3;
static constexpr double A   = angstrom;
static constexpr double fm  = fermi;
static constexpr double fm2 = fermi2;
static constexpr double fm3 = fermi3;
static constexpr double b   = barn;
static constexpr double mb  = millibarn;
static constexpr double ub  = microbarn;
static constexpr double nb  = nanobarn;
static constexpr double pb  = picobarn;

//-- [T: time]

static constexpr double millisecond   = 1.e-03 *second;
static constexpr double microsecond   = 1.e-06 *second;
static constexpr double nanosecond    = 1.e-09 *second;
static constexpr double picosecond    = 1.e-12 *second;
static constexpr double femptosecond  = 1.e-15 *second;
static constexpr double attosecond    = 1.e-18 *second;
static constexpr double zeptosecond   = 1.e-21 *second;
static constexpr double yoctosecond   = 1.e-24 *second;

static constexpr double s  = second;
static constexpr double ms = millisecond;
static constexpr double us = microsecond;
static constexpr double ns = nanosecond;
static constexpr double ps = picosecond;
static constexpr double fs = femptosecond;
static constexpr double as = attosecond;
static constexpr double zs = zeptosecond;
static constexpr double ys = yoctosecond; 

static constexpr double hertz     = 1./second;
static constexpr double kilohertz = 1.e+3*hertz;
static constexpr double megahertz = 1.e+6*hertz;
static constexpr double gigahertz = 1.e+9*hertz;

static constexpr double  Hz  = hertz;
static constexpr double  kHz = kilohertz;
static constexpr double  MHz = megahertz;
static constexpr double  GHz = gigahertz;

//-- [Q: Charge]

static constexpr double qe          = 1.;

//-- [E: Energy]

static constexpr double     electronvolt = 1.e-09 *GeV;
static constexpr double kiloelectronvolt = 1.e+03 *electronvolt;
static constexpr double megaelectronvolt = 1.e+06 *electronvolt ;
static constexpr double teraelectronvolt = 1.e+12 *electronvolt;
static constexpr double petaelectronvolt = 1.e+15 *electronvolt;

static constexpr double  eV = electronvolt;
static constexpr double keV = kiloelectronvolt;
static constexpr double MeV = megaelectronvolt;
static constexpr double TeV = teraelectronvolt;
static constexpr double PeV = petaelectronvolt;

static constexpr double GeV2 = GeV * GeV;
static constexpr double GeV3 = GeV * GeV2;
static constexpr double GeV4 = GeV * GeV3;
static constexpr double GeV5 = GeV * GeV4;

//-- [M: Mass]

static constexpr double      gram = 1.e-3 *kilogram;
static constexpr double milligram = 1.e-3 *gram;

static constexpr double  kg = kilogram;
static constexpr double   g = gram;
static constexpr double  mg = milligram;

//-- [Density]

static constexpr double  kilogram_meter3  = kilogram / meter3;
static constexpr double  gram_centimeter3 = gram     / centimeter3;

static constexpr double kg_m3 = kilogram_meter3;
static constexpr double g_cm3 = gram_centimeter3;

//-- [Dimensionless quantities]

// Angle

static constexpr double radian      = 1.;
static constexpr double milliradian = 1.e-3*radian;
static constexpr double degree      = (3.14159265358979323846/180.0)*radian;
static constexpr double steradian   = 1.;

static constexpr double rad  = radian;
static constexpr double mrad = milliradian;
static constexpr double sr   = steradian;
static constexpr double deg  = degree;

//-- [Etc]

static constexpr double clhep_def_density_unit = g_cm3/(0.62415185185E+19);

} // namespace units
} // namespace genie

#endif // _UNITS_H_
