//____________________________________________________________________________
/*!

\namespace  genie::utils::unit

\brief      Simple unit system utilities

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/
//____________________________________________________________________________

#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Utils/UnitUtils.h"

//____________________________________________________________________________
double genie::utils::units::UnitFromString(string u)
{
// Returns the appropriate unit based on the input string
// The GENIE units are defined in $GENIE/src/Conventions/Units.h

  if      (u == "gigaelectronvolt") return  genie::units::gigaelectronvolt;
  else if (u == "GeV"             ) return  genie::units::GeV;

  else if (u == "meter"           ) return  genie::units::meter;
  else if (u == "kilogram"        ) return  genie::units::kilogram;
  else if (u == "second"          ) return  genie::units::second;

  else if (u == "millimeter"      ) return  genie::units::millimeter;
  else if (u == "millimeter2"     ) return  genie::units::millimeter2;
  else if (u == "millimeter3"     ) return  genie::units::millimeter3;
  else if (u == "centimeter"      ) return  genie::units::centimeter;
  else if (u == "centimeter2"     ) return  genie::units::centimeter2;
  else if (u == "centimeter3"     ) return  genie::units::centimeter3;
  else if (u == "decimeter"       ) return  genie::units::decimeter;
  else if (u == "decimeter2"      ) return  genie::units::decimeter2;
  else if (u == "decimeter3"      ) return  genie::units::decimeter3;
  else if (u == "meter2"          ) return  genie::units::meter2;
  else if (u == "meter3"          ) return  genie::units::meter3;
  else if (u == "micrometer"      ) return  genie::units::micrometer;
  else if (u == "nanometer"       ) return  genie::units::nanometer;
  else if (u == "angstrom"        ) return  genie::units::angstrom;
  else if (u == "fermi"           ) return  genie::units::fermi;
  else if (u == "barn"            ) return  genie::units::barn;
  else if (u == "millibarn"       ) return  genie::units::millibarn;
  else if (u == "microbarn"       ) return  genie::units::microbarn;
  else if (u == "nanobarn"        ) return  genie::units::nanobarn;
  else if (u == "picobarn"        ) return  genie::units::picobarn;

  else if (u == "millisecond"     ) return  genie::units::millisecond;
  else if (u == "microsecond"     ) return  genie::units::microsecond;
  else if (u == "nanosecond"      ) return  genie::units::nanosecond;
  else if (u == "picosecond"      ) return  genie::units::picosecond;
  else if (u == "s"               ) return  genie::units::s;
  else if (u == "ms"              ) return  genie::units::ms;
  else if (u == "us"              ) return  genie::units::us;
  else if (u == "ns"              ) return  genie::units::ns;
  else if (u == "ps"              ) return  genie::units::ps;
  else if (u == "hertz"           ) return  genie::units::hertz;
  else if (u == "kilohertz"       ) return  genie::units::kilohertz;
  else if (u == "megahertz"       ) return  genie::units::megahertz;
  else if (u == "gigahertz"       ) return  genie::units::gigahertz;
  else if (u == "Hz"              ) return  genie::units::Hz;
  else if (u == "kHz"             ) return  genie::units::kHz;
  else if (u == "MHz"             ) return  genie::units::MHz;
  else if (u == "GHz"             ) return  genie::units::GHz;

  else if (u == "qe"              ) return  genie::units::qe;
  else if (u == "qe_coulomb"      ) return  genie::units::qe_coulomb;

  else if (u == "electronvolt"    ) return  genie::units::electronvolt;
  else if (u == "kiloelectronvolt") return  genie::units::kiloelectronvolt;
  else if (u == "megaelectronvolt") return  genie::units::megaelectronvolt;
  else if (u == "teraelectronvolt") return  genie::units::teraelectronvolt;
  else if (u == "petaelectronvolt") return  genie::units::petaelectronvolt;
  else if (u == "eV"              ) return  genie::units::eV;
  else if (u == "keV"             ) return  genie::units::keV;
  else if (u == "MeV"             ) return  genie::units::MeV;
  else if (u == "TeV"             ) return  genie::units::TeV;
  else if (u == "PeV"             ) return  genie::units::PeV;

  else if (u == "gram"            ) return  genie::units::gram;
  else if (u == "milligram"       ) return  genie::units::milligram;
  else if (u == "kg"              ) return  genie::units::kg;
  else if (u == "g"               ) return  genie::units::g;
  else if (u == "mg"              ) return  genie::units::mg;

  else if (u == "radian"          ) return  genie::units::radian;
  else if (u == "milliradian"     ) return  genie::units::milliradian;
  else if (u == "degree"          ) return  genie::units::degree;
  else if (u == "steradian"       ) return  genie::units::steradian;
  else if (u == "rad"             ) return  genie::units::rad;
  else if (u == "mrad"            ) return  genie::units::mrad;
  else if (u == "sr"              ) return  genie::units::sr;
  else if (u == "deg"             ) return  genie::units::deg;

  else {
    LOG("Units", pWARN) << "Unknown units: " << u;
    return 1.;
  }
  return 1.;
}
//____________________________________________________________________________


