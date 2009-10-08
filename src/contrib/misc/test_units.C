//
// Test Units.h
//
// Costas Andreopoulos
//
//

#include "../..//Conventions/Units.h"

using namespace genie;

void test_units(void)
{
  cout << "1 m in GeV^-1 = " << units::m  << endl;
  cout << "1 F in GeV^-1 = " << units::fm << endl;

  double cm2m       = units::cm  / units::m;
  double mev2gev    = units::MeV / units::GeV;
  double to1E_38cm2 = 1. / (1E-38 * units::cm2);

  cout << "cm     -> m          ? => x" << cm2m       << endl;
  cout << "MeV    -> GeV        ? => x" << mev2gev    << endl;
  cout << "GeV^-2 -> 1E-38 cm^2 ? => x" << to1E_38cm2 << endl;
}
