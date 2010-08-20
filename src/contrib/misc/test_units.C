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

  const double cm2m       = units::cm  / units::m;
  const double mev2gev    = units::MeV / units::GeV;
  const double to1E_38cm2 = 1. / (1E-38 * units::cm2);
  const double fm2tomb    = units::fm2 / units::mb;

  cout << "cm     -> m          ? => x" << cm2m       << endl;
  cout << "MeV    -> GeV        ? => x" << mev2gev    << endl;
  cout << "GeV^-2 -> 1E-38 cm^2 ? => x" << to1E_38cm2 << endl;
  cout << "fm^2   -> mbarn      ? => x" << fm2tomb    << endl;
  cout << "MeV    -> GeV        ? => x" << units::MeV   << endl;
}
