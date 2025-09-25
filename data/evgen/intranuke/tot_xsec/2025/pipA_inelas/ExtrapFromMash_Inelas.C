{
#include <iomanip>
  double xs1[21] = {
    12.54,
    573,
    715,
    522.8,
    359.2,
    430.3,
    408.1,
    401.4,
    399.9,
    393.9,
    384,
    371,
    356.5,
    351.5,
    345.5,
    340.5,
    338,
    339,
    342,
    346,
    349.5
  }
  double xs2;
  double ke;
  double A1 = 56;
  double A2[9] = { 1, 7, 14, 40, 63, 93, 120, 208, 209 };
  for (int j=0; j<9; j++) {
    cout << endl;
    cout << setw(4) << "A" << setw(5) << "KE" << "   " << "XSec" << endl;
    for (int i=0; i<21; i++) {
      xs2 = TMath::Power((A2[j]/A1),(.6667)) * xs1[i];
      ke = 50 * i;
      cout << setw(4) << A2[j] << setw(5) << ke << ".  " << setw(5) << xs2 << endl;
    }
  }

}
