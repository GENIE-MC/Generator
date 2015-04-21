{
#include <iomanip>
  double xs1[21] = {
    0.0000E+00,
    0.5400E+02,
    0.9700E+02,
    0.1105E+03,
    0.1005E+03,
    0.9200E+02,
    0.8700E+02,
    0.8350E+02,
    0.2498E+02,
    0.2167E+02,
    0.2031E+02,
    0.2187E+02,
    0.2038E+02,
    0.1464E+02,
    0.1111E+02,
    0.1160E+02,
    0.1468E+02,
    0.1735E+02,
    0.1615E+02,
    0.1184E+02,
    0.8496E+01
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
