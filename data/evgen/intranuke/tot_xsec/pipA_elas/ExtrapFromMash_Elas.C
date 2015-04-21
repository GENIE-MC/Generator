{
#include <iomanip>
  double xs1[21] = {
    6.215,
    310.8,
    672.6,
    687.4,
    694,
    667.7,
    610.6,
    479.7,
    407,
    346,
    295.6,
    254.6,
    222.1,
    191.8,
    168.1,
    150.2,
    137.4,
    129,
    124.4,
    123.1,
    124.5
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
