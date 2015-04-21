{
  #include <iomanip>
  double xs1[21] = {
    0.2762E+02,
    0.1381E+04,
    0.1945E+04,
    0.1879E+04,
    0.1757E+04,
    0.1595E+04,
    0.1445E+04,
    0.1288E+04,
    0.1218E+04,
    0.1159E+04,
    0.1111E+04,
    0.1072E+04,
    0.1042E+04,
    0.1019E+04,
    0.1003E+04,
    0.9927E+03,
    0.9874E+03,
    0.9865E+03,
    0.9894E+03,
    0.9956E+03,
    0.1004E+04
  }
  double xs2;
  double ke;
  double A1 = 56;
  double A2[9] = { 1, 7, 14, 40, 63, 93, 120, 208, 209 };
  for (int j=0; j<9; j++) {
    cout << endl;
    cout << setw(4) << "A" << setw(5) << "KE" << "   " << "TotXS" << endl;
    for (int i=0; i<21; i++) {
      xs2 = TMath::Power((A2[j]/A1),(.6667)) * xs1[i];
      ke = 50 * i;
      cout << setw(4) << A2[j] << setw(5) << ke << ".  " << setw(5) << xs2 << endl;
    }
  }
}
