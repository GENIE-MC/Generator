void glgr(double theta_weinberg)
{
 double sin8w  = TMath::Sin(theta_weinberg);
 double sin28w = TMath::Power(sin8w,2);

 cout << "GL  = " << (1./2.)  - (2./3.)* sin28w << endl;
 cout << "GR  = " <<          - (2./3.)* sin28w << endl;
 cout << "GL' = " << -(1./2.) + (1./3.)* sin28w << endl;
 cout << "GR' = " <<            (1./3.)* sin28w << endl;
}
