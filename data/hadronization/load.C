{
  TTree bfc;
  bfc.ReadFile("./hBkwFwdCorrelation.data","nu/I:nuc/I:Wmin/F:Wmax/F:nF/F:nB/F:dnB/F");

  TTree mult;
  mult.ReadFile("./hMultiplicityVsW2.data","nu/I:nuc/I:tgt/I:W2/F:n/F:dn/F:m/I:xfh/I");

  TTree pi0c;
  pi0c.ReadFile("./hPi0HCorrelation.data","nu/I:nuc/I:Wmin/F:Wmax/F:n/I:npi0/F:dnpi0/F:nt/I");
}
