{
  double p32iso   = 0.50;
  double p12iso   = 0.50;
  double p32rs    = 0.75;
  double p12rs    = 0.25;
 
  const int ndial = 11;
  const double dial_step = 2./(ndial-1);

  double epsilon = 1E-6;

  const int ncosthetabins=1000;

  TH1D * h[ndial];

  for(int i=0; i<ndial; i++) {
    double dial = -1. + i*dial_step;
    h[i] = new TH1D("a","",ncosthetabins,-1.,1.);
    h[i]->SetDirectory(0);
    if( TMath::Abs(dial)    < epsilon )  { h[i]->SetLineWidth(3); }
    if( TMath::Abs(dial-1.) < epsilon )  { h[i]->SetLineWidth(3); h[i]->SetLineStyle(kDashed); }
    for(int j=1; j<=h[i]->GetNbinsX(); j++) {
      double costheta = h[i]->GetBinCenter(j);
      double P2       = 0.5 * (3.*costheta*costheta - 1.);
      double Wiso     = 1 - p32iso * P2 + p12iso * P2; // = 1.0
      double Wrs      = 1 - p32rs  * P2 + p12rs  * P2;
      double Wdef     = Wiso;
      double Wtwk     = dial*Wrs + (1-dial)*Wiso;
      double wght = 1.;
      if(Wdef>0. && Wtwk>0.) {
        wght = Wtwk/Wdef;
      }
      h[i]->SetBinContent(j,wght);
    }//j
  }//i

  for(int i=0; i<ndial; i++) {
    if(i==0) {
       h[i]->Draw();
       h[i]->GetXaxis()->SetTitle("cos#theta_{#pi}");
       h[i]->GetYaxis()->SetTitle("W(cos#theta_{#pi})");
    }
    else {     
       h[i]->Draw("same");
    }
  }

}

