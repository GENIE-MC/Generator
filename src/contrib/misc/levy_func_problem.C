//
// testing problems with Levy function used at AGKY/KNO
//
// INPUTS:
// c      -> Levy function parameter
// method -> 0 : compute P(n) from <n>P(n) = KNO(n/<n>)
//           1 : compute P(n) from a Poisson transformation of the
//               asymptotic scaling function (KNO/Levy) as in 
//               hep-ph/9608479 v1 29-Aug-1996, Eq.5
//        -> 2 : like 1 but psi in Eq.5 is not the Levy function
//               but a Gamma function
//
// Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
// University of Liverpool & STFC Rutherford Appleton Lab
//

double integrand              (double * x, double * par);
double prob_method0           (double n, double nav, double c);
double prob_method1           (double n, double nav, double c);
double prob_method2           (double n, double nav, double c);
double kno_func               (double z, double c);
void   plot_multiplicity_prob (double c);

//.................................................................
void levy_func_problem(double c, int method)
{
  cout << "Using method = " << method << ", Levy(c) = " << c << endl;
  
  const int nnav = 40;

  double nav_min = 0.25;
  double nav_max = 3.00;
  double dnav    = (nav_max-nav_min)/(nnav-1);

  // canvas for plotting <n>_out = f(<n>_in)
  TCanvas * canvas = new TCanvas("canvas","",20,20,700,700);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);  
  TH1F * frame = (TH1F*) canvas->DrawFrame(nav_min,nav_min,nav_max,nav_max);
  frame->Draw();

  double nav_in [nnav];
  double nav_out[nnav];

  for(int i=0; i<nnav; i++) {

    double nav_input = nav_min + i*dnav;

    cout << "** nav(in) = " << nav_input << endl;

    double minmult = 0;
    double maxmult = 10;
    int    nbins   = TMath::Nint(maxmult-minmult+1);

    TH1D * mult_prob = new TH1D("mult_prob",
      "hadronic multiplicity distribution", nbins, minmult-0.5, maxmult+0.5);
    mult_prob->SetDirectory(0);

    int nbins = mult_prob->FindBin(maxmult);
    for(int j=1; j<=nbins; j++) {
       double n = mult_prob->GetBinCenter(j);  // bin centre
              
       double P = 0;
       if      (method==0) P = prob_method0(n,nav_input,c); 
       else if (method==1) P = prob_method1(n,nav_input,c); 
       else if (method==2) P = prob_method2(n,nav_input,c); 

       cout << "   n = " << n << " -> P(n) = " << P << endl;       
       mult_prob->Fill(n,P);
    }    
    double nav_output = mult_prob->GetMean();
    delete mult_prob;

    cout << "   ** nav(out) = " << nav_output << endl;
    
    nav_in[i]  = nav_input;
    nav_out[i] = nav_output;
  }

  TGraph * gr = new TGraph(nnav,nav_in,nav_out);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1);

  TGraph * id = new TGraph(nnav,nav_in,nav_in);
  id->SetLineStyle(1);
  id->SetLineColor(2);

  gr->Draw("P");
  id->Draw("L");

  canvas->Update();
}
//.................................................................
void plot_multiplicity_prob(double c)
{
  const int nnav = 2;
  double nav[nnav] = { 0.5, 4. };
  TH1D * mult_prob[nnav][2];

  double minmult = 0;
  double maxmult = 10;
  int    nbins   = TMath::Nint(maxmult-minmult+1);

  // canvas for plotting sample multiplicity probability distributions
  TCanvas * canvas = new TCanvas("canvas","",120,120,800,800);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);  

  for(int i=0; i<nnav; i++) {

    double nav_input = nav[i];

    cout << "** nav = " << nav_input << endl;

    for(int method=0; method<=1; method++) {
      mult_prob[i][method] = new TH1D("mult_prob",
        "hadronic multiplicity distribution", nbins, minmult-0.5, maxmult+0.5);
      mult_prob[i][method]->SetDirectory(0);

      int nbins = mult_prob[i][method]->FindBin(maxmult);
      for(int j=1; j<=nbins; j++) {
        double n = mult_prob[i][method]->GetBinCenter(j);  // bin centre
              
        double P = 0;
        if      (method==0) P = prob_method0(n,nav_input,c); 
        else if (method==1) P = prob_method1(n,nav_input,c); 

        cout << "   n = " << n << " -> P(n) = " << P << endl;       
        mult_prob[i][method]->Fill(n,P);
      }//bins
      mult_prob[i][method]->Scale(1/mult_prob[i][method]->Integral("width"));
   }//meth    
  }//nnav

  for(int i=0; i<nnav; i++) {
    for(int method=0; method<=1; method++) {
         mult_prob[i][method]->SetLineStyle(1);
         mult_prob[i][method]->SetLineWidth(1);
         mult_prob[i][method]->SetLineColor(method+1);
    }
  }

  canvas->cd();
  canvas->Divide(1,nnav);

  for(int i=0; i<nnav; i++) {
    canvas->cd(1+i);
    mult_prob[i][0]->Draw();
    mult_prob[i][1]->Draw("SAME");
  }
  canvas->Update();
}
//.................................................................
double prob_method2(double n, double nav, double c)
{
// Compute P(n) as in hep-ph/9608479 v1 29-Aug-1996, eq.5
//
  TF1 * intg = new TF1("integrand",integrand,0,10,3);

  intg->SetParameter(0,n);
  intg->SetParameter(1,nav);
  intg->SetParameter(2,-999999);

  double P = intg->Integral(0,10);
  
  delete intg;
  
  return P;
}
//.................................................................
double prob_method1(double n, double nav, double c)
{
// Compute P(n) as in hep-ph/9608479 v1 29-Aug-1996, eq.5
//
  TF1 * intg = new TF1("integrand",integrand,0,10,3);

  intg->SetParameter(0,n);
  intg->SetParameter(1,nav);
  intg->SetParameter(2,c);

  double P = intg->Integral(0,10);
  
  delete intg;
  
  return P;
}
//.................................................................
double prob_method0(double n, double nav, double c)
{
// compute P(n) from <n>P(n) = KNO(n/<n>)
//
  double z   = n/nav;
  double kno = kno_func(z,c);
  double P   = kno/nav;
  return P;
}
//.................................................................
double integrand(double * x, double * par)
{
// preasymptotic multiplicity distribution reconstructed from the
// asymptotic scaling function via Poisson transform

  // inputs
  double xx = x[0];
      
  // parameters
  double n   = par[0];
  double nav = par[1];
  double c   = par[2];

  // compute integrand

  double psi  = 0;
  if(c<0) {
     if(n==0) psi = 0;
     else     psi = TMath::Gamma(xx);
  }
  else    psi = kno_func(xx,c);

  double fct  = TMath::Factorial((int)n);
  double pwr  = TMath::Power(nav*xx,n);
  double expn = TMath::Exp(-nav*xx);
  
  double f = psi * (pwr/fct) * expn;
  return f;
}
//.................................................................
double kno_func(double z, double c)
{
// z reduced multiplicity:
// c KNO param  

  double x   = c*z+1;
  double psi = 2*TMath::Exp(-c)*TMath::Power(c,x)/TMath::Gamma(x);

  return psi;
}
//.................................................................
