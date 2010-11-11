{

  double fXFmin = -1.0;
  double fXFmax =  0.5;
  double fPT2min = 0.001;
  double fPT2max = 0.6;

  double fracerr_xF1pi = 0.20;
  double fracerr_pT1pi = 0.03;


  TF1 * fBaryonXFpdf  = new TF1("fBaryonXFpdf",
                   "0.083*exp(-0.5*pow(x+0.385,2.)/0.131)", fXFmin, fXFmax);  
  TF1 * fBaryonPT2pdf = new TF1("fBaryonPT2pdf", 
                   "exp(-0.214-6.625*x)", fPT2min, fPT2max);  

  double fI0XFpdf  = fBaryonXFpdf ->Integral(fXFmin, fXFmax);
  double fI0PT2pdf = fBaryonPT2pdf->Integral(fPT2min,fPT2max);

  //
  // Now, define same function as above, but insert tweaking dials.
  // - For the xF  PDF add a parameter to shift the peak of the distribution.
  // - For the pT2 PDF add a parameter to shift the <pT2>
  // Include parameter to allow re-normalizing tweaked PDF to 1.

  double fDefPeakBaryonXF = -0.385;
  TF1 * fBaryonXFpdfTwk1  = new TF1("fBaryonXFpdfTwk1",
                   "[0]*0.083*exp(-0.5*pow(x-[1],2.)/0.131)", fXFmin, fXFmax);
  fBaryonXFpdfTwk1->SetParameter(0, 1.); // norm
  fBaryonXFpdfTwk1->SetParameter(1, fDefPeakBaryonXF*(1+fracerr_xF1pi));
  fBaryonXFpdfTwk1->SetParameter(0, fI0XFpdf/fBaryonXFpdfTwk1->Integral(fXFmin, fXFmax) );
  TF1 * fBaryonXFpdfTwk2  = new TF1("fBaryonXFpdfTwk2",
                   "[0]*0.083*exp(-0.5*pow(x-[1],2.)/0.131)", fXFmin, fXFmax);
  fBaryonXFpdfTwk2->SetParameter(0, 1.); // norm
  fBaryonXFpdfTwk2->SetParameter(1, fDefPeakBaryonXF*(1-fracerr_xF1pi));
  fBaryonXFpdfTwk2->SetParameter(0, fI0XFpdf/fBaryonXFpdfTwk2->Integral(fXFmin, fXFmax) );

  double fDefAvgPT2 = 1./6.625;
  TF1 * fBaryonPT2pdfTwk1 = new TF1("fBaryonPT2pdfTwk1", 
                   "[0]*exp(-0.214-x/[1])", fPT2min, fPT2max);  
  fBaryonPT2pdfTwk1->SetParameter(0, 1.); // norm
  fBaryonPT2pdfTwk1->SetParameter(1,fDefAvgPT2*(1+fracerr_pT1pi));
  fBaryonPT2pdfTwk1->SetParameter(0,fI0PT2pdf/fBaryonPT2pdfTwk1->Integral(fPT2min,fPT2max));
  TF1 * fBaryonPT2pdfTwk2 = new TF1("fBaryonPT2pdfTwk2", 
                   "[0]*exp(-0.214-x/[1])", fPT2min, fPT2max);  
  fBaryonPT2pdfTwk2->SetParameter(0, 1.); // norm
  fBaryonPT2pdfTwk2->SetParameter(1,fDefAvgPT2*(1-fracerr_pT1pi));
  fBaryonPT2pdfTwk2->SetParameter(0,fI0PT2pdf/fBaryonPT2pdfTwk2->Integral(fPT2min,fPT2max));

  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  fBaryonXFpdf->Draw();
  fBaryonXFpdfTwk1->SetLineStyle(kDashed);
  fBaryonXFpdfTwk2->SetLineStyle(kDashed);
  fBaryonXFpdfTwk1->Draw("same");
  fBaryonXFpdfTwk2->Draw("same");
  c1->Update();

  TCanvas * c2 = new TCanvas("c2","",120,120,600,600);
  fBaryonPT2pdf->Draw();
  fBaryonPT2pdfTwk1->SetLineStyle(kDashed);
  fBaryonPT2pdfTwk2->SetLineStyle(kDashed);
  fBaryonPT2pdfTwk1->Draw("same");
  fBaryonPT2pdfTwk2->Draw("same");
  c2->Update();
}
