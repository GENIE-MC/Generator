//
// Data vs GENIE neutrino F2 comparisons
//
// Contributed to GENIE/ValidationTools by Debdatta Bhattacharya (Pittsburgh Univ.)
// Jan 15h, 2009
//
// Revisions:
//


{
gROOT->Reset();
gStyle->SetOptStat(0000);

const int nfile = 3672;

// original f2,xf3 from neugen
ifstream base_neugen("data/extract_sf.txt");

Double_t row1_base[nfile];
Double_t row2_base[nfile];
Double_t row3_base[nfile];
Double_t row4_base[nfile];
Double_t row5_base[nfile];
Double_t row6_base[nfile];
Double_t row7_base[nfile];
Double_t row8_base[nfile];
Double_t row9_base[nfile];
Double_t row10_base[nfile];

for(int i=0;i<nfile;i++){
  base_neugen>>row1_base[i]>>row2_base[i]>>row3_base[i]>>row4_base[i]>>row5_base[i]>>row6_base[i]>>row7_base[i]>>row8_base[i]>>row9_base[i]>>row10_base[i];
  //cout<<"  "<<row1_base[i]<<"  "<<row2_base[i]<<"  "<<row3_base[i]<<"  "<<row4_base[i]<<"  "<<row5_base[i]<<endl;
}

ifstream lo_neugen("data/lo_extract_sf.txt");

Double_t lo_row1_base[nfile];
Double_t lo_row2_base[nfile];
Double_t lo_row3_base[nfile];
Double_t lo_row4_base[nfile];
Double_t lo_row5_base[nfile];

for(int i=0;i<nfile;i++){
  lo_neugen>>lo_row1_base[i]>>lo_row2_base[i]>>lo_row3_base[i]>>lo_row4_base[i]>>lo_row5_base[i];
  //cout<<"  "<<row1_base[i]<<"  "<<row2_base[i]<<"  "<<row3_base[i]<<"  "<<row4_base[i]<<"  "<<row5_base[i]<<endl;
}


ifstream newlo_neugen("data/newlo_extract_sf.txt");

Double_t newlo_row1_base[nfile];
Double_t newlo_row2_base[nfile];
Double_t newlo_row3_base[nfile];
Double_t newlo_row4_base[nfile];
Double_t newlo_row5_base[nfile];

for(int i=0;i<nfile;i++){
  newlo_neugen>>newlo_row1_base[i]>>newlo_row2_base[i]>>newlo_row3_base[i]>>newlo_row4_base[i]>>newlo_row5_base[i];
  //cout<<"  "<<row1_base[i]<<"  "<<row2_base[i]<<"  "<<row3_base[i]<<"  "<<row4_base[i]<<"  "<<row5_base[i]<<endl;
}


//
// nutev
//

ifstream nutev("data/final_nutev_xf3.txt",ios::in);

int nnutev = 75;
Double_t nutev1[nnutev];
Double_t nutev2[nnutev];
Double_t nutev3[nnutev];
Double_t nutev4[nnutev];

for(int i=0;i<nnutev;i++){
  nutev>>nutev1[i]>>nutev2[i]>>nutev3[i]>>nutev4[i];
  //cout<<"  "<<nutev1[i]<<"  "<<nutev2[i]<<"  "<<nutev3[i]<<"  "<<nutev4[i]<<endl;
}

//
// cdhsw
//

ifstream cdhsw("data/final_cdhsw_xf3.txt",ios::in);

int ncdhsw = 143;
Double_t cdhsw1[ncdhsw];
Double_t cdhsw2[ncdhsw];
Double_t cdhsw3[ncdhsw];
Double_t cdhsw4[ncdhsw];

for(int i=0;i<ncdhsw;i++){
  cdhsw>>cdhsw1[i]>>cdhsw2[i]>>cdhsw3[i]>>cdhsw4[i];
  //cout<<"  "<<cdhsw1[i]<<"  "<<cdhsw2[i]<<"  "<<cdhsw3[i]<<"  "<<cdhsw4[i]<<endl;
}

//
// ccfr
//

ifstream ccfr("data/ccfr_f2xf3fe.txt",ios::in);

int nccfr = 116;
Double_t ccfr1[nccfr];
Double_t ccfr2[nccfr];
Double_t ccfr3[nccfr];
Double_t ccfr4[nccfr];
Double_t ccfr5[nccfr];
Double_t ccfr6[nccfr];
Double_t ccfr7[nccfr];

for(int i=0;i<nccfr;i++){
  ccfr>>ccfr1[i]>>ccfr2[i]>>ccfr3[i]>>ccfr4[i]>>ccfr5[i]>>ccfr6[i]>>ccfr7[i];
  //cout<<"  "<<ccfr1[i]<<"  "<<ccfr2[i]<<"  "<<ccfr3[i]<<"  "<<ccfr4[i]<<endl;
}


int num = 306;
float x[12];

x[0]  = 0.015;
x[1]  = 0.045;
x[2]  = 0.080;
x[3]  = 0.125;
x[4]  = 0.175;
x[5]  = 0.225;
x[6]  = 0.275;
x[7]  = 0.350;
x[8]  = 0.450;
x[9]  = 0.550;
x[10] = 0.650;
x[11] = 0.750;

// base neugen sf

Double_t q2[12][num];
Double_t f2[12][num];

Int_t count[12];

Int_t num1[12]={5,3,2,5,3,2,5,3,2,5,3,2};

for(int k=0;k<12;k++){
  count[k]=0;
}

for(int k=0;k<nfile;k++){
  for(int l=0;l<12;l++){
    if(fabs(row1_base[k]-x[l])<.001){

      q2[l][count[l]] = row2_base[k];
      f2[l][count[l]] = num1[l]*row4_base[k];

      count[l]++;
    }
  }
}


// lo neugen sf

Double_t lo_q2[12][num];
Double_t lo_f2[12][num];

Int_t lo_count[12];

for(int k=0;k<12;k++){
  lo_count[k]=0;
}

for(int k=0;k<nfile;k++){
  for(int l=0;l<12;l++){
    if(fabs(lo_row1_base[k]-x[l])<.001){

      lo_q2[l][lo_count[l]] = lo_row2_base[k];
      lo_f2[l][lo_count[l]] = num1[l]*lo_row4_base[k];

      lo_count[l]++;
    }
  }
}


// newlo neugen sf

Double_t newlo_q2[12][num];
Double_t newlo_f2[12][num];

Int_t newlo_count[12];

for(int k=0;k<12;k++){
  newlo_count[k]=0;
}

for(int k=0;k<nfile;k++){
  for(int l=0;l<12;l++){
    if(fabs(newlo_row1_base[k]-x[l])<.001){

      newlo_q2[l][newlo_count[l]] = newlo_row2_base[k];
      newlo_f2[l][newlo_count[l]] = num1[l]*newlo_row4_base[k];

      newlo_count[l]++;
    }
  }
}



Double_t nutev_q2[12][num];
Double_t nutev_xf3[12][num];
Double_t nutev_q2err[12][num];
Double_t nutev_err[12][num];
Int_t nutev_count[12];

for(int k=0;k<12;k++){
  nutev_count[k]=0;
}

for(int k=0;k<nnutev;k++){
  for(int l=0;l<12;l++){
    if(fabs(nutev1[k]-x[l])<.001){

      nutev_q2[l][nutev_count[l]] = nutev2[k];
      nutev_q2err[l][nutev_count[l]] = 0;
      nutev_xf3[l][nutev_count[l]] = num1[l]*nutev3[k];
      nutev_err[l][nutev_count[l]] = num1[l]*nutev4[k];

      nutev_count[l]++;
    }
  }
}

//cdhsw
Double_t cdhsw_q2[12][num];
Double_t cdhsw_xf3[12][num];
Double_t cdhsw_q2err[12][num];
Double_t cdhsw_err[12][num];
Int_t cdhsw_count[12];

for(int k=0;k<12;k++){
  cdhsw_count[k]=0;
}

for(int k=0;k<ncdhsw;k++){
  for(int l=0;l<12;l++){
    if(fabs(cdhsw1[k]-x[l])<.001){

      cdhsw_q2[l][cdhsw_count[l]] = cdhsw2[k];
      cdhsw_q2err[l][cdhsw_count[l]] = 0;
      cdhsw_xf3[l][cdhsw_count[l]] = num1[l]*cdhsw3[k];
      cdhsw_err[l][cdhsw_count[l]] = num1[l]*cdhsw4[k];
      cdhsw_count[l]++;
    }
  }
}


//ccfr
Double_t ccfr_q2[12][num];
Double_t ccfr_xf3[12][num];
Double_t ccfr_q2err[12][num];
Double_t ccfr_err[12][num];
Int_t ccfr_count[12];

for(int k=0;k<12;k++){
  ccfr_count[k]=0;
}

for(int k=0;k<nccfr;k++){
  for(int l=0;l<12;l++){
    if(fabs(ccfr1[k]-x[l])<.001){

      ccfr_q2[l][ccfr_count[l]] = ccfr2[k];
      ccfr_q2err[l][ccfr_count[l]] = 0;
      ccfr_xf3[l][ccfr_count[l]] = num1[l]*ccfr5[k];
      ccfr_err[l][ccfr_count[l]] = num1[l]*ccfr6[k];
      ccfr_count[l]++;
    }
  }
}

Int_t count_sp1=0;
Int_t count_sp2=0;

Double_t temp_q21[num];
Double_t temp_sp1[num];
Double_t temp_errsp1[num];

Double_t temp_q22[num];
Double_t temp_sp2[num];
Double_t temp_errsp2[num];

//special ccfr case - x=0.015
for(int k=0;k<nccfr;k++){

  if(fabs(ccfr1[k]-0.0125)<.001){

    temp_q21[count_sp1]= ccfr2[k];
    temp_sp1[count_sp1]= ccfr5[k];
    temp_errsp1[count_sp1]= ccfr6[k];

    //cout<<temp_sp1[count_sp1]<<endl;
    count_sp1++;
  }

  if(fabs(ccfr1[k]-0.0175)<.001){

    temp_q22[count_sp2]= ccfr2[k];
    temp_sp2[count_sp2]= ccfr5[k];
    temp_errsp2[count_sp2]= ccfr6[k];

    //cout<<temp_sp2[count_sp2]<<endl;
    count_sp2++;
  }
}

ccfr_count[0] = count_sp1;

for(int k=0;k<ccfr_count[0];k++){

  ccfr_q2[0][k] = temp_q21[k];
  ccfr_q2err[0][k] = 0;

  ccfr_xf3[0][k] = (20*(temp_sp1[k]+temp_sp2[k]))/2;
  ccfr_err[0][k] = (20*(temp_errsp1[k]+temp_errsp2[k]))/2;
}


count_sp1=0;
count_sp2=0;

//special ccfr case - x=0.0425
for(int k=0;k<nccfr;k++){

  if(fabs(ccfr1[k]-0.035)<.001){

    temp_q21[count_sp1]= ccfr2[k];
    temp_sp1[count_sp1]= ccfr5[k];
    temp_errsp1[count_sp1]= ccfr6[k];

    //cout<<temp_sp1[count_sp1]<<endl;
    count_sp1++;
  }

  if(fabs(ccfr1[k]-0.050)<.001){

    temp_q22[count_sp2]= ccfr2[k];
    temp_sp2[count_sp2]= ccfr5[k];
    temp_errsp2[count_sp2]= ccfr6[k];

    //cout<<temp_sp2[count_sp2]<<endl;
    count_sp2++;
  }
}

ccfr_count[1] = count_sp1;

for(int k=0;k<ccfr_count[1];k++){

  ccfr_q2[1][k] = temp_q21[k];
  ccfr_q2err[1][k] = 0;

  ccfr_xf3[1][k] = (5*(temp_sp1[k]+temp_sp2[k]))/2;
  ccfr_err[1][k] = (5*(temp_errsp1[k]+temp_errsp2[k]))/2;

  //cout<<ccfr_q2[1][k]<<"  "<<ccfr_xf3[1][k]<<endl;
}

count_sp1=0;
count_sp2=0;

//special ccfr case - x=0.0425
for(int k=0;k<nccfr;k++){

  if(fabs(ccfr1[k]-0.070)<.001){

    temp_q21[count_sp1]= ccfr2[k];
    temp_sp1[count_sp1]= ccfr5[k];
    temp_errsp1[count_sp1]= ccfr6[k];

    //cout<<temp_sp1[count_sp1]<<endl;
    count_sp1++;
  }

  if(fabs(ccfr1[k]-0.090)<.001){

    temp_q22[count_sp2]= ccfr2[k];
    temp_sp2[count_sp2]= ccfr5[k];
    temp_errsp2[count_sp2]= ccfr6[k];

    //cout<<temp_sp2[count_sp2]<<endl;
    count_sp2++;
  }
}

ccfr_count[2] = count_sp1;

for(int k=0;k<ccfr_count[2];k++){

  ccfr_q2[2][k] = temp_q21[k];
  ccfr_q2err[2][k] = 0;

  ccfr_xf3[2][k] = ((temp_sp1[k]+temp_sp2[k]))/2;
  ccfr_err[2][k] = ((temp_errsp1[k]+temp_errsp2[k]))/2;

  //cout<<ccfr_q2[2][k]<<"  "<<ccfr_xf3[2][k]<<endl;
}


count_sp1=0;
count_sp2=0;

//special ccfr case - x=0.125
for(int k=0;k<nccfr;k++){

  if(fabs(ccfr1[k]-0.110)<.001){

    temp_q21[count_sp1]= ccfr2[k];
    temp_sp1[count_sp1]= ccfr5[k];
    temp_errsp1[count_sp1]= ccfr6[k];

    //cout<<temp_sp1[count_sp1]<<endl;
    count_sp1++;
  }

  if(fabs(ccfr1[k]-0.140)<.001){

    temp_q22[count_sp2]= ccfr2[k];
    temp_sp2[count_sp2]= ccfr5[k];
    temp_errsp2[count_sp2]= ccfr6[k];

    //cout<<temp_sp2[count_sp2]<<endl;
    count_sp2++;
  }
}

ccfr_count[3] = count_sp1;

for(int k=0;k<ccfr_count[3];k++){

  ccfr_q2[3][k] = temp_q21[k];
  ccfr_q2err[3][k] = 0;

  ccfr_xf3[3][k] = (5*(temp_sp1[k]+temp_sp2[k]))/2;
  ccfr_err[3][k] = (5*(temp_errsp1[k]+temp_errsp2[k]))/2;

  //cout<<ccfr_q2[2][k]<<"  "<<ccfr_xf3[2][k]<<endl;
}


count_sp1=0;
count_sp2=0;

//special ccfr case - x=0.180
for(int k=0;k<nccfr;k++){

  if(fabs(ccfr1[k]-0.180)<.001){

    ccfr_q2[4][ccfr_count[4]] = ccfr2[k];
    ccfr_q2err[4][ccfr_count[4]] = 0;

    ccfr_xf3[4][ccfr_count[4]] = 2*ccfr5[k];
    ccfr_err[4][ccfr_count[4]] = 2*ccfr6[k];

    ccfr_count[4]++;
    }
}



TCanvas c2("c2","",0,0,1300,950);
c2->Divide(2,2);

c2->cd(1);
TH1F *a1 = new TH1F("a1","",10,0.005,2000);
a1->SetAxisRange(0.05,20,"y");
a1->Draw();
a1->GetXaxis()->SetTitle("Q2");
a1->GetYaxis()->SetTitle("F2");
gPad->SetLogy();
gPad->SetLogx();

TGraph *gr_015 = new TGraph(count[0],q2[0],f2[0]);
gr_015->SetLineWidth(2);
gr_015->Draw("C");

TGraph *lo_gr_015 = new TGraph(lo_count[0],lo_q2[0],lo_f2[0]);
lo_gr_015->SetLineWidth(2);
lo_gr_015->SetLineColor(kRed);
lo_gr_015->SetLineStyle(2);
lo_gr_015->Draw("C");

TGraph *newlo_gr_015 = new TGraph(newlo_count[0],newlo_q2[0],newlo_f2[0]);
newlo_gr_015->SetLineWidth(2);
newlo_gr_015->SetLineColor(kRed);
newlo_gr_015->Draw("C");

TLatex *tex = new TLatex();
tex->SetTextAlign(22); //centers the text on the chosen coordinates
tex->SetTextSize(0.045);//a reasonably large size
tex->DrawLatex(28.4, 12.8, "x = 0.015(5*F2)");

TGraph *gr_045 = new TGraph(count[1],q2[1],f2[1]);
gr_045->SetLineWidth(2);
gr_045->Draw("C");

TGraph *lo_gr_045 = new TGraph(lo_count[1],lo_q2[1],lo_f2[1]);
lo_gr_045->SetLineWidth(2);
lo_gr_045->SetLineColor(kRed);
lo_gr_045->SetLineStyle(2);
lo_gr_045->Draw("C");

TGraph *newlo_gr_045 = new TGraph(newlo_count[1],newlo_q2[1],newlo_f2[1]);
newlo_gr_045->SetLineWidth(2);
newlo_gr_045->SetLineColor(kRed);
newlo_gr_045->Draw("C");

tex->DrawLatex(171.55, 5.825, "x = 0.045(3*F2)");

TGraph *gr_080 = new TGraph(count[2],q2[2],f2[2]);
gr_080->SetLineWidth(2);
gr_080->Draw("C");

TGraph *lo_gr_080 = new TGraph(lo_count[2],lo_q2[2],lo_f2[2]);
lo_gr_080->SetLineWidth(2);
lo_gr_080->SetLineColor(kRed);
lo_gr_080->SetLineStyle(2);
lo_gr_080->Draw("C");

TGraph *newlo_gr_080 = new TGraph(newlo_count[2],newlo_q2[2],newlo_f2[2]);
newlo_gr_080->SetLineWidth(2);
newlo_gr_080->SetLineColor(kRed);
newlo_gr_080->Draw("C");

tex->DrawLatex(74.85, 2.076, "x = 0.080(2*F2)");

TLegend *leg1 = new TLegend(0.356,0.262,0.81668,0.541984);
leg1->AddEntry(gr_015,"standard neugen F2","l");
leg1->AddEntry(lo_gr_015,"turn off B-Y correction","l");
leg1->AddEntry(newlo_gr_015,"+ turn off q2 freezing","l");
leg1->Draw();

c2->cd(2);
TH1F *a2 = new TH1F("a2","",10,0.005,2000);
a2->SetAxisRange(0.05,20,"y");
a2->Draw();
a2->GetXaxis()->SetTitle("Q2");
a2->GetYaxis()->SetTitle("F2");
gPad->SetLogy();
gPad->SetLogx();

TGraph *gr_125 = new TGraph(count[3],q2[3],f2[3]);
gr_125->SetLineWidth(2);
gr_125->Draw("C");

TGraph *lo_gr_125 = new TGraph(lo_count[3],lo_q2[3],lo_f2[3]);
lo_gr_125->SetLineWidth(2);
lo_gr_125->SetLineColor(kRed);
lo_gr_125->SetLineStyle(2);
lo_gr_125->Draw("C");

TGraph *newlo_gr_125 = new TGraph(newlo_count[3],newlo_q2[3],newlo_f2[3]);
newlo_gr_125->SetLineWidth(2);
newlo_gr_125->SetLineColor(kRed);
newlo_gr_125->Draw("C");

tex->DrawLatex(80, 9, "x = 0.125(5*F2)");

TGraph *gr_175 = new TGraph(count[4],q2[4],f2[4]);
gr_175->SetLineWidth(2);
gr_175->Draw("C");

TGraph *lo_gr_175 = new TGraph(lo_count[4],lo_q2[4],lo_f2[4]);
lo_gr_175->SetLineWidth(2);
lo_gr_175->SetLineColor(kRed);
lo_gr_175->SetLineStyle(2);
lo_gr_175->Draw("C");

TGraph *newlo_gr_175 = new TGraph(newlo_count[4],newlo_q2[4],newlo_f2[4]);
newlo_gr_175->SetLineWidth(2);
newlo_gr_175->SetLineColor(kRed);
newlo_gr_175->Draw("C");

tex->DrawLatex(119.76, 4.23, "x = 0.175(3*F2)");

TGraph *gr_225 = new TGraph(count[5],q2[5],f2[5]);
gr_225->SetLineWidth(2);
gr_225->Draw("C");

TGraph *lo_gr_225 = new TGraph(lo_count[5],lo_q2[5],lo_f2[5]);
lo_gr_225->SetLineWidth(2);
lo_gr_225->SetLineColor(kRed);
lo_gr_225->SetLineStyle(2);
lo_gr_225->Draw("C");

TGraph *newlo_gr_225 = new TGraph(newlo_count[5],newlo_q2[5],newlo_f2[5]);
newlo_gr_225->SetLineWidth(2);
newlo_gr_225->SetLineColor(kRed);
newlo_gr_225->Draw("C");

tex->DrawLatex(95.99, 1.32, "x = 0.225(2*F2)");


c2->cd(3);
TH1F *a3 = new TH1F("a3","",10,0.005,2000);
a3->SetAxisRange(0.05,20,"y");
a3->Draw();
a3->GetXaxis()->SetTitle("Q2");
a3->GetYaxis()->SetTitle("F2");
gPad->SetLogy();
gPad->SetLogx();

TGraph *gr_275 = new TGraph(count[6],q2[6],f2[6]);
gr_275->SetLineWidth(2);
gr_275->Draw("C");

TGraph *lo_gr_275 = new TGraph(lo_count[6],lo_q2[6],lo_f2[6]);
lo_gr_275->SetLineWidth(2);
lo_gr_275->SetLineColor(kRed);
lo_gr_275->SetLineStyle(2);
lo_gr_275->Draw("C");

TGraph *newlo_gr_275 = new TGraph(newlo_count[6],newlo_q2[6],newlo_f2[6]);
newlo_gr_275->SetLineWidth(2);
newlo_gr_275->SetLineColor(kRed);
newlo_gr_275->Draw("C");

tex->DrawLatex(79.106,6.16, "x = 0.275(5*F2)");

TGraph *gr_350 = new TGraph(count[7],q2[7],f2[7]);
gr_350->SetLineWidth(2);
gr_350->Draw("C");

TGraph *lo_gr_350 = new TGraph(lo_count[7],lo_q2[7],lo_f2[7]);
lo_gr_350->SetLineWidth(2);
lo_gr_350->SetLineColor(kRed);
lo_gr_350->SetLineStyle(2);
lo_gr_350->Draw("C");

TGraph *newlo_gr_350 = new TGraph(newlo_count[7],newlo_q2[7],newlo_f2[7]);
newlo_gr_350->SetLineWidth(2);
newlo_gr_350->SetLineColor(kRed);
newlo_gr_350->Draw("C");

tex->DrawLatex(202.5, 2.45, "x = 0.350(3*F2)");

TGraph *gr_450 = new TGraph(count[8],q2[8],f2[8]);
gr_450->SetLineWidth(2);
gr_450->Draw("C");

TGraph *lo_gr_450 = new TGraph(lo_count[8],lo_q2[8],lo_f2[8]);
lo_gr_450->SetLineWidth(2);
lo_gr_450->SetLineColor(kRed);
lo_gr_450->SetLineStyle(2);
lo_gr_450->Draw("C");

TGraph *newlo_gr_450 = new TGraph(newlo_count[8],newlo_q2[8],newlo_f2[8]);
newlo_gr_450->SetLineWidth(2);
newlo_gr_450->SetLineColor(kRed);
newlo_gr_450->Draw("C");

tex->DrawLatex(80, 0.4545, "x = 0.450(2*F2)");

c2->cd(4);
TH1F *a4 = new TH1F("a4","",10,0.005,2000);
a4->SetAxisRange(0.01,20,"y");
a4->Draw();
a4->GetXaxis()->SetTitle("Q2");
a4->GetYaxis()->SetTitle("F2");
gPad->SetLogy();
gPad->SetLogx();

TGraph *gr_550 = new TGraph(count[9],q2[9],f2[9]);
gr_550->SetLineWidth(2);
gr_550->Draw("C");

TGraph *lo_gr_550 = new TGraph(lo_count[9],lo_q2[9],lo_f2[9]);
lo_gr_550->SetLineWidth(2);
lo_gr_550->SetLineColor(kRed);
lo_gr_550->SetLineStyle(2);
lo_gr_550->Draw("C");

TGraph *newlo_gr_550 = new TGraph(newlo_count[9],newlo_q2[9],newlo_f2[9]);
newlo_gr_550->SetLineWidth(2);
newlo_gr_550->SetLineColor(kRed);
newlo_gr_550->Draw("C");

tex->DrawLatex(107.223, 2.000, "x = 0.550(5*F2)");

TGraph *gr_650 = new TGraph(count[10],q2[10],f2[10]);
gr_650->SetLineWidth(2);
gr_650->Draw("C");

TGraph *lo_gr_650 = new TGraph(lo_count[10],lo_q2[10],lo_f2[10]);
lo_gr_650->SetLineWidth(2);
lo_gr_650->SetLineColor(kRed);
lo_gr_650->SetLineStyle(2);
lo_gr_650->Draw("C");

TGraph *newlo_gr_650 = new TGraph(newlo_count[10],newlo_q2[10],newlo_f2[10]);
newlo_gr_650->SetLineWidth(2);
newlo_gr_650->SetLineColor(kRed);
newlo_gr_650->Draw("C");

tex->DrawLatex(252.645, 0.4062, "x = 0.650(3*F2)");

TGraph *gr_750 = new TGraph(count[11],q2[11],f2[11]);
gr_750->SetLineWidth(2);
gr_750->Draw("C");

TGraph *lo_gr_750 = new TGraph(lo_count[11],lo_q2[11],lo_f2[11]);
lo_gr_750->SetLineWidth(2);
lo_gr_750->SetLineColor(kRed);
lo_gr_750->SetLineStyle(2);
lo_gr_750->Draw("C");

TGraph *newlo_gr_750 = new TGraph(newlo_count[11],newlo_q2[11],newlo_f2[11]);
newlo_gr_750->SetLineWidth(2);
newlo_gr_750->SetLineColor(kRed);
newlo_gr_750->Draw("C");

tex->DrawLatex(83.603, 0.0303, "x = 0.750(2*F2)");


c2->Print("f2_nu.eps");


}



/*TGraphErrors *grnutev_015 = new TGraphErrors(nutev_count[0],nutev_q2[0],nutev_xf3[0],nutev_q2err[0],nutev_err[0]);
grnutev_015->Draw("P*");
grnutev_015->SetMarkerColor(kMagenta);
grnutev_015->SetMarkerStyle(20);
grnutev_015->SetMarkerSize(0.9);
grnutev_015->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_015 = new TGraphErrors(cdhsw_count[0],cdhsw_q2[0],cdhsw_xf3[0],cdhsw_q2err[0],cdhsw_err[0]);
grcdhsw_015->Draw("P*");
grcdhsw_015->SetMarkerColor(kGreen);
grcdhsw_015->SetLineColor(kGreen);
TGraphErrors *grccfr_015 = new TGraphErrors(ccfr_count[0],ccfr_q2[0],ccfr_xf3[0],ccfr_q2err[0],ccfr_err[0]);
grccfr_015->Draw("P*");
grccfr_015->SetMarkerColor(kMagenta);
grccfr_015->SetLineColor(kMagenta);*/


/*TGraphErrors *grnutev_045 = new TGraphErrors(nutev_count[1],nutev_q2[1],nutev_xf3[1],nutev_q2err[1],nutev_q2err[1]);
grnutev_045->Draw("P*");
grnutev_045->SetMarkerColor(kMagenta);
grnutev_045->SetMarkerStyle(20);
grnutev_045->SetMarkerSize(0.9);
grnutev_045->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_045 = new TGraphErrors(cdhsw_count[1],cdhsw_q2[1],cdhsw_xf3[1],cdhsw_q2err[1],cdhsw_q2err[1]);
grcdhsw_045->Draw("P*");
grcdhsw_045->SetMarkerColor(kGreen);
grcdhsw_045->SetLineColor(kGreen);
TGraphErrors *grccfr_045 = new TGraphErrors(ccfr_count[1],ccfr_q2[1],ccfr_xf3[1],ccfr_q2err[1],ccfr_err[1]);
grccfr_045->Draw("P*");
grccfr_045->SetMarkerColor(kMagenta);
grccfr_045->SetLineColor(kMagenta);*/

/*TGraphErrors *grnutev_080 = new TGraphErrors(nutev_count[2],nutev_q2[2],nutev_xf3[2],nutev_q2err[2],nutev_q2err[2]);
grnutev_080->Draw("P*");
grnutev_080->SetMarkerColor(kMagenta);
grnutev_080->SetMarkerStyle(20);
grnutev_080->SetMarkerSize(0.9);
grnutev_080->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_080 = new TGraphErrors(cdhsw_count[2],cdhsw_q2[2],cdhsw_xf3[2],cdhsw_q2err[2],cdhsw_q2err[2]);
grcdhsw_080->Draw("P*");
grcdhsw_080->SetMarkerColor(kGreen);
grcdhsw_080->SetLineColor(kGreen);
TGraphErrors *grccfr_080 = new TGraphErrors(ccfr_count[2],ccfr_q2[2],ccfr_xf3[2],ccfr_q2err[2],ccfr_err[2]);
grccfr_080->Draw("P*");
grccfr_080->SetMarkerColor(kMagenta);
grccfr_080->SetLineColor(kMagenta);*/

/*TGraphErrors *grnutev_125 = new TGraphErrors(nutev_count[3],nutev_q2[3],nutev_xf3[3],nutev_q2err[3],nutev_q2err[3]);
grnutev_125->Draw("P*");
grnutev_125->SetMarkerColor(kMagenta);
grnutev_125->SetMarkerStyle(20);
grnutev_125->SetMarkerSize(0.9);
grnutev_125->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_125 = new TGraphErrors(cdhsw_count[3],cdhsw_q2[3],cdhsw_xf3[3],cdhsw_q2err[3],cdhsw_q2err[3]);
grcdhsw_125->Draw("P*");
grcdhsw_125->SetMarkerColor(kGreen);
grcdhsw_125->SetLineColor(kGreen);
TGraphErrors *grccfr_125 = new TGraphErrors(ccfr_count[3],ccfr_q2[3],ccfr_xf3[3],ccfr_q2err[3],ccfr_q2err[3]);
grccfr_125->Draw("P*");
grccfr_125->SetMarkerColor(kMagenta);
grccfr_125->SetLineColor(kMagenta);*/

/*TGraphErrors *grnutev_175 = new TGraphErrors(nutev_count[4],nutev_q2[4],nutev_xf3[4],nutev_q2err[4],nutev_q2err[4]);
grnutev_175->Draw("P*");
grnutev_175->SetMarkerColor(kMagenta);
grnutev_175->SetMarkerStyle(20);
grnutev_175->SetMarkerSize(0.9);
grnutev_175->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_175 = new TGraphErrors(cdhsw_count[4],cdhsw_q2[4],cdhsw_xf3[4],cdhsw_q2err[4],cdhsw_q2err[4]);
grcdhsw_175->Draw("P*");
grcdhsw_175->SetMarkerColor(kGreen);
grcdhsw_175->SetLineColor(kGreen);
TGraphErrors *grccfr_175 = new TGraphErrors(ccfr_count[4],ccfr_q2[4],ccfr_xf3[4],ccfr_q2err[4],ccfr_q2err[4]);
grccfr_175->Draw("P*");
grccfr_175->SetMarkerColor(kMagenta);
grccfr_175->SetLineColor(kMagenta);*/


/*TGraphErrors *grnutev_225 = new TGraphErrors(nutev_count[5],nutev_q2[5],nutev_xf3[5],nutev_q2err[5],nutev_q2err[5]);
grnutev_225->Draw("P*");
grnutev_225->SetMarkerColor(kMagenta);
grnutev_225->SetMarkerStyle(20);
grnutev_225->SetMarkerSize(0.9);
grnutev_225->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_225 = new TGraphErrors(cdhsw_count[5],cdhsw_q2[5],cdhsw_xf3[5],cdhsw_q2err[5],cdhsw_q2err[5]);
grcdhsw_225->Draw("P*");
grcdhsw_225->SetMarkerColor(kGreen);
grcdhsw_225->SetLineColor(kGreen);
TGraphErrors *grccfr_225 = new TGraphErrors(ccfr_count[5],ccfr_q2[5],ccfr_xf3[5],ccfr_q2err[5],ccfr_q2err[5]);
grccfr_225->Draw("P*");
grccfr_225->SetMarkerColor(kMagenta);
grccfr_225->SetLineColor(kMagenta);*/


/*TGraphErrors *grnutev_275 = new TGraphErrors(nutev_count[6],nutev_q2[6],nutev_xf3[6],nutev_q2err[6],nutev_q2err[6]);
grnutev_275->Draw("P*");
grnutev_275->SetMarkerColor(kMagenta);
grnutev_275->SetMarkerStyle(20);
grnutev_275->SetMarkerSize(0.9);
grnutev_275->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_275 = new TGraphErrors(cdhsw_count[6],cdhsw_q2[6],cdhsw_xf3[6],cdhsw_q2err[6],cdhsw_q2err[6]);
grcdhsw_275->Draw("P*");
grcdhsw_275->SetMarkerColor(kGreen);
grcdhsw_275->SetLineColor(kGreen);
TGraphErrors *grccfr_275 = new TGraphErrors(ccfr_count[6],ccfr_q2[6],ccfr_xf3[6],ccfr_q2err[6],ccfr_q2err[6]);
grccfr_275->Draw("P*");
grccfr_275->SetMarkerColor(kMagenta);
grccfr_275->SetLineColor(kMagenta);*/

/*TGraphErrors *grnutev_350 = new TGraphErrors(nutev_count[7],nutev_q2[7],nutev_xf3[7],nutev_q2err[7],nutev_q2err[7]);
grnutev_350->Draw("P*");
grnutev_350->SetMarkerColor(kMagenta);
grnutev_350->SetMarkerStyle(20);
grnutev_350->SetMarkerSize(0.9);
grnutev_350->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_350 = new TGraphErrors(cdhsw_count[7],cdhsw_q2[7],cdhsw_xf3[7],cdhsw_q2err[7],cdhsw_q2err[7]);
grcdhsw_350->Draw("P*");
grcdhsw_350->SetMarkerColor(kGreen);
grcdhsw_350->SetLineColor(kGreen);
TGraphErrors *grccfr_350 = new TGraphErrors(ccfr_count[7],ccfr_q2[7],ccfr_xf3[7],ccfr_q2err[7],ccfr_q2err[7]);
grccfr_350->Draw("P*");
grccfr_350->SetMarkerColor(kMagenta);
grccfr_350->SetLineColor(kMagenta);*/

/*TGraphErrors *grnutev_450 = new TGraphErrors(nutev_count[8],nutev_q2[8],nutev_xf3[8],nutev_q2err[8],nutev_q2err[8]);
grnutev_450->Draw("P*");
grnutev_450->SetMarkerColor(kMagenta);
grnutev_450->SetMarkerStyle(20);
grnutev_450->SetMarkerSize(0.9);
grnutev_450->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_450 = new TGraphErrors(cdhsw_count[8],cdhsw_q2[8],cdhsw_xf3[8],cdhsw_q2err[8],cdhsw_q2err[8]);
grcdhsw_450->Draw("P*");
grcdhsw_450->SetMarkerColor(kGreen);
grcdhsw_450->SetLineColor(kGreen);
TGraphErrors *grccfr_450 = new TGraphErrors(ccfr_count[8],ccfr_q2[8],ccfr_xf3[8],ccfr_q2err[8],ccfr_q2err[8]);
grccfr_450->Draw("P*");
grccfr_450->SetMarkerColor(kMagenta);
grccfr_450->SetLineColor(kMagenta);*/

/*TGraphErrors *grnutev_550 = new TGraphErrors(nutev_count[9],nutev_q2[9],nutev_xf3[9],nutev_q2err[9],nutev_q2err[9]);
grnutev_550->Draw("P*");
grnutev_550->SetMarkerColor(kMagenta);
grnutev_550->SetMarkerStyle(20);
grnutev_550->SetMarkerSize(0.9);
grnutev_550->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_550 = new TGraphErrors(cdhsw_count[9],cdhsw_q2[9],cdhsw_xf3[9],cdhsw_q2err[9],cdhsw_q2err[9]);
grcdhsw_550->Draw("P*");
grcdhsw_550->SetMarkerColor(kGreen);
grcdhsw_550->SetLineColor(kGreen);
TGraphErrors *grccfr_550 = new TGraphErrors(ccfr_count[9],ccfr_q2[9],ccfr_xf3[9],ccfr_q2err[9],ccfr_q2err[9]);
grccfr_550->Draw("P*");
grccfr_550->SetMarkerColor(kMagenta);
grccfr_550->SetLineColor(kMagenta);*/


/*TGraphErrors *grnutev_650 = new TGraphErrors(nutev_count[10],nutev_q2[10],nutev_xf3[10],nutev_q2err[10],nutev_q2err[10]);
grnutev_650->Draw("P*");
grnutev_650->SetMarkerColor(kMagenta);
grnutev_650->SetMarkerStyle(20);
grnutev_650->SetMarkerSize(0.9);
grnutev_650->SetLineColor(kMagenta);
TGraphErrors *grcdhsw_650 = new TGraphErrors(cdhsw_count[10],cdhsw_q2[10],cdhsw_xf3[10],cdhsw_q2err[10],cdhsw_q2err[10]);
grcdhsw_650->Draw("P*");
grcdhsw_650->SetMarkerColor(kGreen);
grcdhsw_650->SetLineColor(kGreen);
TGraphErrors *grccfr_650 = new TGraphErrors(ccfr_count[10],ccfr_q2[10],ccfr_xf3[10],ccfr_q2err[10],ccfr_q2err[10]);
grccfr_650->Draw("P*");
grccfr_650->SetMarkerColor(kMagenta);
grccfr_650->SetLineColor(kMagenta);*/


/*TGraphErrors *grnutev_750 = new TGraphErrors(nutev_count[11],nutev_q2[11],nutev_xf3[11],nutev_q2err[11],nutev_q2err[11]);
grnutev_750->Draw("P*");
grnutev_750->SetMarkerColor(kMagenta);
grnutev_750->SetMarkerStyle(20);
grnutev_750->SetMarkerSize(0.9);
grnutev_750->SetLineColor(kMagenta);
TGraphErrors *grccfr_750 = new TGraphErrors(ccfr_count[11],ccfr_q2[11],ccfr_xf3[11],ccfr_q2err[11],ccfr_q2err[11]);
grccfr_750->Draw("P*");
grccfr_750->SetMarkerColor(kMagenta);
grccfr_750->SetLineColor(kMagenta);*/
