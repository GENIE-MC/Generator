{
TFile favg    ("./ccqe_avg_defkF.root",     "read");
TFile fnoavg  ("./ccqe_noavg_defkF.root",   "read");
TFile favgkFlo("./ccqe_avg_kF200.root",     "read");
TFile favgkFhi("./ccqe_avg_kF250.root",     "read");
TFile favgMalo("./ccqe_avg_Mam15pct.root",  "read");
TFile favgMahi("./ccqe_avg_Map15pct.root",  "read");

TDirectory * davg     = (TDirectory *) favg     -> Get("nu_mu_O16");
TDirectory * dnoavg   = (TDirectory *) fnoavg   -> Get("nu_mu_O16");
TDirectory * davgkFlo = (TDirectory *) favgkFlo -> Get("nu_mu_O16");
TDirectory * davgkFhi = (TDirectory *) favgkFhi -> Get("nu_mu_O16");
TDirectory * davgMalo = (TDirectory *) favgMalo -> Get("nu_mu_O16");
TDirectory * davgMahi = (TDirectory *) favgMahi -> Get("nu_mu_O16");

TGraph * xsavg     = (TGraph*) davg     -> Get("qel_cc_n");
TGraph * xsnoavg   = (TGraph*) dnoavg   -> Get("qel_cc_n");
TGraph * xsavgkFlo = (TGraph*) davgkFlo -> Get("qel_cc_n");
TGraph * xsavgkFhi = (TGraph*) davgkFhi -> Get("qel_cc_n");
TGraph * xsavgMalo = (TGraph*) davgMalo -> Get("qel_cc_n");
TGraph * xsavgMahi = (TGraph*) davgMahi -> Get("qel_cc_n");

TCanvas * c1 = new TCanvas();
TH1F * hframe1 = (TH1F*)c1->DrawFrame(0,0,1.8,10);

xsavg->SetLineWidth(3);
xsavgkFlo->SetLineColor(kGreen);
xsavgkFhi->SetLineColor(kRed);

xsavg->Draw("L");  
xsnoavg->Draw("LSAME");  
xsavgkFlo->Draw("LSAME");  
xsavgkFhi->Draw("LSAME");  

const int    nsteps = 100;
const double emin   = 0.20;
const double emax   = 0.70;
const double step   = (emax-emin)/(nsteps-1);

double errp   [nsteps];
double errn   [nsteps];
double errp_pF[nsteps];
double errn_pF[nsteps];
double errp_kF[nsteps];
double errn_kF[nsteps];
double errp_Ma[nsteps];
double errn_Ma[nsteps];
double enu    [nsteps];

for(int i=0; i<nsteps; i++) {

        double x          = emin + i*step;
        double yavg       = xsavg->Eval(x);
        double ynoavg     = xsnoavg->Eval(x);
        double yavgnoavgc = 0.5*(yavg+ynoavg);
        double yavgkFlo   = xsavgkFlo->Eval(x);
        double yavgkFhi   = xsavgkFhi->Eval(x);
        double yavgMalo   = xsavgMalo->Eval(x);
        double yavgMahi   = xsavgMahi->Eval(x);

        enu [i] = x;
        
        errp_pF[i] =  100*( TMath::Max(yavg,ynoavg) - yavgnoavgc)/yavgnoavgc;
        errn_pF[i] =  100*( TMath::Min(yavg,ynoavg) - yavgnoavgc)/yavgnoavgc;
        errp_kF[i] =  100*( TMath::Max(yavgkFhi,yavgkFlo) - yavg)/yavg;
        errn_kF[i] =  100*( TMath::Min(yavgkFhi,yavgkFlo) - yavg)/yavg;
        errp_Ma[i] =  100*( TMath::Max(yavgMahi,yavgMalo) - yavg)/yavg;
        errn_Ma[i] =  100*( TMath::Min(yavgMahi,yavgMalo) - yavg)/yavg;
        errp   [i] =  TMath::Sqrt(errp_pF[i]*errp_pF[i] + errp_kF[i]*errp_kF[i] + errp_Ma[i]*errp_Ma[i]);
        errn   [i] = -TMath::Sqrt(errn_pF[i]*errn_pF[i] + errn_kF[i]*errn_kF[i] + errn_Ma[i]*errn_Ma[i]);
}

TGraph * grerrp_pF = new TGraph(nsteps, enu, errp_pF);
TGraph * grerrn_pF = new TGraph(nsteps, enu, errn_pF);
TGraph * grerrp_kF = new TGraph(nsteps, enu, errp_kF);
TGraph * grerrn_kF = new TGraph(nsteps, enu, errn_kF);
TGraph * grerrp_Ma = new TGraph(nsteps, enu, errp_Ma);
TGraph * grerrn_Ma = new TGraph(nsteps, enu, errn_Ma);
TGraph * grerrp    = new TGraph(nsteps, enu, errp);
TGraph * grerrn    = new TGraph(nsteps, enu, errn);

TCanvas * c2 = new TCanvas();
TH1F * hframe2 = (TH1F*)c2->DrawFrame(0.2,-50,0.7,50);
hframe2->GetXaxis()->SetTitle("E_{#nu} (GeV)");
hframe2->GetYaxis()->SetTitle("#sigma spread %");
hframe2->Draw();

grerrp   ->SetLineWidth(3);
grerrn   ->SetLineWidth(3);
grerrp_pF->SetLineStyle(kDotted);
grerrn_pF->SetLineStyle(kDotted);
grerrp_kF->SetLineStyle(kDashed);
grerrn_kF->SetLineStyle(kDashed);

grerrp   ->Draw("L");
grerrn   ->Draw("L");
grerrp_pF->Draw("L");
grerrn_pF->Draw("L");
grerrp_kF->Draw("L");
grerrn_kF->Draw("L");
grerrp_Ma->Draw("L");
grerrn_Ma->Draw("L");

TLegend * legend = new TLegend(0.6, 0.6, 0.9, 0.9);
legend->SetFillColor(0);
legend->SetBorderSize(0);
legend->AddEntry( grerrp, "Combined error for #nu_{#mu} O^{16} CCQE", "L");
legend->AddEntry( grerrp_pF, "Effect of E_{#nu} smearing from the nucleon momentum distribution", "L");
legend->AddEntry( grerrp_kF, "k_{F} #pm 25 MeV (~10%) (effect of modifying Pauli-blocking term)", "L");
legend->AddEntry( grerrp_Ma, "M_{A} #pm 15%", "L");
legend->Draw();
}
