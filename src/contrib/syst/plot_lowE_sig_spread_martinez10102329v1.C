{

double xlo[12] = {
0.18142,
0.21958,
0.25929,
0.30509,
0.35853,
0.39975,
0.46386,
0.51727,
0.57220,
0.63170,
0.68813,
0.74914
};

double ylo[12] = {
0.16782,
0.63281,
1.30182,
1.99975,
2.78485,
3.39549,
4.18018,
4.70288,
5.19638,
5.60223,
5.92074,
6.20991
};

double xhi[18] = {
0.17535,
0.20894,
0.23644,
0.26394,
0.29449,
0.32961,
0.36474,
0.41209,
0.44112,
0.47626,
0.51139,
0.54651,
0.58315,
0.61672,
0.65487,
0.68996,
0.71894,
0.75097
};

double yhi[18] = {
0.40129,
0.92477,
1.50679,
2.08881,
2.61241,
3.22329,
3.80502,
4.70697,
5.37640,
6.13306,
6.71478,
7.29651,
7.76156,
8.11011,
8.48764,
8.77782,
8.95162,
9.12531
};


double xlo_no2p2h[17] = {
0.15701,
0.18445,
0.20884,
0.24238,
0.28049,
0.31707,
0.35213,
0.39482,
0.42835,
0.46646,
0.51677,
0.56555,
0.60823,
0.65091,
0.69207,
0.72409,
0.74695,
};
double ylo_no2p2h[17] = {
0.02907,
0.26163,
0.55233,
1.04651,
1.62791,
2.20930,
2.70349,
3.28488,
3.77907,
4.15698,
4.70930,
5.17442,
5.49419,
5.72674,
5.95930,
6.10465,
6.25000,
};
double xhi_no2p2h[20] = {
0.14024,
0.17683,
0.19817,
0.21646,
0.24238,
0.26677,
0.29421,
0.32622,
0.35823,
0.38720,
0.42073,
0.45427,
0.49238,
0.52287,
0.55945,
0.60671,
0.64329,
0.68445,
0.71646,
0.74848,
};
double yhi_no2p2h[20] = {
0.08721,
0.46512,
0.75581,
1.19186,
1.65698,
2.09302,
2.67442,
3.28488,
3.72093,
4.21512,
4.65116,
5.11628,
5.52326,
5.78488,
6.13372,
6.45349,
6.68605,
6.88953,
7.03488,
7.12209,
};

TSpline3 lo    ("lo",xlo,ylo,12,"0");
TSpline3 hi    ("hi",xhi,yhi,18,"0");
TSpline3 lo_no2("lo",xlo_no2p2h,ylo_no2p2h,17,"0");
TSpline3 hi_no2("hi",xhi_no2p2h,yhi_no2p2h,20,"0");

const int    nsteps = 100;
const double emin   = 0.20;
const double emax   = 0.70;
const double step   = (emax-emin)/(nsteps-1);

double errp    [nsteps];
double errn    [nsteps];
double errp_no2[nsteps];
double errn_no2[nsteps];
double enu     [nsteps];

for(int i=0; i<nsteps; i++) {

	double x      = emin + i*step;
	double yh     = hi->Eval(x);
	double yl     = lo->Eval(x);
	double yc     = 0.5*(yh+yl);
	double yh_no2 = hi_no2->Eval(x);
	double yl_no2 = lo_no2->Eval(x);
	double yc_no2 = 0.5*(yh_no2+yl_no2);

	errp    [i] = 100*(yh-yc)/yc;
	errn    [i] = 100*(yl-yc)/yc;
	errp_no2[i] = 100*(yh_no2-yc_no2)/yc_no2;
	errn_no2[i] = 100*(yl_no2-yc_no2)/yc_no2;
	enu [i] = x;
}

TGraph * grerrp     = new TGraph(nsteps, enu, errp);
TGraph * grerrn     = new TGraph(nsteps, enu, errn);
TGraph * grerrp_no2 = new TGraph(nsteps, enu, errp_no2);
TGraph * grerrn_no2 = new TGraph(nsteps, enu, errn_no2);

grerrp->SetLineStyle(kDashed);
grerrn->SetLineStyle(kDashed);

TCanvas * c = new TCanvas();
TH1F * hframe = (TH1F*)c->DrawFrame(0.2,-50,0.7,50);
hframe->GetXaxis()->SetTitle("E_{#nu} (GeV)");
hframe->GetYaxis()->SetTitle("#sigma spread %");
hframe->Draw();

grerrp->Draw("L");
grerrn->Draw("L");
grerrp_no2->Draw("L");
grerrn_no2->Draw("L");

TLegend * legend = new TLegend(0.6, 0.6, 0.9, 0.9);
legend->SetFillColor(0);
legend->SetBorderSize(0);
legend->AddEntry( grerrp,     "Model spread (from Fig.1 in Martini arXiv:1010.2329v1", "L");
legend->AddEntry( grerrp_no2, "Model spread (as above but excluding RPA-2p2h)", "L");
legend->Draw();

}
