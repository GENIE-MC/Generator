#include <TCanvas.h>
#include <TF2.h>
#include <TGraph.h>
#include <TH1D.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gst.h"

using namespace std;

namespace {


const double q2_error = 0.05;

void fill_delnu_plot_genie(const gst& gtree, const double q2,  TH1D* delnu_histo) {
	const double nucleon_mass = 0.939;
	if (fabs(q2 - gtree.Q2) < q2_error) {
		const double nu = gtree.Ev - gtree.El;
		const double delnu = nu - gtree.Q2 / (2 * nucleon_mass);
		delnu_histo->Fill(delnu);
	}
}

const char* get_nucleus_name(const int nucZ) {
	switch(nucZ) {
		case 1:
			return "H";
		case 2:
			return "He";
		case 6:
			return "C";
		case 10:
			return "Ne";
		case 13:
			return "Al";
		case 18:
			return "Ar";
		case 26:
			return "Fe";
		case 82:
			return "Pb";
	}
	return "ERROR";
}

TH1D* get_delnu_plot_genie(const int nucZ, const int nucA, const double q2) {
	char histo_title[500];
	sprintf(histo_title, "^{%d}%s, Q^{2} = (%0.1f #pm %0.2f) GeV",
	        nucA, get_nucleus_name(nucZ), q2, q2_error);
	TH1D* delnu = new TH1D("genie", histo_title, 100, -1.0, 2.0);
	delnu->GetXaxis()->SetTitle("#nu - Q^{2} / (2 * M_{p}) (GeV)");
	delnu->GetYaxis()->SetTitle("(1 / #sigma) d#sigma / d#nu (GeV^{-1})");
	delnu->SetLineColor(kRed);
	const int start_run = 500000;
	const int end_run = 500010;
	for (int run = start_run; run <= end_run; run++) {
		char data_file[500];
		sprintf(data_file, "./data/EFF/nu_14/%d_%d/gntp.%d.gst.root", nucZ, nucA, run);
		gst gtree(data_file);
		const int num_events_per_file = 50000;
		for (int event = 0; event < num_events_per_file; event++) {
			gtree.GetEntry(event);
			fill_delnu_plot_genie(gtree, q2, delnu);
		}
	}
	delnu->Scale(1.0 / delnu->Integral("width"));
	return delnu;
}

bool fill_superscaling_vectors(const int nucZ, const int nucA, const double q2,
                               vector<float>* delnu_vals, vector<float>* relative_xs) {
	char superscaling_data_file[500];
	sprintf(superscaling_data_file, "superscaling/%d_%d_%0.1f.txt", nucZ, nucA, q2);
	ifstream infile(superscaling_data_file);
	if (!infile.is_open()) {
		return false;
	}
	string line;
	while (getline(infile, line)) {
		float delnu, xs, garbage;
		if (!(istringstream(line) >> delnu >> garbage >> xs)) {
			cerr << "Couldn't process line: " << line << endl;
			return false;
		}
		delnu_vals->push_back(delnu);
		relative_xs->push_back(xs);
	}
	return true;
}

TH1D* get_delnu_plot_superscaling(const int nucZ, const int nucA, const double q2) {
	vector<float> delnu_vals;
	vector<float> relative_xs;
	if(!fill_superscaling_vectors(nucZ, nucA, q2, &delnu_vals, &relative_xs)) {
		return new TH1D("", "", 100, -1, 1);
	}
	TH1D* delnu = new TH1D("SuSc", "SuSc", delnu_vals.size() - 1, &delnu_vals[0]);
	for (int i = 0; i < (int)delnu_vals.size(); i++) {
		delnu->Fill(delnu_vals[i], relative_xs[i]);
	}
	delnu->SetLineColor(kBlack);
	delnu->Scale(1.0 / delnu->Integral("width"));
	return delnu;
}

void validation_plot_delnu(const int nucZ, const int nucA, const double q2) {
	TCanvas canvas;
	TH1D* genie_delnu = get_delnu_plot_genie(nucZ, nucA, q2);
	TH1D* susc_delnu = get_delnu_plot_superscaling(nucZ, nucA, q2);
	genie_delnu->Draw();
	susc_delnu->Draw("SAME");
	char outfile[500];
	sprintf(outfile, "./plots/delnu_%d_%d_%0.1f.pdf", nucZ, nucA, q2);
	for (int i = 1; i < (int)strlen(outfile) - 4; i++) {
		outfile[i] = outfile[i] == '.' ? 'p' : outfile[i];
	}
	canvas.Print(outfile);
	delete genie_delnu;
	delete susc_delnu;
}

void validation_plot_energy_xs_nu_nubar(bool is_bar) {
	TCanvas canvas;
	TFile def_file("./splines/DEF_QE_splines.root", "READ");
	TFile eff_file("./splines/TE_QE_splines.root", "READ");
	char hist_name[500];
	sprintf(hist_name, "nu_mu_%sC12/qel_cc_%s",
	        (is_bar ? "bar_" : ""), (is_bar ? "p" : "n"));
	TGraph* def_energy_xs = (TGraph*)def_file.Get(hist_name);
	def_energy_xs->SetLineColor(kOrange);
	TGraph* eff_energy_xs = (TGraph*)eff_file.Get(hist_name);
	eff_energy_xs->SetLineColor(kRed);
	def_energy_xs->Apply(&TF2("name,", "y / 6.0"));
	eff_energy_xs->Apply(&TF2("name", "y / 6.0"));
	char title[5000];
	sprintf(title, "%sNeutrino Cross Section - %s + %s #rightarrow %s + #mu^%s",
	        (is_bar ? "Anti-" : ""), (is_bar ? "#bar{#nu}" : "#nu"),
	        (is_bar ? "p" : "n"), (is_bar ? "n" : "p"), (is_bar ? "{+}" : "{-}"));
	TH1F* frame = canvas.DrawFrame(1e-1, 0, 1e2, (is_bar ? 1.2 : 1.6));
	frame->SetTitle(title);
	frame->GetYaxis()->SetTitle("#sigma (10^{-38} cm^{2})");
	frame->GetXaxis()->SetTitle("E_{#nu} (GeV)");
	frame->Draw();
	eff_energy_xs->Draw("SAME");
	def_energy_xs->Draw("SAME");
	canvas.SetLogx();
	char file_name[500];
	sprintf(file_name, "./plots/nu_%sxs_energy.pdf", (is_bar ? "bar_" : ""));
	canvas.Print(file_name);
}

void validation_plot_energy_xs() {
	validation_plot_energy_xs_nu_nubar(true);
	validation_plot_energy_xs_nu_nubar(false);
}

void validation_plot_energy_ratio_nu_nubar(bool is_bar) {
	TCanvas canvas;
	TFile def_file("./splines/DEF_QE_splines.root", "READ");
	TFile eff_file("./splines/TE_QE_splines.root", "READ");
	char hist_name[500];
	sprintf(hist_name, "nu_mu_%sC12/qel_cc_%s",
	        (is_bar ? "bar_" : ""), (is_bar ? "p" : "n"));
	TGraph* def_energy_xs = (TGraph*)def_file.Get(hist_name);
	TGraph* eff_energy_xs = (TGraph*)eff_file.Get(hist_name);
	TGraph ratio;
	ratio.Set(def_energy_xs->GetN());
	for (int i = 0; i < def_energy_xs->GetN(); i++) {
		double x_val, def_point, eff_point;
		def_energy_xs->GetPoint(i, x_val, def_point);
		eff_energy_xs->GetPoint(i, x_val, eff_point);
		if (def_point > 0) {
			ratio.SetPoint(i, x_val,  eff_point / def_point);
		} else {
			ratio.SetPoint(i, x_val, 0);
		}
	}
	ratio.SetLineColor(kRed);
	ratio.Draw();
	char title[5000];
	sprintf(title, "%s + %s #rightarrow %s + #mu^%s, (Trans Enh in G^{#nu}_{M}, M_{A} = 1.014) / (M_{A} = 1.014)",
	        (is_bar ? "#bar{#nu}" : "#nu"), (is_bar ? "p" : "n"), (is_bar ? "n" : "p"), (is_bar ? "{+}" : "{-}"));
	TH1F* frame = canvas.DrawFrame(0, (is_bar ? 0.85 : 0.9), 10, (is_bar ? 1.4 : 1.5));
	frame->SetTitle(title);
	frame->GetYaxis()->SetTitle("ratio #sigma");
	frame->GetXaxis()->SetTitle("E^{#nu} (GeV)");
	frame->Draw();
	ratio.Draw("SAME");
	char file_name[500];
	sprintf(file_name, "./plots/nu_%sratio_energy.pdf", (is_bar ? "bar_" : ""));
	canvas.Print(file_name);
}

void validation_plot_energy_ratio() {
	validation_plot_energy_ratio_nu_nubar(true);
	validation_plot_energy_ratio_nu_nubar(false);
}

double get_qe_increase_factor(double is_bar, double energy = 3.0) {
	TFile def_file("./splines/DEF_QE_splines.root", "READ");
	TFile eff_file("./splines/TE_QE_splines.root", "READ");
	char hist_name[500];
	sprintf(hist_name, "nu_mu_%sC12/qel_cc_%s",
	        (is_bar ? "bar_" : ""), (is_bar ? "p" : "n"));
	TGraph* def_energy_xs = (TGraph*)def_file.Get(hist_name);
	TGraph* eff_energy_xs = (TGraph*)eff_file.Get(hist_name);
	for (int i = 0; i < def_energy_xs->GetN(); i++) {
		double point_energy, def_xs, eff_xs;
		def_energy_xs->GetPoint(i, point_energy, def_xs);
		if (point_energy > energy) {
			eff_energy_xs->GetPoint(i, point_energy, eff_xs);
			return eff_xs / def_xs;
		}
	}
	return -1;
}

void populate_q2_plot(bool is_bar, const char* tag, TH1D* q2) {
	const int start_run = 500000;
	const int end_run = 500010;
	for (int run = start_run; run <= end_run; run++) {
		char data_file[500];
		sprintf(data_file, "./data/%s/nu_%d/6_12/gntp.%d.gst.root", tag,
		        (is_bar ? -14 : 14), run);
		gst gtree(data_file);
		const int num_events_per_file = 50000;
		for (int event = 0; event < num_events_per_file; event++) {
			gtree.GetEntry(event);
			q2->Fill(gtree.Q2);
		}
	}
}


void validation_plot_q2_ratio_nu_nubar(bool is_bar) {
	TCanvas canvas;
	const double increase_factor = get_qe_increase_factor(is_bar);
	char title[500];
	sprintf(title, "E_{#nu} = 3 GeV, %s + %s #rightarrow %s + #mu^%s, (Trans Enh in G^{#nu}_{M}, M_{A} = 1.014) / (M_{A} = 1.014)",
	        (is_bar ? "#bar{#nu}" : "#nu"), (is_bar ? "p" : "n"), (is_bar ? "n" : "p"), (is_bar ? "{+}" : "{-}"));
	TH1D te_q2("", title, 50, 0, 2.5);
	te_q2.GetYaxis()->SetTitle("ratio d#sigma / dQ^{2}");
	te_q2.GetXaxis()->SetTitle("Q^{2} (GeV / c)^{2}");
	TH1D def_q2("", "",50, 0, 2.5);
	populate_q2_plot(is_bar, "TE", &te_q2);
	populate_q2_plot(is_bar, "DEF", &def_q2);
	te_q2.Scale(increase_factor);
	te_q2.Divide(&def_q2);
	te_q2.SetLineColor(kRed);
	te_q2.Draw();
	char filename[500];
	sprintf(filename, "./plots/nu%s_ratio_q2.pdf", (is_bar ? "_bar" : ""));
	canvas.Print(filename);
}


void validation_plot_q2_ratio() {
	validation_plot_q2_ratio_nu_nubar(true);
	validation_plot_q2_ratio_nu_nubar(false);
}

}  // namespace

void Validation() {
	const int nucZs[] = {1, 2, 6,  10, 13, 18, 26, 82 };
	const int nucAs[] = {2, 4, 12, 20, 27, 40, 56, 208};
	const int num_nucs = 8;

	const double q2s[] = {0.1, 0.3, 0.5, 0.7, 1.0, 1.2, 1.5, 2.0};
	const int num_q2s = 8;
	validation_plot_energy_xs();
	validation_plot_energy_ratio();
	validation_plot_q2_ratio();
	
	for (int i = 0; i < num_nucs; i++) {
		for (int j = 0; j < num_q2s; j++) {
			validation_plot_delnu(nucZs[i], nucAs[i], q2s[j]);
		}
	}
}

int main() {
	Validation();
	return 0;
}
