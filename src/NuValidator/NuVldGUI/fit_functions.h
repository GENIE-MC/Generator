// Simple fit functions
//

#ifndef _FITFUNC_H_
#define _FITFUNC_H_

//--- using the new neugen-wrapper

double xsec_fitfunc             (double * x, double * par); // incl. neutrino-xsec
double xsec_fitfunc_e           (double * x, double * par); // incl. neutrino-xsec/E
double exclusive_xsec_fitfunc   (double * x, double * par); // excl. neutrino-xsec
double exclusive_xsec_fitfunc_e (double * x, double * par); // excl. neutrino-xsec/E

//--- util method: fit-params -> neugen-config

void update_neugen_config_with_fit_params(double * par);
void update_neugen_config(double par, int iparam);

void fcn_mgr_floating_norm(
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
void fcn_mgr_floating_norm_e(
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

#endif
