//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenConstants

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_CONSTANTS_H_
#define _NEUGEN_CONSTANTS_H_

namespace genie {
namespace nuvld {

const char * k_neugen_plot_type[] = {
   "total xsec", "diff. xsec", "struc. func", 0
};

const char * k_scaling_flux[] = {
   "none", "ANL", "GGM", "BNL", "BEBC", 0
};

const char * k_neugen_plot_variable[] = {
   "none", "|q^2|", "W", "x", "y", 0
};

const char * k_neugen_sf_plot_variable[] = {
   "|q^2|", "x" , 0
};

const char * k_neugen_plot_range_option[] = {
   "automatic", "custom", 0
};

const char * k_neugen_sf[] = { "F2", "xF3", 0 };

const char * k_neugen_nu[] = {
   "nu_e", "nu_e_bar", "nu_mu", "nu_mu_bar", "nu_tau", "nu_tau_bar", 0
};
const char * k_neugen_wcurr[] = { "CC", "NC", "CC+NC", 0 };

const char * k_neugen_initstate[] = {
   "nu + p", "nu + n", "nu_bar + p", "nu_bar + n", "nu + N (isoscalar)", "nu_bar + N (isoscalar)",
   "l + p", "l + n", "l + N (isoscalar)", 0
};

const char * k_neugen_finstate[] = {
   /*QEL*/   
   "l(-) + p ", "l(+) + n ", "nu + p ", "nu + n ", "nu_bar + p ", "nu_bar + n ",
   
   /*RES - nu CC*/
   "l(-) + p + pi(+)", "l(-) + p + pi(0) ", "l(-) + n + pi(+)",              

   /*RES - nu_bar CC*/   
   "l(+) + n + pi(-)", "l(+) + n + pi(0)", 

   /*RES - nu NC*/
   "nu + p + pi(0)", "nu + n + pi(+)", "nu + n + pi(0)", "nu + p + pi(-)",  

   /*RES - nu_bar NC*/
   "nu_bar + p + pi(0)", "nu_bar + n + pi(+)", "nu_bar + n + pi(0)", "nu_bar + p + pi(-)",  

    /*DIS*/
   "l(-) + X", "l(+) + X", "nu + X", "nu_bar + X", 0                       
};

const char * k_neugen_cut_variable[] = {
   "none", "|q^2|", "W", "x", "y", 0
};

} // nuvld namespace
} // genie namespace

#endif


