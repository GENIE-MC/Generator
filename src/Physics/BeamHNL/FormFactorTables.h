//----------------------------------------------------------------------------
/*!

  Lookup tables provided for the functions `f_{1}(x)` and `f_{2}(x)`
  that enter the N --> \nu \ell_{\beta} \ell_{\beta} channel

\namespace  genie::HNL::HNLSelector

\brief      Form factor lookup tables

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    January 25th, 2022

\cpright    ??? - TBD

*/
//----------------------------------------------------------------------------

#ifndef _HNL_JFORMFACTORTABLES_H_
#define _HNL_JFORMFACTORTABLES_H_

namespace genie {
namespace HNL {

    namespace HNLSelector {

      // lookup tables, 50x10 + 1
      static const double PARTWIDTH = 1e-3;
      
      static const double FormfactorX[] = { 0.000, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 
					0.020, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 
					0.040, 0.041, 0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.050, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 0.058, 0.059, 
					0.060, 0.061, 0.062, 0.063, 0.064, 0.065, 0.066, 0.067, 0.068, 0.069, 0.070, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 
					0.080, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.090, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099,
					0.100, 0.101, 0.102, 0.103, 0.104, 0.105, 0.106, 0.107, 0.108, 0.109, 0.110, 0.111, 0.112, 0.113, 0.114, 0.115, 0.116, 0.117, 0.118, 0.119,
					0.120, 0.121, 0.122, 0.123, 0.124, 0.125, 0.126, 0.127, 0.128, 0.129, 0.130, 0.131, 0.132, 0.133, 0.134, 0.135, 0.136, 0.137, 0.138, 0.139,
					0.140, 0.141, 0.142, 0.143, 0.144, 0.145, 0.146, 0.147, 0.148, 0.149, 0.150, 0.151, 0.152, 0.153, 0.154, 0.155, 0.156, 0.157, 0.158, 0.159,
					0.160, 0.161, 0.162, 0.163, 0.164, 0.165, 0.166, 0.167, 0.168, 0.169, 0.170, 0.171, 0.172, 0.173, 0.174, 0.175, 0.176, 0.177, 0.178, 0.179,
					0.180, 0.181, 0.182, 0.183, 0.184, 0.185, 0.186, 0.187, 0.188, 0.189, 0.190, 0.191, 0.192, 0.193, 0.194, 0.195, 0.196, 0.197, 0.198, 0.199,
					0.200, 0.201, 0.202, 0.203, 0.204, 0.205, 0.206, 0.207, 0.208, 0.209, 0.210, 0.211, 0.212, 0.213, 0.214, 0.215, 0.216, 0.217, 0.218, 0.219,
					0.220, 0.221, 0.222, 0.223, 0.224, 0.225, 0.226, 0.227, 0.228, 0.229, 0.230, 0.231, 0.232, 0.233, 0.234, 0.235, 0.236, 0.237, 0.238, 0.239,
					0.240, 0.241, 0.242, 0.243, 0.244, 0.245, 0.246, 0.247, 0.248, 0.249, 0.250, 0.251, 0.252, 0.253, 0.254, 0.255, 0.256, 0.257, 0.258, 0.259,
					0.260, 0.261, 0.262, 0.263, 0.264, 0.265, 0.266, 0.267, 0.268, 0.269, 0.270, 0.271, 0.272, 0.273, 0.274, 0.275, 0.276, 0.277, 0.278, 0.279,
					0.280, 0.281, 0.282, 0.283, 0.284, 0.285, 0.286, 0.287, 0.288, 0.289, 0.290, 0.291, 0.292, 0.293, 0.294, 0.295, 0.296, 0.297, 0.298, 0.299,
					0.300, 0.301, 0.302, 0.303, 0.304, 0.305, 0.306, 0.307, 0.308, 0.309, 0.310, 0.311, 0.312, 0.313, 0.314, 0.315, 0.316, 0.317, 0.318, 0.319,
					0.320, 0.321, 0.322, 0.323, 0.324, 0.325, 0.326, 0.327, 0.328, 0.329, 0.330, 0.331, 0.332, 0.333, 0.334, 0.335, 0.336, 0.337, 0.338, 0.339,
					0.340, 0.341, 0.342, 0.343, 0.344, 0.345, 0.346, 0.347, 0.348, 0.349, 0.350, 0.351, 0.352, 0.353, 0.354, 0.355, 0.356, 0.357, 0.358, 0.359,
					0.360, 0.361, 0.362, 0.363, 0.364, 0.365, 0.366, 0.367, 0.368, 0.369, 0.370, 0.371, 0.372, 0.373, 0.374, 0.375, 0.376, 0.377, 0.378, 0.379,
					0.380, 0.381, 0.382, 0.383, 0.384, 0.385, 0.386, 0.387, 0.388, 0.389, 0.390, 0.391, 0.392, 0.393, 0.394, 0.395, 0.396, 0.397, 0.398, 0.399,
					0.400, 0.401, 0.402, 0.403, 0.404, 0.405, 0.406, 0.407, 0.408, 0.409, 0.410, 0.411, 0.412, 0.413, 0.414, 0.415, 0.416, 0.417, 0.418, 0.419,
					0.420, 0.421, 0.422, 0.423, 0.424, 0.425, 0.426, 0.427, 0.428, 0.429, 0.430, 0.431, 0.432, 0.433, 0.434, 0.435, 0.436, 0.437, 0.438, 0.439,
					0.440, 0.441, 0.442, 0.443, 0.444, 0.445, 0.446, 0.447, 0.448, 0.449, 0.450, 0.451, 0.452, 0.453, 0.454, 0.455, 0.456, 0.457, 0.458, 0.459,
					0.460, 0.461, 0.462, 0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.470, 0.471, 0.472, 0.473, 0.474, 0.475, 0.476, 0.477, 0.478, 0.479,
					0.480, 0.481, 0.482, 0.483, 0.484, 0.485, 0.486, 0.487, 0.488, 0.489, 0.490, 0.491, 0.492, 0.493, 0.494, 0.495, 0.496, 0.497, 0.498, 0.499, 0.500 }; //where F1, F2 are evaluated

      static const double FormfactorF1[] = { 1.000e+00, 1.000e+00, 9.999e-01, 9.999e-01, 9.997e-01, 9.996e-01, 9.994e-01, 9.992e-01, 9.990e-01, 9.987e-01,
					 9.984e-01, 9.981e-01, 9.977e-01, 9.973e-01, 9.969e-01, 9.964e-01, 9.959e-01, 9.954e-01, 9.948e-01, 9.943e-01,
					 9.936e-01, 9.930e-01, 9.923e-01, 9.916e-01, 9.909e-01, 9.901e-01, 9.893e-01, 9.884e-01, 9.876e-01, 9.867e-01,
					 9.858e-01, 9.848e-01, 9.838e-01, 9.828e-01, 9.818e-01, 9.807e-01, 9.796e-01, 9.784e-01, 9.773e-01, 9.761e-01,
					 9.749e-01, 9.736e-01, 9.723e-01, 9.710e-01, 9.697e-01, 9.683e-01, 9.669e-01, 9.655e-01, 9.640e-01, 9.626e-01,
					 9.610e-01, 9.595e-01, 9.579e-01, 9.564e-01, 9.547e-01, 9.531e-01, 9.514e-01, 9.497e-01, 9.480e-01, 9.462e-01,
					 9.445e-01, 9.427e-01, 9.408e-01, 9.390e-01, 9.371e-01, 9.352e-01, 9.332e-01, 9.313e-01, 9.293e-01, 9.273e-01,
					 9.252e-01, 9.232e-01, 9.211e-01, 9.190e-01, 9.168e-01, 9.147e-01, 9.125e-01, 9.103e-01, 9.081e-01, 9.058e-01,
					 9.035e-01, 9.012e-01, 8.989e-01, 8.966e-01, 8.942e-01, 8.918e-01, 8.894e-01, 8.870e-01, 8.845e-01, 8.820e-01,
					 8.795e-01, 8.770e-01, 8.745e-01, 8.719e-01, 8.693e-01, 8.667e-01, 8.641e-01, 8.615e-01, 8.588e-01, 8.561e-01,
					 8.534e-01, 8.507e-01, 8.480e-01, 8.452e-01, 8.424e-01, 8.396e-01, 8.368e-01, 8.340e-01, 8.311e-01, 8.283e-01,
					 8.254e-01, 8.225e-01, 8.195e-01, 8.166e-01, 8.136e-01, 8.107e-01, 8.077e-01, 8.047e-01, 8.017e-01, 7.986e-01,
					 7.956e-01, 7.925e-01, 7.894e-01, 7.863e-01, 7.832e-01, 7.801e-01, 7.770e-01, 7.738e-01, 7.706e-01, 7.675e-01,
					 7.643e-01, 7.611e-01, 7.578e-01, 7.546e-01, 7.514e-01, 7.481e-01, 7.448e-01, 7.415e-01, 7.382e-01, 7.349e-01,
					 7.316e-01, 7.283e-01, 7.250e-01, 7.216e-01, 7.182e-01, 7.149e-01, 7.115e-01, 7.081e-01, 7.047e-01, 7.013e-01,
					 6.979e-01, 6.944e-01, 6.910e-01, 6.876e-01, 6.841e-01, 6.806e-01, 6.772e-01, 6.737e-01, 6.702e-01, 6.667e-01,
					 6.632e-01, 6.597e-01, 6.562e-01, 6.527e-01, 6.491e-01, 6.456e-01, 6.421e-01, 6.385e-01, 6.350e-01, 6.314e-01,
					 6.279e-01, 6.243e-01, 6.207e-01, 6.171e-01, 6.136e-01, 6.100e-01, 6.064e-01, 6.028e-01, 5.992e-01, 5.956e-01,
					 5.920e-01, 5.884e-01, 5.848e-01, 5.812e-01, 5.776e-01, 5.740e-01, 5.704e-01, 5.668e-01, 5.632e-01, 5.595e-01,
					 5.559e-01, 5.523e-01, 5.487e-01, 5.451e-01, 5.415e-01, 5.378e-01, 5.342e-01, 5.306e-01, 5.270e-01, 5.234e-01,
					 5.198e-01, 5.162e-01, 5.125e-01, 5.089e-01, 5.053e-01, 5.017e-01, 4.981e-01, 4.945e-01, 4.909e-01, 4.873e-01,
					 4.837e-01, 4.802e-01, 4.766e-01, 4.730e-01, 4.694e-01, 4.659e-01, 4.623e-01, 4.587e-01, 4.552e-01, 4.516e-01,
					 4.481e-01, 4.445e-01, 4.410e-01, 4.375e-01, 4.339e-01, 4.304e-01, 4.269e-01, 4.234e-01, 4.199e-01, 4.164e-01,
					 4.129e-01, 4.094e-01, 4.060e-01, 4.025e-01, 3.990e-01, 3.956e-01, 3.922e-01, 3.887e-01, 3.853e-01, 3.819e-01,
					 3.785e-01, 3.751e-01, 3.717e-01, 3.683e-01, 3.649e-01, 3.616e-01, 3.582e-01, 3.549e-01, 3.515e-01, 3.482e-01,
					 3.449e-01, 3.416e-01, 3.383e-01, 3.350e-01, 3.318e-01, 3.285e-01, 3.253e-01, 3.220e-01, 3.188e-01, 3.156e-01,
					 3.124e-01, 3.092e-01, 3.060e-01, 3.029e-01, 2.997e-01, 2.966e-01, 2.935e-01, 2.904e-01, 2.873e-01, 2.842e-01,
					 2.811e-01, 2.780e-01, 2.750e-01, 2.720e-01, 2.689e-01, 2.659e-01, 2.630e-01, 2.600e-01, 2.570e-01, 2.541e-01,
					 2.511e-01, 2.482e-01, 2.453e-01, 2.424e-01, 2.396e-01, 2.367e-01, 2.339e-01, 2.310e-01, 2.282e-01, 2.254e-01,
					 2.226e-01, 2.199e-01, 2.171e-01, 2.144e-01, 2.117e-01, 2.090e-01, 2.063e-01, 2.036e-01, 2.010e-01, 1.983e-01,
					 1.957e-01, 1.931e-01, 1.905e-01, 1.880e-01, 1.854e-01, 1.829e-01, 1.804e-01, 1.779e-01, 1.754e-01, 1.729e-01,
					 1.705e-01, 1.681e-01, 1.656e-01, 1.632e-01, 1.609e-01, 1.585e-01, 1.562e-01, 1.539e-01, 1.516e-01, 1.493e-01,
					 1.470e-01, 1.447e-01, 1.425e-01, 1.403e-01, 1.381e-01, 1.359e-01, 1.338e-01, 1.316e-01, 1.295e-01, 1.274e-01,
					 1.253e-01, 1.233e-01, 1.212e-01, 1.192e-01, 1.172e-01, 1.152e-01, 1.132e-01, 1.113e-01, 1.093e-01, 1.074e-01,
					 1.055e-01, 1.037e-01, 1.018e-01, 9.996e-02, 9.814e-02, 9.634e-02, 9.455e-02, 9.279e-02, 9.105e-02, 8.932e-02,
					 8.761e-02, 8.593e-02, 8.426e-02, 8.261e-02, 8.098e-02, 7.936e-02, 7.777e-02, 7.620e-02, 7.464e-02, 7.310e-02,
					 7.159e-02, 7.009e-02, 6.861e-02, 6.714e-02, 6.570e-02, 6.428e-02, 6.287e-02, 6.148e-02, 6.011e-02, 5.876e-02,
					 5.743e-02, 5.612e-02, 5.482e-02, 5.354e-02, 5.229e-02, 5.104e-02, 4.982e-02, 4.862e-02, 4.743e-02, 4.626e-02,
					 4.511e-02, 4.398e-02, 4.286e-02, 4.177e-02, 4.069e-02, 3.962e-02, 3.858e-02, 3.755e-02, 3.654e-02, 3.555e-02,
					 3.457e-02, 3.361e-02, 3.267e-02, 3.174e-02, 3.083e-02, 2.994e-02, 2.906e-02, 2.820e-02, 2.736e-02, 2.653e-02,
					 2.572e-02, 2.493e-02, 2.415e-02, 2.339e-02, 2.264e-02, 2.191e-02, 2.119e-02, 2.049e-02, 1.980e-02, 1.913e-02,
					 1.847e-02, 1.783e-02, 1.720e-02, 1.659e-02, 1.599e-02, 1.541e-02, 1.484e-02, 1.428e-02, 1.374e-02, 1.321e-02,
					 1.270e-02, 1.220e-02, 1.171e-02, 1.123e-02, 1.077e-02, 1.032e-02, 9.884e-03, 9.459e-03, 9.047e-03, 8.647e-03,
					 8.258e-03, 7.882e-03, 7.517e-03, 7.163e-03, 6.821e-03, 6.490e-03, 6.170e-03, 5.860e-03, 5.562e-03, 5.273e-03,
					 4.995e-03, 4.727e-03, 4.469e-03, 4.220e-03, 3.981e-03, 3.752e-03, 3.531e-03, 3.320e-03, 3.117e-03, 2.923e-03,
					 2.737e-03, 2.559e-03, 2.390e-03, 2.228e-03, 2.074e-03, 1.927e-03, 1.788e-03, 1.656e-03, 1.530e-03, 1.412e-03,
					 1.300e-03, 1.194e-03, 1.094e-03, 1.000e-03, 9.118e-04, 8.292e-04, 7.519e-04, 6.797e-04, 6.125e-04, 5.500e-04,
					 4.922e-04, 4.387e-04, 3.894e-04, 3.440e-04, 3.026e-04, 2.647e-04, 2.303e-04, 1.991e-04, 1.710e-04, 1.458e-04,
					 1.234e-04, 1.035e-04, 8.594e-05, 7.061e-05, 5.731e-05, 4.589e-05, 3.617e-05, 2.801e-05, 2.124e-05, 1.572e-05,
					 1.130e-05, 7.842e-06, 5.211e-06, 3.277e-06, 1.917e-06, 1.016e-06, 4.671e-07, 1.713e-07, 4.158e-08, 3.688e-09, 0.000e+00 };

      static const double FormfactorF2[] = { 0.000e+00, 8.000e-06, 3.199e-05, 7.196e-05, 1.279e-04, 1.997e-04, 2.874e-04, 3.909e-04, 5.102e-04, 6.452e-04,
					 7.958e-04, 9.620e-04, 1.144e-03, 1.341e-03, 1.553e-03, 1.781e-03, 2.024e-03, 2.281e-03, 2.554e-03, 2.842e-03,
					 3.144e-03, 3.461e-03, 3.792e-03, 4.137e-03, 4.497e-03, 4.871e-03, 5.259e-03, 5.661e-03, 6.076e-03, 6.505e-03,
					 6.947e-03, 7.403e-03, 7.871e-03, 8.353e-03, 8.847e-03, 9.354e-03, 9.874e-03, 1.041e-02, 1.095e-02, 1.151e-02,
					 1.207e-02, 1.265e-02, 1.324e-02, 1.384e-02, 1.446e-02, 1.508e-02, 1.572e-02, 1.636e-02, 1.702e-02, 1.768e-02,
					 1.836e-02, 1.905e-02, 1.974e-02, 2.045e-02, 2.116e-02, 2.189e-02, 2.262e-02, 2.336e-02, 2.411e-02, 2.487e-02,
					 2.563e-02, 2.641e-02, 2.719e-02, 2.798e-02, 2.878e-02, 2.958e-02, 3.039e-02, 3.121e-02, 3.203e-02, 3.286e-02,
					 3.370e-02, 3.454e-02, 3.539e-02, 3.625e-02, 3.711e-02, 3.797e-02, 3.884e-02, 3.972e-02, 4.059e-02, 4.148e-02,
					 4.237e-02, 4.326e-02, 4.415e-02, 4.505e-02, 4.595e-02, 4.686e-02, 4.777e-02, 4.868e-02, 4.959e-02, 5.051e-02,
					 5.143e-02, 5.235e-02, 5.327e-02, 5.419e-02, 5.512e-02, 5.604e-02, 5.697e-02, 5.790e-02, 5.883e-02, 5.976e-02,
					 6.068e-02, 6.161e-02, 6.254e-02, 6.347e-02, 6.440e-02, 6.533e-02, 6.625e-02, 6.718e-02, 6.810e-02, 6.902e-02,
					 6.995e-02, 7.086e-02, 7.178e-02, 7.270e-02, 7.361e-02, 7.452e-02, 7.542e-02, 7.633e-02, 7.723e-02, 7.813e-02,
					 7.902e-02, 7.991e-02, 8.080e-02, 8.168e-02, 8.256e-02, 8.343e-02, 8.430e-02, 8.517e-02, 8.603e-02, 8.689e-02,
					 8.774e-02, 8.858e-02, 8.942e-02, 9.026e-02, 9.108e-02, 9.191e-02, 9.272e-02, 9.353e-02, 9.434e-02, 9.514e-02,
					 9.593e-02, 9.671e-02, 9.749e-02, 9.826e-02, 9.902e-02, 9.978e-02, 1.005e-01, 1.013e-01, 1.020e-01, 1.027e-01,
					 1.034e-01, 1.042e-01, 1.049e-01, 1.055e-01, 1.062e-01, 1.069e-01, 1.076e-01, 1.082e-01, 1.089e-01, 1.095e-01,
					 1.102e-01, 1.108e-01, 1.114e-01, 1.120e-01, 1.126e-01, 1.132e-01, 1.137e-01, 1.143e-01, 1.149e-01, 1.154e-01,
					 1.159e-01, 1.165e-01, 1.170e-01, 1.175e-01, 1.180e-01, 1.185e-01, 1.189e-01, 1.194e-01, 1.198e-01, 1.203e-01,
					 1.207e-01, 1.211e-01, 1.215e-01, 1.219e-01, 1.223e-01, 1.227e-01, 1.231e-01, 1.234e-01, 1.238e-01, 1.241e-01,
					 1.244e-01, 1.247e-01, 1.250e-01, 1.253e-01, 1.256e-01, 1.258e-01, 1.261e-01, 1.263e-01, 1.265e-01, 1.267e-01,
					 1.269e-01, 1.271e-01, 1.273e-01, 1.275e-01, 1.276e-01, 1.278e-01, 1.279e-01, 1.280e-01, 1.281e-01, 1.282e-01,
					 1.283e-01, 1.284e-01, 1.285e-01, 1.285e-01, 1.285e-01, 1.286e-01, 1.286e-01, 1.286e-01, 1.286e-01, 1.285e-01,
					 1.285e-01, 1.285e-01, 1.284e-01, 1.283e-01, 1.283e-01, 1.282e-01, 1.281e-01, 1.279e-01, 1.278e-01, 1.277e-01,
					 1.275e-01, 1.274e-01, 1.272e-01, 1.270e-01, 1.268e-01, 1.266e-01, 1.264e-01, 1.262e-01, 1.259e-01, 1.257e-01,
					 1.254e-01, 1.251e-01, 1.248e-01, 1.246e-01, 1.242e-01, 1.239e-01, 1.236e-01, 1.233e-01, 1.229e-01, 1.226e-01,
					 1.222e-01, 1.218e-01, 1.214e-01, 1.210e-01, 1.206e-01, 1.202e-01, 1.198e-01, 1.193e-01, 1.189e-01, 1.184e-01,
					 1.180e-01, 1.175e-01, 1.170e-01, 1.165e-01, 1.160e-01, 1.155e-01, 1.150e-01, 1.144e-01, 1.139e-01, 1.134e-01,
					 1.128e-01, 1.122e-01, 1.117e-01, 1.111e-01, 1.105e-01, 1.099e-01, 1.093e-01, 1.087e-01, 1.081e-01, 1.075e-01,
					 1.068e-01, 1.062e-01, 1.055e-01, 1.049e-01, 1.042e-01, 1.036e-01, 1.029e-01, 1.022e-01, 1.015e-01, 1.008e-01,
					 1.001e-01, 9.943e-02, 9.872e-02, 9.800e-02, 9.728e-02, 9.656e-02, 9.583e-02, 9.509e-02, 9.435e-02, 9.361e-02,
					 9.286e-02, 9.210e-02, 9.134e-02, 9.058e-02, 8.981e-02, 8.904e-02, 8.827e-02, 8.749e-02, 8.671e-02, 8.592e-02,
					 8.514e-02, 8.435e-02, 8.355e-02, 8.275e-02, 8.196e-02, 8.115e-02, 8.035e-02, 7.954e-02, 7.874e-02, 7.793e-02,
					 7.711e-02, 7.630e-02, 7.549e-02, 7.467e-02, 7.385e-02, 7.304e-02, 7.222e-02, 7.140e-02, 7.058e-02, 6.976e-02,
					 6.894e-02, 6.812e-02, 6.730e-02, 6.648e-02, 6.566e-02, 6.484e-02, 6.402e-02, 6.320e-02, 6.238e-02, 6.157e-02,
					 6.075e-02, 5.994e-02, 5.913e-02, 5.832e-02, 5.751e-02, 5.670e-02, 5.590e-02, 5.509e-02, 5.429e-02, 5.349e-02,
					 5.270e-02, 5.191e-02, 5.112e-02, 5.033e-02, 4.955e-02, 4.877e-02, 4.799e-02, 4.722e-02, 4.645e-02, 4.568e-02,
					 4.492e-02, 4.416e-02, 4.341e-02, 4.266e-02, 4.191e-02, 4.117e-02, 4.044e-02, 3.971e-02, 3.898e-02, 3.826e-02,
					 3.754e-02, 3.683e-02, 3.612e-02, 3.542e-02, 3.473e-02, 3.404e-02, 3.335e-02, 3.267e-02, 3.200e-02, 3.133e-02,
					 3.067e-02, 3.002e-02, 2.937e-02, 2.873e-02, 2.809e-02, 2.746e-02, 2.684e-02, 2.623e-02, 2.562e-02, 2.501e-02,
					 2.442e-02, 2.383e-02, 2.325e-02, 2.267e-02, 2.211e-02, 2.154e-02, 2.099e-02, 2.045e-02, 1.991e-02, 1.938e-02,
					 1.885e-02, 1.834e-02, 1.783e-02, 1.733e-02, 1.683e-02, 1.635e-02, 1.587e-02, 1.540e-02, 1.493e-02, 1.448e-02,
					 1.403e-02, 1.359e-02, 1.316e-02, 1.273e-02, 1.232e-02, 1.191e-02, 1.151e-02, 1.111e-02, 1.073e-02, 1.035e-02,
					 9.982e-03, 9.620e-03, 9.266e-03, 8.920e-03, 8.582e-03, 8.251e-03, 7.928e-03, 7.613e-03, 7.305e-03, 7.005e-03,
					 6.712e-03, 6.427e-03, 6.150e-03, 5.880e-03, 5.617e-03, 5.362e-03, 5.114e-03, 4.873e-03, 4.639e-03, 4.413e-03,
					 4.194e-03, 3.981e-03, 3.776e-03, 3.577e-03, 3.385e-03, 3.200e-03, 3.021e-03, 2.849e-03, 2.683e-03, 2.524e-03,
					 2.371e-03, 2.224e-03, 2.083e-03, 1.948e-03, 1.819e-03, 1.696e-03, 1.578e-03, 1.465e-03, 1.359e-03, 1.257e-03,
					 1.161e-03, 1.069e-03, 9.827e-04, 9.011e-04, 8.241e-04, 7.516e-04, 6.835e-04, 6.197e-04, 5.601e-04, 5.045e-04,
					 4.527e-04, 4.047e-04, 3.602e-04, 3.192e-04, 2.815e-04, 2.470e-04, 2.155e-04, 1.869e-04, 1.610e-04, 1.377e-04,
					 1.168e-04, 9.823e-05, 8.181e-05, 6.740e-05, 5.486e-05, 4.405e-05, 3.482e-05, 2.703e-05, 2.056e-05, 1.526e-05,
					 1.100e-05, 7.653e-06, 5.100e-06, 3.216e-06, 1.887e-06, 1.003e-06, 4.621e-07, 1.699e-07, 4.136e-08, 3.678e-09, 0.000e+00 };

    } // namespace HNLSelector

} // namespace HNL
} // namespace genie

#endif // #ifndef _HNL_JFORMFACTORTABLES_H_