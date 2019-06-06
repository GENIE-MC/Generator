//____________________________________________________________________________
/*!

\class    genie::mueloss::BetheBlochMaterialParams

\brief    Bethe Bloch parameters for various materials

\ref      W.Lohmann, R.Kopp and R.Voss,
          Energy Loss of Muons in the Energy Range 1-10000 GeV, CERN 85-03

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 10, 2003

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BETHE_BLOCH_MATERIAL_PARAMS_H_
#define _BETHE_BLOCH_MATERIAL_PARAMS_H_

namespace genie   {
namespace mueloss {

#include "Physics/MuonEnergyLoss/MuELMaterial.h"

class BetheBlochMaterialParams
{
public:
  //____________________________________________________________________
  static double IonizationPotential(MuELMaterial_t material) 
  {
   // returns the ionization potential for the input material (in eV)

     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return 375.9; break;
      case eMuBismuthGermanate:  return 534.1; break;
      case eMuPyrex:             return 134.0; break;
      case eMuCalciumCarbonate:  return 136.4; break;
      case eMuConcrete:          return 135.2; break;
      case eMuFreon12:           return 143.0; break;
      case eMuFreon13B1:         return 210.5; break;
      case eMuLeadOxide:         return 766.7; break;
      case eMuLithiumFluoride:   return  94.0; break;
      case eMuLucite:            return  74.0; break;
      case eMuPolyethylene:      return  57.4; break;
      case eMuPolystyrene:       return  68.7; break;
      case eMuLiquidPropane:     return  52.0; break;
      case eMuSiliconDioxide:    return 139.2; break;
      case eMuSodiumIodide:      return 452.0; break;
      case eMuStandardRock:      return 136.4; break;
      case eMuUraniumOxide:      return 720.6; break;
      case eMuWater:             return  75.0; break;

      //- elements
      case eMuHydrogen:          return  21.8; break;
      case eMuDeuterium:         return  21.8; break;
      case eMuHelium:            return  41.8; break;
      case eMuLithium:           return  40.0; break;
      case eMuBeryllium:         return  63.7; break;
      case eMuBoron:             return  76.0; break;
      case eMuCarbon:            return  78.0; break;
      case eMuNitrogen:          return  82.0; break;
      case eMuOxygen:            return  95.0; break;
      case eMuFluorine:          return 115.0; break;
      case eMuNeon:              return 137.0; break;
      case eMuSodium:            return 149.0; break;
      case eMuMagnesium:         return 156.0; break;
      case eMuAluminium:         return 166.0; break;
      case eMuSilicon:           return 173.0; break;
      case eMuSulphur:           return 180.0; break;
      case eMuChlorine:          return 174.0; break;
      case eMuArgon:             return 188.0; break;
      case eMuPotassium:         return 190.0; break;
      case eMuCalcium:           return 191.0; break;
      case eMuChromium:          return 257.0; break;
      case eMuManganese:         return 272.0; break;
      case eMuIron:              return 286.0; break;
      case eMuNickel:            return 311.0; break;
      case eMuCopper:            return 322.0; break;
      case eMuZinc:              return 330.0; break;
      case eMuGermanium:         return 350.0; break;
      case eMuBromine:           return 343.0; break;
      case eMuTin:               return 488.0; break;
      case eMuIodine:            return 491.0; break;
      case eMuBarium:            return 491.0; break;
      case eMuTungsten:          return 727.0; break;
      case eMuLead:              return 823.0; break;
      case eMuBismuth:           return 823.0; break;
      case eMuUranium:           return 890.0; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
  static double DensityCorrection_C(MuELMaterial_t material)
  {
   // returns the density correction factor C

     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return -5.412; break;
      case eMuBismuthGermanate:  return -5.741; break;
      case eMuPyrex:             return -3.971; break;
      case eMuCalciumCarbonate:  return -3.774; break;
      case eMuConcrete:          return -3.946; break;
      case eMuFreon12:           return -4.825; break;
      case eMuFreon13B1:         return -5.356; break;
      case eMuLeadOxide:         return -6.216; break;
      case eMuLithiumFluoride:   return -3.167; break;
      case eMuLucite:            return -3.330; break;
      case eMuPolyethylene:      return -3.002; break;
      case eMuPolystyrene:       return -3.300; break;
      case eMuLiquidPropane:     return -3.553; break;
      case eMuSiliconDioxide:    return -4.003; break;
      case eMuSodiumIodide:      return -6.057; break;
      case eMuStandardRock:      return -3.774; break;
      case eMuUraniumOxide:      return -5.961; break;
      case eMuWater:             return -3.502; break;

      //- elements
      case eMuHydrogen:          return -3.263; break;
      case eMuDeuterium:         return -2.942; break;
      case eMuHelium:            return -4.517; break;
      case eMuLithium:           return -3.122; break;
      case eMuBeryllium:         return -2.785; break;
      case eMuBoron:             return -2.848; break;
      case eMuCarbon:            return -2.868; break;
      case eMuNitrogen:          return -3.998; break;
      case eMuOxygen:            return -3.948; break;
      case eMuFluorine:          return -4.413; break;
      case eMuNeon:              return -4.632; break;
      case eMuSodium:            return -5.053; break;
      case eMuMagnesium:         return -4.530; break;
      case eMuAluminium:         return -4.240; break;
      case eMuSilicon:           return -4.435; break;
      case eMuSulphur:           return -4.666; break;
      case eMuChlorine:          return -4.887; break;
      case eMuArgon:             return -5.217; break;
      case eMuPotassium:         return -5.642; break;
      case eMuCalcium:           return -5.040; break;
      case eMuChromium:          return -4.178; break;
      case eMuManganese:         return -4.270; break;
      case eMuIron:              return -4.291; break;
      case eMuNickel:            return -4.312; break;
      case eMuCopper:            return -4.419; break;
      case eMuZinc:              return -4.691; break;
      case eMuGermanium:         return -5.141; break;
      case eMuBromine:           return -5.641; break;
      case eMuTin:               return -5.534; break;
      case eMuIodine:            return -5.949; break;
      case eMuBarium:            return -6.315; break;
      case eMuTungsten:          return -5.406; break;
      case eMuLead:              return -6.202; break;
      case eMuBismuth:           return -6.351; break;
      case eMuUranium:           return -5.869; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
  static double DensityCorrection_X0(MuELMaterial_t material)
  {
   // returns the density correction constant X0

     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return -0.010; break;
      case eMuBismuthGermanate:  return  0.046; break;
      case eMuPyrex:             return  0.148; break;
      case eMuCalciumCarbonate:  return  0.049; break;
      case eMuConcrete:          return  0.130; break;
      case eMuFreon12:           return  0.304; break;
      case eMuFreon13B1:         return  0.352; break;
      case eMuLeadOxide:         return  0.036; break;
      case eMuLithiumFluoride:   return  0.017; break;
      case eMuLucite:            return  0.182; break;
      case eMuPolyethylene:      return  0.137; break;
      case eMuPolystyrene:       return  0.165; break;
      case eMuLiquidPropane:     return  0.286; break;
      case eMuSiliconDioxide:    return  0.139; break;
      case eMuSodiumIodide:      return  0.120; break;
      case eMuStandardRock:      return  0.049; break;
      case eMuUraniumOxide:      return -0.194; break;
      case eMuWater:             return  0.240; break;

      //- elements
      case eMuHydrogen:          return  0.476; break;
      case eMuDeuterium:         return  0.200; break;
      case eMuHelium:            return  0.473; break;
      case eMuLithium:           return  0.130; break;
      case eMuBeryllium:         return  0.059; break;
      case eMuBoron:             return  0.031; break;
      case eMuCarbon:            return -0.018; break;
      case eMuNitrogen:          return  0.304; break;
      case eMuOxygen:            return  0.287; break;
      case eMuFluorine:          return  0.200; break;
      case eMuNeon:              return  0.200; break;
      case eMuSodium:            return  0.288; break;
      case eMuMagnesium:         return  0.150; break;
      case eMuAluminium:         return  0.171; break;
      case eMuSilicon:           return  0.201; break;
      case eMuSulphur:           return  0.158; break;
      case eMuChlorine:          return  0.200; break;
      case eMuArgon:             return  0.201; break;
      case eMuPotassium:         return  0.385; break;
      case eMuCalcium:           return  0.323; break;
      case eMuChromium:          return  0.034; break;
      case eMuManganese:         return  0.045; break;
      case eMuIron:              return -0.001; break;
      case eMuNickel:            return -0.057; break;
      case eMuCopper:            return -0.025; break;
      case eMuZinc:              return  0.005; break;
      case eMuGermanium:         return  0.338; break;
      case eMuBromine:           return  0.339; break;
      case eMuTin:               return  0.288; break;
      case eMuIodine:            return  0.055; break;
      case eMuBarium:            return  0.419; break;
      case eMuTungsten:          return  0.217; break;
      case eMuLead:              return  0.378; break;
      case eMuBismuth:           return  0.415; break;
      case eMuUranium:           return  0.226; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
  static double DensityCorrection_X1(MuELMaterial_t material)
  {
   // returns the density correction constant X1

     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return 3.387; break;
      case eMuBismuthGermanate:  return 3.782; break;
      case eMuPyrex:             return 2.993; break;
      case eMuCalciumCarbonate:  return 3.055; break;
      case eMuConcrete:          return 3.047; break;
      case eMuFreon12:           return 3.266; break;
      case eMuFreon13B1:         return 3.755; break;
      case eMuLeadOxide:         return 3.546; break;
      case eMuLithiumFluoride:   return 2.705; break;
      case eMuLucite:            return 2.668; break;
      case eMuPolyethylene:      return 2.518; break;
      case eMuPolystyrene:       return 2.503; break;
      case eMuLiquidPropane:     return 2.657; break;
      case eMuSiliconDioxide:    return 3.003; break;
      case eMuSodiumIodide:      return 3.592; break;
      case eMuStandardRock:      return 3.055; break;
      case eMuUraniumOxide:      return 3.529; break;
      case eMuWater:             return 2.800; break;

      //- elements
      case eMuHydrogen:          return 1.922; break;
      case eMuDeuterium:         return 2.000; break;
      case eMuHelium:            return 2.000; break;
      case eMuLithium:           return 1.640; break;
      case eMuBeryllium:         return 1.692; break;
      case eMuBoron:             return 1.969; break;
      case eMuCarbon:            return 2.342; break;
      case eMuNitrogen:          return 2.000; break;
      case eMuOxygen:            return 2.000; break;
      case eMuFluorine:          return 3.000; break;
      case eMuNeon:              return 3.000; break;
      case eMuSodium:            return 3.196; break;
      case eMuMagnesium:         return 3.067; break;
      case eMuAluminium:         return 3.013; break;
      case eMuSilicon:           return 2.872; break;
      case eMuSulphur:           return 2.716; break;
      case eMuChlorine:          return 3.000; break;
      case eMuArgon:             return 3.000; break;
      case eMuPotassium:         return 3.172; break;
      case eMuCalcium:           return 3.119; break;
      case eMuChromium:          return 3.045; break;
      case eMuManganese:         return 3.107; break;
      case eMuIron:              return 3.153; break;
      case eMuNickel:            return 3.185; break;
      case eMuCopper:            return 3.279; break;
      case eMuZinc:              return 3.367; break;
      case eMuGermanium:         return 3.610; break;
      case eMuBromine:           return 3.000; break;
      case eMuTin:               return 3.296; break;
      case eMuIodine:            return 3.260; break;
      case eMuBarium:            return 3.455; break;
      case eMuTungsten:          return 3.496; break;
      case eMuLead:              return 3.807; break;
      case eMuBismuth:           return 3.825; break;
      case eMuUranium:           return 3.372; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
  static double DensityCorrection_a(MuELMaterial_t material)
  {
   // returns the density correction constant a

     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return 0.160; break;
      case eMuBismuthGermanate:  return 0.096; break;
      case eMuPyrex:             return 0.083; break;
      case eMuCalciumCarbonate:  return 0.083; break;
      case eMuConcrete:          return 0.075; break;
      case eMuFreon12:           return 0.080; break;
      case eMuFreon13B1:         return 0.039; break;
      case eMuLeadOxide:         return 0.196; break;
      case eMuLithiumFluoride:   return 0.076; break;
      case eMuLucite:            return 0.114; break;
      case eMuPolyethylene:      return 0.121; break;
      case eMuPolystyrene:       return 0.165; break;
      case eMuLiquidPropane:     return 0.103; break;
      case eMuSiliconDioxide:    return 0.084; break;
      case eMuSodiumIodide:      return 0.125; break;
      case eMuStandardRock:      return 0.083; break;
      case eMuUraniumOxide:      return 0.205; break;
      case eMuWater:             return 0.091; break;

      //- elements
      case eMuHydrogen:          return 0.135; break;
      case eMuDeuterium:         return 0.347; break;
      case eMuHelium:            return 0.657; break;
      case eMuLithium:           return 0.951; break;
      case eMuBeryllium:         return 0.804; break;
      case eMuBoron:             return 0.562; break;
      case eMuCarbon:            return 0.261; break;
      case eMuNitrogen:          return 0.533; break;
      case eMuOxygen:            return 0.523; break;
      case eMuFluorine:          return 0.159; break;
      case eMuNeon:              return 0.169; break;
      case eMuSodium:            return 0.078; break;
      case eMuMagnesium:         return 0.082; break;
      case eMuAluminium:         return 0.080; break;
      case eMuSilicon:           return 0.146; break;
      case eMuSulphur:           return 0.340; break;
      case eMuChlorine:          return 0.181; break;
      case eMuArgon:             return 0.196; break;
      case eMuPotassium:         return 0.198; break;
      case eMuCalcium:           return 0.156; break;
      case eMuChromium:          return 0.154; break;
      case eMuManganese:         return 0.150; break;
      case eMuIron:              return 0.147; break;
      case eMuNickel:            return 0.165; break;
      case eMuCopper:            return 0.143; break;
      case eMuZinc:              return 0.147; break;
      case eMuGermanium:         return 0.072; break;
      case eMuBromine:           return 0.217; break;
      case eMuTin:               return 0.187; break;
      case eMuIodine:            return 0.238; break;
      case eMuBarium:            return 0.183; break;
      case eMuTungsten:          return 0.155; break;
      case eMuLead:              return 0.094; break;
      case eMuBismuth:           return 0.094; break;
      case eMuUranium:           return 0.197; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
  static double DensityCorrection_m(MuELMaterial_t material)
  {
   // returns the density correction constant m

     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return 2.887; break;
      case eMuBismuthGermanate:  return 3.078; break;
      case eMuPyrex:             return 3.522; break;
      case eMuCalciumCarbonate:  return 3.412; break;
      case eMuConcrete:          return 3.547; break;
      case eMuFreon12:           return 3.463; break;
      case eMuFreon13B1:         return 3.719; break;
      case eMuLeadOxide:         return 2.730; break;
      case eMuLithiumFluoride:   return 3.748; break;
      case eMuLucite:            return 3.384; break;
      case eMuPolyethylene:      return 3.429; break;
      case eMuPolystyrene:       return 3.222; break;
      case eMuLiquidPropane:     return 3.562; break;
      case eMuSiliconDioxide:    return 3.506; break;
      case eMuSodiumIodide:      return 3.040; break;
      case eMuStandardRock:      return 3.412; break;
      case eMuUraniumOxide:      return 2.671; break;
      case eMuWater:             return 3.477; break;

      //- elements
      case eMuHydrogen:          return 5.625; break;
      case eMuDeuterium:         return 3.000; break;
      case eMuHelium:            return 3.000; break;
      case eMuLithium:           return 2.499; break;
      case eMuBeryllium:         return 2.434; break;
      case eMuBoron:             return 2.451; break;
      case eMuCarbon:            return 2.870; break;
      case eMuNitrogen:          return 3.000; break;
      case eMuOxygen:            return 3.000; break;
      case eMuFluorine:          return 3.000; break;
      case eMuNeon:              return 3.000; break;
      case eMuSodium:            return 3.645; break;
      case eMuMagnesium:         return 3.617; break;
      case eMuAluminium:         return 3.635; break;
      case eMuSilicon:           return 3.255; break;
      case eMuSulphur:           return 2.646; break;
      case eMuChlorine:          return 3.000; break;
      case eMuArgon:             return 3.000; break;
      case eMuPotassium:         return 2.923; break;
      case eMuCalcium:           return 3.075; break;
      case eMuChromium:          return 2.990; break;
      case eMuManganese:         return 2.980; break;
      case eMuIron:              return 2.963; break;
      case eMuNickel:            return 2.843; break;
      case eMuCopper:            return 2.904; break;
      case eMuZinc:              return 2.865; break;
      case eMuGermanium:         return 3.331; break;
      case eMuBromine:           return 3.000; break;
      case eMuTin:               return 2.858; break;
      case eMuIodine:            return 2.728; break;
      case eMuBarium:            return 2.891; break;
      case eMuTungsten:          return 2.845; break;
      case eMuLead:              return 3.161; break;
      case eMuBismuth:           return 3.167; break;
      case eMuUranium:           return 2.817; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
};

}      // mueloss namespace
}      // genie   namespace

#endif // _BETHE_BLOCH_MATERIAL_PARAMS_H_

