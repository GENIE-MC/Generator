//____________________________________________________________________________
/*!

\class    genie::mueloss::MuELMaterial

\brief    Enumeration of materials for which the MuELoss package knows how to
          calculate muon energy losses.

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

#ifndef _MUELOSS_MATERIAL_H_
#define _MUELOSS_MATERIAL_H_

namespace genie   {
namespace mueloss {

typedef enum EMuELMaterial {

  eMuUndefined = 0,

  // ---- elements

  eMuHydrogen         = 101,
  eMuDeuterium        = 102,
  eMuHelium           = 103,
  eMuLithium          = 104,
  eMuBeryllium        = 105,
  eMuBoron            = 106,
  eMuCarbon           = 107,
  eMuNitrogen         = 108,
  eMuOxygen           = 109,
  eMuFluorine         = 110,
  eMuNeon             = 111,
  eMuSodium           = 112,
  eMuMagnesium        = 113,
  eMuAluminium        = 114,
  eMuSilicon          = 115,
  eMuSulphur          = 116,
  eMuChlorine         = 117,
  eMuArgon            = 118,
  eMuPotassium        = 119,
  eMuCalcium          = 120,
  eMuChromium         = 121,
  eMuManganese        = 122,
  eMuIron             = 123,
  eMuNickel           = 124,
  eMuCopper           = 125,
  eMuZinc             = 126,
  eMuGermanium        = 127,
  eMuBromine          = 128,
  eMuTin              = 129,
  eMuIodine           = 130,
  eMuBarium           = 131,
  eMuTungsten         = 132,
  eMuLead             = 133,
  eMuBismuth          = 134,
  eMuUranium          = 135,

  // ---- compound materials 

  eMuBariumFluoride   = 201,
  eMuBismuthGermanate = 202,
  eMuPyrex            = 203, /* <-- SiO2(80%),B2O3(12%),Na2O(5%) */
  eMuCalciumCarbonate = 204,
  eMuConcrete         = 205, /* <-- O2(52.9%),Si(33.7%),Ca(4.4%),Al(3.4%),Na(1.6%),Fe(1.4%),K(1.3%),H2(1%) */
  eMuFreon12          = 206,
  eMuFreon13B1        = 207,
  eMuLeadGlassSF5     = 208, /* <-- PbO(55%),SiO2 (38%),K2O(5%),Na2O(1%) */
  eMuLeadOxide        = 209,
  eMuLithiumFluoride  = 210,
  eMuLucite           = 211,
  eMuPolyethylene     = 212,
  eMuPolystyrene      = 213,
  eMuLiquidPropane    = 214,
  eMuSiliconDioxide   = 215,
  eMuSodiumIodide     = 216,
  eMuStandardRock     = 217,
  eMuUraniumOxide     = 218,
  eMuWater            = 219

} MuELMaterial_t;

class MuELMaterial
{
public:
  //____________________________________________________________________
  static const char * AsString(MuELMaterial_t material)
  {
     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return "Barium Fluoride";   break;
      case eMuBismuthGermanate:  return "Bismuth Germanate"; break;
      case eMuPyrex:             return "Pyrex";             break;
      case eMuCalciumCarbonate:  return "Calcium Carbonate"; break;
      case eMuConcrete:          return "Concrete";          break;
      case eMuFreon12:           return "Freon 12";          break;
      case eMuFreon13B1:         return "Freon 13B1";        break;
      case eMuLeadOxide:         return "Lead Oxide";        break;
      case eMuLithiumFluoride:   return "Lithium Fluoride";  break;
      case eMuLucite:            return "Lucite";            break;
      case eMuPolyethylene:      return "Polyethylene";      break;
      case eMuPolystyrene:       return "Polystyrene";       break;
      case eMuLiquidPropane:     return "Liquid Propane";    break;
      case eMuSiliconDioxide:    return "Silicon Dioxide";   break;
      case eMuSodiumIodide:      return "Sodium Iodide";     break;
      case eMuStandardRock:      return "Standard Rock";     break;
      case eMuUraniumOxide:      return "Uranium Oxide";     break;
      case eMuWater:             return "Water";             break;

      //- elements
      case eMuHydrogen:          return "Hydrogen";          break;
      case eMuDeuterium:         return "Deuterium";         break;
      case eMuHelium:            return "Helium";            break;
      case eMuLithium:           return "Lithium";           break;
      case eMuBeryllium:         return "Beryllium";         break;
      case eMuBoron:             return "Boron";             break;
      case eMuCarbon:            return "Carbon";            break;
      case eMuNitrogen:          return "Nitrogen";          break;
      case eMuOxygen:            return "Oxygen";            break;
      case eMuFluorine:          return "Fluorine";          break;
      case eMuNeon:              return "Neon";              break;
      case eMuSodium:            return "Sodium";            break;
      case eMuMagnesium:         return "Magnesium";         break;
      case eMuAluminium:         return "Aluminium";         break;
      case eMuSilicon:           return "Silicon";           break;
      case eMuSulphur:           return "Sulphur";           break;
      case eMuChlorine:          return "Chlorine";          break;
      case eMuArgon:             return "Argon";             break;
      case eMuPotassium:         return "Potassium";         break;
      case eMuCalcium:           return "Calcium";           break;
      case eMuChromium:          return "Chromium";          break;
      case eMuManganese:         return "Manganese";         break;
      case eMuIron:              return "Iron";              break;
      case eMuNickel:            return "Nickel";            break;
      case eMuCopper:            return "Copper";            break;
      case eMuZinc:              return "Zinc";              break;
      case eMuGermanium:         return "Germanium";         break;
      case eMuBromine:           return "Bromine";           break;
      case eMuTin:               return "Tin";               break;
      case eMuIodine:            return "Iodine";            break;
      case eMuBarium:            return "Barium";            break;
      case eMuTungsten:          return "Tungsten";          break;
      case eMuLead:              return "Lead";              break;
      case eMuBismuth:           return "Bismuth";           break;
      case eMuUranium:           return "Uranium";           break;

      case eMuUndefined:
      default:
         return "*** unknown material ***";
     }
     return "*** unknown material ***";
  }
  //____________________________________________________________________
  // material density in gr/cm^3
  static double Density(MuELMaterial_t material)
  {
     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return  4.830; break;
      case eMuBismuthGermanate:  return  7.100; break;
      case eMuPyrex:             return  2.230; break;
      case eMuCalciumCarbonate:  return  2.800; break;
      case eMuConcrete:          return  2.500; break;
      case eMuFreon12:           return  1.120; break;
      case eMuFreon13B1:         return  1.500; break;
      case eMuLeadOxide:         return  9.530; break;
      case eMuLithiumFluoride:   return  2.635; break;
      case eMuLucite:            return  1.190; break;
      case eMuPolyethylene:      return  0.940; break;
      case eMuPolystyrene:       return  1.060; break;
      case eMuLiquidPropane:     return  0.430; break;
      case eMuSiliconDioxide:    return  2.320; break;
      case eMuSodiumIodide:      return  3.367; break;
      case eMuStandardRock:      return  2.650; break;
      case eMuUraniumOxide:      return 10.960; break;
      case eMuWater:             return  1.000; break;

      //- elements
      case eMuHydrogen:          return  0.063; break;
      case eMuDeuterium:         return  0.140; break;
      case eMuHelium:            return  0.125; break;
      case eMuLithium:           return  0.534; break;
      case eMuBeryllium:         return  1.848; break;
      case eMuBoron:             return  2.370; break;
      case eMuCarbon:            return  2.265; break;
      case eMuNitrogen:          return  0.808; break;
      case eMuOxygen:            return  1.140; break;
      case eMuFluorine:          return  1.108; break;
      case eMuNeon:              return  1.207; break;
      case eMuSodium:            return  0.971; break;
      case eMuMagnesium:         return  1.740; break;
      case eMuAluminium:         return  2.699; break;
      case eMuSilicon:           return  2.330; break;
      case eMuSulphur:           return  2.000; break;
      case eMuChlorine:          return  1.560; break;
      case eMuArgon:             return  1.393; break;
      case eMuPotassium:         return  0.862; break;
      case eMuCalcium:           return  1.550; break;
      case eMuChromium:          return  7.180; break;
      case eMuManganese:         return  7.440; break;
      case eMuIron:              return  7.874; break;
      case eMuNickel:            return  8.902; break;
      case eMuCopper:            return  8.960; break;
      case eMuZinc:              return  7.133; break;
      case eMuGermanium:         return  5.323; break;
      case eMuBromine:           return  3.120; break;
      case eMuTin:               return  7.310; break;
      case eMuIodine:            return  4.930; break;
      case eMuBarium:            return  3.500; break;
      case eMuTungsten:          return 19.300; break;
      case eMuLead:              return 11.350; break;
      case eMuBismuth:           return  9.747; break;
      case eMuUranium:           return 18.950; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
  static double Z(MuELMaterial_t material)
  {
     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return 0.4221; break;   // For compound materials this is
      case eMuBismuthGermanate:  return 0.4207; break;   // actually Z/A, and A will be set to 1
      case eMuPyrex:             return 0.4971; break;   // OK for now because it is Z/A we
      case eMuCalciumCarbonate:  return 0.4996; break;   // use... but *change that*
      case eMuConcrete:          return 0.5027; break;
      case eMuFreon12:           return 0.4797; break;
      case eMuFreon13B1:         return 0.4567; break;
      case eMuLeadOxide:         return 0.4032; break;
      case eMuLithiumFluoride:   return 0.4626; break;
      case eMuLucite:            return 0.5394; break;
      case eMuPolyethylene:      return 0.5703; break;
      case eMuPolystyrene:       return 0.5377; break;
      case eMuLiquidPropane:     return 0.5896; break;
      case eMuSiliconDioxide:    return 0.4993; break;
      case eMuSodiumIodide:      return 0.4270; break;
      case eMuStandardRock:      return 0.5000; break;
      case eMuUraniumOxide:      return 0.4000; break;
      case eMuWater:             return 0.5551; break;

      //- elements
      case eMuHydrogen:          return  1.; break;
      case eMuDeuterium:         return  1.; break;
      case eMuHelium:            return  2.; break;
      case eMuLithium:           return  3.; break;
      case eMuBeryllium:         return  4.; break;
      case eMuBoron:             return  5.; break;
      case eMuCarbon:            return  6.; break;
      case eMuNitrogen:          return  7.; break;
      case eMuOxygen:            return  8.; break;
      case eMuFluorine:          return  9.; break;
      case eMuNeon:              return 10.; break;
      case eMuSodium:            return 11.; break;
      case eMuMagnesium:         return 12.; break;
      case eMuAluminium:         return 13.; break;
      case eMuSilicon:           return 14.; break;
      case eMuSulphur:           return 16.; break;
      case eMuChlorine:          return 17.; break;
      case eMuArgon:             return 18.; break;
      case eMuPotassium:         return 19.; break;
      case eMuCalcium:           return 20.; break;
      case eMuChromium:          return 24.; break;
      case eMuManganese:         return 25.; break;
      case eMuIron:              return 26.; break;
      case eMuNickel:            return 28.; break;
      case eMuCopper:            return 29.; break;
      case eMuZinc:              return 30.; break;
      case eMuGermanium:         return 32.; break;
      case eMuBromine:           return 35.; break;
      case eMuTin:               return 50.; break;
      case eMuIodine:            return 53.; break;
      case eMuBarium:            return 56.; break;
      case eMuTungsten:          return 74.; break;
      case eMuLead:              return 82.; break;
      case eMuBismuth:           return 83.; break;
      case eMuUranium:           return 92.; break;

      case eMuUndefined:
      default:
         return 0;
     }
     return 0;
  }
  //____________________________________________________________________
  static double A(MuELMaterial_t material)
  {
     switch(material) {

      //- compound materials
      case eMuBariumFluoride:    return 1.0; break;    // A for compound materials is set to 1
      case eMuBismuthGermanate:  return 1.0; break;    // because Z was set to Z/A.
      case eMuPyrex:             return 1.0; break;    // OK for now because it is Z/A we
      case eMuCalciumCarbonate:  return 1.0; break;    // use... but *change that*
      case eMuConcrete:          return 1.0; break;
      case eMuFreon12:           return 1.0; break;
      case eMuFreon13B1:         return 1.0; break;
      case eMuLeadOxide:         return 1.0; break;
      case eMuLithiumFluoride:   return 1.0; break;
      case eMuLucite:            return 1.0; break;
      case eMuPolyethylene:      return 1.0; break;
      case eMuPolystyrene:       return 1.0; break;
      case eMuLiquidPropane:     return 1.0; break;
      case eMuSiliconDioxide:    return 1.0; break;
      case eMuSodiumIodide:      return 1.0; break;
      case eMuStandardRock:      return 1.0; break;
      case eMuUraniumOxide:      return 1.0; break;
      case eMuWater:             return 1.0; break;

      //- elements
      case eMuHydrogen:          return   1.008; break;
      case eMuDeuterium:         return   2.014; break;
      case eMuHelium:            return   4.003; break;
      case eMuLithium:           return   6.940; break;
      case eMuBeryllium:         return   9.012; break;
      case eMuBoron:             return  10.810; break;
      case eMuCarbon:            return  12.011; break;
      case eMuNitrogen:          return  14.007; break;
      case eMuOxygen:            return  15.999; break;
      case eMuFluorine:          return  18.998; break;
      case eMuNeon:              return  20.170; break;
      case eMuSodium:            return  22.990; break;
      case eMuMagnesium:         return  24.305; break;
      case eMuAluminium:         return  26.982; break;
      case eMuSilicon:           return  28.086; break;
      case eMuSulphur:           return  32.060; break;
      case eMuChlorine:          return  35.453; break;
      case eMuArgon:             return  39.948; break;
      case eMuPotassium:         return  39.098; break;
      case eMuCalcium:           return  40.080; break;
      case eMuChromium:          return  51.996; break;
      case eMuManganese:         return  54.938; break;
      case eMuIron:              return  55.847; break;
      case eMuNickel:            return  58.710; break;
      case eMuCopper:            return  63.546; break;
      case eMuZinc:              return  65.380; break;
      case eMuGermanium:         return  72.590; break;
      case eMuBromine:           return  79.904; break;
      case eMuTin:               return 118.690; break;
      case eMuIodine:            return 126.905; break;
      case eMuBarium:            return 137.330; break;
      case eMuTungsten:          return 183.850; break;
      case eMuLead:              return 207.200; break;
      case eMuBismuth:           return 208.980; break;
      case eMuUranium:           return 238.029; break;

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

#endif // _MUELOSS_MATERIAL_H_

