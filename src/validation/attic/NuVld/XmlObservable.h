//____________________________________________________________________________
/*!

\class    genie::nuvld::XmlObservable

\brief                                                           

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug, 2003

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _OBSERVABLE_H_
#define _OBSERVABLE_H_

#include <string>

using std::string;

namespace genie {
namespace nuvld {
  
typedef enum EXmlObservable {
  
   e_undefined = -1,
   e_qes_xsec,
   e_spp_xsec,
   e_mpp_xsec,
   e_tot_xsec,
   e_coh_xsec,
   e_qes_diff_xsec,
   e_spp_diff_xsec,
   e_tot_diff_xsec,
   e_electron_diff_xsec,
   e_f2,
   e_xf3

} XmlObservable_t;


class XmlObservable {

public:

    //_____________________________________________________________________________________
    static const char * AsString(XmlObservable_t obs) {
       
        switch(obs) {
          
        case e_qes_xsec:            return "v QES cross section"; 
                                    break;
        case e_spp_xsec:            return "v SPP cross section"; 
                                    break;
        case e_mpp_xsec:            return "v MPP cross section"; 
                                    break;
        case e_tot_xsec:            return "v TOT cross section"; 
                                    break;
        case e_coh_xsec:            return "v COH pi production cross section"; 
                                    break;
        case e_qes_diff_xsec:       return "v QES differential cross section"; 
                                    break;
        case e_spp_diff_xsec:       return "v SPP differential cross section"; 
                                    break;
        case e_tot_diff_xsec:       return "v TOT differential cross sections"; 
                                    break;
        case e_electron_diff_xsec:  return "electron differential cross sections"; 
                                    break;
        case e_f2:                  return "structure function F2";
                                    break;
        case e_xf3:                 return "structure function xF3";
                                    break;
        case e_undefined:           return "unknown XmlObservable";      
                                    break;  	
        default:                    return "unknown XmlObservable";      
                                    break;  
        }
        return "unknown XmlObservable";      
    }
    //_____________________________________________________________________________________
    static const char * Code(XmlObservable_t obs) {
       
        switch(obs) {
        case e_qes_xsec:            return "QES_XSEC"; 
                                    break;
        case e_spp_xsec:            return "SPP_XSEC"; 
                                    break;
        case e_mpp_xsec:            return "MPP_XSEC"; 
                                    break;
        case e_tot_xsec:            return "TOT_XSEC"; 
                                    break;
        case e_coh_xsec:            return "COH_XSEC"; 
                                    break;
        case e_qes_diff_xsec:       return "QES_PXSEC"; 
                                    break;
        case e_spp_diff_xsec:       return "SPP_PXSEC"; 
                                    break;
        case e_tot_diff_xsec:       return "TOT_PXSEC"; 
                                    break;
        case e_electron_diff_xsec:  return "ELEC_PXSEC"; 
                                    break;
        case e_f2:                  return "F2"; 
                                    break;
        case e_xf3:                 return "xF3"; 
                                    break;
        case e_undefined:        
        default:                    return "unknown XmlObservable";      
                                    break;  
        }
        return "unknown XmlObservable";      
    }
    //_____________________________________________________________________________________
    static XmlObservable_t GetXmlObservable(string str_observable) {

       if     ( str_observable.compare("qes_xsec")           == 0 ) return e_qes_xsec;
       else if( str_observable.compare("spp_xsec")           == 0 ) return e_spp_xsec;
       else if( str_observable.compare("mpp_xsec")           == 0 ) return e_mpp_xsec;
       else if( str_observable.compare("tot_xsec")           == 0 ) return e_tot_xsec;
       else if( str_observable.compare("coh_xsec")           == 0 ) return e_coh_xsec;
       else if( str_observable.compare("qes_diff_xsec")      == 0 ) return e_qes_diff_xsec;
       else if( str_observable.compare("spp_diff_xsec")      == 0 ) return e_spp_diff_xsec;
       else if( str_observable.compare("tot_diff_xsec")      == 0 ) return e_tot_diff_xsec;
       else if( str_observable.compare("electron_diff_xsec") == 0 ) return e_electron_diff_xsec;
       else if( str_observable.compare("f2")                 == 0 ) return e_f2;
       else if( str_observable.compare("xf3")                == 0 ) return e_xf3;
       else                                                         return e_undefined;
    }
    //_____________________________________________________________________________________

};

} // nuvld namespace
} // genie namespace

#endif // _OBSERVABLE_H_
