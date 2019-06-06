//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 06, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jul 14, 2009 - RH
   Tweak PrintBanner() to fix problem with random garbage trailing the banner.
*/
//____________________________________________________________________________

#include <iostream>
#include <fstream>
#include <sstream>

#include <TSystem.h>

#include "Framework/Conventions/GVersion.h"
#include "Framework/Utils/PrintUtils.h"

using std::ostringstream;
using std::cout;
using std::endl;
using std::ios;
using std::ifstream;

//____________________________________________________________________________
string genie::utils::print::P4AsString(const TLorentzVector * p)
{
  ostringstream fmt;

  fmt << "(E = " << p->Energy()
      << ", px = " << p->Px()
      << ", py = " << p->Py()
      << ", pz = " << p->Pz() << ")";

  double m2 = p->Mag2();
  if(m2>0.) fmt << " / M = " << TMath::Sqrt(m2);
  else      fmt << " / M^2 = " << m2;

  fmt << " / P = " << p->P();

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::P4AsShortString(const TLorentzVector * p)
{
  ostringstream fmt;

  fmt << "(E = " << p->Energy()
      << ", px = " << p->Px()
      << ", py = " << p->Py()
      << ", pz = " << p->Pz() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::X4AsString(const TLorentzVector * vec4)
{
  ostringstream fmt;

  fmt << "(t = "  << vec4->T()
      << ", x = " << vec4->X()
      << ", y = " << vec4->Y()
      << ", z = " << vec4->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::P3AsString(const TVector3 * vec)
{
  ostringstream fmt;

  fmt << "(px = "  << vec->X()
      << ", py = " << vec->Y()
      << ", pz = " << vec->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::Vec3AsString(const TVector3 * vec)
{
  ostringstream fmt;

  fmt << "(x = "  << vec->X()
      << ", y = " << vec->Y()
      << ", z = " << vec->Z() << ")";

  return fmt.str();
}
//____________________________________________________________________________
string genie::utils::print::BoolAsString(bool b)
{
  return BoolAsTFString(b);
}
//____________________________________________________________________________
string genie::utils::print::BoolAsTFString(bool b)
{
  if(b) return "true";
  else  return "false";
}
//____________________________________________________________________________
string genie::utils::print::BoolAsIOString(bool b)
{
  if(b) return "ON";
  else  return "OFF";
}
//____________________________________________________________________________
string genie::utils::print::BoolAsYNString(bool b)
{
  if(b) return "YES";
  else  return "NO";
}
//____________________________________________________________________________
void genie::utils::print::PrintBanner(void)
{
// loads & prints the GENIE banner

  string base_dir = string(gSystem->Getenv("GENIE"));

#ifdef __GENIE_DEVEL_VERSION__
  string warn_dev_banner = 
      base_dir + 
      string("/data/logo/warning_development_version.txt");
  PrintBanner(warn_dev_banner, 0);
#endif

#ifdef __GENIE_RELEASE_CANDIDATE__ 
  string warn_rc_banner = 
      base_dir + 
      string("/data/logo/warning_release_candidate.txt");
  PrintBanner(warn_rc_banner, 0);
#endif

  string main_banner = 
      base_dir + 
      string("/data/logo/genie_banner_long.txt");
  PrintBanner(main_banner, 0);
}
//___________________________________________________________________________
void genie::utils::print::PrintBanner(string filename, UInt_t wait_msec)
{
  ifstream banner(filename.c_str(), ios::in);

  if( banner.is_open() ) {
      banner.seekg(0, ios::end);

      int    length = banner.tellg();
      char * buffer = new char[length];

      banner.seekg(0, ios::beg);
      banner.read(buffer, length);

      //cout << "\n\n" << buffer << "\n" << endl;
      cout << "\n\n";
      cout.write(buffer,length);
      cout << "\n" << endl;

      delete [] buffer;

      gSystem->Sleep(wait_msec); // watch the banner for a while
  }
}
//___________________________________________________________________________
string genie::utils::print::PrintFramedMesg(
                                  string mesg, unsigned int nl, const char f)
{
  string frame(4+mesg.size(),f);

  string framed_mesg = string("\n") + 
                       frame + string("\n") + 
                       string("  ") + mesg + string("  ") + string("\n") + 
                       frame;

  for(unsigned il=0; il<nl; il++) { framed_mesg += string("\n"); }

  return framed_mesg;
}
//___________________________________________________________________________
