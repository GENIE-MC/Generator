//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Oct 08, 2009 - CA
   That file was added in 2.5.1
 @ Apr 20, 2012 - CA
   Added LocalTimeAsString(string format) to tag validation program outputs.

*/
//____________________________________________________________________________

#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>

#include <dirent.h>
#include <ctime>

#include <TSystem.h>


#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"

using std::atoi;

//___________________________________________________________________________
vector<string>
  genie::utils::system::GetAllFilesInPath(string path, string extension)
{
  vector<string> files;

  DIR *dp = NULL;
  struct dirent *dirp = NULL;
  if((dp = opendir(path.c_str())) == NULL) {
    LOG("System", pERROR) << "Can not open directory: " << path;
  }
  else {
    while ((dirp = readdir(dp)) != NULL) {
      string filename = path + "/" + string(dirp->d_name);
      bool match = false;
      if(extension.size()==0) match=true; // get all files
      else {
        // match extension
        vector<string> filenamev = genie::utils::str::Split(filename,".");
        if(filenamev.size()==0) break;
        string file_extension = filenamev[filenamev.size()-1];
        match = (file_extension.find(extension) != string::npos);
      }
      if(match) {
        files.push_back(filename);
      }
    }
    closedir(dp);
  }
  return files;
}
//___________________________________________________________________________
int genie::utils::system::GenieMajorVrsNum(string tag)
{
  vector<string> vrs = utils::str::Split(tag,".");
  assert(vrs.size() == 3);
  return atoi(vrs[0].c_str());
}
//___________________________________________________________________________
int genie::utils::system::GenieMinorVrsNum(string tag)
{
  vector<string> vrs = utils::str::Split(tag,".");
  assert(vrs.size() == 3);
  return atoi(vrs[1].c_str());
}
//___________________________________________________________________________
int genie::utils::system::GenieRevisVrsNum(string tag)
{
  vector<string> vrs = utils::str::Split(tag,".");
  assert(vrs.size() == 3);
  return atoi(vrs[2].c_str());
}
//___________________________________________________________________________
bool genie::utils::system::FileExists(string filename)
{
  if(filename.size() == 0) return false;

  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (is_accessible) return true;

  return false;
}

//___________________________________________________________________________

bool genie::utils::system::DirectoryExists( const char * path ) {

  struct stat info;

  if( stat( path, &info ) != 0 ) {
    return false ;
  } else if(info.st_mode & S_IFDIR) {
    return true ;
  } else {
    return false ;
  }
}


//___________________________________________________________________________
string genie::utils::system::LocalTimeAsString(string format)
{
  time_t now = time(0);
  tm* local_time = localtime(&now);

  int yr  = local_time->tm_year + 1900;
  int mon = local_time->tm_mon + 1;
  int day = local_time->tm_mday;
  int hr  = local_time->tm_hour + 1;
  int min = local_time->tm_min;
  int sec = local_time->tm_sec;

  // daylight saving
  if(local_time->tm_isdst > 0)
  {
    if(hr > 0) {
       hr--;
    }
    else
    if(hr == 0) {
       hr = 23;
       if(day > 1) {
          day--;
       }
       else {
          mon--;
          if(mon == 1 || mon == 3 || mon ==  5 ||
             mon == 7 || mon == 8 || mon == 10 || mon == 12)
          {
             day = 31;
          }
          else
          if(mon == 2)
          {
             day = 28;
          }
          else
          {
             day = 30;
          }
       }
    }
  }

  string local_time_as_string =
        Form(format.c_str(),yr,mon,day,hr,min,sec);

  return local_time_as_string;
}
//___________________________________________________________________________
