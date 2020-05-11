#include "Tools/EvtLib/Utils.h"

#include "Framework/Messenger/Messenger.h"

#include <wordexp.h>

//___________________________________________________________________________
void genie::evtlib::Expand(std::string& s)
{
  wordexp_t p;
  const int status = wordexp(s.c_str(), &p, WRDE_SHOWERR | WRDE_UNDEF);
  if(status != 0){
    LOG("EvtLib", pFATAL) << "String '" << s
                          << "' returned error " << status << " from wordexp().";
    exit(1);
  }

  if(p.we_wordc == 0){
    LOG("EvtLib", pFATAL) << "String '" << s
                          << "' didn't expand to anything.";
    exit(1);
  }

  if(p.we_wordc > 1){
    LOG("EvtLib", pFATAL) << "String '" << s
                          << "' expanded to " << p.we_wordc << " locations.";
    exit(1);
  }

  s = p.we_wordv[0];

  wordfree(&p);
}
