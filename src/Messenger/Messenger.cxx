//____________________________________________________________________________
/*!

\class    genie::Messenger

\brief    A more convenient interface to the log4cpp Message Service

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 16, 2004

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"

using namespace genie;

using std::cout;

//____________________________________________________________________________
Messenger * Messenger::fInstance = 0;
//____________________________________________________________________________
Messenger::Messenger()
{
  fInstance =  0;
}
//____________________________________________________________________________
Messenger::~Messenger()
{
  fInstance = 0;
}
//____________________________________________________________________________
Messenger * Messenger::Instance()
{
  if(fInstance == 0) {

    static Messenger::Cleaner cleaner;

    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new Messenger;

    log4cpp::Appender * appender;

    appender = new log4cpp::OstreamAppender("default", &cout);
    appender->setLayout(new log4cpp::BasicLayout());

    log4cpp::Category & MSG = log4cpp::Category::getRoot();

    MSG.setAdditivity(false);
    MSG.addAppender(appender);
  }
  return fInstance;
}
//____________________________________________________________________________
log4cpp::Category & Messenger::operator () (const char * stream)
{
  log4cpp::Category & MSG = log4cpp::Category::getInstance(stream);

  return MSG;
}
//____________________________________________________________________________
void Messenger::SetPriorityLevel(
                       const char * stream, log4cpp::Priority::Value priority)
{
  log4cpp::Category & MSG = log4cpp::Category::getInstance(stream);

  MSG.setPriority(priority);
}
//____________________________________________________________________________

