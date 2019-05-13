//____________________________________________________________________________
/*!

\class    genie::Messenger

\brief    A more convenient interface to the log4cpp Message Service

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  June 16, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MESSENGER_H_
#define _MESSENGER_H_

#include <iostream>
#include <cstring>
#include <string>
#include <map>

// ROOT5 has difficulty with parsing log4cpp headers
#if !defined(__CINT__) && !defined(__MAKECINT__)
  #include "log4cpp/Category.hh"
  #include "log4cpp/Appender.hh"
  #include "log4cpp/OstreamAppender.hh"
  #include "log4cpp/BasicLayout.hh"
  #include "log4cpp/Priority.hh"
#else
  namespace log4cpp {
    class Priority {
      typedef int Value;
    };
    class Category;
  }
#endif

#include "Framework/Conventions/GBuild.h"

using std::string;

// comment defined priority levels for the document generator
/*! \def pFATAL  \brief Defines the FATAL  priority level */
/*! \def pALERT  \brief Defines the ALERT  priority level */
/*! \def pCRIT   \brief Defines the CRIT   priority level */
/*! \def pERROR  \brief Defines the ERROR  priority level */
/*! \def pWARN   \brief Defines the WARN   priority level */
/*! \def pNOTICE \brief Defines the NOTICE priority level */
/*! \def pINFO   \brief Defines the INFO   priority level */
/*! \def pDEBUG  \brief Defines the DEBUG  priority level */

#define pFATAL  log4cpp::Priority::FATAL
#define pALERT  log4cpp::Priority::ALERT
#define pCRIT   log4cpp::Priority::CRIT
#define pERROR  log4cpp::Priority::ERROR
#define pWARN   log4cpp::Priority::WARN
#define pNOTICE log4cpp::Priority::NOTICE
#define pINFO   log4cpp::Priority::INFO
#define pDEBUG  log4cpp::Priority::DEBUG

/*! \def ENDL  \brief A shortcut for log4cpp's CategoryStream::ENDLINE or std manipulators*/

#ifdef __GENIE_USES_LOG4CPP_VERSION__
  #if __GENIE_USES_LOG4CPP_VERSION__==0
    #define ENDL log4cpp::CategoryStream::ENDLINE
  #else
    #define ENDL std::endl
  #endif
#else
  #define ENDL std::endl
#endif

/*!
  \def   SLOG(stream, priority)
  \brief A macro that returns the requested log4cpp::Category
         appending a short string (using the __FUNCTION__ and __LINE__ macros)
         with information for the calling method [produces short message].
*/

#define SLOG(stream, priority) \
           (*Messenger::Instance())(stream) \
               << priority << "[s] <" \
               << __FUNCTION__ << " (" << __LINE__ << ")> : "

/*!
  \def   LOG(stream, priority)
  \brief A macro that returns the requested log4cpp::Category
         appending a string (using the __FILE__, __FUNCTION__ and __LINE__ macros)
         with information for the calling method [produces normal messages].
*/

#define LOG(stream, priority) \
           (*Messenger::Instance())(stream) \
               << priority << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

/*!
  \def   HIDE_GENIE_MSG_LOG_MACROS
  \brief Use this cpp flag variable to ensure LOG_DEBUG ... LOG_FATAL macros are not exposed.
         This allows this header to be used in conjunction with the art framework's
         conflicting MessengeFacility's macros of the same name.  The two argument
         LOG macro (see above) is still available for use.
         Currently this comes up only via Algorithm.h's inclusion of Algorithm.icc
         which included Messenger.h.
*/
#ifndef HIDE_GENIE_MSG_LOG_MACROS

#define LOG_FATAL(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::FATAL << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_ALERT(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::ALERT << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_CRIT(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::CRIT << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_ERROR(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::ERROR << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_WARN(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::WARN << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_NOTICE(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::NOTICE << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_INFO(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::INFO << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_DEBUG(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::DEBUG << "[n] <" \
               << __FILE__ << "::" << __FUNCTION__ << " (" << __LINE__ << ")> : "

#endif // HIDE_GENIE_MSG_LOG_MACROS

/*!
  \def   LLOG(stream, priority)
  \brief A macro that returns the requested log4cpp::Category
         appending a string (using the __PRETTY_FUNCTION__ and __LINE__ macros)
         with information for the calling method [produces long messages].
*/

#define LLOG(stream, priority) \
           (*Messenger::Instance())(stream) \
               << priority << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_FATAL(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::FATAL << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_ALERT(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::ALERT << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_CRIT(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::CRIT << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_ERROR(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::ERROR << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_WARN(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::WARN << "'[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_NOTICE(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::NOTICE << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_INFO(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::INFO << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LLOG_DEBUG(stream) \
          (*Messenger::Instance())(stream) \
               << log4cpp::Priority::DEBUG << "[l] <" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

/*!
  \def   BLOG(stream, priority)
  \brief A macro that returns the requested log4cpp::Category appending no
         additional information
*/

#define BLOG(stream, priority) \
          (*Messenger::Instance())(stream) << priority

/*!
  \def   MAXSLOG(stream, priority, maxcount)
  \brief Similar to SLOG(stream,priority) but quits after "maxcount" messages

  \def   MAXLOG(stream, priority, maxcount)
  \brief Similar to LOG(stream,priority) but quits after "maxcount" messages

  \def   MAXLLOG(stream, priority, maxcount)
  \brief Similar to LLOG(stream,priority) but quits after "maxcount" messages


*/

// Macro to concatenate two symbols:
#define TOKCAT(x,y) x##y
// Macro to expand preprocessor variables and concatenate:
#define TOKCAT2(x,y) TOKCAT(x,y)
// Macro to concatenate source line with a symbol:
#define LINECAT(x) TOKCAT2(x, __LINE__ )

#define MAXSLOG(s,l,c)  \
  static int  LINECAT(MSGCNT) = 0; \
  const char* LINECAT(MSGADD) = (++LINECAT(MSGCNT)==c) ? "..Last Message .. " : ""; \
  if (LINECAT(MSGCNT) > c || LINECAT(MSGCNT) < 0) \
     {;} else SLOG(s,l) << LINECAT(MSGADD)

#define MAXLOG(s,l,c)  \
  static int  LINECAT(MSGCNT) = 0; \
  const char* LINECAT(MSGADD) = (++LINECAT(MSGCNT)==c) ? "..Last Message .. " : ""; \
  if (LINECAT(MSGCNT) > c || LINECAT(MSGCNT) < 0) \
     {;} else LOG(s,l) << LINECAT(MSGADD)

#define MAXLLOG(s,l,c)  \
  static int  LINECAT(MSGCNT) = 0; \
  const char* LINECAT(MSGADD) = (++LINECAT(MSGCNT)==c) ? "..Last Message .. " : ""; \
  if (LINECAT(MSGCNT) > c || LINECAT(MSGCNT) < 0) \
     {;} else LLOG(s,l) << LINECAT(MSGADD)


namespace genie {

extern bool gAbortingInErr;

class Messenger
{
public:
  static Messenger * Instance(void);

  log4cpp::Category & operator () (const char * stream);
  void SetPriorityLevel(const char * stream, log4cpp::Priority::Value p);

  bool SetPrioritiesFromXmlFile(string filename);

private:
  Messenger();
  Messenger(const Messenger & config_pool);
  virtual ~Messenger();

  static Messenger * fInstance;

  void Configure(void);

  log4cpp::Priority::Value PriorityFromString(string priority);

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (Messenger::fInstance !=0) {
            delete Messenger::fInstance;
            Messenger::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;
};

}      // genie namespace
#endif // _MESSENGER_H_
