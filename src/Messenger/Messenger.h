//____________________________________________________________________________
/*!

\class    genie::Messenger

\brief    A more convenient interface to the log4cpp Message Service

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 16, 2004
 
*/
//____________________________________________________________________________

#ifndef _MESSENGER_H_
#define _MESSENGER_H_

#include <iostream>
#include <string>
#include <map>

#include "log4cpp/Category.hh"
#include "log4cpp/Appender.hh"
#include "log4cpp/OstreamAppender.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/Priority.hh"

// comment defined priority levels for the document generator
/*! \def pFATAL  \brief Defines the FATAL priority level */
/*! \def pALERT  \brief Defines the ALERT priority level */
/*! \def pCRIT   \brief Defines the ALERT priority level */
/*! \def pERROR  \brief Defines the ALERT priority level */
/*! \def pWARN   \brief Defines the ALERT priority level */
/*! \def pNOTICE \brief Defines the ALERT priority level */
/*! \def pINFO   \brief Defines the ALERT priority level */
/*! \def pDEBUG  \brief Defines the ALERT priority level */

#define pFATAL  log4cpp::Priority::FATAL
#define pALERT  log4cpp::Priority::ALERT
#define pCRIT   log4cpp::Priority::CRIT
#define pERROR  log4cpp::Priority::ERROR
#define pWARN   log4cpp::Priority::WARN
#define pNOTICE log4cpp::Priority::NOTICE
#define pINFO   log4cpp::Priority::INFO
#define pDEBUG  log4cpp::Priority::DEBUG

/*! \def ENDL  \brief A shortcut for log4cpp's CategoryStream::ENDLINE */

#define ENDL log4cpp::CategoryStream::ENDLINE

/*!
  \def   SLOG(stream, priority)
  \brief A macro that returns the requested log4cpp::Category
         appending a short string (using the __FUNCTION__ and __LINE__ macros)
         with information for the calling method.
*/

#define SLOG(stream, priority) \
	   (*Messenger::Instance())(stream) \
               << priority << "<" \
               << __FUNCTION__ << " (" << __LINE__ << ")> : "

/*!
  \def   SLOG(stream, priority)
  \brief A macro that returns the requested log4cpp::Category
         appending a string (using the __PRETTY_FUNCTION__ and __LINE__ macros)
         with information for the calling method.
*/
               
#define LOG(stream, priority) \
	   (*Messenger::Instance())(stream) \
               << priority << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_FATAL(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::FATAL << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_ALERT(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::ALERT << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_CRIT(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::CRIT << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_ERROR(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::ERROR << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_WARN(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::WARN << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_NOTICE(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::NOTICE << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_INFO(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::INFO << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

#define LOG_DEBUG(stream) \
	  (*Messenger::Instance())(stream) \
               << log4cpp::Priority::DEBUG << "<" \
               << __PRETTY_FUNCTION__ << " (" << __LINE__ << ")> : "

namespace genie {

class Messenger 
{
public:

  static Messenger * Instance(void);

  log4cpp::Category & operator () (const char * stream);

  void SetPriorityLevel(const char * stream, 
                                 log4cpp::Priority::Value priority);

private:

  Messenger();
  Messenger(const Messenger & config_pool);
  virtual ~Messenger();

  static Messenger * fInstance;
  
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
