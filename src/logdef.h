#ifndef __LOGDEF_H__
#define __LOGDEF_H__

#include "logtypes.h"

#define __LOG(s,lvl) (s);
/* original: std::cout << __FILE__ << "#" << __LINE__ << " [" << __FUNCTION__ << "] " << (lvl) << ": " << s; */

#ifndef LOG_LEVEL
#define LOG_LEVEL LOGLEVEL_INFO
#endif

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOGLEVEL_ERROR)
#define LOG_ERROR(s) __LOG(s, "error")
#else
#define LOG_ERROR(s) // empty
#endif // LOG_LEVEL

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOGLEVEL_WARN)
#define LOG_WARN(s) __LOG(s, "warning")
#else
#define LOG_WARN(s) // empty
#endif // LOG_LEVEL

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOGLEVEL_INFO)
#define LOG_INFO(s) __LOG(s, "info")
#else
#define LOG_INFO(s) // empty
#endif // LOG_LEVEL

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOGLEVEL_DEBUG)
#define LOG_DEBUG(s) __LOG(s, "debug")
#else
#define LOG_DEBUG(s) // empty
#endif // LOG_LEVEL

#endif // __LOGDEF_H__
