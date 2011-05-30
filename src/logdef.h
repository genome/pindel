#ifndef __LOGDEF_H__
#define __LOGDEF_H__

#include "logtypes.h"

#define __LOG(s,lvl) \
	std::cout << __FILE__ << "#" << __LINE__ << " [" << __FUNCTION__ << "] " << (lvl) << ": "; \
	s;

#ifndef LOG_LEVEL
#warning "No log level found; logging is disabled"
#define LOG_LEVEL LOG_DISABLED
#endif

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOG_ERROR)
#define LOG_ERR(s) __LOG(s, "error")
#else
#define LOG_ERR(s) // empty
#endif // LOG_LEVEL

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOG_WARN)
#define LOG_WARN(s) __LOG(s, "warning")
#else
#define LOG_WARN(s) // empty
#endif // LOG_LEVEL

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOG_INFO)
#define LOG_INFO(s) __LOG(s, "info")
#else
#define LOG_INFO(s) // empty
#endif // LOG_LEVEL

#if defined(LOG_LEVEL) && (LOG_LEVEL >= LOG_DEBUG)
#define LOG_DBG(s) __LOG(s, "debug")
#else
#define LOG_DBG(s) // empty
#endif // LOG_LEVEL

#endif // __LOGDEF_H__
