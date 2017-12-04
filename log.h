#ifndef _LOG_H_SBWT
#define _LOG_H_SBWT
#include <string>
#include <iostream>


#ifndef LOGDEBUG
#define LOGDEBUG(msg)  std::cout << __FILE__ << ":" << __LINE__ << ": " << msg << std::endl
#endif /* LOGDEBUG */

namespace logger {

using std::string;

void LogDebug(const string &msg);

} /* logger */

#endif /* _LOG_H_SBWT */

