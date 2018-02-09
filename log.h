#ifndef _LOG_H_SBWT
#define _LOG_H_SBWT
#include <string>
#include <iostream>


#ifndef LOGDEBUG
#define LOGDEBUG(msg) std::cout << __FILE__ << ":" << __LINE__ << ": " << msg << std::endl
#endif /* LOGDEBUG */

#define LOGERROR(msg) std::cerr << __FILE__ << ":" << __LINE__ << ":\033[0;31m\t" << msg << "\033[0m" << std::endl;

#define LOGEINFO(msg) std::cerr << "[Info]\t\033[0;33m" << msg << "\033[0m" << std::endl;

namespace logger {

using std::string;

void LogDebug(const string &);
void LogError(const string &);

} /* logger */

#endif /* _LOG_H_SBWT */

