#ifndef _LOG_H_SBWT
#define _LOG_H_SBWT
#include <string>
#include <iostream>


#ifndef LOGDEBUG
#define LOGDEBUG(msg) std::cout << __FILE__ << ":" << __LINE__ << ":[Debug]\t" << msg << std::endl
#endif /* LOGDEBUG */

#define LOGERROR(msg) std::cerr << __FILE__ << ":" << __LINE__ << ":[Error]\t" << msg << std::endl;

#define LOGINFO(msg) std::cerr << "[Info]\t" << msg;

#define LOGPUT(msg) std::cerr << msg;

namespace logger {

using std::string;

void LogDebug(const string &);
void LogError(const string &);

} /* logger */

#endif /* _LOG_H_SBWT */

