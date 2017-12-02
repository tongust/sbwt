#include "log.h"

#include <iostream>
#include <string>
#include <ctime>
#include <vector>

namespace logger {
using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::time_t;
void LogDebug(const string &msg) {
	time_t result = std::time(nullptr);
	cout << "[Debug]\t"
		<< std::asctime(std::localtime(&result))
		<< msg << endl;
	return;
}
} /* logger */
