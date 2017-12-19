#ifndef SBWT_IO_BUILD_INDEX_H
#define SBWT_IO_BUILD_INDEX_H

#include <string>

#include "sbwt.h"

namespace sbwt{
using std::string;
void WriteIntoDiskBuildIndex(sbwt::BuildIndexRawData&, const string&);

} /* namespace sbwt */


#endif /* SBWT_IO_BUILD_INDEX_H */