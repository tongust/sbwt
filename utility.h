#ifndef SBWT_UTILITY_H
#define SBWT_UTILITY_H

#include <stdlib.h>  // for strtol uint32_t etc.

#include "sbwt.h"

namespace utility {
char *Extract(char *, const uint32_t &, uint32_t &);
char *ReadFasta(char *, uint32_t &);
bool IsDNA(char);
bool IsN(char);
uint32_t GetUint(int, char *);
void CountSeedOccurrence(sbwt::BuildIndexRawData &, uint32_t);
void PrintHelp_BuildIndex(int, char**);
void PrintHelp_CountOcc(int, char**);
void PrintHelp_SbwtAligner(int, char**);
void PrintHelp_SbwtTestBitset(int, char **);
void PrintHelp_SbwtTestExactMatch(int, char**);
} /* namespace utility */
#endif /* SBWT_UTILITY_H */
