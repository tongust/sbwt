#ifndef _SBWT_H_
#define _SBWT_H_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <time.h>
#include <tuple>
#include <vector>
#include <string>

#include "const.h"
//#include "patternShadow.h"
//#include "verification.h"
namespace sbwt {

using std::string;

struct BuildIndexRawData{
	BuildIndexRawData();
	BuildIndexRawData(char*, size_t, const uint32_t &, const uint32_t&);
        /// Build index from reading files
        BuildIndexRawData(const string&);
	~BuildIndexRawData();
	char *seq_raw;			/* Reference sequence */
	char *seq_transformed;		/* Transformed reference sequence */
	uint32_t **occurrence;		/* occurrence of A/C/G/T */
	uint32_t *suffix_array;		/* Suffix array */
	uint32_t first_column[4];	/* C in the formula, the first column of sbwt matrix*/

	uint32_t length_ref;		/* Length of reference sequence including $s*/
	uint32_t num_block_sort;	/* Number of blocks, 4 for 256 */
	uint32_t num_dollar;		/* The number of $s those are appended to the tail of reference sequence */
	uint32_t period;		/* The period for sbwt. 1 is used for normal bwt */
};

void SortSbwt(char*, uint32_t*, uint32_t, uint32_t, uint32_t, const uint32_t&, const uint32_t&);
inline void VectorSwap(uint32_t, uint32_t, uint32_t, uint32_t*);
void BuildIndex(BuildIndexRawData&);
void BuildIndexBlockwise(BuildIndexRawData&);
void CountOccurrence(BuildIndexRawData&);
void PrintFullSearchMatrix(BuildIndexRawData&);
void SortSbwt(BuildIndexRawData&);
void SortSbwtBlockwise( char*, uint32_t*, uint32_t, uint32_t, uint32_t, const uint32_t&, const uint32_t&, const uint32_t&);
void SortSbwtBlockwise(BuildIndexRawData&);
void Transform(BuildIndexRawData&);

#ifndef SNIPPET_COUNTOCC
#define SNIPPET_COUNTOCC()  \
        O[0][i] = O[0][i-1];\
        O[1][i] = O[1][i-1];\
        O[2][i] = O[2][i-1];\
        O[3][i] = O[3][i-1]
#endif /* SNIPPET_COUNTOCC */

} /* namespace sbwt */
#endif /* _SBWT_H_ */
