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

#include "const.h"
//#include "patternShadow.h"
//#include "verification.h"
namespace sbwt {

struct BuidlIndexRawData{
	BuidlIndexRawData();
	BuidlIndexRawData(char*, size_t, const uint32_t &, const uint32_t&);
	~BuidlIndexRawData();
	char *seq_raw;			/* Reference sequence */
	char *seq_transformed;		/* Transformed reference sequence */
	uint32_t **occurency;		/* Occurency of A/C/G/T */
	uint32_t *suffix_array;		/* Suffix array */
	uint32_t first_coloumn[4];	/* C in the formula, the first column of sbwt matrix*/
	uint32_t length_ref;		/* Length of reference sequence */
	uint32_t num_block_sort;	/* Number of blocks, 4 for 256 */
	uint32_t num_dollar;		/* The number of $s those are appended to the tail of reference sequence */
	uint32_t period;		/* The period for sbwt. 1 is used for normal bwt */
};

void SbwtSort(char*, uint32_t*, uint32_t, uint32_t, uint32_t, const uint32_t&, const uint32_t&);


} /* namespace sbwt */
#endif /* _SBWT_H_ */
