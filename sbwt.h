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
	/// Build index block-wise (with read reference function)
        BuildIndexRawData(char *, const uint32_t &/*period*/, const uint32_t&/*# of blocks*/);
        /// Build index block-wise
	BuildIndexRawData(char*, size_t, const uint32_t &, const uint32_t&);
        /// Build index from index files
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

        uint8_t *bin_8bit;              /* 8-bit-packed binary sequence */
        uint32_t size_bin_8bit;         /* Length of packed binary sequence */
};


class SecondIndex{
public:
	SecondIndex();
	/// build second index from files stream
        SecondIndex(const string &);

        /// array consists of header and sequence
	///     1. header: the position in reference;
        ///     2. size of single segment;
	///     3. sequence: sorted index (position in SSA).
	uint16_t *array_ptr;
	/// size of array_ptr
	uint64_t size;
        /// min size of array to sort
        static uint32_t size_min;
        /// size of seed
        uint32_t size_seed;

	~SecondIndex();

public:
        void RebuildIndexInit(BuildIndexRawData&, uint32_t);
        void RebuildIndex(BuildIndexRawData&);
	void PrintSecondIndex(BuildIndexRawData&);
        bool Empty();
};


void SortSbwt(char*, uint32_t*, uint32_t, uint32_t, uint32_t, const uint32_t&, const uint32_t&);
void SortSbwt(char*, uint32_t*, uint32_t, uint32_t, uint32_t, const uint32_t&, const uint32_t&, const uint32_t&/*limited length*/);
inline void VectorSwap(uint32_t, uint32_t, uint32_t, uint32_t*);
void BuildIndex(BuildIndexRawData&);
void BuildIndexBlockwise(BuildIndexRawData&);
void BuildSortedIndexBlockwise(BuildIndexRawData&);
void BuildSortedIndexTransCountOcc(BuildIndexRawData&);
void CountOccurrence(BuildIndexRawData&);
void PrintFullSearchMatrix(BuildIndexRawData&);
void PrintFullSearchMatrix(uint32_t *SA, char *X, uint32_t N, uint32_t period);
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
