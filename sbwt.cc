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
#include "sbwt.h"
namespace sbwt {

BuidlIndexRawData::BuidlIndexRawData():
	length_ref(0),
	num_block_sort(4),
	num_dollar(2),
	occurency(nullptr),
	period(2),
	seq_raw(nullptr),
	seq_transformed(nullptr),
	suffix_array(nullptr)
{
	for (int i = 0; i < 4; ++i) {
		first_coloumn[i] = 0;
	}
}

BuidlIndexRawData::BuidlIndexRawData (char *seq_dna, size_t n, const uint32_t &per = 2, const uint32_t &nb = 4):
	seq_raw(seq_dna), length_ref(n), period(per), num_block_sort(nb) {  

	num_dollar = period - (length_ref % period);
	length_ref += num_dollar;
	seq_raw = (char *)realloc(seq_dna, length_ref*sizeof(char));

	for (uint32_t i = 0; i < num_dollar; ++i) {
		seq_raw[length_ref-i-1] = '$';
	}		 

	for (int i = 0; i < 4; ++i) { first_coloumn[i] = 0; }

	occurency = new uint32_t*[4];
	for (int i = 0; i != 4; ++i)
	{
		occurency[i] = new uint32_t[length_ref]();
		occurency[i][0] = 0;
	}

	seq_transformed = new char[length_ref];

	suffix_array = new uint32_t[length_ref];
	/* initialize suffix_array with 0,1,...,N-1 */
	for (size_t i = 0; i != length_ref; ++i) suffix_array[i] = i;

}
BuidlIndexRawData::~BuidlIndexRawData()
{

        if (seq_raw) { delete[] seq_raw; }

        if (suffix_array) { delete[] suffix_array; }

        if (occurency) {
		for (int i = 0; i != 4; ++i) { if (occurency[i]) delete[] occurency[i]; }
		delete[] occurency;
        }

	if (seq_transformed) { delete[] seq_transformed; }
}


void SbwtSort(
	char *seq, uint32_t *seq_index,
	uint32_t begin, uint32_t end, uint32_t depth,
	const uint32_t &length_ref, const uint32_t &step) {
	/* Condition of end */
	if (begin+1 >= end || depth >= length_ref) return;

	int64_t a = 0, b = 0, c = 0,
		d = 0, r = 0, v = 0,
		distance = 0;
	uint32_t tmpval = 0, tmpval1 = 0;
	distance = end - begin;
	a = (rand() % distance) + begin;

	/* swap begin with a */
	tmpval = seq_index[begin];
	seq_index[begin] = seq_index[a];
	seq_index[a] = tmpval;

	tmpval1 = seq_index[begin] + depth;
	/* Guarantee that: seq[index] or '$' if index > length(seq) */
	v = tmpval1 >= length_ref ? '$' : seq[tmpval1];
	a = b = begin + 1;
	c = d = end - 1;

	for (;;) {
		while ( b <= c &&
			(tmpval1 = seq_index[b]+depth, r = (tmpval1 >= length_ref ? '$' : seq[tmpval1]) - v) <= 0
			) {
			if (r == 0) {
				/* swap a with b */
				tmpval = seq_index[a];
				seq_index[a] = seq_index[b];
				seq_index[b] = tmpval;

				a++;
			}
			b++;
		}

		while (
			b <= c &&
			(tmpval1 = seq_index[c]+depth, r = (tmpval1 >= length_ref ? '$' : seq[tmpval1]) - v) >= 0
			) {
			if (r == 0) {
				/* swap c with d */
				tmpval = seq_index[c];
				seq_index[c] = seq_index[d];
				seq_index[d] = tmpval;
				
				d--;
			}
			c--;
		}
		if (b > c) break;

		/* swap b and c */
		tmpval = seq_index[b];
		seq_index[b] = seq_index[c];
		seq_index[c] = tmpval;
		b++;
		c--;
	}

	
	






	return;
}

} /* namespace sbwt */
