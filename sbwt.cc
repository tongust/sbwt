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
#include "log.h"
namespace sbwt {

BuildIndexRawData::BuildIndexRawData():
	seq_raw(nullptr),
	seq_transformed(nullptr),
	occurrence(nullptr),
	suffix_array(nullptr),
	length_ref(0),
	num_block_sort(4),
	num_dollar(2),
	period(2)
{
	for (int i = 0; i < 4; ++i) {
		first_coloumn[i] = 0;
	}
}

BuildIndexRawData::BuildIndexRawData (char *seq_dna, size_t n, const uint32_t &per = 2, const uint32_t &nb = 4):
	seq_raw(seq_dna),
        length_ref(n),
        num_block_sort(nb),
        period(per)
{  

	num_dollar = period - (length_ref % period);
	length_ref += num_dollar;
	seq_raw = (char *)realloc(seq_dna, length_ref*sizeof(char));

	for (uint32_t i = 0; i < num_dollar; ++i) {
		seq_raw[length_ref-i-1] = '$';
	}		 

	for (int i = 0; i < 4; ++i) { first_coloumn[i] = 0; }

	occurrence = new uint32_t*[4];
	for (int i = 0; i != 4; ++i)
	{
		occurrence[i] = new uint32_t[length_ref]();
		occurrence[i][0] = 0;
	}

	seq_transformed = new char[length_ref];

	suffix_array = new uint32_t[length_ref];
	/* initialize suffix_array with 0,1,...,N-1 */
	for (size_t i = 0; i != length_ref; ++i) suffix_array[i] = i;

}
BuildIndexRawData::~BuildIndexRawData()
{

        if (seq_raw) { delete[] seq_raw; }

        if (suffix_array) { delete[] suffix_array; }

        if (occurrence) {
		for (int i = 0; i != 4; ++i) { if (occurrence[i]) delete[] occurrence[i]; }
		delete[] occurrence;
        }

	if (seq_transformed) { delete[] seq_transformed; }
}


/* Swap within array seq_index*/
void VectorSwap(uint32_t i, uint32_t j, uint32_t n, uint32_t* seq_index)
{
        static uint32_t tmpval = 0;
        while (n-- > 0) {
                /* swap */
                tmpval = seq_index[i];
                seq_index[i] = seq_index[j];
                seq_index[j] = tmpval;

                ++i;
                ++j;
        }
        return;
}

void SortSbwt(
	char *seq,
        uint32_t *seq_index,
	uint32_t begin,
        uint32_t end,
        uint32_t depth,
	const uint32_t &length_ref,
        const uint32_t &step
        ) {
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

				++a;
			}
			++b;
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
				
				--d;
			}
			--c;
		}
		if (b > c) break;

		/* swap b and c */
		tmpval = seq_index[b];
		seq_index[b] = seq_index[c];
		seq_index[c] = tmpval;
		++b;
		--c;
	}

        int64_t t0 = a - begin;
        int64_t t1 = b - a;
        r = t0 < t1 ? t0 : t1;
        VectorSwap(begin, b-r, r, seq_index);

        t0 = d - c;
        t1 = end - 1 - d;
        r = t0 < t1 ? t0 : t1;
        VectorSwap(b, end-r, r, seq_index);
        r = b - a + begin;

        SortSbwt(seq, seq_index, begin, r, depth, length_ref, step);

        tmpval = seq_index[r] + depth;
        if (tmpval < length_ref) {
                SortSbwt(seq, seq_index, r, end-d+c, depth+step, length_ref, step);
        }

        r = d - c;
        SortSbwt(seq, seq_index, end-r, end, depth, length_ref, step);
        return;
}
void SortSbwt(BuildIndexRawData &build_index) {
        return SortSbwt(build_index.seq_raw,
                        build_index.suffix_array,
                        0, build_index.length_ref,
                        0, build_index.length_ref,
                        build_index.period);
}

void Transform(BuildIndexRawData &build_index) {
        if (!build_index.seq_transformed) {
                logger::LogDebug("seq_transformed is empty.");
                return;
        }
        uint32_t tmpval = 0;
        auto sa = build_index.suffix_array;
        auto seq = build_index.seq_transformed;
        for (size_t i = 0; i != build_index.length_ref; ++i) {
                tmpval = *sa;
                if (tmpval < build_index.period) {
                        *seq = '$';
                } else {
                        *seq = build_index.seq_raw[tmpval - build_index.period];
                }
                ++sa;
                ++seq;
        }
        return;
}

void CountOccourrence(BuildIndexRawData &build_index) {
        auto N = build_index.length_ref;
        auto B = build_index.seq_transformed;
        auto O = build_index.occurrence;
        auto C = build_index.first_coloumn;
        
        if (!B) return;
        switch (B[0]) {
                case 'A':
                        O[0][0] = 1;
                        ++C[0];
                        break;
                case 'C':
                        O[1][0] = 1;
                        ++C[1];
                        break;
                case 'G':
                        O[2][0] = 1;
                        ++C[2];
                        break;
                case 'T':
                        O[3][0] = 1;
                        break;
                default:
                        break;
        }

        for (size_t i = 1; i < N; ++i) {
        switch (B[i]) {
                case 'A':
                        SNIPPET_COUNTOCC();
                        ++C[0];
                        O[0][i] = O[0][i-1] + 1;
                        break;
                case 'C':
                        SNIPPET_COUNTOCC();
                        ++C[1];
                        O[1][i] = O[1][i-1] + 1;
                        break;
                case 'G':
                        SNIPPET_COUNTOCC();
                        ++C[2];
                        O[2][i] = O[2][i-1] + 1;
                        break;
                case 'T':
                        SNIPPET_COUNTOCC();
                        O[3][i] = O[3][i-1] + 1;
                        break;
                default:
                        SNIPPET_COUNTOCC();
                        break;
                }
        }
        C[3] = C[0] + C[1] + C[2];
        C[2] = C[0] + C[1];
        C[1] = C[0];
        C[0] = 0;

        return;
}

void BuildIndex(BuildIndexRawData &build_index) {
        if (build_index.suffix_array
                        && build_index.seq_raw) {
                SortSbwt(build_index);
                Transform(build_index);
                CountOccourrence(build_index);
        }
}

void PrintFullSearchMatrix(BuildIndexRawData &build_index) {
        using std::endl;
        using std::cout;

        auto SA = build_index.suffix_array;
        auto X = build_index.seq_raw;
        auto N = build_index.length_ref;
        auto B = build_index.seq_transformed;
        auto O = build_index.occurrence;
        auto C = build_index.first_coloumn;

        auto count10 = [](uint32_t i)->uint32_t{
                uint32_t ret = 0;
                for (;;) {
                        i /= 10;
                        if (i==0) break;
                        ++ret;
                }
                return ret;
        };

        if (!X || !SA) return;

        cout << "\nRef: " << endl;

        for (size_t i = 0; i != N; ++i) cout << X[i];

        cout << endl << "i\tSA\tFM\tBWT\n \t ";

        int tcc = 0;
        for (size_t i = 0; i != N; ++i) {
                if (i%5 == 0) {
                        cout << "\t" <<  i;
                        tcc = count10(i);
                }
                else {
                        if (tcc<=0) cout << " ";
                        else --tcc;
                }
        }
        cout << endl;

        for (size_t i = 0; i != N; ++i) {
                cout << i << "\t" << SA[i] <<"\t";
                uint32_t j = 0;
                for (j = 0; j < N-1; ++j) {
                        if (j && j%5 == 0) cout << "\t";
                        cout << X[(j+SA[i])%N];
                }
                j = N - 1;
                if (j%5 == 0)cout << "\t";
                cout << "" << X[(j+SA[i])%N];
                cout << endl;
        }
        cout << "\nspaced BWT:\n";
        for (uint32_t i = 0; i != N-1; ++i) {
                cout << B[i] << ",";
        }
        cout << B[N-1] << "\n";
        cout << "Occ\nA\tC\tG\tT\n";
        for (int i = 0; i != 4; ++i)
                cout << C[i] << "\t";
        cout << "\nOcc\nindex\tA\tC\tG\tT\n";
        for (uint32_t i = 0; i != N; ++i) {
                cout << i << "\t";
                for (int j = 0; j != 4; ++j)
                        cout << O[j][i] << "\t";
        cout << "\n";
        }
        return;
}




} /* namespace sbwt */
