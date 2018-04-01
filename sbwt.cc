#include <stdlib.h>
#include <stdint.h>

#include <algorithm>
#include <bitset>
#include <chrono>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <iomanip>
#include <vector>
#include <string>
#include <unordered_map>

#include "sbwt.h"
#include "log.h"
#include "alphabet.h"
#include "ref_read.h"
#include "word_io.h"
#include "sequence_pack.h"
#include "utility.h"

namespace sbwt {
using std::vector;
using std::unordered_map;
using namespace utility;

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
		first_column[i] = 0;
	}
}


BuildIndexRawData::BuildIndexRawData(char *file_name, const uint32_t &per, const uint32_t &nb):
        seq_raw(nullptr),
        bin_8bit(nullptr),
        size_bin_8bit(0)
{
        period = per;
        num_block_sort = nb;

        uint32_t read_length = 0;

        /// Read raw sequence file in fasta format
        char *buffer = nullptr;
        {
                uint32_t string_length;
		FILE *handler = fopen(file_name, "r");

		if (!handler) {
                        logger::LogError("failure to open");
		}
                fseek(handler, 0, SEEK_END);    /* Seek the last byte */
                string_length = ftell(handler); /* Calculate the offset from head to tail */
                rewind(handler);                /* Back to the head */

                //buffer = (char *) malloc(sizeof(char) * (string_length + 1));
                buffer = new char[sizeof(char) * (string_length + 1)];
                read_length = fread(buffer, sizeof(char), string_length, handler);

                buffer[string_length] = '\0';

                if (read_length != string_length) {
                        delete []buffer;
                        logger::LogDebug("Errors in fread");
                        return;
                }
        }

        /// Extract pure DNA nucleotides (A/C/G/T)
        {
                uint32_t total_num = 0;
                char *ptr = buffer;
                for (uint32_t i = 0; i < read_length; ++i) {
                        if (IsDNA(*ptr)) {
                                ++total_num;
                        }
                        ++ptr;
                }

                /// enough memory for $s and "boarder case"
                /// And the period must be less than 1024
                /// In case the boarder of sequence array will be reached.
                uint32_t n_alloc = ((total_num/1024)+4)*1024;
                seq_raw = new char[n_alloc]();
                //for (uint32_t i = 0; i != n_alloc; ++i) { seq_raw[i] = 0; }

                char *ptr0 = seq_raw;
                ptr = buffer;
                for (uint32_t i = 0; i < read_length; ++i) {
                        if (IsDNA(*ptr)) {
                                *ptr0 = *ptr;
                                ++ptr0;
                        }
                        ++ptr;
                }

                length_ref = total_num;
        }

        /// Init
        {
                if (period > 1024) {
                        period = 1024;
                        logger::LogError("The period must be less than 1024.");
                }

		num_dollar = period - (length_ref % period);
                if (length_ref % period) {
                        num_dollar += period;
                }
		length_ref += num_dollar;

		for (uint32_t i = 0; i < num_dollar; ++i) {
			seq_raw[length_ref-i-1] = '$';
		}

		for (int i = 0; i < 4; ++i) { first_column[i] = 0; }

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

        delete[] buffer;
}


BuildIndexRawData::BuildIndexRawData (char *seq_dna, size_t n, const uint32_t &per = 2, const uint32_t &nb = 4):
	seq_raw(seq_dna),
        length_ref(n),
        num_block_sort(nb),
        period(per),
        bin_8bit(nullptr),
        size_bin_8bit(0)
{

        if (period > 1024) {
                period = 1024;
                logger::LogError("The period must be less than 1024.");
        }

	num_dollar = period - (length_ref % period);
        if (length_ref % period) {
                num_dollar += period;
        }
	length_ref += num_dollar;
        /// TODO Bad idea to use realloc()
        /// How to free memory created by realloc(), free() or delete[]?
        /// Suggestion: use alloc and free
        /// Update: no need to use malloc, there is enough space for $s, however the upper-bound of period is
        /// set to 1024.
	/// seq_raw = (char *)realloc(seq_dna, length_ref*sizeof(char));

	for (uint32_t i = 0; i < num_dollar; ++i) {
		seq_raw[length_ref-i-1] = '$';
	}

	for (int i = 0; i < 4; ++i) { first_column[i] = 0; }

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

BuildIndexRawData::BuildIndexRawData(const string &prefix_filename)
{

        string file_array_filename = prefix_filename + ".array.sbwt";
        string file_meta_filename = prefix_filename + ".meta.sbwt";

        std::ifstream array_fin(file_array_filename.c_str(), std::ios_base::in | ios::binary);
        std::ifstream meta_fin(file_meta_filename.c_str(), std::ios_base::in | ios::binary);

        if ((!array_fin.is_open()) || (!meta_fin.is_open())) {
                std::cerr << "Cannot open index files: " << prefix_filename << std::endl;
                return;
        }

        uint32_t tmp32 = readU32(meta_fin, true);
        bool is_big_endian = tmp32 != 0;

        /**
         * Meta information
         */
        length_ref          = readU32(meta_fin, is_big_endian);/* Length of reference sequence including $s*/
        num_block_sort      = readU32(meta_fin, is_big_endian);/* Number of blocks, 4 for 256 */
        num_dollar          = readU32(meta_fin, is_big_endian);/* The # of $s those are appended */
        period              = readU32(meta_fin, is_big_endian);/* The period of sbwt */
        size_bin_8bit       = readU32(meta_fin, is_big_endian);/* The size of packed sequence */

        /// read first column
        for (int i = 0; i != 4; ++i) {
                first_column[i] = readU32(meta_fin, is_big_endian);
        }

        /**
         * Array and sequence
         */
/// 64-bit packed sequence version
#if 0
        {
                std::shared_ptr<uint64_t > binary_seq64_ptr(new uint64_t[size_packed_seq]);

                uint64_t *beg = binary_seq64_ptr.get();
                uint64_t *end = beg + (size_packed_seq-1);
                for (;beg != end; ++beg) {
                        *beg = readU64(array_fin, is_big_endian);
                }
        }
#endif

        /// 8-bit version
        /// Watch out for the boarder
        bin_8bit = new uint8_t[(size_bin_8bit * 4) + 1024];
        array_fin.read((char*) bin_8bit, size_bin_8bit*4);

        /// The raw sequence
        /// TODO map directly the memory to files
        seq_raw = new char[length_ref];
        array_fin.read(seq_raw, length_ref);

        seq_transformed = nullptr;
        //seq_transformed = new char[length_ref];
        //array_fin.read(seq_transformed, length_ref);

        /// Occurrence
        occurrence = new uint32_t*[4];
        for (int i = 0; i != 4; ++i) {
                occurrence[i] = new uint32_t[length_ref];
                uint32_t *beg = occurrence[i];
                uint32_t *end = beg+length_ref;
                while (beg != end) {
                        *beg = readU32(array_fin, is_big_endian);
                        ++beg;
                }
        }

        /// suffix array
        suffix_array = new uint32_t[length_ref];
        uint32_t *beg = suffix_array;
        uint32_t *end = suffix_array + length_ref;
        while (beg != end) {
                *beg = readU32(array_fin, is_big_endian);
                ++beg;
        }

        array_fin.close();
        meta_fin.close();

}


BuildIndexRawData::~BuildIndexRawData()
{

        delete[] seq_raw;

        delete[] suffix_array;

        if (occurrence) {
		for (int i = 0; i != 4; ++i) { if (occurrence[i]) delete[] occurrence[i]; }
		delete[] occurrence;
        }

	if (seq_transformed) {
                delete[] seq_transformed;
        }

        if (bin_8bit) {
                delete[] bin_8bit;
        }
}


/* Swap within array seq_index*/
void VectorSwap(uint32_t i, uint32_t j, uint32_t n, uint32_t* seq_index)
{
        uint32_t tmpval = 0;
        while (n-- > 0) {
                /* swap */
                tmpval = seq_index[i];
                seq_index[i] = seq_index[j];
                seq_index[j] = tmpval;

                ++i;
                ++j;
        }
}

void SortSbwt(
	char *seq,
        uint32_t *seq_index,
	uint32_t begin,
        uint32_t end,
        uint32_t depth,
	const uint32_t &length_ref,
        const uint32_t &step
        )
{
	/* Condition of end */
	if (begin+1 >= end || depth >= length_ref) return;

        if (end - begin == 2) {
                bool smaller = true;
                uint32_t i1 = seq_index[begin],
                         i2 = seq_index[begin+1];
                for (uint32_t i = depth; i < length_ref; i+=step) {
                        char c1 = i1+i >= length_ref ? '$' : seq[i1+i];
                        char c2 = i2+i >= length_ref ? '$' : seq[i2+i];
                        if (c1 == '$' || c2 == '$') {
                                if (c1 != '$') {
                                        smaller = false;
                                }
                                break;
                        } else if (c1 != c2) {
                                smaller = c1 < c2;
                                break;
                        }
                }
                if (!smaller) {
                        seq_index[begin] = i2;
                        seq_index[begin+1] = i1;
                }
                return;
        }

	int64_t a = 0, b = 0, c = 0,
		d = 0, r = 0, v = 0,
		distance = 0;
	uint64_t tmpval = 0, tmpval1 = 0;
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

        uint64_t size_ref64 = length_ref;

	for (;;) {
#if 1
                while (true) {
                        if (b <= c)  {
                                tmpval1 = seq_index[b];
                                tmpval1 += depth;
                                r = (tmpval1 >= size_ref64) ? '$' : seq[tmpval1];
                                r -= v;
                                if (r > 0) {
                                        break;
                                }
                        } else {
                                break;
                        }

			if (r == 0) {
				/* swap a with b */
				tmpval = seq_index[a];
				seq_index[a] = seq_index[b];
				seq_index[b] = tmpval;

				++a;
			}
			++b;
		}

                while(true) {
                        if (b <= c) {
                                tmpval1 = seq_index[c];
                                tmpval1 += depth;
                                r = (tmpval1 >= size_ref64) ? '$' : seq[tmpval1];
                                r -= v;
                                if (r < 0) {
                                        break;
                                }
                        } else {
                                break;
                        }

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
#endif

#if 0
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
#endif
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
}

        void SortSbwt(
                        char *seq,
                        uint32_t *seq_index,
                        uint32_t begin,
                        uint32_t end,
                        uint32_t depth,
                        const uint32_t &length_ref,
                        const uint32_t &step,
                        const uint32_t &length_limited
        )
        {
                /* Condition of end */
                if (begin+1 >= end || depth >= length_ref || depth >= length_limited) return;

                if (end - begin == 2) {
                        bool smaller = true;
                        uint32_t i1 = seq_index[begin],
                                        i2 = seq_index[begin+1];
                        for (uint32_t i = depth; i < length_ref; i+=step) {
                                char c1 = i1+i >= length_ref ? '$' : seq[i1+i];
                                char c2 = i2+i >= length_ref ? '$' : seq[i2+i];
                                if (c1 == '$' || c2 == '$') {
                                        if (c1 != '$') {
                                                smaller = false;
                                        }
                                        break;
                                } else if (c1 != c2) {
                                        smaller = c1 < c2;
                                        break;
                                }
                        }
                        if (!smaller) {
                                seq_index[begin] = i2;
                                seq_index[begin+1] = i1;
                        }
                        return;
                }

                int64_t a = 0, b = 0, c = 0,
                                d = 0, r = 0, v = 0,
                                distance = 0;
                uint64_t tmpval = 0, tmpval1 = 0;
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

                uint64_t size_ref64 = length_ref;

                for (;;) {
#if 1
                        while (true) {
                                if (b <= c)  {
                                        tmpval1 = seq_index[b];
                                        tmpval1 += depth;
                                        r = (tmpval1 >= size_ref64) ? '$' : seq[tmpval1];
                                        r -= v;
                                        if (r > 0) {
                                                break;
                                        }
                                } else {
                                        break;
                                }

                                if (r == 0) {
                                        /* swap a with b */
                                        tmpval = seq_index[a];
                                        seq_index[a] = seq_index[b];
                                        seq_index[b] = tmpval;

                                        ++a;
                                }
                                ++b;
                        }

                        while(true) {
                                if (b <= c) {
                                        tmpval1 = seq_index[c];
                                        tmpval1 += depth;
                                        r = (tmpval1 >= size_ref64) ? '$' : seq[tmpval1];
                                        r -= v;
                                        if (r < 0) {
                                                break;
                                        }
                                } else {
                                        break;
                                }

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
#endif

#if 0
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
#endif
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
        }
void SortSbwt(BuildIndexRawData &build_index)
{
        return SortSbwt(build_index.seq_raw,
                        build_index.suffix_array,
                        0, build_index.length_ref,
                        0, build_index.length_ref,
                        build_index.period);
}

void Transform(BuildIndexRawData &build_index)
{
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
                        /// Bug here. Feb-02
                        //*seq = '$';
                        *seq = build_index.seq_raw[tmpval + build_index.length_ref - build_index.period];
                } else {
                        *seq = build_index.seq_raw[tmpval - build_index.period];
                }
                ++sa;
                ++seq;
        }
}

void CountOccurrence(BuildIndexRawData &build_index)
{
        auto N = build_index.length_ref;
        auto B = build_index.seq_transformed;
        auto O = build_index.occurrence;
        auto C = &build_index.first_column[0];
        
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
        C[0] += build_index.num_dollar;
        C[3] = C[0] + C[1] + C[2];
        C[2] = C[0] + C[1];
        C[1] = C[0];
        C[0] = build_index.num_dollar;

}

void BuildIndex(BuildIndexRawData &build_index) {
        if (build_index.suffix_array && build_index.seq_raw) {
                SortSbwt(build_index);
                Transform(build_index);
                CountOccurrence(build_index);
        }
}


void BuildSortedIndexBlockwise(BuildIndexRawData &build_index)
{
         if (build_index.suffix_array && build_index.seq_raw) {
                LOGINFO("Sort sbwt block-wise...\n")
                SortSbwtBlockwise(build_index);
                LOGINFO("SortSbwtBlockwise done\n");
        }
}
void BuildSortedIndexTransCountOcc(BuildIndexRawData &build_index)
{         if (build_index.suffix_array && build_index.seq_raw) {
                LOGINFO("Transform...\t");
                Transform(build_index);
                LOGPUT("Done\n");
                LOGINFO("CountOccurrence...\t");
                CountOccurrence(build_index);
                LOGPUT("Done\n");
        }

}


void BuildIndexBlockwise(BuildIndexRawData &build_index) {
        if (build_index.suffix_array && build_index.seq_raw) {
                LOGINFO("Sort sbwt block-wise...\n")
                SortSbwtBlockwise(build_index);
                LOGINFO("SortSbwtBlockwise done\n");
                LOGINFO("Transform...\t");
                Transform(build_index);
                LOGPUT("Done\n");
                LOGINFO("CountOccurrence...\t");
                CountOccurrence(build_index);
                LOGPUT("Done\n");
        }
}

void PrintFullSearchMatrix(BuildIndexRawData &build_index)
{
        using std::endl;
        using std::cout;
        using std::bitset;

        auto SA = build_index.suffix_array;
        auto X = build_index.seq_raw;
        auto N = build_index.length_ref;
        auto B = build_index.seq_transformed;
        auto O = build_index.occurrence;
        auto C = &build_index.first_column[0];
        auto T = build_index.period;

        if (!X || !SA) return;

        cout << "\033[0;33mLength of reference:\033[0m\t"
             << build_index.length_ref
             << endl
             << "\033[0;33mPeriod:\033[0m\t"
             << build_index.period
             << endl;

        auto count10 = [](uint32_t i)->uint32_t{
                uint32_t ret = 0;
                for (;;) {
                        i /= 10;
                        if (i==0) break;
                        ++ret;
                }
                return ret;
        };

        cout << "\033[0;33mRef:\033[0m" << endl;
        for (size_t i = 0; i != N; ++i) {
                cout << X[i];
                if (i % 4 == 3 && i + 1 != N) cout << "     ";
        } cout << endl;

        if (build_index.bin_8bit) {
                for (int k = 0; k != 4; ++k) {
                        uint32_t lmd_i = (k+1)*build_index.size_bin_8bit;
                        for (uint32_t i = k*build_index.size_bin_8bit; i != lmd_i; ++i) {
                                cout << bitset<8>(build_index.bin_8bit[i]) << " ";
                        } cout << endl;
                }
        }

        cout << endl << "\033[0;33mi\tSA\tFM\n\033[0m";
        const int print_width = 10;
        cout << "\t\t";
        std::cout.fill(' ');
        for (uint32_t i = 0; i != N; ++i) {
                if (i % print_width == 0) {
                        std::cout.width(print_width);
                        cout << std::left << i;
                }
        } cout << endl;

        for (size_t i = 0; i != N; ++i) {
                cout << i << "\t" << SA[i] <<"\t";
                uint32_t j = 0;
                for (j = 0; j <= N-1; ++j) {
                        //if (j && j%print_width == 0) cout << "\t";
                        if (j + SA[i] >= N) {
                                cout << "\033[0;30m" << X[(j+SA[i])%N] << "\033[0m";
                                continue;
                        }
                        if (j % T == 0) {
                                cout << "\033[1;31m" << X[(j+SA[i])%N] << "\033[0m";
                        }
                        else {
                                cout << "\033[0;40m" << X[(j+SA[i])%N] << "\033[0m";
                        }
                } cout << endl;
        }
        cout << "\nspaced BWT:\n";
        if (B!= nullptr) {
                for (uint32_t i = 0; i != N-1; ++i) { cout << B[i] << ","; }
                cout << B[N-1] << "\n";
        }

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
}


void PrintFullSearchMatrix(uint32_t *SA, char *X, uint32_t N, uint32_t period)
{
        using std::endl;
        using std::cout;
        using std::bitset;

        if (!X || !SA) return;

        cout << "Length of reference:\t"
             << N
             << endl
             << "Period:\t"
             << period
             << endl;

        auto count10 = [](uint32_t i)->uint32_t{
                uint32_t ret = 0;
                for (;;) {
                        i /= 10;
                        if (i==0) break;
                        ++ret;
                }
                return ret;
        };

        cout << "\nRef: " << endl;
        for (size_t i = 0; i != N; ++i) {
                cout << X[i];
                if (i % 4 == 3 && i + 1 != N) cout << "     ";
        } cout << endl;

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
}


/**Build sbwt index blockwise for large genomes such homo, however the max length
 * should be less than 4G. Otherwise, you should split the reference into muliple
 * partation.
 */
/// TODO The way to build index of huge reference. Distribution system?
void SortSbwtBlockwise(BuildIndexRawData &build_index)
{
	char *seq = build_index.seq_raw;
        uint32_t *seq_index = build_index.suffix_array;
        const uint32_t &length_ref = build_index.length_ref;
        const uint32_t &num_block_sort = build_index.num_block_sort;
        uint32_t N = length_ref;
        const uint32_t &period = build_index.period;

        /// Firstly, split the sequence rotation matrix into 4^num_block blocks
        LOGINFO("Firstly, split the sequence rotation matrix into 4^"<< num_block_sort << " blocks\n");
        try {
                SortSbwtBlockwise( seq, seq_index, 0, N, 0, length_ref, period, num_block_sort );
        } catch (...) {
                LOGERROR("SortSbwtBlockwise");
                throw;
        }

        {
                /**
                 * Splitting blocks:
                 * The segments in range [beg0, end0) share with same
                 * spaced prefix (length: num_block_sort * period).
                 * Those blocks sharing with same prefix are sorted
                 * again according to their suffixes.
                 */
                LOGINFO("Splitting blocks...\n");
                vector<char> tmpcv(num_block_sort, '\0');
                uint32_t j = 0;
                uint32_t t0 = 0;
                bool flg;
                uint32_t beg0 = 0, end0 = 1;
                for (j = 0; j < num_block_sort; ++j) {
                       tmpcv[j] = seq[(j*period+seq_index[0])%N];
                }
                for (uint32_t i = 1; i < N; ++i) {
                        flg = true;
                        for (j = 0; j < num_block_sort; ++j) {
                                t0 = (j*period+seq_index[i]) % N;
                                if (seq[t0] != tmpcv[j] && flg) {
                                        flg = false;
                                }
                                tmpcv[j] = seq[t0];
                        }
                        if (!flg) {
                                end0 = i;
                                if (end0 > 1 + beg0) {
                                        try {
                                                LOGINFO("Sort block in [" << beg0 << ",\t" << end0 << ")...\t");
                                                SortSbwt(seq, seq_index, beg0, end0, num_block_sort*period, N, period);
                                                LOGPUT("Done\n");
                                        } catch (...) {
                                                LOGERROR("[" << beg0 << ",\t" << end0 << ")\n");
                                                throw;
                                        }

                                }
                                beg0 = i;
                        }
                }
                /// tail
                LOGINFO("Tail case...\t");
                if (N > 1 + beg0) {
                        try {
                                SortSbwt(seq, seq_index, beg0, N, num_block_sort*period, N, period);
                        } catch (...) {
                                LOGERROR("SortSbwt");
                                throw;
                        }

                }
                LOGPUT("Done\n");
        }
}

void SortSbwtBlockwise(
	char *seq,                      /* sequence */
        uint32_t *seq_index,            /* suffix array */
	uint32_t begin,
        uint32_t end,
        uint32_t depth,
	const uint32_t &length_ref,
        const uint32_t &step,           /* a.k.a period*/
        const uint32_t &num_block
        )
{
	/* Condition of end */
        ///Bug: depth >= num_block
	if (depth >= step*num_block /* Condition for blockwise-building */|| begin+1 >= end || depth >= length_ref) return;

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

        SortSbwtBlockwise(seq, seq_index, begin, r, depth, length_ref, step, num_block);

        tmpval = seq_index[r] + depth;
        if (tmpval < length_ref) {
                SortSbwtBlockwise(seq, seq_index, r, end-d+c, depth+step, length_ref, step, num_block);
        }

        r = d - c;
        SortSbwtBlockwise(seq, seq_index, end-r, end, depth, length_ref, step, num_block);
}


SecondIndex::SecondIndex(): array_ptr(nullptr), size(0) { }

SecondIndex::SecondIndex(const string &prefix_filename)
                : array_ptr(nullptr)
{
        string file_second_filename = prefix_filename + ".second.sbwt";
        std::ifstream second_fin(file_second_filename.c_str(), std::ios_base::in | ios::binary);

        if (!second_fin.is_open()) {
                std::cerr << "Cannot open index files: " << prefix_filename << std::endl;
                return;
        }

        uint16_t tmp16 = readU16(second_fin);
        bool is_big_endian = tmp16 != 0;

        if (is_big_endian) {
                LOGERROR("Current platform is big endian."
                         << " SBWT will not work on big-endian platform.");
        } else {
                this->size = readU64(second_fin, is_big_endian);
                this->size_min = readU32(second_fin, is_big_endian);
                this->size_seed = readU32(second_fin, is_big_endian);

                if (size > 0) {
                        this->array_ptr = new uint16_t[size];
                }

                if (this->array_ptr != nullptr) {
                        uint16_t *ptr = array_ptr;

                        for (uint32_t i = 0; i < size; ++i) {
                                *ptr = readU16(second_fin);
                                ++ptr;
                        }
                }
        }

        second_fin.close();
}

uint32_t SecondIndex::size_min = 200;

SecondIndex::~SecondIndex()
{
        delete []array_ptr;
}


void SecondIndex::RebuildIndex(BuildIndexRawData &build_index) {
        std::string str_iter(size_seed, 0);

        auto seq = build_index.seq_raw;
        auto SA = build_index.suffix_array;
        auto X = build_index.seq_raw;
        auto N = build_index.length_ref;
        auto period = build_index.period;

        if (N <= size_seed) {
                std::cerr << "The length is less than " << size_seed << std::endl;
                return;
        }

        uint32_t count;
        uint32_t loop_end = size_seed * period;
        /// Build second index
        count = 0;
        uint32_t index_array = 0;
        uint32_t beg0 = 0;
        uint32_t end0 = 0;
        uint32_t depth0 = 0;
        for (uint32_t i = 0; i < N; ++i) {
                for (uint32_t j = 0; j < loop_end; j+=period) {
                        if (str_iter[j/period] != X[(j+SA[i])%N]) {
                                if (count != 0) {
                                        /// Sort blockwise
                                        if (count >= size_min && count < 65536) {
                                                beg0 = i - count;
                                                end0 = i;
                                                uint32_t pos_saved = SA[beg0];
#if DEBUG_SECONDINDEX
                                                cout << "beg0: " << beg0
                                                     << ",end0: " << end0
                                                     <<  ",count: " << count
                                                     << ",SA[beg0]: " << SA[beg0]
                                                     << ",index_array: " << index_array
                                                     << endl;
#endif
                                                unordered_map<uint32_t, uint16_t > mmap;
                                                /// to restore SA[beg0, end0)
                                                vector<uint32_t > index_backup(count,0);
                                                uint16_t tmp0 = 0;
                                                for (uint32_t u = beg0; u < end0; ++u) {
                                                        mmap[SA[u]] = tmp0;
                                                        index_backup[tmp0] = SA[u];
                                                        ++tmp0;
                                                }

                                                for (uint32_t tau = 0; tau != period; ++tau) {
                                                        depth0 = tau * size_seed;
                                                        SortSbwt(seq, SA, beg0, end0, depth0, N, 1, depth0 + size_seed);
#if DEBUG_SECONDINDEX
                                                        for (uint32_t u = beg0; u < end0; ++u) {
                                                                for (uint32_t k0 = 0; k0 < loop_end; ++k0) {
                                                                        if (k0 % period == 0) {
                                                                                cout << "\033[1;31m" << X[(SA[u] + k0) % N] << "\033[0m";
                                                                        }
                                                                        else {
                                                                                cout << X[(SA[u] + k0) % N];
                                                                        }
                                                                        if (k0 % size_seed == (size_seed-1)) {
                                                                                cout << " ";
                                                                        }
                                                                }
                                                                cout << endl;
                                                        }
#endif
                                                        for (uint32_t u = beg0; u < end0; ++u) {
                                                                uint32_t i_shift = (u-beg0) + index_array + 3 + tau*count;
                                                                *(this->array_ptr + i_shift) = mmap[SA[u]];
#if DEBUG_SECONDINDEX
                                                                cout << *(array_ptr +i_shift) << " ";
#endif
                                                        }
#if DEBUG_SECONDINDEX
                                                        cout << endl;
#endif
                                                }

                                                /// exchange
                                                uint16_t *p16 = this->array_ptr + index_array;
                                                uint32_t *p32 = (uint32_t*) p16;
                                                *p32 = pos_saved;
                                                /// restore SA
                                                for (uint32_t u = beg0; u < end0; ++u) {
                                                        SA[u] = index_backup[u - beg0];
                                                }
                                                /// change the first element of SA
                                                SA[beg0] = index_array;
                                                *(this->array_ptr+index_array+2) = (uint16_t)count;

                                                index_array += count*period + 3;
                                        }
                                }

                                count = 0;
                                for (uint32_t k = 0; k != size_seed; ++k) {
                                        str_iter[k] = X[(k*period+SA[i])%N];
                                }
                                break;
                        }
                } /* j */
                ++count;
        } /* i */

        /// tail case
        if (count != 0) {
                if (count >= size_min && count < 65536) {
                        beg0 = N-count;
                        end0 = N;
                        uint32_t pos_saved = SA[beg0];
#if DEBUG_SECONDINDEX
                        cout << beg0 << ", " << end0 <<  "   " << count << endl;
#endif
                        unordered_map<uint32_t, uint16_t > mmap;
                        /// to restore SA[beg0, end0)
                        vector<uint32_t > index_backup(count,0);
                        uint16_t tmp0 = 0;
                        for (uint32_t u = beg0; u < end0; ++u) {
                                mmap[SA[u]] = tmp0;
                                index_backup[tmp0] = SA[u];
                                ++tmp0;
                        }

                        for (uint32_t tau = 0; tau != period; ++tau) {
                                depth0 = tau * size_seed;
                                SortSbwt(seq, SA, beg0, end0, depth0, N, 1, depth0 + size_seed);
                                for (uint32_t u = beg0; u < end0; ++u) {
                                        uint32_t i_shift = (u-beg0) + index_array + 3 + tau*count;
                                        *(this->array_ptr + i_shift) = mmap[SA[u]];
                                }
                        }

                        /// exchange
                        uint16_t *p16 = this->array_ptr + index_array;
                        uint32_t *p32 = (uint32_t*) p16;
                        *p32 = pos_saved;
                        /// restore SA
                        for (uint32_t u = beg0; u < end0; ++u) {
                                SA[u] = index_backup[u - beg0];
                        }
                        /// change the first element of SA
                        SA[beg0] = index_array;
                        *(this->array_ptr + index_array + 2) = (uint16_t)count;
                }
        }

        /*
         * Cannot show it. Because the BuildIndexRawData is disordered.
         */
#if 0
        /// visualization in array
        this->PrintSecondIndex(build_index);
#endif
}


void SecondIndex::PrintSecondIndex(BuildIndexRawData &build_index)
{
        if (this->array_ptr == nullptr) {
                return;
        }

        std::string str_iter(size_seed, 0);

        auto SA = build_index.suffix_array;
        auto X = build_index.seq_raw;
        auto N = build_index.length_ref;
        auto period = build_index.period;

        if (N <= size_seed) {
                std::cerr << "The length is less than " << size_seed << std::endl;
                return;
        }

        uint32_t count;
        uint32_t loop_end = size_seed * period;

        count = 0;
        for (uint32_t i = 0; i < N; ++i) {
                for (uint32_t j = 0; j < loop_end; j+=period) {
                        if (str_iter[j/period] != X[(j+SA[i])%N]) {
                                if (count != 0) {
                                        /// Sort blockwise
                                        if (count >= size_min && count < 65536) {
                                                /// Caution: may bring into SEGERROR
                                                /// because SA[i] in "if (str_iter[j/period] != X[(j+SA[i])%N])" is
                                                /// changed into index of array.
                                                uint32_t beg0 = i - count - 1;
                                                uint32_t end0 = i;
                                                uint32_t index_array = SA[beg0];
                                                uint16_t *p16 = this->array_ptr + index_array;
                                                uint32_t *p32 = (uint32_t*) p16;
                                                uint32_t pos = *p32;
                                                cout << "beg0: "        << beg0
                                                     << ",end0: "       << end0
                                                     << ",SA[beg0]: "   << SA[beg0]
                                                     << ",seg_size: "   << p16[2]
                                                     << ",pos: "        << pos
                                                     << ",period: "     << period
                                                     << endl;
                                                uint16_t size_seg = p16[2];
                                                cout << "Before second index" << endl;
                                                for (uint32_t i0 = beg0; i0 < end0; ++i0) {
                                                        uint32_t p0 = SA[i0];
                                                        if (i0 == beg0) {
                                                                p0 = pos;
                                                        }
                                                        for (uint32_t j0 = 0; j0 < loop_end; ++j0) {
                                                                if (j0 % period == 0) {
                                                                        cout << "\033[1;31m" << X[(j0 + p0) % N] << "\033[0m";
                                                                }
                                                                else {
                                                                        cout << X[(j0 + p0) % N];
                                                                }
                                                                if (j0 % size_seed == (size_seed-1)) {
                                                                        cout << " ";
                                                                }
                                                        }
                                                        cout << endl;
                                                }
                                                cout << "Second index" << endl;
                                                p16 += 3;
                                                for (uint32_t tau = 0; tau < period; ++tau) {
                                                        cout << endl;
                                                        for (uint32_t i0 = 0; i0 < size_seg; ++i0) {
                                                                cout << p16[i0] << " ";
                                                        }
                                                        cout << endl;

                                                        for (uint32_t i0 = 0; i0 < size_seg; ++i0) {
                                                                uint32_t p0 = SA[p16[i0] + beg0];
                                                                if (p0 == 0) {
                                                                        p0 = pos;
                                                                }

                                                                for (uint32_t j0 = 0; j0 < loop_end; ++j0) {
                                                                        if (j0 % period == 0) {
                                                                                cout << "\033[1;31m" << X[(j0 + p0) % N] << "\033[0m";
                                                                        }
                                                                        else {
                                                                                cout << X[(j0 + p0) % N];
                                                                        }
                                                                        if (j0 % size_seed == (size_seed-1)) {
                                                                                cout << " ";
                                                                        }
                                                                }
                                                                cout << endl;
                                                        }
                                                        p16 += size_seg;
                                                }
                                        }
                                }

                                count = 0;
                                for (uint32_t k = 0; k != size_seed; ++k) {
                                        str_iter[k] = X[(k*period+SA[i])%N];
                                }
                                break;
                        }
                } /* j */
                ++count;
        } /* i */

        /// tail case
        if (count != 0) {
                if (count >= size_min && count < 65536) {
                        /// Ignore
                }
        }

}



bool SecondIndex::Empty()
{
        return this->array_ptr == nullptr;
}


/// Init SecondIndex data structure
void SecondIndex::RebuildIndexInit(BuildIndexRawData &build_index, uint32_t size_seed)
{
        this->size_seed = size_seed;
        std::string str_iter(size_seed, 0);
        size = 0;

        auto SA = build_index.suffix_array;
        auto X = build_index.seq_raw;
        auto N = build_index.length_ref;
        auto period = build_index.period;

        if (N <= size_seed) {
                std::cerr << "The length is less than " << size_seed << std::endl;
                return;
        }

        uint32_t count;
        uint32_t loop_end = size_seed * period;

        count = 0;
        for (uint32_t i = 0; i < N; ++i) {
                for (uint32_t j = 0; j < loop_end; j+=period) {
                        if (str_iter[j/period] != X[(j+SA[i])%N]) {
                                if (count != 0) {
                                        /// Sort blockwise
                                        //cout << (i - count) << "," << (i) << endl;
                                        if (count >= size_min && count < 65536/* TODO */) {
                                                size += 3 + count*period;
                                        }
                                }

                                count = 0;
                                for (uint32_t k = 0; k != size_seed; ++k) {
                                        str_iter[k] = X[(k*period+SA[i])%N];
                                }
                                break;
                        }
                } /* j */
                ++count;
        } /* i */

        /// tail case
        if (count != 0) {
                if (count >= size_min && count < 65536) {
                        size += 3 + count*period;
                }
        }

        array_ptr = new uint16_t[size];

}


} /* namespace sbwt */
