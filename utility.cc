#include <limits.h>  // for INT_MAX
#include <stdio.h>
#include <stdlib.h>  // for strtol uint32_t etc.

#include <iostream>
#include <string>
#include <vector>

#include "const.h"
#include "log.h"
#include "sbwt.h"
#include "utility.h"

namespace utility
{

using std::string;
using std::cout;
using std::endl;

/* Extract pure DNA nucleotides (A/C/G/T)
 */
char *Extract(char *seq, const uint32_t &seq_size, uint32_t &trimmed_length)
{
        uint32_t total_num = 0;
        char *ptr = seq;
        for (uint32_t i = 0; i < seq_size; ++i) {
                if (IsDNA(*ptr)) {
                        ++total_num;
                }
                ++ptr;
        }

        /// enough memory for $s and "boarder case"
        /// And the period must be less than 1024
        /// In case the boarder of sequence array will be reached.
        char *ret_ptr = new char[((total_num/1024)+4)*1024];

        char *ptr0 = ret_ptr;
        ptr = seq;
        for (uint32_t i = 0; i < seq_size; ++i) {
                if (IsDNA(*ptr)) {
                        *ptr0 = *ptr;
                        ++ptr0;
                }
                ++ptr;
        }

        trimmed_length = total_num;

        return ret_ptr;

}

/* Read raw sequence file in fasta format
 * TODO
 * Raw fasta file to DNA string.
 * Better to read and extract DNA string blockwise.
 */
char *ReadFasta(char *file_name, uint32_t &read_length)
{
        char *buffer = nullptr;
        uint32_t string_length;
        FILE *handler = fopen(file_name, "r");

        if (handler) {
                fseek(handler, 0, SEEK_END);    /* Seek the last byte */
                string_length = ftell(handler); /* Calculate the offset from head to tail */
                rewind(handler);                /* Back to the head */

                //buffer = (char *) malloc(sizeof(char) * (string_length + 1));
                buffer = new char[sizeof(char) * (string_length + 1)];
                read_length = 0;
                read_length = fread(buffer, sizeof(char), string_length, handler);

                buffer[string_length] = '\0';

                if (read_length != string_length) {
                        delete []buffer;
                        buffer = nullptr;
                        logger::LogDebug("Errors in fread");
                        return buffer;
                }
        }
        return buffer;
}

bool IsDNA(char c)
{
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

bool IsN(char c)
{
        return c == 'N';
}

/** Get integers from argv.
 */
uint32_t GetUint(int argc, char *parameter)
{
        errno = 0;
        uint32_t ret;
        char *tmpptr;
        auto tmpval= strtoul(parameter, &tmpptr, 10);

        if (errno != 0 || *tmpptr != '\0' || tmpval > UINT_MAX) {
                LOGDEBUG("Error: invalid inputs");
                exit(1);
        } else {
                ret = tmpval;
        }
        return ret;
}

/* Count the occurrence of seeds (25 bp) in reference.
 * */
void CountSeedOccurrence(sbwt::BuildIndexRawData &build_index, uint32_t seed_length)
{
        using std::vector;
        std::string  str_iter(seed_length, 0);

        auto SA = build_index.suffix_array;
        auto X = build_index.seq_raw;
        auto N = build_index.length_ref;
        auto period = build_index.period;

        if (N <= seed_length) {
                logger::LogDebug("The length is less than 25");
                return;
        }

        uint32_t count = 0;  /* count the # of same DNA sequence */

        vector<uint32_t> counts((1<<20), 0);

        count = 0;
        uint32_t loop_end = seed_length * period;
        for ( size_t i = 0; i < N; ++i) {
                for (uint32_t j = 0; j < loop_end; j+=period) {
                        if (str_iter[j/period] != X[(j+SA[i])%N]) {
                                ++counts[count];
                                count = 0;
                                for (uint32_t k = 0; k != seed_length; ++k) {
                                        str_iter[k] = X[(k*period+SA[i])%N];
                                }
                                break;
                        }
                } /* j */
                ++count;

        } /* i */

        /// tail case
        if (count != 0) {

                ++counts[count];
        }
        for (uint32_t i = 1; i < counts.size(); ++i) {
                if (counts[i] != 0) {
                         cout << i << "\t" << counts[i] << "\n";
                }
        }
}

void PrintHelp_BuildIndex(int argc, char **argv)
{
        cout << "usage: build_index [fa] [period] <size_seed>"
             << endl;
}

void PrintHelp_CountOcc(int argc, char **argv)
{
        cout << "usage: count_occ "
             << "[prefix name of index files] [seed_seed]"
             << endl;
}
        
void PrintHelp_SbwtAligner(int argc, char **argv)
{
        cout << "usage: sbwt "
             << "[fasta file] [prefix of index files] <size of replica, default: "
             << sbwt::SecondIndex::size_min
             << ">"
             << endl;

}


void PrintHelp_SbwtTestBitset(int argc, char **argv)
{
        cout << "usage: sbwt_test(test bit set)[fasta file]"
             << "[prefix of index files] [option<1: 64bit; "
             << "2: 8bit; otherwise: naive>]"
             << endl;
}


void PrintHelp_SbwtTestExactMatch(int argc, char **argv)
{
        cout << "usage: sbwt_test(test exact match) "
             << "[reads.fa] [index]"
             << endl;
}

} /* namespace utility */
