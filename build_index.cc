#include <stdio.h>   /* printf, scanf, puts, NULL */
#include <stdlib.h>  // for strtol uint32_t etc.
#include <time.h>    /* time */

#include <string>
#include <iostream>
#include <vector>
#include <memory> // for shared_ptr
#include <fstream> // for ofstream, ifstream

#include "sbwt.h"
#include "utility.h"
#include "io_build_index.h"
#include "log.h"

using std::shared_ptr;
using std::string;
using namespace utility;

int main(int argc, char **argv)
{
        if (argc != 3) {
                PrintHelp_BuildIndex(argc, argv);
                return 1;
        }

        /* initialize random seed: */
        //srand (time(NULL));

        uint32_t period = GetUint(argc, argv[1]);
        uint32_t num_block_sort = 5;
        char *file_name = argv[2];
#if 0
        uint32_t read_length = 0;
        char *seq = nullptr;

        /* read reference file */
        {
                LOGDEBUG("Read reference...");
                uint32_t file_size = 0;
                shared_ptr<char> seq_shared(ReadFasta(file_name, file_size));
                seq = Extract(seq_shared.get(), file_size, read_length);
                LOGDEBUG("Total length: ");
                LOGDEBUG(read_length);

        }

        LOGDEBUG("Init build index...");
        sbwt::BuildIndexRawData build_index(seq, read_length, period, num_block_sort);
#endif
#if 1
        LOGINFO("Read reference and init index...\n");
        sbwt::BuildIndexRawData build_index(file_name, period, num_block_sort);
        LOGINFO("Total length: " << build_index.length_ref << "\n");
#endif

#define SECOND_INDEX

#ifdef SECOND_INDEX
        LOGINFO("Building index...\n");
        sbwt::BuildSortedIndexBlockwise(build_index);

        /*
         * Prerequisite in sorted index:
         *     The length of reads must be less than length of seed by period.
         *     For example, if processing Illumina's reads (150 bp),
         *     it is recommended that pair of period and seed: 2 with 75, 3 with 50,
         *     5 with 30, 6 with 25.
         * Proof:
         *     Because of resorting, the suffix array range in same-sampled-prefix matrix
         *     is disordered. And the Occ and C will be changed and the spaced suffix array
         *     ranging search equation will not hold. If the length of reads is larger than
         *     period by length of seed, the left and right range of reads may be covered.
         * Mine Field:
         *     The range of sorted index.
         */
        sbwt::SecondIndex secondIndex;
        /// Init SecondIndex
        secondIndex.RebuildIndexInit(build_index, 10);
        /// Sorting
        secondIndex.RebuildIndex(build_index);

        /// must add transform and count occurrence
        sbwt::BuildSortedIndexTransCountOcc(build_index);

#if DEBUG_SECONDINDEX
        sbwt::PrintFullSearchMatrix(build_index);
#endif
        /// write into disk
        LOGINFO("Write into disk...\n");
        sbwt::WriteIntoDiskBuildIndex(build_index, string(file_name));

        LOGINFO("Write second index...\n");
        sbwt::WriteIntoDiskBuildSecondIndex(build_index, string(file_name), secondIndex);
        LOGINFO("Done\n");

        LOGINFO("Build index done\n");


#else
        LOGINFO("Building index...\n");
        sbwt::BuildIndexBlockwise(build_index);

        //sbwt::PrintFullSearchMatrix(build_index);

        /// write into disk
        LOGINFO("Write into disk...\n");
        sbwt::WriteIntoDiskBuildIndex(build_index, string(file_name));
        LOGINFO("Build index done\n");

#endif

        return 0;
}
