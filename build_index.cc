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
        if (argc < 3) {
                PrintHelp_BuildIndex(argc, argv);
                return 1;
        }

        /* initialize random seed: */
        //srand (time(NULL));

        char *file_name = argv[1];
        uint32_t period = GetUint(argc, argv[2]);
        uint32_t size_seed = 0;
        if (argc >= 4) {
                size_seed = GetUint(argc, argv[3]);
        }
        uint32_t num_block_sort = 5;

        LOGINFO("Read reference and init index...\n");
        sbwt::BuildIndexRawData build_index(file_name, period, num_block_sort);
        LOGINFO("Total length: " << build_index.length_ref << "\n");

        if (size_seed > 0) {
                LOGINFO("Building index...\n");

                /// Build SA, Occ, B, and C
                sbwt::BuildIndexBlockwise(build_index);

#if DEBUG_SECONDINDEX
                //sbwt::PrintFullSearchMatrix(build_index);
#endif

                LOGINFO("Build second index...\n");
                sbwt::SecondIndex secondIndex;
                /// Init SecondIndex
                secondIndex.RebuildIndexInit(build_index, size_seed);
                /// Sorting
                /// Change the first index in SA
                secondIndex.RebuildIndex(build_index);
                LOGINFO("Size of second index: " << secondIndex.size << "\n");

                /// write SA(changed), Occ, B and C into disk
                LOGINFO("Write into disk...\n");
                sbwt::WriteIntoDiskBuildIndex(build_index, string(file_name));


                sbwt::WriteIntoDiskBuildSecondIndex(build_index, string(file_name), secondIndex);
                LOGINFO("Done\n");

                LOGINFO("Build index done\n");


        }
        else {
                LOGINFO("Building index...\n");
                sbwt::BuildIndexBlockwise(build_index);

                //sbwt::PrintFullSearchMatrix(build_index);

                /// write into disk
                LOGINFO("Write into disk...\n");
                sbwt::WriteIntoDiskBuildIndex(build_index, string(file_name));
                LOGINFO("Build index done\n");

        }

        return 0;
}
