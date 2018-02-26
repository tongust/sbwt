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
        LOGINFO("Building index...\n");
        sbwt::BuildIndexBlockwise(build_index);

        //sbwt::PrintFullSearchMatrix(build_index);

        /// write into disk
        LOGINFO("Write into disk...\n");
        sbwt::WriteIntoDiskBuildIndex(build_index, string(file_name));
        LOGINFO("Build index done\n");

        return 0;
}
