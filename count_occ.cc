/*
 * Count the occurrence of the short read in reference read
 */
#include <errno.h>   // for errno
#include <limits.h>  // for INT_MAX
#include <stdlib.h>  // for strtol uint32_t etc.

#include <string>
#include <iostream>
#include <vector>
#include <memory> // for shared_ptr

#include "sbwt.h"
#include "utility.h"

using std::string;
using std::cout;
using std::endl;
using std::shared_ptr;

int main(int argc, char **argv)
{
        using namespace utility;

        if (argc != 4) {
                PrintHelp_CountOcc(argc, argv);
                return 1;
        }

        uint32_t period = GetUint(argc, argv[1]);
        uint32_t num_block_sort = 4;
        uint32_t seed_length = GetUint(argc, argv[3]);
        uint32_t read_length = 0;
        char *file_name = argv[2];
        char *seq = nullptr;
        {
                uint32_t file_size = 0;
                shared_ptr<char> seq_shared(ReadFasta(file_name, file_size));
                seq = Extract(seq_shared.get(), file_size, read_length);
        }

        sbwt::BuildIndexRawData build_index(seq, read_length, period, num_block_sort);
        sbwt::BuildIndex(build_index);

        CountSeedOccurrence(build_index, seed_length);

        return 0;
}
