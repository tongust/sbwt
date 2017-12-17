#include <stdlib.h>  // for strtol uint32_t etc.
#include <string.h>
#include <stdio.h>

#include <string>
#include <iostream>
#include <vector>
#include <memory> // for shared_ptr

#include "log.h"
#include "sbwt.h"
#include "const.h"
#include "utility.h"

using std::shared_ptr;

int main(int argc, char **argv)
{
        using namespace utility;

        if (argc != 3) {
                PrintHelp_BuildIndex(argc, argv);
                return 1;
        }

        uint32_t period = GetUint(argc, argv[1]);
        uint32_t num_block_sort = 4;
        uint32_t read_length = 0;
        char *file_name = argv[2];

        char *seq = nullptr;
        /* read reference file */
        {
                uint32_t file_size = 0;
                shared_ptr<char> seq_shared(ReadFasta(file_name, file_size));
                seq = Extract(seq_shared.get(), file_size, read_length);
        }

        sbwt::BuildIndexRawData build_index(seq, read_length, period, num_block_sort);
        sbwt::BuildIndex(build_index);

        sbwt::PrintFullSearchMatrix(build_index);

        return 0;
}
