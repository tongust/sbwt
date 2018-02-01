#include <stdlib.h>  // for strtol uint32_t etc.

#include <string>
#include <iostream>
#include <vector>
#include <memory> // for shared_ptr
#include <fstream> // for ofstream, ifstream

#include "sbwt.h"
#include "utility.h"
#include "io_build_index.h"

using std::shared_ptr;
using std::string;
using namespace utility;

int main(int argc, char **argv)
{
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
        sbwt::BuildIndexBlockwise(build_index);

        //sbwt::PrintFullSearchMatrix(build_index);

        /// write into disk
        sbwt::WriteIntoDiskBuildIndex(build_index, string(file_name));

        return 0;
}
