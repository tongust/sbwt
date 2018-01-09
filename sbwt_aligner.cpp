#include <algorithm>
#include <iostream>
#include <string>
#include <bitset>

#include "io_build_index.h"
#include "utility.h"
#include "sbwt_search.h"

using std::shared_ptr;
using std::string;
using namespace utility;
using std::bitset;

int main(int argc, char **argv)
{
        if (argc != 3) {
                PrintHelp_SbwtAligner(argc, argv);
                return 1;
        }
#ifdef SBWT_VERBOSE
        using std::string;
        string file_fasta(argv[1]);
        string file_index(argv[2]);
        sbwt::Test_reads_buffer(file_fasta, file_index);
#endif /* SBWT_VERBOSE */

#if 0
        string reads_filename = string(argv[1]);
        string prefix_filename = string(argv[2]);
        sbwt::BuildIndexRawData build_index(prefix_filename);
#endif
        return 0;
}
