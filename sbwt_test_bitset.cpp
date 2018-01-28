#include <algorithm>
#include <iostream>
#include <string>
#include <bitset>
#include <memory>

#include "io_build_index.h"
#include "utility.h"
#include "sbwt_search.h"
#include "sequence_pack.h"

using std::shared_ptr;
using std::string;
using namespace utility;
using std::bitset;


int main(int argc, char **argv)
{
        if (argc != 4) {
                PrintHelp_SbwtTestBitset(argc, argv);
                return 1;
        }

#ifdef SBWT_VERBOSE
        using std::string;
        string file_fasta(argv[1]);
        string file_index(argv[2]);
        char option = argv[3][0];
        if (option == '1') {
                sbwt::TestBit64Map(file_fasta, file_index);
        } else if (option == '2') {
                sbwt::TestBit8Map(file_fasta, file_index);
        } else {
                sbwt::TestNaiveMap(file_fasta, file_index);
        }
#endif /* SBWT_VERBOSE */
        return 0;
}