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

        if (argc != 3) {
                PrintHelp_CountOcc(argc, argv);
                return 1;
        }

        string prefix_filename = string(argv[1]);
        uint32_t seed_length = GetUint(argc, argv[2]);

        sbwt::BuildIndexRawData build_index(prefix_filename);

        CountSeedOccurrence(build_index, seed_length);

        return 0;
}
