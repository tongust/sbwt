#include <algorithm>
#include <iostream>
#include <memory>

#include "io_build_index.h"
#include "utility.h"

using std::shared_ptr;
using std::string;
using namespace utility;

int main(int argc, char **argv)
{
        if (argc != 2) {
                PrintHelp_SbwtAligner(argc, argv);
                return 1;
        }

        string prefix_filename = string(argv[1]);
        sbwt::BuildIndexRawData build_index(prefix_filename);
        //sbwt::PrintFullSearchMatrix(build_index);

        return 0;
}