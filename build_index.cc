#include <stdlib.h>

#include <string>

#include "log.h"
#include "sbwt.h"
#include "const.h"

using std::string;
using std::cout;
using std::endl;

int main(int argc, char **argv)
{
        logger::LogDebug("Test of Build-Index\n");

        string dna = "ACGCATGATAGCAATGATAGTAGCTA";

        uint32_t length_seq = 1<<4;
        uint32_t period = 1;
        std::cin >> period;
        uint32_t num_block_sort = 4;


        char *seq = new char[length_seq];
        for (uint32_t i = 0; i < length_seq; ++i) {
                seq[i] = dna[i];
        }
        sbwt::BuildIndexRawData build_index(seq, length_seq, period, num_block_sort);
        sbwt::BuildIndex(build_index);
        sbwt::PrintFullSearchMatrix(build_index);
	return 0;
}
