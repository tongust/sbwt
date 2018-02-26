#include <algorithm>
#include <iostream>
#include <bitset>
#include <memory>

#include "alphabet.h"
#include "io_build_index.h"
#include "sbwt_search.h"
#include "sequence_pack.h"
#include "utility.h"

using std::shared_ptr;
using std::string;
using namespace utility;
using std::bitset;

namespace sbwt_test
{
void TestBitset(int argc, char **argv);
void TestSbwtExactMatch(int argc, char **argv);
} /* namespace sbwt_test */

int main(int argc, char **argv)
{
        using namespace sbwt_test;

        TestSbwtExactMatch(argc, argv);
        return 0;
}

namespace sbwt_test
{
void TestBitset(int argc, char **argv)
{
        if (argc != 4) {
                PrintHelp_SbwtTestBitset(argc, argv);
                return;
        }

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
}


void TestSbwtExactMatch(int argc, char **argv)
{
        using namespace sbwt;
        using std::cout;
        using std::endl;

        if (argc != 3) {
                PrintHelp_SbwtTestExactMatch(argc, argv);
                return;
        }

        /// Build index from files
        string reads_filename = string(argv[1]);
        string prefix_filename = string(argv[2]);
        BuildIndexRawData build_index(prefix_filename);
        //PrintFullSearchMatrix(build_index);

        uint32_t period = build_index.period;
        uint32_t N = build_index.length_ref;
        uint32_t *C = build_index.first_column;
        uint32_t **Occ = build_index.occurrence;
        uint32_t *SA = build_index.suffix_array;

        /// reads
        /// Read first segment.
        reads_buffer rb_reads(reads_filename);
        for (;;) {
                rb_reads.ReadNext();/// header
                if (rb_reads.length_read < 1) { break; }
                rb_reads.ReadNext();/// sequence
                if (rb_reads.length_read < 1) { break; }

                uint32_t size_read_char = rb_reads.length_read;
                uint32_t L, R;

                //L = 0;
                R = N - 1;

                /// Assuming the sequence has been down-sampled.
                char *ptr = rb_reads.buffer + (size_read_char - 1);
                static char key = 0;
                static uint32_t a = 0;
                int lmd = size_read_char;

                /// Init
                key = *ptr;
                a = charToDna5[key];
                L = C[a];
                R = C[a] + Occ[a][R] - 1;
                //cout << key << " " << L << " " << R << endl;

                for (int i = 1; i != lmd; ++i) {
                        if (L > R) {
                                break;
                        }
                        --ptr;
                        key = *ptr;
                        a = charToDna5[key];
                        L = C[a] + Occ[a][L-1];
                        R = C[a] + Occ[a][R] - 1;
                        //cout << key << " " << L << " " << R << endl;
                }

                if (L > R) {
                        cout << -1 << "\n";
                        continue;
                } else {
                        for (uint32_t l = L; l <= R; ++l) {
                                cout << SA[l];
                                if (l != R) {
                                        cout << ",";
                                }
                        } cout << endl;
                }
        }
}
} /* namespace sbwt_test */
