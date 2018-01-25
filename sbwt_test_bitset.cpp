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


void test_bit8_rc(char **argv);

int main(int argc, char **argv)
{
        if (argc != 4) {
                PrintHelp_SbwtTestBitset(argc, argv);
                return 1;
        }

        test_bit8_rc(argv);
        return 0;

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


void test_bit8_rc(char **argv)
{
        using std::string;
        using std::shared_ptr;
        using sbwt::reads_buffer;
        using std::cout;
        using std::endl;
        using std::bitset;

        string file_fasta(argv[1]);
        string file_index(argv[2]);

        reads_buffer rb_fasta(file_fasta);
        reads_buffer rb_index(file_index);

        rb_index.ReadNext();
        rb_index.ReadNext();

        uint32_t size_ref_8bit = rb_index.length_read;
        uint32_t size_ref_char = rb_index.length_read;

        size_ref_8bit = size_ref_8bit % 4 == 0 ? size_ref_8bit >> 2 : (size_ref_8bit>>2) + 1;

        shared_ptr<uint8_t > ref_bin_sptr(new uint8_t[size_ref_char+1024]);
        shared_ptr<uint8_t > ref_bin_rc_sptr(new uint8_t[size_ref_char+1024]);
        shared_ptr<uint8_t > ref_bin_rc_sptr1(new uint8_t[size_ref_char+1024]);

        uint8_t *ref_bin_ptr = ref_bin_sptr.get();
        uint8_t *ref_bin_rc_ptr = ref_bin_rc_sptr.get();
        uint8_t *ref_bin_rc_ptr1 = ref_bin_rc_sptr1.get();


        sbwtio::BaseChar2Binary8B(rb_index.buffer, size_ref_8bit, ref_bin_ptr);
        for (uint32_t i = 1<<24; i != 0; --i) {
                sbwtio::BaseChar2Binary8B_RC(ref_bin_ptr, size_ref_char, ref_bin_rc_ptr);
                sbwtio::BaseChar2Binary8B_RC(rb_index.buffer, size_ref_char, ref_bin_rc_ptr1);
        }

}
