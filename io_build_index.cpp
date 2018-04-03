#include <ctype.h>

#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cassert>
#include <memory> // for shared_ptr

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>

#include "io_build_index.h"
#include "alphabet.h"
#include "ref_read.h"
#include "endian_swap.h"
#include "word_io.h"
#include "sequence_pack.h"
#include "sbwt.h"
#include "log.h"

using namespace seqan;
using namespace std;

namespace sbwt
{

void WriteIntoDiskBuildIndex(sbwt::BuildIndexRawData &build_index, const string &prefix_filename)
{
        vector<uint32_t > nums(build_index.length_ref);
        for (int i = 0; i != build_index.length_ref; ++i) {nums[i] = build_index.suffix_array[i];}
        string file_array_filename = prefix_filename
                                     + "."
                                     + std::to_string(build_index.period)
                                     + ".array.sbwt";
        string file_meta_filename = prefix_filename
                                    + "."
                                    + std::to_string(build_index.period)
                                    + ".meta.sbwt";

        std::ofstream array_fout(file_array_filename.c_str(), std::ios::binary);
        std::ofstream meta_fout(file_meta_filename.c_str(), std::ios::binary);

        uint32_t size_packed_seq_8bit = build_index.length_ref / 4;
        size_packed_seq_8bit += build_index.length_ref % 4 == 0 ? 0 : 1;
        ///uint64_t size_packed_seq = build_index.length_ref / 32;


        /**
         * Meta information
         */
        LOGINFO("Write meta information...\n");
        bool is_bigendian = currentlyBigEndian();
        writeU32(meta_fout, (uint32_t)is_bigendian, is_bigendian);      /* The endian flag */
        writeU32(meta_fout, build_index.length_ref, is_bigendian);      /* Length of reference sequence including $s*/
        writeU32(meta_fout, build_index.num_block_sort, is_bigendian);  /* Number of blocks, 4 for 256 */
        writeU32(meta_fout, build_index.num_dollar, is_bigendian);      /* The # of $s those are appended */
        writeU32(meta_fout, build_index.period, is_bigendian);          /* The period of sbwt */
        writeU32(meta_fout, size_packed_seq_8bit, is_bigendian);        /* The size of packed sequence */

        /// write first column
        for (int i = 0; i != 4; ++i) {
                writeU32(meta_fout, build_index.first_column[i], is_bigendian);
        }

        /**
         * Array and sequence
         */
/// 64-bit version
# if 0
        {
                std::shared_ptr<uint64_t > binary_seq64_ptr(new uint64_t[size_packed_seq]);
                sbwtio::BaseChar2Binary64B(build_index.seq_raw, build_index.length_ref, binary_seq64_ptr.get());

                uint64_t *beg = binary_seq64_ptr.get();
                uint64_t *end = beg + (size_packed_seq-1);
                for (;beg != end; ++beg) {
                        writeU64(array_fout, *beg);
                }
        }
#endif
        /// 8-bit version
        {
                /// Watch out for the boarder
                LOGINFO("Packed binary sequence...\n");
                std::shared_ptr<uint8_t > binary_8bit_sptr(new uint8_t[(size_packed_seq_8bit * 4) + 1024]);

                for (int i = 0; i != 4; ++i) {
                        sbwtio::BaseChar2Binary8B(build_index.seq_raw + i,
                                                  size_packed_seq_8bit,
                                                  binary_8bit_sptr.get() + i*size_packed_seq_8bit);

                        array_fout.write((const char*)(binary_8bit_sptr.get() + i*size_packed_seq_8bit),
                                         size_packed_seq_8bit);
                }
        }

        /// The raw sequence
        /// TODO map directly the memory to files
        LOGINFO("Write raw sequence...\n");
        array_fout.write(build_index.seq_raw, build_index.length_ref);

        //array_fout.write(build_index.seq_transformed, build_index.length_ref);

        /// Occurrence
        LOGINFO("Write raw occurrence...\n");
        for (int i = 0; i != 4; ++i) {
                uint32_t *beg = build_index.occurrence[i];
                uint32_t *end = beg+build_index.length_ref;
                while (beg != end) {
                        writeU32(array_fout, *beg, is_bigendian);
                        ++beg;
                }
        }

        /// suffix array
        LOGINFO("Write raw suffix array...\n");
        uint32_t *beg = build_index.suffix_array;
        uint32_t *end = build_index.suffix_array + build_index.length_ref;
        while (beg != end) {
                writeU32(array_fout, *beg, is_bigendian);
                ++beg;
        }

        array_fout.flush();
        array_fout.close();
        meta_fout.flush();
        meta_fout.close();
}


void WriteIntoDiskBuildSecondIndex(
                sbwt::BuildIndexRawData &build_index,
                const string &prefix_filename,
                SecondIndex &second_index)
{
        string file_second_filename = prefix_filename
                                      + "."
                                      + std::to_string(build_index.period)
                                      + "."
                                      + std::to_string(second_index.size_seed)
                                      + ".second.sbwt";
        std::ofstream second_fout(file_second_filename.c_str(), std::ios::binary);
        LOGINFO("Write second index...\n");
        bool is_bigendian = currentlyBigEndian();
        if (is_bigendian) {
                LOGERROR("Current platform is big endian."
                         << " SBWT will not work on big-endian platform.");
                writeU16(second_fout, (uint16_t)is_bigendian);
                second_fout.flush();
                second_fout.close();
                return;
        }

        /// Write flag on endianness
        writeU16(second_fout, (uint16_t)is_bigendian);
        /// meta information
        writeU64(second_fout, second_index.size);
        writeU32(second_fout, second_index.size_min);
        writeU32(second_fout, second_index.size_seed);

        /// Write array
        auto ptr = second_index.array_ptr;
        for (uint32_t i = 0; i < second_index.size; ++i) {
                writeU16(second_fout, *ptr);
                ++ptr;
        }

        second_fout.flush();
        second_fout.close();
}


} /* namespace sbwt */

void TestWriteBinary(vector<int> &nums, const string &file1) {
        ofstream fout1(file1.c_str(), ios::binary);
        static bool is_bigendian = currentlyBigEndian();
        /* write into the size of vector */
        writeI32(fout1, (uint32_t) is_bigendian/**/);
        writeU32(fout1, (uint32_t)nums.size());
        for (auto i : nums) {
                writeI32(fout1, i, is_bigendian);
        }

        fout1.flush();
        fout1.close();

        return;
}

void TestReadBinary(const string &file_name)
{
        ifstream fin(file_name.c_str(), ios_base::in | ios::binary);
        if (!fin.is_open()) {
                cout << "Cannot open " << file_name << endl;
                return;
        }
        static bool is_bigendian = currentlyBigEndian();
        size_t size_vec = readU32(fin, is_bigendian);

        vector<int> nums(size_vec);
        for (auto &i : nums) {
                i = readI32(fin, is_bigendian);
                cout << i << " ";
        }
        fin.close();
        return;
}

uint32_t CountDna(char *file_name, const bool &count_N = true) {
        vector<istream* > istream_vec;

        RefReadInParams ref_readin_params(-1, -1, false, true);

        istream_vec.push_back(new ifstream(file_name, ios::in));

        vector<RefRecord> ref_record_vec;

        uint32_t size_total_seq = 0;
        size_total_seq = fastaRefReadSizes(istream_vec, ref_record_vec, ref_readin_params);

        return size_total_seq;
}

int TestReadWriteBinaryData(int argc, char **argv)
{
        vector<int> nums = {1,2,3,4,1};
        string file_name = "vec.sbwt.bin";

        TestWriteBinary(nums, file_name);
        TestReadBinary(file_name);

        return 0;
}

void ReadFasta(char **argv)
{
        fstream file_stream;
        file_stream.open(argv[1], ios_base::in/* un-compressed */);

        String<char> fasta_tag;
        String<Dna> fasta_seq;

        readMeta(file_stream, fasta_tag, Fasta());
        {

                cout << fasta_tag << "\n";

                read(file_stream, fasta_seq, Fasta());

                cout << fasta_seq << "\n";

        }
        file_stream.close();

        return;
}
