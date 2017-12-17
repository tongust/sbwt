#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cassert>

#include <ctype.h>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>

#include "test_io.h"
#include "alphabet.h"
#include "ref_read.h"
#include "endian_swap.h"
#include "word_io.h"

using namespace seqan;
using namespace std;

void PrintMsg();

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

int main(int argc, char **argv)
{
        /*
        if (argc != 2) {
                PrintMsg();
                return 1;
        }
        cout << "total sizes: " << CountDna(argv[1], true) << endl;
        */
        vector<int> nums = {1,2,3,4,1};
/*        for (auto &i : nums) {
                cin >> i;
        }
*/
        string file_name = "vec.sbwt.bin";
        TestWriteBinary(nums, file_name);
        TestReadBinary(file_name);

        return 0;
}
void PrintMsg()
{
        cout << "Usage: [fasta file]\n";
        return;
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
