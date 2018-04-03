#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cstdint>
#include <iostream>
#include <bitset>
#include <memory>
#include <chrono>

#include "sbwt_search.h"
#include "sbwt.h"
#include "sequence_pack.h"
#include "log.h"
#include "alphabet.h"
#include "utility.h"

namespace sbwt
{
        const uint32_t FLAG32[32] = {
                        1,      1<<1,   1<<2,   1<<3,   1<<4,   1<<5,   1<<6,   1<<7,
                        1<<8,   1<<9,   1<<10,  1<<11,  1<<12,  1<<13,  1<<14,  1<<15,
                        1<<16,  1<<17,  1<<18,  1<<19,  1<<20,  1<<21,  1<<22,  1<<23,
                        1<<24,  1<<25,  1<<26,  1<<27,  1<<28,  1<<29,  1<<30,  1<<31
        };

        /**
         * Construction of bitset64
         * @param rhs
         */
        bitset64::bitset64(const bitset64 &rhs): flag(rhs.flag)
        {
                uint64_t *p = array;
                for (auto i : rhs.array) {
                        *p = i;
                        ++p;
                }
        }

        bitset64& bitset64::operator=(bitset64& rhs)
        {
                if (&array[0] == &rhs.array[0]) {
                        return *this;
                }
                flag = rhs.flag;

                uint64_t *p = array;
                for (auto i : rhs.array) {
                        *p = i;
                        ++p;
                }
                return *this;
        }

        /// Shift bits
        /// right must be not greater than 64
        void bitset64::RightShift(uint32_t right, uint32_t size)
        {
                right >>= 1;/// divided by 2
                if (flag & (1<<right)) {
                        return;
                } else {
                        flag |= 1<<right;
                }

                uint64_t *p = array + (right<<5);

                right <<= 1;
                uint32_t left = 64 - right;

                uint64_t *q     = array;
                uint64_t *q_end = array + size;

                *p = (*q)>>right;
                ++q;
                ++p;
                for (; q != q_end; ++q) {
                        *p = (*(q-1))<<left | (*q)>>right;
                        ++p;
                }
                *p = (*(q-1))<<left;
        }

        /// Shift DNA string
        void bitset64::RightShiftDna(uint32_t right, uint32_t size)
        {
                if (flag & (1<<right)) {
                        return;
                } else {
                        flag |= 1<<right;
                }

                uint64_t *p = array + (right<<5);

                right <<= 1;
                uint32_t left = 64 - right;

                uint64_t *q     = array;
                /// TODO SEG ERROR
                uint64_t *q_end = array + size;

                *p = (*q)<<right;
                ++q;
                ++p;
                for (; q != q_end; ++q) {
                        *p = (*(q-1))>>left | (*q)<<right;
                        ++p;
                }
                *p = (*(q-1))>>left;
        }

        void bitset64::clear()
        {
                flag = 1;
        }

        uint32_t bitset64::size()
        {
                return MAX_READS_BINARY_LENGTH * 32;
        }

#ifdef SBWT_VERBOSE
        using std::cout;
        using std::endl;
        using std::shared_ptr;

        void bitset64::Print()
        {
                using std::bitset;
                using std::cout;
                using std::endl;
                using std::string;

                cout << "[-]\t";
                for (int i = 0; i != 32; ++i) {
                        cout << "[" << i << "]" << (i<10 ? "  " : " ");
                        for (int j = 1; j != 16; ++j) {
                                cout << j << (j < 10 ? "    " : "   ");
                        }
                        if (i != 31)
                                cout << " |  ";
                } cout << endl;
                auto print_str = [](const uint64_t &t)->void {
                        auto s = bitset<64>(t).to_string();
                        for (int i = 0; i != 64; ++i) {
                                if (i && i % 4 == 0 && i != 63) {
                                        cout << " ";
                                }
                                cout << s[i];
                        }
                };
                for (int i = 0; i != size(); ++i) {
                        if ( (i % MAX_READS_BINARY_LENGTH) == 0) {
                                cout << "[" << (i>>5) << "]\t";
                                print_str(array[i]);
                        } else if ((i+1) % MAX_READS_BINARY_LENGTH == 0){
                                cout << "  |  ";
                                print_str(array[i]);
                                cout << endl;
                        } else {
                                cout << "  |  ";
                                print_str(array[i]);
                        }
                }
        }

        static void PrintBit64Index()
        {
                for (int i = 0; i != 64; ++i) {
                        if (i % 4) continue;
                        if (i < 10)
                                cout << i << "    ";
                        else
                                cout << i << "   ";
                } cout << "\n";
        }

        static void PrintBit64Binary(uint64_t val, const string &msg = "")
        {
                string s = bitset<64>(val).to_string();
                int i = 0;
                for (auto c : s) {
                        cout << c;
                        if (i%4 == 3 && i != 63) {
                                cout << " ";
                        }
                        ++i;
                }
                if (!msg.empty())
                        cout << "------" << msg << endl;
        }

        static void PrintBit64Char(uint64_t val, const string &msg = "")
        {
                static char bit2char[2][2] = {{'A', 'C'}, {'T', 'G'}};
                static string dnas('A', 32);
                string s = bitset<64>(val).to_string();
                for (int i = 0; i != 32; ++i) {
                        dnas[i] = bit2char[s[i*2] == '1'][s[i*2+1] == '1'];
                }
                for (int i = 0; i != 16; ++i) {
                        cout << dnas[i*2] << " " << dnas[i*2 + 1] << "  ";
                }
                if (!msg.empty())
                        cout << "-----" << msg << endl;
        }

        static void PrintBit64WithIndex(uint64_t val, const string &msg = "")
        {
                PrintBit64Index();
                PrintBit64Binary(val, msg);
        }

        static void PrintBit64CharWithIndex(uint64_t val, const string &msg = "")
        {
                PrintBit64Index();
                PrintBit64Char(val, msg);
        }


        void bitset64::Print(uint64_t* bin_seq, uint32_t size)
        {
                using std::bitset;
                using std::cout;
                using std::endl;
                using std::string;

                uint64_t *end = bin_seq + size;
                for (uint64_t *it = bin_seq; it != end; ++it) {
                        cout << bitset<64>(*it);
                        if (it+1 != end) {
                                cout << " | ";
                        } else {
                                cout << endl;
                        }
                }
        }

        /// Test bitset64::RightShift(uint32_t, uint32_t);
        void TestBitShift()
        {
                using sbwt::bitset64;
                using std::bitset;
                using std::cout;
                using std::endl;
                using std::string;

                auto PrintBitStr = [](const string &str)->void {
                        auto size = str.size();
                        for (size_t i = 0; i != size; ++i) {
                                if (i % 8 == 0 && i && i+1 != size) {
                                        cout << " - ";
                                }
                                cout << str[i];
                        }
                        cout << endl;
                };

                string bs = "01010101001010101000100001010100101011111111111111100100111111111111110001010011111111110001000000011110000000111111100011110011";
                string bs0(bs.begin(), bs.begin() + 64),
                                bs1(bs.begin()+64, bs.end());
                uint64_t array[2];
                array[0] = bitset<64>(bs0).to_ullong();
                array[1] = bitset<64>(bs1).to_ullong();

                bitset<128> b2(bs);
                PrintBitStr(b2.to_string());

                bitset64 b64(array, 2);

                for (uint32_t i = 0; i != 10000; ++i) {
                        uint32_t t = i % 64;
                        b64.RightShift(t, 2);
                }

                b64.Print();
        }

        static const uint64_t DnaStringRightShiftMask[32] = {
                ~0ULL,                   ~0x3ULL,                 ~0xfULL,                 ~0x3fULL,
                ~0xffULL,                ~0x3ffULL,               ~0xfffULL,               ~0x3fffULL,
                ~0xffffULL,              ~0x3ffffULL,             ~0xfffffULL,             ~0x3fffffULL,
                ~0xffffffULL,            ~0x3ffffffULL,           ~0xfffffffULL,           ~0x3fffffffULL,
                ~0xffffffffULL,          ~0x3ffffffffULL,         ~0xfffffffffULL,         ~0x3fffffffffULL,
                ~0xffffffffffULL,        ~0x3ffffffffffULL,       ~0xfffffffffffULL,       ~0x3fffffffffffULL,
                ~0xffffffffffffULL,      ~0x3ffffffffffffULL,     ~0xfffffffffffffULL,     ~0x3fffffffffffffULL,
                ~0xffffffffffffffULL,    ~0x3ffffffffffffffULL,   ~0xfffffffffffffffULL,   ~0x3fffffffffffffffULL
        };
        const uint64_t DnaStringRightShiftMaskReverse[32] = {
                0ULL,                   0x3ULL,                 0xfULL,                 0x3fULL,
                0xffULL,                0x3ffULL,               0xfffULL,               0x3fffULL,
                0xffffULL,              0x3ffffULL,             0xfffffULL,             0x3fffffULL,
                0xffffffULL,            0x3ffffffULL,           0xfffffffULL,           0x3fffffffULL,
                0xffffffffULL,          0x3ffffffffULL,         0xfffffffffULL,         0x3fffffffffULL,
                0xffffffffffULL,        0x3ffffffffffULL,       0xfffffffffffULL,       0x3fffffffffffULL,
                0xffffffffffffULL,      0x3ffffffffffffULL,     0xfffffffffffffULL,     0x3fffffffffffffULL,
                0xffffffffffffffULL,    0x3ffffffffffffffULL,   0xfffffffffffffffULL,   0x3fffffffffffffffULL
        };

        void TestNaiveMap(const std::string &file_fasta, const std::string &file_index)
        {
                reads_buffer rb_fasta(file_fasta);
                reads_buffer rb_index(file_index);

                /// reads reference file
                rb_index.ReadNext();
                rb_index.ReadNext();

                const uint32_t size_read = 150;

                uint32_t diff_count = 0;

                char *p = nullptr;
                char *p_end = nullptr;
                char *q = nullptr;

                for (uint32_t i = 0;;++i) {

                        q = rb_index.buffer + i;
                        /**
                         * read head
                         */
                        rb_fasta.ReadNext();
                        if (rb_fasta.length_read < 1) break;

                        /**
                         * read DNAs
                         */
                        rb_fasta.ReadNext();
                        diff_count = 0;

                        for (int k = 0; k != 20; ++k) {
                                p = rb_fasta.buffer;
                                p_end = p + size_read;

                                for (; p != p_end; ++p, q++) {
                                        if (*p != *q) {
                                                ++diff_count;
                                        }
                                }
                        }
                        cout << diff_count;
                }
        }

        void TestBit8Map(const std::string &file_fasta, const std::string &file_index)
        {
                reads_buffer rb_fasta(file_fasta);
                reads_buffer rb_index(file_index);

                /// read reference file
                rb_index.ReadNext();
                rb_index.ReadNext();

                /// the size of reference sequence will be discarded to be 4^n.
                uint32_t size_ref_8bit = rb_index.length_read >> 2;
                --size_ref_8bit;
                uint32_t size_ref_char = size_ref_8bit << 2;

                shared_ptr<uint8_t > ref_bin_sptr(new uint8_t[size_ref_char]);
                shared_ptr<uint8_t *> ref_bin_sptr_array(new uint8_t*[4]);

                uint8_t *ref_bin_ptr = ref_bin_sptr.get();
                uint8_t **ref_bin_ptr_array = ref_bin_sptr_array.get();

                ref_bin_ptr_array[0] = ref_bin_ptr + 0*size_ref_8bit;
                ref_bin_ptr_array[1] = ref_bin_ptr + 1*size_ref_8bit;
                ref_bin_ptr_array[2] = ref_bin_ptr + 2*size_ref_8bit;
                ref_bin_ptr_array[3] = ref_bin_ptr + 3*size_ref_8bit;

                sbwtio::BaseChar2Binary8B(rb_index.buffer+0, size_ref_8bit, ref_bin_ptr_array[0]);
                sbwtio::BaseChar2Binary8B(rb_index.buffer+1, size_ref_8bit, ref_bin_ptr_array[1]);
                sbwtio::BaseChar2Binary8B(rb_index.buffer+2, size_ref_8bit, ref_bin_ptr_array[2]);
                sbwtio::BaseChar2Binary8B(rb_index.buffer+3, size_ref_8bit, ref_bin_ptr_array[3]);

                uint64_t read_bin_buffer[256] = {0};
                uint8_t *read_bin = (uint8_t*)read_bin_buffer;
                uint32_t size_read_char = 150;
                uint32_t size_read_mod4 = size_read_char % 4;
                uint32_t size_read_bit8 = (size_read_mod4) ? (size_read_char >> 2) + 1 : size_read_char >> 2;

                uint32_t size_read_mod32 = size_read_char % 32;
                uint32_t size_read_bit32 = size_read_mod32 ?
                                           (size_read_char >> 5) + 1:
                                           size_read_char >> 5;
                uint32_t size_read_bit32_1 = size_read_bit32 - 1;

                uint64_t popcount_mask = DnaStringRightShiftMaskReverse[size_read_mod32];


                uint64_t *p = nullptr;
                uint64_t *p_end = nullptr;
                uint64_t *q = read_bin_buffer;
                uint64_t *q_end = nullptr;
                if (size_read_mod32) {
                        q_end = q + size_read_bit32_1;
                } else {
                        q_end = q + size_read_bit32;
                }
                uint8_t *b = nullptr;
                uint64_t v = 0;
                for (uint32_t i = 0; ; ++i) {
                        rb_fasta.ReadNext();
                        if (rb_fasta.length_read < 1) {
                                break;
                        }
                        rb_fasta.ReadNext();

                        static int diff_count = 0;
                        diff_count = 0;
                        sbwtio::BaseChar2Binary8B(rb_fasta.buffer, size_read_bit8, read_bin);

                        for (int k = 0; k != 20; ++k) {
                                b = ref_bin_ptr_array[i % 4] + (i >> 2);
                                p = (uint64_t *) b;
                                q = read_bin_buffer;/// Bug Here

                                if (size_read_mod32) {
                                        p_end = p + size_read_bit32_1;
                                        v = (*p_end) ^ (*q_end);
                                        v &= popcount_mask;
                                        diff_count += HammingWeightDna64(v);
                                } else {
                                        p_end = p + size_read_bit32;
                                }
                                for (; p != p_end;) {
                                        v = (*p++) ^ (*q++);
                                        diff_count += HammingWeightDna64(v);
                                }
                        }
                        cout << diff_count;
                }
        }

        void TestBit64Map(const std::string &file_fasta, const std::string &file_index)
        {
                reads_buffer rb_fasta(file_fasta);
                reads_buffer rb_index(file_index);

                /// reads reference file
                rb_index.ReadNext();
                rb_index.ReadNext();

                uint32_t size_ref_char = rb_index.length_read;
                uint32_t size_ref_64bit = (size_ref_char >> 5) + 1;
                uint32_t size_read_char = 150;
                uint32_t size_read_char_mod32 = size_read_char & 31;
                uint32_t size_read_64bit = (150/32)+1;

                shared_ptr<uint64_t > ref_bin_sptr(new uint64_t[size_ref_64bit]);
                uint64_t *ref_bin_ptr = ref_bin_sptr.get();

                sbwtio::BaseChar2Binary64B(rb_index.buffer, size_ref_char, ref_bin_ptr);

                /// print reference sequence
                //bitset64::Print(ref_bin_ptr, size_read_64bit);

                bitset64 bt;

                uint32_t index = 0;

                uint64_t v2 = 0;

                int diff_count = 0;

                uint64_t *p = nullptr;
                uint64_t *p_end = nullptr;
                uint64_t *q = nullptr;

                if (size_read_char < 32) {
                        logger::LogError("Size of read is less than 32");
                        return;
                }

                for (uint32_t i = 0;;++i) {
                        /**
                         * read head
                         */
                        rb_fasta.ReadNext();
                        if (rb_fasta.length_read < 1) break;

                        /**
                         * read DNAs
                         */
                        rb_fasta.ReadNext();

#if 1
                        sbwtio::BaseChar2Binary64B(rb_fasta.buffer, rb_fasta.length_read, bt.array);
                        bt.clear();

                        /// Profile?
                        //index = i % 32;
                        index = i & 31;

                        diff_count = 0;

                        if (index != 0) {
                                p = bt.array + (index<<5);
                                bt.RightShiftDna(index, size_read_char);

                                uint32_t tmp = size_read_char + index;
                                tmp = (tmp & 31 == 0) ? tmp >> 5 : (tmp>>5)+1;
                                p_end = p + tmp;////// Here. Bugs!
                                q = ref_bin_ptr + (i >> 5);

                                /// the first segment
                                v2 = (*p) ^ (*q);
                                v2 &= DnaStringRightShiftMask[index];
                                diff_count += HammingWeightDna64(v2);
                                ++p;
                                ++q;

                                /// the last segment
                                index = (index + size_read_char) & 31;

                                if (index != 0) {
                                        --p_end;
                                        for (; p != p_end;) {
                                                v2 = (*p) ^(*q);
                                                diff_count += HammingWeightDna64(v2);
                                                ++p;
                                                ++q;
                                        }

                                        ++p_end;
                                        v2 = (*p) ^ (*q);
                                        v2 &= DnaStringRightShiftMaskReverse[index];
                                        diff_count += HammingWeightDna64(v2);
                                } else {
                                        --p_end;/// Bug here
                                        for (; p != p_end;) {
                                                v2 = (*p) ^(*q);
                                                diff_count += HammingWeightDna64(v2);
                                                ++p;
                                                ++q;
                                        }
                                }
                        } else {
                                p = bt.array;

                                p_end = p + size_read_64bit;
                                q = ref_bin_ptr + (i >> 5);

                                /// the last segment
                                if (size_read_char_mod32 != 0) {
                                        --p_end;
                                        for (; p != p_end;) {
                                                v2 = (*p) ^(*q);
                                                diff_count += HammingWeightDna64(v2);
                                                ++p;
                                                ++q;
                                        }

                                        ++p_end;
                                        v2 = (*p) ^ (*q);
                                        v2 &= DnaStringRightShiftMaskReverse[size_read_char_mod32];
                                        diff_count += HammingWeightDna64(v2);
                                } else {
                                        for (; p != p_end;) {
                                                v2 = (*p) ^(*q);
                                                diff_count += HammingWeightDna64(v2);
                                                ++p;
                                                ++q;
                                        }
                                }
                        }
#endif
                        cout << diff_count;
                        //bt.Print();
                }

        }

        /// Test: read the dnas string and transform them into bits reverse-complemently.
        void TestBaseChar2BinaryReverseComplement(char **argv)
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

                /// method 1
                sbwtio::BaseChar2Binary8B(rb_index.buffer, size_ref_8bit, ref_bin_ptr);
                //sbwtio::BaseChar2Binary8B_RC(ref_bin_ptr, size_ref_char, ref_bin_rc_ptr);

                /// method 2
                //sbwtio::BaseChar2Binary8B_RC(rb_index.buffer, size_ref_char, ref_bin_rc_ptr1);

                uint32_t size_char = rb_index.length_read;
                uint32_t mod_4_by2 = (size_char & 3) << 1;
                uint32_t mod_4_r_by2 = 8 - mod_4_by2;
                uint32_t size_8bit = mod_4_by2 ? (size_char>>2) + 1 : size_char >> 2;
                //sbwtio::BaseChar2Binary8B_RC_Exter(rb_index.buffer, ref_bin_rc_ptr1, size_char, size_8bit, mod_4_by2, mod_4_r_by2);


                //for (uint32_t i = 0; i != size_8bit; ++i) { cout << bitset<8>(ref_bin_rc_ptr[i]) << ", "; } cout << endl;
                //for (uint32_t i = 0; i != size_8bit; ++i) { cout << bitset<8>(ref_bin_rc_ptr1[i]) << ", "; } cout << endl;

                for (uint32_t i = 1<<24; i != 0; --i) {
                        /// method 1
                        ///sbwtio::BaseChar2Binary8B_RC(ref_bin_ptr, size_ref_char, ref_bin_rc_ptr);
                        /// method 2
                        //sbwtio::BaseChar2Binary8B_RC(rb_index.buffer, size_ref_char, ref_bin_rc_ptr1);
                        /// method 3
                        sbwtio::BaseChar2Binary8B_RC_Exter(rb_index.buffer, ref_bin_rc_ptr1, size_char, size_8bit, mod_4_by2, mod_4_r_by2);
                }
        }

#endif /* SBWT_VERBOSE */

        bitset64::bitset64(uint64_t *nums, uint32_t size): flag(1)
        {
                uint64_t *p = array,
                         *nums_end = nums + size;
                for (; nums != nums_end; ++nums, ++p) {
                        *p = *nums;
                }
        }

        reads_buffer::reads_buffer():
                length_buffer(MAX_READS_LENGTH+1),
                length_read(0),
                file_stream(nullptr)
        {
                Init();
        }

        reads_buffer::reads_buffer(char *fn, uint32_t size):
                length_buffer(MAX_READS_LENGTH+1),
                length_read(0)
        {
                file_name.resize(size);
                std::copy(fn, fn+size, file_name.begin());
                Init();
        }

        reads_buffer::reads_buffer(const string &fn):
                length_buffer(MAX_READS_LENGTH+1),
                length_read(0)
        {
                file_name = fn;
                Init();
        }

        reads_buffer::reads_buffer(const reads_buffer& rb)
        {
                if (buffer == rb.buffer) return;
                this->Dispose();
                /// TODO same file description?
                file_name = rb.file_name;
                Init();
        }

        void reads_buffer::Init()
        {
                try {
                        buffer = (char*) malloc(length_buffer);
                        if (buffer == nullptr) {
                                logger::LogError("malloc buffer");
                        }
                        ptr_buffer              = &buffer;
                        ptr_length_buffer       = &length_buffer;

                        if (!file_name.empty()) {
                                file_stream = fopen(file_name.c_str(), "r");

                                if (!file_stream) {
                                        logger::LogError("fopen");
                                }
                        }
                }
                catch (...) {
                        logger::LogError("Construct reads_buffer");
                        throw;
                }
        }

        void reads_buffer::Dispose()
        {
                try {
                        free(buffer);
                        ptr_length_buffer       = nullptr;
                        ptr_buffer              = nullptr;
                        buffer                  = nullptr;

                        if (file_stream) {
                                if (fclose(file_stream) != 0) {
                                        logger::LogError("close file");
                                }
                        }
                }
                catch (...) {
                        logger::LogError("Dispose reads_buffer");
                }
        }

        reads_buffer::~reads_buffer()
        {
                Dispose();
        }

        ssize_t reads_buffer::Size() { return length_read; }

        /**
         * ssize_t getline(char **lineptr, size_t *n, FILE *stream);
         *
         * "If *lineptr is set to NULL and *n is set 0 before the call, then
         * getline() will allocate a buffer for storing the line.  This buffer
         * should be freed by the user program even if getline() failed."
         *
         * Reference: http://man7.org/linux/man-pages/man3/getline.3.html
         * @param filename
         */
        void reads_buffer::ReadNext()
        {
                length_read = getline(ptr_buffer, ptr_length_buffer, file_stream);
                if (length_read == -1)
                        return;
                /*
                if (length_read == -1) {
                        logger::LogError("ReadNext");
                        return;
                }
                */
                --length_read;
        }


        void UnsortedUnpackedSearch(char **argv)
        {
                /// Build index from files
                string reads_filename = string(argv[1]);
                string prefix_filename = string(argv[2]);
                BuildIndexRawData build_index(prefix_filename);
                //sbwt::PrintFullSearchMatrix(build_index);

                uint32_t period = build_index.period;
                uint32_t N = build_index.length_ref;
                uint32_t *C = build_index.first_column;
                uint32_t **Occ = build_index.occurrence;
                uint32_t *SA = build_index.suffix_array;

                static uint8_t *ref_bin_ptr_array[4] = {nullptr};
                ref_bin_ptr_array[0] = build_index.bin_8bit;
                ref_bin_ptr_array[1] = build_index.bin_8bit + 1*build_index.size_bin_8bit;
                ref_bin_ptr_array[2] = build_index.bin_8bit + 2*build_index.size_bin_8bit;
                ref_bin_ptr_array[3] = build_index.bin_8bit + 3*build_index.size_bin_8bit;

                /// reads
                /// Read first segment.
                reads_buffer rb_reads(reads_filename);
                rb_reads.ReadNext();/// header
                if (rb_reads.length_read < 1) {
                        LOGERROR("Empty reads");
                        return;
                }
                char head_str[1024*1024] = {'\0'};
                strcpy(head_str, rb_reads.buffer);
                /// strip the new line char.
                head_str[rb_reads.length_read] = '\0';

                rb_reads.ReadNext();/// sequence
                if (rb_reads.length_read < period) {
                        LOGERROR("The length of single segment is less than the period");
                        return;
                }

                // Assuming the length of single sequence segment is kept same.
                uint32_t size_read_char = rb_reads.length_read;
                uint32_t size_read_mod4 = size_read_char % 4;
                uint32_t size_read_bit8 = (size_read_mod4) ? (size_read_char >> 2) + 1 : size_read_char >> 2;
                uint32_t size_read_mod32 = size_read_char % 32;
                uint32_t size_read_bit32 = size_read_mod32 ?
                                           (size_read_char >> 5) + 1:
                                           size_read_char >> 5;
                uint32_t size_read_bit32_1 = size_read_bit32 - 1;
                uint64_t popcount_mask = DnaStringRightShiftMaskReverse[size_read_mod32];

                /// binary buffer for segment
                uint64_t read_bin_buffer[256] = {0};
                uint8_t *read_bin = (uint8_t*)read_bin_buffer;
                uint64_t *p = nullptr;
                uint64_t *p_end = nullptr;
                uint64_t *q = read_bin_buffer;
                uint64_t *q_end = nullptr;
                if (size_read_mod32) {
                        q_end = q + size_read_bit32_1;
                } else {
                        q_end = q + size_read_bit32;
                }

                uint32_t L = 0, R = 0;
                uint32_t N_1 = N - 1;

                // Not sure the addr. of buffer will be changed.
                // ---
                // Alternatively, before calling getline(), *lineptr can contain a
                // pointer to a malloc(3)-allocated buffer *n bytes in size.  If the
                // buffer is not large enough to hold the line, getline() resizes it
                // with realloc(3), updating *lineptr and *n as necessary.
                // ---
                // Reference: http://man7.org/linux/man-pages/man3/getline.3.html

                static char key = 0;
                static uint32_t a = 0;
                uint32_t lmd = 0;
                uint32_t *end_array = new uint32_t[period];
                char *ptr = rb_reads.buffer + (size_read_char - 1);
                char *ptr_reads = rb_reads.buffer;
                char *ptr_reads_end = ptr_reads + size_read_char;
                char *ptr_ref = build_index.seq_raw;
                char **ptr_array = new char*[period];
                uint32_t *begin_index = new uint32_t[period];
                uint32_t index_tmp, index;

                uint32_t *psa, *psa_end;

                {
                        uint32_t tmp0 = size_read_char / period;
                        uint32_t tmp1 = size_read_char % period;

                        for (uint32_t i = 0; i != period; ++i) {
                                end_array[i] = tmp0 + (i < tmp1);
                                ptr_array[i] = ptr - i;
                                begin_index[i] = size_read_char - ((end_array[i] - 1)*period + i + 1);
                        }
                }

                uint32_t count  = 0;
                char *align_p = nullptr;
                char *align_q = nullptr;

/// break-version
#if 0
                bool flag_map = false;
                for(;;) {
                        flag_map = false;

                        for (uint32_t i = 0; i != period; ++i)
                        {
                                lmd = end_array[i];
                                ptr = ptr_array[i];

                                // Init
                                key = *ptr;
                                a = charToDna5_32bit[key];
                                L = C[a];
                                R = C[a] + Occ[a][N_1] - 1;

                                for (uint32_t j = 1; j != lmd; ++j) {
                                        ptr -= period;
                                        key = *ptr;
                                        a = charToDna5_32bit[key];
                                        L = C[a] + Occ[a][L-1];
                                        R = C[a] + Occ[a][R] - 1;
                                        if (L > R) {
                                                break;
                                        }
                                }

                                if (L > R) {
                                        continue;
                                }

                                psa = SA + L;
                                psa_end = SA + R + 1;
                                index_tmp = begin_index[i];

                                for (; psa != psa_end; ++psa) {
                                        index = *psa - index_tmp;
                                        count  = 0;

                                        // Brute force
                                        align_p = ptr_reads;
                                        align_q = ptr_ref + index;
                                        flag_map = true;

                                        for (; align_p != ptr_reads_end;) {
                                                if (*align_p++ != *align_q++) {
                                                        if (++count == period) {
                                                                flag_map = false;
                                                                break;
                                                        }
                                                }
                                        }
                                        if (flag_map) {
                                                break;
                                        }
                                }
                                if (flag_map) {
                                        //cout << "matched: " << index << "\n";
                                        break;
                                }
                        }/* i */

                        if (flag_map) {
                                cout << "matched\t" << index << "\n";
                        } else {
                                cout << "*\n";
                        }

                        rb_reads.ReadNext();/// header
                        if (rb_reads.length_read < 1) {
                                LOGEINFO("Empty reads");
                                break;
                        }
                        rb_reads.ReadNext();/// sequence

                } /* Readings */
#endif /* break-version */

/// goto-version
#if 1
                char *p_end_rc = nullptr;
                char *p_rc = nullptr;
                char *p0_rc = nullptr;
                char tmp_rc;
                char flag_minus_plus = '+';
                {
                        uint32_t i;

                        for(;;) {
                                flag_minus_plus = '+';
                                i = 0;
                                for (;;)
                                {
                                        lmd = end_array[i];
                                        ptr = ptr_array[i];

                                        /// Init
                                        key = *ptr;
                                        a = charToDna5_32bit[key];
                                        L = C[a];
                                        R = C[a] + Occ[a][N_1] - 1;

                                        for (uint32_t j = 1; j != lmd; ++j) {
                                                ptr -= period;
                                                key = *ptr;
                                                a = charToDna5_32bit[key];
                                                L = C[a] + Occ[a][L-1];
                                                R = C[a] + Occ[a][R] - 1;
                                                if (L > R) {
                                                        goto loop_verification;
                                                }
                                        }

                                        psa = SA + L;
                                        psa_end = SA + R + 1;
                                        index_tmp = begin_index[i];

                                        for (;;) {
                                                psa_loop:
                                                if (psa == psa_end) {
                                                        goto loop_verification;
                                                }

                                                index = *psa - index_tmp;
                                                count  = 0;

                                                // Brute force
                                                align_p = ptr_reads;
                                                align_q = ptr_ref + index;

                                                for (; align_p != ptr_reads_end;) {
                                                        if (*align_p++ != *align_q++) {
                                                                if (++count == period) {
                                                                        ++psa;// Careful about its position.
                                                                        goto psa_loop;
                                                                }
                                                        }
                                                } // align p and q

                                                goto match_success;

                                        } // psa

                                        /// loop verification
                                        loop_verification:
                                        if (++i == period) {
                                                /// Still not matched.
                                                goto match_failure;
                                        }
                                }/* i */

                                /// failed in forward searching
                                match_failure:
                                /// Reverse Complement
                                p0_rc = rb_reads.buffer;
                                p_end_rc = p0_rc + ( (size_read_char+1) >> 1 );
                                p_rc = p0_rc + (size_read_char - 1);
                                for (; p0_rc != p_end_rc; ) {
                                        tmp_rc = sbwtio::DnaCharMapReverseComplement[*p0_rc];
                                        *p0_rc++ = sbwtio::DnaCharMapReverseComplement[*p_rc];
                                        *p_rc-- = tmp_rc;
                                }

                                flag_minus_plus = '-';
                                i = 0;
                                for (;;)
                                {
                                        lmd = end_array[i];
                                        ptr = ptr_array[i];

                                        // Init
                                        key = *ptr;
                                        a = charToDna5_32bit[key];
                                        L = C[a];
                                        R = C[a] + Occ[a][N_1] - 1;

                                        for (uint32_t j = 1; j != lmd; ++j) {
                                                ptr -= period;
                                                key = *ptr;
                                                a = charToDna5_32bit[key];
                                                L = C[a] + Occ[a][L-1];
                                                R = C[a] + Occ[a][R] - 1;
                                                if (L > R) {
                                                        goto loop_verification_rc;
                                                }
                                        }

                                        psa = SA + L;
                                        psa_end = SA + R + 1;
                                        index_tmp = begin_index[i];

                                        for (;;) {
                                                psa_loop_rc:
                                                if (psa == psa_end) {
                                                        goto loop_verification_rc;
                                                }

                                                index = *psa - index_tmp;
                                                count  = 0;

                                                // Brute force
                                                align_p = ptr_reads;
                                                align_q = ptr_ref + index;

                                                for (; align_p != ptr_reads_end;) {
                                                        if (*align_p++ != *align_q++) {
                                                                if (++count == period) {
                                                                        ++psa;// Careful about its position.
                                                                        goto psa_loop_rc;
                                                                }
                                                        }
                                                } // align p and q

                                                goto match_success;

                                        } // psa

                                        /// loop verification
                                        loop_verification_rc:
                                        if (++i == period) {
                                                /// Still not matched.
                                                goto match_failure_rc;
                                        }
                                }/* i */

                                match_failure_rc:
                                rb_reads.ReadNext();/// header
                                if (rb_reads.length_read < 1) {
                                        LOGERROR("Empty reads");
                                        break;
                                }
                                strcpy(head_str, rb_reads.buffer);
                                /// strip the new line char.
                                head_str[rb_reads.length_read] = '\0';

                                rb_reads.ReadNext();/// sequence
                                continue;

                                match_success:
                                cout << head_str << "\t"
                                     << flag_minus_plus << "\tchr1\t"
                                     << index << "\t" << rb_reads.buffer;
                                rb_reads.ReadNext();/// header
                                if (rb_reads.length_read < 1) {
                                        LOGERROR("Empty reads");
                                        break;
                                }
                                strcpy(head_str, rb_reads.buffer);
                                /// strip the new line char.
                                head_str[rb_reads.length_read] = '\0';

                                rb_reads.ReadNext();/// sequence
                        } /* Readings */
                }
#endif/* goto-version */

                delete[] end_array;
                delete[] ptr_array;
                delete[] begin_index;
        }


        void UnsortedPackedSearch(char **argv)
        {
                /// Build index from files
                auto begin_time_index = std::chrono::high_resolution_clock::now();

                string reads_filename = string(argv[1]);
                string prefix_filename = string(argv[2]);
                BuildIndexRawData build_index(prefix_filename);

                auto end_time_index = std::chrono::high_resolution_clock::now();
                auto duration_index = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_index-begin_time_index).count();
                LOGINFO( "Load index elapsed time:\t"
                                << duration_index
                                << " ms"
                                << std::endl);
                //sbwt::PrintFullSearchMatrix(build_index);

                auto begin_time_search = std::chrono::high_resolution_clock::now();

                uint32_t period = build_index.period;
                uint32_t N = build_index.length_ref;
                uint32_t *C = build_index.first_column;
                uint32_t **Occ = build_index.occurrence;
                uint32_t *SA = build_index.suffix_array;

                static uint8_t *ref_bin_ptr_array[4] = {nullptr};
                ref_bin_ptr_array[0] = build_index.bin_8bit;
                ref_bin_ptr_array[1] = build_index.bin_8bit + 1*build_index.size_bin_8bit;
                ref_bin_ptr_array[2] = build_index.bin_8bit + 2*build_index.size_bin_8bit;
                ref_bin_ptr_array[3] = build_index.bin_8bit + 3*build_index.size_bin_8bit;

                /// reads
                /// Read first segment.
                reads_buffer rb_reads(reads_filename);
                rb_reads.ReadNext();/// header
                if (rb_reads.length_read < 1) {
                        LOGINFO("Empty reads\n");
                        return;
                }
                char head_str[1024*1024] = {'\0'};
                strcpy(head_str, rb_reads.buffer);
                /// strip the new line char.
                head_str[rb_reads.length_read] = '\0';

                /// count matched segment
                uint64_t total_reads = 1;
                uint64_t matched_reads = 0;

                rb_reads.ReadNext();/// sequence
                if (rb_reads.length_read < period) {
                        LOGERROR("The length of single segment is less than the period");
                        return;
                }

                // Assuming the length of single sequence segment is kept same.
                uint32_t size_read_char = rb_reads.length_read;
                uint32_t size_read_mod4 = size_read_char % 4;
                uint32_t size_read_bit8 = (size_read_mod4) ? (size_read_char >> 2) + 1 : size_read_char >> 2;
                uint32_t size_read_mod32 = size_read_char % 32;
                uint32_t size_read_bit32 = size_read_mod32 ?
                                           (size_read_char >> 5) + 1:
                                           size_read_char >> 5;
                uint32_t size_read_bit32_1 = size_read_bit32 - 1;
                uint64_t popcount_mask = DnaStringRightShiftMaskReverse[size_read_mod32];

                /// binary buffer for segment
                uint64_t read_bin_buffer[256] = {0};
                uint64_t read_bin_buffer_rc[256] = {0};
                uint8_t *read_bin = (uint8_t*)read_bin_buffer;
                uint8_t *read_bin_rc = (uint8_t*)read_bin_buffer_rc;
                uint64_t *p = nullptr;
                uint64_t *p_end = nullptr;
                uint64_t *q = read_bin_buffer;
                uint64_t *q_end = nullptr;
                uint64_t *q_end_rc = nullptr;
                if (size_read_mod32) {
                        q_end = q + size_read_bit32_1;
                        q_end_rc = read_bin_buffer_rc + size_read_bit32_1;
                } else {
                        q_end = q + size_read_bit32;
                        q_end_rc = read_bin_buffer_rc + size_read_bit32;
                }

                uint32_t L = 0, R = 0;
                uint32_t N_1 = N - 1;

                // Not sure the addr. of buffer will be changed.
                // ---
                // Alternatively, before calling getline(), *lineptr can contain a
                // pointer to a malloc(3)-allocated buffer *n bytes in size.  If the
                // buffer is not large enough to hold the line, getline() resizes it
                // with realloc(3), updating *lineptr and *n as necessary.
                // ---
                // Reference: http://man7.org/linux/man-pages/man3/getline.3.html

                char key = 0;
                static uint32_t a = 0;
                uint32_t lmd = 0;
                uint32_t *end_array     = new uint32_t[period];
                char *ptr               = rb_reads.buffer + (size_read_char - 1);
                char *ptr_reads         = rb_reads.buffer;
                char **ptr_array        = new char*[period];
                char **ptr_array_rc     = new char*[period];
                uint32_t *begin_index   = new uint32_t[period];

                uint32_t index_tmp, index;

                uint32_t *psa, *psa_end;

                {
                        uint32_t tmp0 = size_read_char / period;
                        uint32_t tmp1 = size_read_char % period;

                        for (uint32_t i = 0; i != period; ++i) {
                                end_array[i]    = tmp0 + (i < tmp1);
                                ptr_array[i]    = ptr - i;
                                ptr_array_rc[i] = ptr_reads + i;
                                begin_index[i]  = size_read_char - ((end_array[i] - 1)*period + i + 1);
                        }
                }

                uint32_t count  = 0;

                uint8_t *b              = nullptr;
                uint64_t v = 0;

                char flag_minus_plus;

                uint32_t i;
                for(;;) {
                        /// forward
                        flag_minus_plus = '+';
                        sbwtio::BaseChar2Binary8B(ptr_reads, size_read_bit8, read_bin);
                        i = 0;
                        for (;;)
                        {
                                lmd = end_array[i];
                                ptr = ptr_array[i];

                                /// Init
                                key = *ptr;
                                a = charToDna5_32bit[key];
                                L = C[a];
                                R = C[a] + Occ[a][N_1] - 1;

                                for (uint32_t j = 1; j != lmd; ++j) {
                                        ptr -= period;
                                        key = *ptr;
                                        a = charToDna5_32bit[key];
                                        L = C[a] + Occ[a][L-1];
                                        R = C[a] + Occ[a][R] - 1;
                                        if (L > R) {
                                                goto loop_verification;
                                        }
                                }

                                psa = SA + L;
                                psa_end = SA + R + 1;
                                index_tmp = begin_index[i];

                                // packed method
                                if (size_read_mod32) {
                                        for (;;) {
                                                psa_loop_if:
                                                if (psa == psa_end) {
                                                        goto loop_verification;
                                                }
                                                /// Bug: it may happen that: *psa < index_tmp
                                                /// For example:
                                                ///     0    5
                                                /// Ref ACCCGTAGATTGAGTTGC...
                                                ///   0  3
                                                /// P TTACCCGTAGATTGAGTTGC...
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -= index_tmp;

                                                b = ref_bin_ptr_array[index & 3/*mod 4*/] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer;/// Bug Here

                                                p_end = p + size_read_bit32_1;
                                                v = (*p_end) ^ (*q_end);
                                                v &= popcount_mask;
                                                count = HammingWeightDna64(v);

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_if;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                } else {
                                        for (;;) {
                                                psa_loop_else:
                                                if (psa == psa_end) {
                                                        goto loop_verification;
                                                }
                                                // Bug here
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -= index_tmp;

                                                count  = 0;

                                                b = ref_bin_ptr_array[index & 3/* mod 4 */] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer;/// Bug Here

                                                p_end = p + size_read_bit32;

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_else;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                }

                                /// loop verification
                                loop_verification:
                                if (++i == period) {
                                        /// Still not matched.
                                        goto match_failure;
                                }
                        }/* i in [0, period) */

                        /// failed in forward searching
                        match_failure:
                        /// Reverse Complement Binary Stream
                        sbwtio::BaseChar2Binary8B_RC(read_bin, size_read_char, read_bin_rc);

                        i = 0;
                        flag_minus_plus = '-';
                        for (;;)
                        {
                                lmd = end_array[i];
                                ptr = ptr_array_rc[i];/// rc

                                /// Init
                                key = *ptr;
                                a = rcCharToDna5_32bit[key];/// rc
                                L = C[a];
                                R = C[a] + Occ[a][N_1] - 1;

                                for (uint32_t j = 1; j != lmd; ++j) {
                                        ptr += period;/// rc
                                        key = *ptr;
                                        /// Bugs here
                                        a = rcCharToDna5_32bit[key];/// rc
                                        L = C[a] + Occ[a][L-1];
                                        R = C[a] + Occ[a][R] - 1;
                                        if (L > R) {
                                                goto loop_verification_rc;
                                        }
                                }

                                psa = SA + L;
                                psa_end = SA + R + 1;
                                index_tmp = begin_index[i];

                                // packed method
                                if (size_read_mod32) {
                                        for (;;) {
                                                psa_loop_if_rc:
                                                if (psa == psa_end) {
                                                        goto loop_verification_rc;
                                                }

                                                //index = *psa - index_tmp;
                                                // Bug here
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -=  index_tmp;

                                                b = ref_bin_ptr_array[index & 3/* mod 4 */] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer_rc;/// Bug Here

                                                p_end = p + size_read_bit32_1;
                                                v = (*p_end) ^ (*q_end_rc);/// rc
                                                v &= popcount_mask;
                                                count = HammingWeightDna64(v);

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_if_rc;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                } else {
                                        for (;;) {
                                                psa_loop_else_rc:
                                                if (psa == psa_end) {
                                                        goto loop_verification_rc;
                                                }

                                                //index = *psa - index_tmp;
                                                // Bug here
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -= index_tmp;

                                                count  = 0;

                                                b = ref_bin_ptr_array[index & 3/* mod 4 */] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer_rc;/// Bug Here

                                                p_end = p + size_read_bit32;

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_else_rc;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                }

                                /// loop verification
                                loop_verification_rc:
                                if (++i == period) {
                                        /// Still not matched.
                                        goto match_failure_rc;
                                }
                        }/* i in [0, period) */

                        match_failure_rc:
                        ++total_reads;
                        rb_reads.ReadNext();/// header
                        if (rb_reads.length_read < 1) {
                                LOGINFO("Search done\n");
                                break;
                        }
                        strcpy(head_str, rb_reads.buffer);
                        /// strip the new line char.
                        head_str[rb_reads.length_read] = '\0';

                        rb_reads.ReadNext();/// sequence
                        continue;

                        match_success:
                        ++matched_reads;
                        ++total_reads;
                        cout << (head_str+1) << "\t"
                             << flag_minus_plus << "\tchr1\t"
                             << index << "\t" << rb_reads.buffer;
                        rb_reads.ReadNext();/// header
                        if (rb_reads.length_read < 1) {
                                LOGINFO("Search done\n");
                                break;
                        }
                        strcpy(head_str, rb_reads.buffer);
                        /// strip the new line char.
                        head_str[rb_reads.length_read] = '\0';

                        rb_reads.ReadNext();/// sequence
                } /* Readings */

                auto end_time_search = std::chrono::high_resolution_clock::now();
                auto duration_search = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_search-begin_time_search).count();
                LOGINFO( "Search phase elapsed time:\t"
                                         << duration_search
                                         << " ms"
                                         << std::endl);
                LOGINFO("Reads processed:\t" << total_reads << "\n");
                LOGINFO("Reads with alignment:\t"
                                << matched_reads
                                << " ("
                                << (matched_reads*100.0/total_reads)
                                << "%)\n");

                delete[] end_array;
                delete[] ptr_array;
                delete[] begin_index;
                delete[] ptr_array_rc;
        }


        void SortedPackedSearch(int argc, char **argv)
        {
                LOGINFO("Sorted searching...\n");
                 /// Build index from files
                auto begin_time_index = std::chrono::high_resolution_clock::now();

                string reads_filename = string(argv[1]);
                string prefix_filename = string(argv[2]);
                uint32_t L_R_min = SecondIndex::size_min - 2;
                if (argc >= 4) {
                        uint32_t get_uint = utility::GetUint(argc, argv[3]);
                        if (get_uint > SecondIndex::size_min) {
                                L_R_min = get_uint;
                        }
                }
                LOGINFO("Size of replica: " << L_R_min << "\n");
                BuildIndexRawData build_index(prefix_filename);
                /// test second index
                SecondIndex second_index(prefix_filename);
                if (second_index.Empty()) {
                        LOGERROR("The second index is empty");
                }

                auto end_time_index = std::chrono::high_resolution_clock::now();
                auto duration_index = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_index-begin_time_index).count();
                LOGINFO( "Load index elapsed time:\t"
                                << duration_index
                                << " ms"
                                << std::endl);
#if 0
                {

                        uint32_t count1 = 0;
                        for (;;) {
                                if (count1 >= second_index.size) {
                                        break;
                                }
                                uint16_t *p16 = second_index.array_ptr + count1;
                                uint32_t pos = *((uint32_t*)p16);
                                uint16_t size_seg = p16[2];
                                count1 += size_seg * build_index.period + 3;
                                p16 += 3;
                                cout << "pos: " << pos << ", size_seg: " << size_seg << endl;
                                for (uint32_t tau = 0; tau < build_index.period; ++tau) {
                                        for (uint32_t i0 = 0; i0 < size_seg; ++i0) {
                                                cout << p16[tau*size_seg+i0] << " ";
                                        }
                                        cout << endl;
                                }

                        }
                        //sbwt::PrintFullSearchMatrix(build_index);
                        second_index.PrintSecondIndex(build_index);
                }
#endif

                auto begin_time_search = std::chrono::high_resolution_clock::now();

                uint32_t period = build_index.period;
                uint32_t N = build_index.length_ref;
                uint32_t *C = build_index.first_column;
                uint32_t **Occ = build_index.occurrence;
                uint32_t *SA = build_index.suffix_array;

                static uint8_t *ref_bin_ptr_array[4] = {nullptr};
                ref_bin_ptr_array[0] = build_index.bin_8bit;
                ref_bin_ptr_array[1] = build_index.bin_8bit + 1*build_index.size_bin_8bit;
                ref_bin_ptr_array[2] = build_index.bin_8bit + 2*build_index.size_bin_8bit;
                ref_bin_ptr_array[3] = build_index.bin_8bit + 3*build_index.size_bin_8bit;

                /// reads
                /// Read first segment.
                reads_buffer rb_reads(reads_filename);
                rb_reads.ReadNext();/// header
                if (rb_reads.length_read < 1) {
                        LOGINFO("Empty reads\n");
                        return;
                }
                char head_str[1024*1024] = {'\0'};
                strcpy(head_str, rb_reads.buffer);
                /// strip the new line char.
                head_str[rb_reads.length_read] = '\0';

                /// count matched segment
                uint64_t total_reads = 1;
                uint64_t matched_reads = 0;

                rb_reads.ReadNext();/// sequence

                /// For second index
                /// the minimum distance between L and R
                uint32_t size_seed = second_index.size_seed;
                uint32_t size_std = size_seed * period;
                uint16_t *ptr16_right = nullptr;
                uint16_t *ptr16_left = nullptr;
                char *ptr_iter = nullptr;
                uint16_t *ptr16_mid = nullptr;
                bool is_2nd_success = false;/* seed */
                bool is_2nd_extend_success = false;/* extend */
                int compare_flag = 0;
                uint32_t current_pos = 0;
                uint32_t distance_2nd = 0;
                uint32_t mid_index_2nd = 0;
                uint32_t power_2nd_count = 0;
                uint32_t mid_pos_2nd = 0;
                uint16_t *ptr16_2nd_begin = second_index.array_ptr;
                uint16_t *ptr16_2nd = nullptr;
                uint32_t index_2nd = 0;
                uint32_t size_range = 0;
                uint32_t pos_original = 0;

                if (rb_reads.length_read < period || size_std < rb_reads.length_read) {
                        LOGERROR("the size of read is shorter than "
                                                 << size_std);
                        return;
                }

                // Assuming the length of single sequence segment is kept same.
                uint32_t size_read_char = rb_reads.length_read;
                uint32_t size_read_mod4 = size_read_char % 4;
                uint32_t size_read_bit8 = (size_read_mod4) ? (size_read_char >> 2) + 1 : size_read_char >> 2;
                uint32_t size_read_mod32 = size_read_char % 32;
                uint32_t size_read_bit32 = size_read_mod32 ?
                                           (size_read_char >> 5) + 1:
                                           size_read_char >> 5;
                uint32_t size_read_bit32_1 = size_read_bit32 - 1;
                uint64_t popcount_mask = DnaStringRightShiftMaskReverse[size_read_mod32];

                /// binary buffer for segment
                uint64_t read_bin_buffer[256] = {0};
                uint64_t read_bin_buffer_rc[256] = {0};
                uint8_t *read_bin = (uint8_t*)read_bin_buffer;
                uint8_t *read_bin_rc = (uint8_t*)read_bin_buffer_rc;
                uint64_t *p = nullptr;
                uint64_t *p_end = nullptr;
                uint64_t *q = read_bin_buffer;
                uint64_t *q_end = nullptr;
                uint64_t *q_end_rc = nullptr;
                if (size_read_mod32) {
                        q_end = q + size_read_bit32_1;
                        q_end_rc = read_bin_buffer_rc + size_read_bit32_1;
                } else {
                        q_end = q + size_read_bit32;
                        q_end_rc = read_bin_buffer_rc + size_read_bit32;
                }

                uint32_t L = 0, R = 0;
                uint32_t N_1 = N - 1;

                // Not sure the addr. of buffer will be changed.
                // ---
                // Alternatively, before calling getline(), *lineptr can contain a
                // pointer to a malloc(3)-allocated buffer *n bytes in size.  If the
                // buffer is not large enough to hold the line, getline() resizes it
                // with realloc(3), updating *lineptr and *n as necessary.
                // ---
                // Reference: http://man7.org/linux/man-pages/man3/getline.3.html

                char key = 0;
                static uint32_t a = 0;
                uint32_t lmd = 0;
                uint32_t *end_array     = new uint32_t[period];
                char *ptr               = rb_reads.buffer + (size_read_char - 1);
                char *ptr_reads         = rb_reads.buffer;
                char *ptr_end           = rb_reads.buffer + size_read_char;
                char **ptr_array        = new char*[period];
                char **ptr_array_rc     = new char*[period];
                char *X = build_index.seq_raw;
                uint32_t *begin_index   = new uint32_t[period];

                uint32_t index_tmp, index;

                uint32_t *psa, *psa_end;

                {
                        uint32_t tmp0 = size_read_char / period;
                        uint32_t tmp1 = size_read_char % period;

                        for (uint32_t i = 0; i != period; ++i) {
                                end_array[i]    = tmp0 + (i < tmp1);
                                ptr_array[i]    = ptr - i;
                                ptr_array_rc[i] = ptr_reads + i;
                                begin_index[i]  = size_read_char - ((end_array[i] - 1)*period + i + 1);
                        }
                }

                uint32_t count  = 0;

                uint8_t *b              = nullptr;
                uint64_t v = 0;

                char flag_minus_plus;

                uint32_t i;
                for(;;) {
                        /// forward
                        flag_minus_plus = '+';
                        sbwtio::BaseChar2Binary8B(ptr_reads, size_read_bit8, read_bin);
                        i = 0;
                        for (;;)
                        {
                                lmd = end_array[i];
                                ptr = ptr_array[i];

                                /// Init
                                key = *ptr;
                                a = charToDna5_32bit[key];
                                L = C[a];
                                R = C[a] + Occ[a][N_1] - 1;

                                for (uint32_t j = 1; j != lmd; ++j) {
                                        ptr -= period;
                                        key = *ptr;
                                        a = charToDna5_32bit[key];
                                        L = C[a] + Occ[a][L-1];
                                        R = C[a] + Occ[a][R] - 1;
                                        if (L > R) {
                                                goto loop_verification;
                                        }
                                }

                                /// Use second index to power searching
                                if (R - L > L_R_min && R-L < 65535/* TODO */) {
                                        ++power_2nd_count;
                                        is_2nd_extend_success = false;
                                        psa = SA + L;
                                        size_range = R - L + 1;
                                        psa_end = psa + size_range;

                                        index_2nd = *psa;
                                        ptr16_2nd = ptr16_2nd_begin + index_2nd;
                                        /// if index_2nd == 0
                                        pos_original = *(uint32_t*)ptr16_2nd;

                                        ptr16_2nd += 3;
#if 0
                                        cout << "L: " << L
                                             << " R: " << R
                                             << " *psa: " << *psa
                                             << " *psa_end: " << *psa_end
                                             << " size_range: " << size_range
                                             << " seed: ";
                                        for (uint32_t u00 = 0; u00 < size_seed; ++u00) {
                                                cout << ptr[u00*period];
                                        }
                                        cout << endl;
#endif

#if 0
                                        /// visualization
                                        cout << "Original SA" << endl;
                                        for (uint32_t u = L; u <= R; ++u) {
                                                uint32_t u0 = *(SA + u);
                                                if (u == L) {
                                                        u0 = pos_original;
                                                }
                                                for (uint32_t j0 = 0; j0 < size_std; ++j0) {
                                                                if (j0 % period == 0) {
                                                                        cout << "\033[1;31m" << X[(j0 + u0) % N] << "\033[0m";
                                                                }
                                                                else {
                                                                        cout << X[(j0 + u0) % N];
                                                                }
                                                                if (j0 % size_seed == (size_seed-1)) {
                                                                        cout << " ";
                                                                }

                                                }
                                                cout << endl;
                                        }
#endif
                                        char *cmp_p = ptr;
                                        char *cmp_r = nullptr;

                                        current_pos = 0;
                                        for (uint32_t tau = 0; tau != period; ++tau) {
                                                if (ptr_end < ptr + size_seed) {
                                                        break;
                                                }
                                                /// array: ptr16_2nd,
                                                /// size:  size_range
                                                ptr16_left = ptr16_2nd;
                                                ptr16_right = ptr16_left + size_range;
                                                uint16_t *ptr16_2nd_head = ptr16_2nd;
                                                uint16_t *ptr16_2nd_tail = ptr16_right - 1;
                                                is_2nd_success = false;
                      #define COMPARE_2ND_INDEX if (mid_index_2nd == 0) \
                                                        mid_pos_2nd = pos_original;\
                                                else\
                                                        mid_pos_2nd = *(psa + mid_index_2nd);\
                                                compare_flag = strncmp(X + mid_pos_2nd + current_pos, ptr, size_seed)

                                                /// binary search
                                                while (ptr16_left < ptr16_right) {
                                                        distance_2nd = ptr16_right - ptr16_left;
                                                        if (distance_2nd == 1) {
                                                                ptr16_mid = ptr16_left;
                                                                mid_index_2nd = *ptr16_mid;
                                                                COMPARE_2ND_INDEX;
                                                                is_2nd_success = compare_flag == 0;
                                                                break;
                                                        } else {
                                                                ptr16_mid = ptr16_left + (distance_2nd / 2);
                                                                mid_index_2nd = *ptr16_mid;
                                                                COMPARE_2ND_INDEX;
                                                                if (compare_flag == 0) {
                                                                        is_2nd_success = true;
                                                                        break;
                                                                } else if (compare_flag > 0)/* mid_key > key */ {
                                                                        ptr16_right = ptr16_mid;
                                                                } else {
                                                                        ptr16_left = ptr16_mid;
                                                                }
                                                        }
                                                }

                                                /// Whether seeding successful
                                                if (is_2nd_success) {
                                                        /// Second index success, need verification.
                                                        /// range search [left, right]
                                                        /// Mid
                                                        /// Caution
                                                        if (mid_pos_2nd < current_pos) {
                                                                continue;
                                                        }

                              #define EXTEND_2nd_INDEX  ptr_iter = ptr_reads;\
                                                        index = mid_pos_2nd - (cmp_p - ptr_reads);\
                                                        cmp_r = X + index;/* Compare from front of read */\
                                                        count = 0;\
                                                        for (;ptr_iter != ptr_end;) {\
                                                                if (*cmp_r++ != *ptr_iter++) {\
                                                                        ++count;\
                                                                        if (count >= period) {\
                                                                                break;\
                                                                        }\
                                                                }\
                                                        }\
                                                        is_2nd_extend_success = count < period;\
                                                        if(is_2nd_extend_success)\
                                                                goto match_success

                                                        EXTEND_2nd_INDEX;
                                                        /// left <==
                                                        ptr16_left = ptr16_mid;
                                                        /// check the head
                                                        if (ptr16_left != ptr16_2nd_head) {
                                                                while (--ptr16_left != ptr16_2nd_head) {
                                                                        mid_index_2nd = *ptr16_left;
                                                                        COMPARE_2ND_INDEX;
                                                                        /// seeding successful
                                                                        if (compare_flag == 0) {
                                                                                /// Caufion
                                                                                if (mid_pos_2nd < current_pos) {
                                                                                        continue;
                                                                                }
                                                                                EXTEND_2nd_INDEX;
                                                                        }
                                                                }
                                                        }
                                                        /// ==> right
                                                        ptr16_right = ptr16_mid;
                                                        /// check head
                                                        if (ptr16_right != ptr16_2nd_tail) {
                                                                while (++ptr16_right != ptr16_2nd_tail) {
                                                                        mid_index_2nd = *ptr16_right;
                                                                        COMPARE_2ND_INDEX;
                                                                        /// seeding successful
                                                                        if (compare_flag == 0) {
                                                                                /// Caution
                                                                                if (mid_pos_2nd < current_pos) {
                                                                                        continue;
                                                                                }
                                                                                EXTEND_2nd_INDEX;
                                                                        }
                                                                }
                                                        }
                                                }

                                                //for (uint32_t u = 0; u != size_seed; ++u) { cout << ptr[u]; } cout << endl;
                                                ptr16_2nd += 3 + size_range;
                                                ptr += size_seed;
                                                current_pos += size_seed;
                                        }
                                        /// second index failed
                                        goto loop_verification;
                                }

                                psa = SA + L;
                                psa_end = SA + R + 1;
                                index_tmp = begin_index[i];

                                // packed method
                                if (size_read_mod32) {
                                        for (;;) {
                                                psa_loop_if:
                                                if (psa == psa_end) {
                                                        goto loop_verification;
                                                }
                                                /// Bug: it may happen that: *psa < index_tmp
                                                /// For example:
                                                ///     0    5
                                                /// Ref ACCCGTAGATTGAGTTGC...
                                                ///   0  3
                                                /// P TTACCCGTAGATTGAGTTGC...
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -= index_tmp;

                                                b = ref_bin_ptr_array[index & 3/*mod 4*/] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer;/// Bug Here

                                                p_end = p + size_read_bit32_1;
                                                v = (*p_end) ^ (*q_end);
                                                v &= popcount_mask;
                                                count = HammingWeightDna64(v);

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_if;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                } else {
                                        for (;;) {
                                                psa_loop_else:
                                                if (psa == psa_end) {
                                                        goto loop_verification;
                                                }
                                                // Bug here
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -= index_tmp;

                                                count  = 0;

                                                b = ref_bin_ptr_array[index & 3/* mod 4 */] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer;/// Bug Here

                                                p_end = p + size_read_bit32;

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_else;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                }

                                /// loop verification
                                loop_verification:
                                if (++i == period) {
                                        /// Still not matched.
                                        goto match_failure;
                                }
                        }/* i in [0, period) */

                        /// failed in forward searching
                        match_failure:
                        /// Reverse Complement Binary Stream
                        sbwtio::BaseChar2Binary8B_RC(read_bin, size_read_char, read_bin_rc);

                        i = 0;
                        flag_minus_plus = '-';
                        for (;;)
                        {
                                lmd = end_array[i];
                                ptr = ptr_array_rc[i];/// rc

                                /// Init
                                key = *ptr;
                                a = rcCharToDna5_32bit[key];/// rc
                                L = C[a];
                                R = C[a] + Occ[a][N_1] - 1;

                                for (uint32_t j = 1; j != lmd; ++j) {
                                        ptr += period;/// rc
                                        key = *ptr;
                                        /// Bugs here
                                        a = rcCharToDna5_32bit[key];/// rc
                                        L = C[a] + Occ[a][L-1];
                                        R = C[a] + Occ[a][R] - 1;
                                        if (L > R) {
                                                goto loop_verification_rc;
                                        }
                                }

                                // TODO Use 2nd index to power reverse-complemently searching
                                psa = SA + L;
                                psa_end = SA + R + 1;
                                index_tmp = begin_index[i];

                                // packed method
                                if (size_read_mod32) {
                                        for (;;) {
                                                psa_loop_if_rc:
                                                if (psa == psa_end) {
                                                        goto loop_verification_rc;
                                                }

                                                //index = *psa - index_tmp;
                                                // Bug here
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -=  index_tmp;

                                                b = ref_bin_ptr_array[index & 3/* mod 4 */] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer_rc;/// Bug Here

                                                p_end = p + size_read_bit32_1;
                                                v = (*p_end) ^ (*q_end_rc);/// rc
                                                v &= popcount_mask;
                                                count = HammingWeightDna64(v);

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_if_rc;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                } else {
                                        for (;;) {
                                                psa_loop_else_rc:
                                                if (psa == psa_end) {
                                                        goto loop_verification_rc;
                                                }

                                                //index = *psa - index_tmp;
                                                // Bug here
                                                index = *psa;
                                                if (index < index_tmp) {
                                                        ++psa;
                                                        continue;
                                                }
                                                index -= index_tmp;

                                                count  = 0;

                                                b = ref_bin_ptr_array[index & 3/* mod 4 */] + (index >> 2);
                                                p = (uint64_t *) b;
                                                q = read_bin_buffer_rc;/// Bug Here

                                                p_end = p + size_read_bit32;

                                                for (; p != p_end;) {
                                                        v = (*p++) ^ (*q++);
                                                        count += HammingWeightDna64(v);
                                                        if (count > period) {
                                                                ++psa;// Careful about its position
                                                                goto psa_loop_else_rc;
                                                        }
                                                }

                                                goto match_success;
                                        } // psa
                                }

                                /// loop verification
                                loop_verification_rc:
                                if (++i == period) {
                                        /// Still not matched.
                                        goto match_failure_rc;
                                }
                        }/* i in [0, period) */

                        match_failure_rc:
                        ++total_reads;

                        rb_reads.ReadNext();/// header
                        if (rb_reads.length_read < 1) {
                                LOGINFO("Search done\n");
                                break;
                        }
                        strcpy(head_str, rb_reads.buffer);
                        /// strip the new line char.
                        head_str[rb_reads.length_read] = '\0';

                        rb_reads.ReadNext();/// sequence
                        continue;

                        match_success:
                        ++matched_reads;
                        ++total_reads;
                        cout << (head_str+1) << "\t"
                             << flag_minus_plus << "\tchr1\t"
                             << index << "\t" << rb_reads.buffer;
                        rb_reads.ReadNext();/// header
                        if (rb_reads.length_read < 1) {
                                LOGINFO("Search done\n");
                                break;
                        }
                        strcpy(head_str, rb_reads.buffer);
                        /// strip the new line char.
                        head_str[rb_reads.length_read] = '\0';

                        rb_reads.ReadNext();/// sequence
                } /* Readings */

                auto end_time_search = std::chrono::high_resolution_clock::now();
                auto duration_search = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_search-begin_time_search).count();
                LOGINFO( "Search phase elapsed time:\t"
                                         << duration_search
                                         << " ms"
                                         << std::endl);
                LOGINFO("Reads processed:\t" << total_reads << "\n");
                LOGINFO("Reads with alignment:\t"
                                        << matched_reads
                                        << " ("
                                        << (matched_reads*100.0/total_reads)
                                        << "%)\n");

                LOGINFO("Second index powering searching: " << power_2nd_count << "\n");
                delete[] end_array;
                delete[] ptr_array;
                delete[] begin_index;
                delete[] ptr_array_rc;

        }


} /* namespace sbwt */
