#include <stdio.h>
#include <stdlib.h>

#include <cstdint>
#include <iostream>
#include <string>
#include <bitset>

#include "sbwt_search.h"
#include "sequence_pack.h"
#include "log.h"

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

        void bitset64::Print()
        {
                using std::bitset;
                using std::cout;
                using std::endl;
                using std::string;

                auto print_str = [](const uint64_t &t)->void {
                        auto s = bitset<64>(t).to_string();
                        for (int i = 0; i != 64; ++i) {
                                if (i && i % 8 == 0 && i != 63) {
        //                                cout << " ";
                                }
                                cout << s[i];
                        }
                };
                for (int i = 0; i != size(); ++i) {
                        if ( (i % MAX_READS_BINARY_LENGTH) == 0) {
                                cout << "[" << (i>>5) << "]\t";
                                print_str(array[i]);
                        } else if ((i+1) % MAX_READS_BINARY_LENGTH == 0){
                                cout << " | ";
                                print_str(array[i]);
                                cout << endl;
                        } else {
                                cout << " | ";
                                print_str(array[i]);
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

        void Test_reads_buffer(const std::string &file_name)
        {

                string bs = "01010101001010101000100001010100101011111111111111100100111111111111110001010011111111110001000000011110000000111111100011110011";
                string bs0(bs.begin(), bs.begin() + 64),
                                bs1(bs.begin()+64, bs.end());
                uint64_t array[2];
                bitset<64> b0(bs0);
                bitset<64> b1(bs1);

                array[0] = bitset<64>(bs0).to_ullong();
                array[1] = bitset<64>(bs1).to_ullong();

                b0 = b1;

                b0[32] = 1;
                b0[33] = 1;
                b0[1] = ~b0[1];
                b0[2] = ~b0[2];
                b0[3] = ~b0[3];
                cout << b0 << endl << b1 << endl;

                uint64_t v0 = b0.to_ullong();
                uint64_t v1 = b1.to_ullong();
                uint64_t v2 = 0;
                int diff_count = 0;

                HammingWeightDna64(v0, v1);
                cout << diff_count << endl;
                /// hammming weight
#if 0
                reads_buffer rb(file_name);
                bitset64 bt;
                rb.ReadNext();
                rb.ReadNext();
                std::cout << rb.length_read << std::endl;
                for (int i = 0; i < rb.length_read; ++i)
                        std::cout << rb.buffer[i]; std::cout << std::endl;
                sbwtio::BaseChar2Binary64B(rb.buffer, rb.length_read, bt.array);
                for (int i = 0; i < 32; ++i)
                        bt.RightShiftDna(i, 2);
                bt.Print();
#endif
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
                length_buffer(MAX_READS_BINARY_LENGTH+1),
                length_read(0),
                file_stream(nullptr)
        {
                Init();
        }

        reads_buffer::reads_buffer(char *fn, uint32_t size):
                length_buffer(MAX_READS_BINARY_LENGTH+1),
                length_read(0)
        {
                file_name.resize(size);
                std::copy(fn, fn+size, file_name.begin());
                Init();
        }

        reads_buffer::reads_buffer(const string &fn):
                length_buffer(MAX_READS_BINARY_LENGTH+1),
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
                if (length_read == -1) {
                        logger::LogError("ReadNext");
                        return;
                }
                --length_read;
        }

} /* namespace sbwt */
