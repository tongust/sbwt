#ifndef SBWT_SEARCH_H
#define SBWT_SEARCH_H

#include <stdlib.h>
#include <stdio.h>

#include <cstdint>
#include <string>

namespace sbwt
{
#define _GNU_SOURCE
#define MAX_PERIOD 10
#define CHUNCK_SIZE 1024
#define MAX_READS_LENGTH 1024
#define MAX_READS_BINARY_LENGTH 32

        using std::string;

        /// Variables
        extern const uint32_t FLAG32[32];
        extern const uint64_t DnaStringRightShiftMaskReverse[32];

        class bitset64
        {
        public:
                bitset64():flag(1) {}
                bitset64(uint64_t *, uint32_t);
                bitset64(const bitset64&);
                bitset64& operator=(bitset64&);

        public:
                void RightShift(uint32_t, uint32_t);
                void RightShiftDna(uint32_t, uint32_t);
                void clear();
                uint32_t size();
#ifdef SBWT_VERBOSE
                void Print();
                static void Print(uint64_t*, uint32_t);
#endif

        public:
                uint64_t array[32*MAX_READS_BINARY_LENGTH+1024];//////////////////////////////////// Why ???
        private:
                uint32_t flag;
        };

#ifdef SBWT_VERBOSE
        void TestBitShift();
        void TestBit64Map(const std::string&, const std::string&);
        void TestBit8Map(const std::string&, const std::string&);
        void TestNaiveMap(const std::string&, const std::string&);
        void TestBaseChar2BinaryReverseComplement(char **argv);
#endif /* SBWT_VERBOSE */

        class reads_buffer
        {
        public:
                reads_buffer();
                reads_buffer(char *, uint32_t);
                reads_buffer(const string&);
                reads_buffer(const reads_buffer&);
                ~reads_buffer();

        public:
                char            *buffer;
                char            **ptr_buffer;
                uint64_t        length_buffer;
                uint64_t        *ptr_length_buffer;
                ssize_t         length_read;

                FILE            *file_stream;
                string          file_name;

        public:
                ssize_t Size();
                void ReadNext();
                void Init();
                void Dispose();
        };


        void UnsortedUnpackedSearch(char **);
        void UnsortedPackedSearch(char **);
        void SortedPackedSearch(int argc, char **);
} /* namespace sbwt */

#endif /* SBWT_SEARCH_H */
