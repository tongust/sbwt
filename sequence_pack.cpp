#include <vector>
#include <string>
#include <bitset>
#include <iostream>

#include "sequence_pack.h"

namespace sbwtio
{
using std::vector;
using std::bitset;
using std::cout;
using std::endl;
using std::string;


/// For converting from ASCII to the Dna5 code where A=0, C=1, G=2,
/// T=3, N=4
uint8_t charToDna5[] = {
                /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
                /*    A     C           G                    N */
                /*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*             T */
                /*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0,
                /*    a     c           g                    n */
                /* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*             t */
                /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

/// For converting from ASCII to the reverse-complement Dna5 code where
/// A=3, C=2, G=1, T=0, N=4
uint8_t rcCharToDna5[] = {
                /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*  64 */ 0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 4, 0,
                /*    A     C           G                    N */
                /*  80 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*             T */
                /*  96 */ 0, 3, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 4, 0,
                /*    a     c           g                    n */
                /* 112 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /*             t */
                /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

        /*
         * Look-up table: reverse complement byte-wise. Such 10-00-11-01 to 11-01-10-00
         */
        const uint8_t Array256Swap2bitTableComplement[256] = {
                        170,234,42,106,186,250,58,122,138,202,10,74,154,218,26,90,
                        174,238,46,110,190,254,62,126,142,206,14,78,158,222,30,94,
                        162,226,34,98,178,242,50,114,130,194,2,66,146,210,18,82,
                        166,230,38,102,182,246,54,118,134,198,6,70,150,214,22,86,
                        171,235,43,107,187,251,59,123,139,203,11,75,155,219,27,91,
                        175,239,47,111,191,255,63,127,143,207,15,79,159,223,31,95,
                        163,227,35,99,179,243,51,115,131,195,3,67,147,211,19,83,
                        167,231,39,103,183,247,55,119,135,199,7,71,151,215,23,87,
                        168,232,40,104,184,248,56,120,136,200,8,72,152,216,24,88,
                        172,236,44,108,188,252,60,124,140,204,12,76,156,220,28,92,
                        160,224,32,96,176,240,48,112,128,192,0,64,144,208,16,80,
                        164,228,36,100,180,244,52,116,132,196,4,68,148,212,20,84,
                        169,233,41,105,185,249,57,121,137,201,9,73,153,217,25,89,
                        173,237,45,109,189,253,61,125,141,205,13,77,157,221,29,93,
                        161,225,33,97,177,241,49,113,129,193,1,65,145,209,17,81,
                        165,229,37,101,181,245,53,117,133,197,5,69,149,213,21,85
        };
        /*
         * The table turns over a bite 2-bit-wise. Such 10-00-11-01 to 01-11-00-10
         * */
        const uint8_t Array256Swap2bitTable[256] = {
                        0,  64, 128, 192,  16,  80, 144, 208,
                        32,  96, 160, 224,  48, 112, 176, 240,
                        4,  68, 132, 196,  20,  84, 148, 212,
                        36, 100, 164, 228,  52, 116, 180, 244,
                        8,  72, 136, 200,  24,  88, 152, 216,
                        40, 104, 168, 232,  56, 120, 184, 248,
                        12,  76, 140, 204,  28,  92, 156, 220,
                        44, 108, 172, 236,  60, 124, 188, 252,
                        1,  65, 129, 193,  17,  81, 145, 209,
                        33,  97, 161, 225,  49, 113, 177, 241,
                        5,  69, 133, 197,  21,  85, 149, 213,
                        37, 101, 165, 229,  53, 117, 181, 245,
                        9,  73, 137, 201,  25,  89, 153, 217,
                        41, 105, 169, 233,  57, 121, 185, 249,
                        13,  77, 141, 205,  29,  93, 157, 221,
                        45, 109, 173, 237,  61, 125, 189, 253,
                        2,  66, 130, 194,  18,  82, 146, 210,
                        34,  98, 162, 226,  50, 114, 178, 242,
                        6,  70, 134, 198,  22,  86, 150, 214,
                        38, 102, 166, 230,  54, 118, 182, 246,
                        10,  74, 138, 202,  26,  90, 154, 218,
                        42, 106, 170, 234,  58, 122, 186, 250,
                        14,  78, 142, 206,  30,  94, 158, 222,
                        46, 110, 174, 238,  62, 126, 190, 254,
                        3,  67, 131, 195,  19,  83, 147, 211,
                        35,  99, 163, 227,  51, 115, 179, 243,
                        7,  71, 135, 199,  23,  87, 151, 215,
                        39, 103, 167, 231,  55, 119, 183, 247,
                        11,  75, 139, 203,  27,  91, 155, 219,
                        43, 107, 171, 235,  59, 123, 187, 251,
                        15,  79, 143, 207,  31,  95, 159, 223,
                        47, 111, 175, 239,  63, 127, 191, 255
        };

template <typename T>
void PrintBinary(T c, const std::string &end_str)
{
        auto s = bitset<sizeof(T)*8>(c).to_string();
        using std::cout;
        for (int i = 0; i < s.size(); ++i)
        {
                if (i != 0 && i % 8 == 0)
                        cout << " ";
                cout << s[i];
        }
        cout << end_str;
        return;
}

#ifdef SBWT_DEBUG
//#define SBWT_DEBUG_BASECHAR2BINARY8B
#else
#undef SBWT_DEBUG_BASECHAR2BINARY8B
#endif
//#define SBWT_DEBUG_BASECHAR2BINARY8B

        /// Turn string of DNAs to packed 8-bit sequence.
        /// @param size: the one of binary sequence.
        void BaseChar2Binary8B(char *buffer, uint32_t size, uint8_t *bin)
        {
#ifdef SBWT_DEBUG_BASECHAR2BINARY8B
                auto print = [](uint32_t v)->void {
                        auto s = bitset<32>(v).to_string();
                        for (int i = 0; i < 32; ++i) {
                                cout << s[i];
                                if (i % 2 == 1)
                                        cout << " ";
                                if (i % 4 == 3)
                                        cout << " ";
                                if (i % 8 == 7)
                                        cout << "  ";
                        } cout << endl;
                };
#endif
                uint32_t *p = (uint32_t*)buffer;
                uint32_t *p_end = p + size;
                static uint8_t array[4] = {0};
                static uint32_t &c = *((uint32_t*)array);
                for (; p != p_end; ++p) {
                        c = ((*p) >> 1) & 0x03030303;
#ifdef SBWT_DEBUG_BASECHAR2BINARY8B
                        static int count = 0;
                        cout << ++count << endl;
                        print(*p);
                        print(c);
#endif
                        c |= (c >> 6);
#ifdef SBWT_DEBUG_BASECHAR2BINARY8B
                        print(c);
#endif
                        c |= c >> 12;
#ifdef SBWT_DEBUG_BASECHAR2BINARY8B
                        print(c);
#endif
                        *bin = array[0];
                        ++bin;
#ifdef SBWT_DEBUG_BASECHAR2BINARY8B
                        cout << bitset<8>(c) << endl;
#endif
                }
        }


        /// Version 1: Streamed from BaseChar2Binary8B (as benchmark)
        /// One-byte-wise
        /// Profile 100%
        void BaseChar2Binary8B_RC(uint8_t *bin, uint32_t size, uint8_t *bin_rc)
        {
                uint32_t size_mod8 = size & 3;// mod 4
                uint32_t size_8bit = size >> 2;
                static uint8_t v2[2] = {0};
                static uint16_t *p = (uint16_t*)v2;
                static uint16_t &v = *p;
                static uint8_t &v2_lo = v2[0];
                if (size_mod8 != 0) {
                        size_mod8 <<= 1;
                        uint8_t *bin_begin = bin + (size_8bit - 1);
                        for (uint32_t i = 0; i != size_8bit; ++i) {
                                v = *((uint16_t *) (bin_begin));
                                v >>= size_mod8;
                                *bin_rc++ = Array256Swap2bitTableComplement[v2_lo];
                                --bin_begin;
                        }
                        *bin_rc = Array256Swap2bitTableComplement[(*bin)] >> (8-size_mod8);
                } else {
                        bin_rc += size_8bit - 1;
                        for (uint32_t i = size_8bit; i != 0; --i) {
                                *bin_rc-- = Array256Swap2bitTableComplement[*bin++];
                        }
                }
        }


        /// N as G --> C
        const char DnaCharMapReverseComplement[256] = {
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 84, 0, 71, 0, 0, 0, 67, 0, 0, 0, 0, 0, 0, 67/*N*/, 0,
                        0, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };
        /// Version 2: reverse the dnas and then use BaseChar2Binary8B
        /// Profile 41%
        void BaseChar2Binary8B_RC(char *buffer, uint32_t size, uint8_t *bin)
        {
                static char *p_end = nullptr;
                static char *p = nullptr;
                static char *p0 = nullptr;

                p0 = buffer;
                p_end = buffer + ( (size+1) >> 1 );
                p = buffer + (size - 1);

                for (; buffer != p_end; ) { static char tmp = 0; tmp = DnaCharMapReverseComplement[*buffer]; *buffer++ = DnaCharMapReverseComplement[*p]; *p-- = tmp; }
                size = size & 3 == 0 ? size >> 2 : (size >> 2) + 1;
                BaseChar2Binary8B(p0, size, bin);
        }


        /// Version 3: turned from dnas sequence as BaseChar2Binary8B does.
        /// size_char must be not less than 4.
        /// Profile 79%
        void BaseChar2Binary8B_RC_Exter(
                        char *buffer,
                        uint8_t *bin,
                        uint32_t size_char,
                        uint32_t size_8bit,
                        uint32_t mod_4_by2,         /* need multiplyied by 2 */
                        uint32_t mod_4_r_by2        /* need multiplyied by 2 */
        )
        {
                static uint8_t array[4] = {0};
                static uint32_t &c = *((uint32_t*)array);
                if (mod_4_by2) {
                        uint32_t *p_end = (uint32_t*)buffer;
                        uint32_t *p = (uint32_t*)(buffer + (size_char - 4));
                        uint8_t *bin_end = bin + (size_8bit-1);
                        for (; bin != bin_end;) {
                                c = ((*p) >> 1) & 0x03030303;
                                c |= (c >> 6);
                                c |= c >> 12;
                                *bin = Array256Swap2bitTableComplement[array[0]];
                                --p;
                                ++bin;
                        }
                        /// the last one
                        c = ((*p_end) >> 1) & 0x03030303;
                        c |= (c >> 6);
                        c |= c >> 12;
                        *bin = Array256Swap2bitTableComplement[array[0]] >> mod_4_r_by2;
                } else {
                        uint32_t *p = (uint32_t*)buffer;
                        uint32_t *p_end = p + size_8bit;
                        bin += (size_8bit-1);
                        for (; p != p_end; ++p) {
                                c = ((*p) >> 1) & 0x03030303;
                                c |= (c >> 6);
                                c |= c >> 12;
                                *bin-- = Array256Swap2bitTableComplement[array[0]];
                        }
                }
        }


/* It seems that 'inline' does not work here */
/* TODO should be included */
void BaseChar2Binary64B(char *buffer, uint64_t size_buffer, uint64_t* binary_seq)
{
        uint64_t *p = (uint64_t*)buffer;
        uint64_t c = 0, v = 0, i = 0;
        for (;;) {
                v = 0;
                if (i >= size_buffer) { *binary_seq = v; break; }
                else { i += 8; }
                c = *p;
                UINT64_TO_8BP(c);
                v |= c       &  0x000000000000FFFFULL;

                if (i >= size_buffer) { *binary_seq = v; break; }
                else { i += 8; }
                ++p;
                c = *p;
                UINT64_TO_8BP(c);
                v |= (c<<16) &  0x00000000FFFF0000ULL;

                if (i >= size_buffer) { *binary_seq = v; break; }
                else { i += 8; }
                ++p;
                c = *p;
                UINT64_TO_8BP(c);
                v |= (c<<32) &  0x0000FFFF00000000ULL;

                if (i >= size_buffer) { *binary_seq = v; break; }
                else { i += 8; }
                ++p;
                c = *p;
                UINT64_TO_8BP(c);
                v |= (c<<48) &  0xFFFF000000000000ULL;

                *binary_seq = v;
                ++p;
                ++binary_seq;
        }
}

/* Reverse Complement version */
/// 0100 0[00]1 - A
/// 0100 0[01]1 - C
/// 0100 0[11]1 - G
/// 0101 0[10]0 - T
/// 0100 1[11]0 - N (G)
void BaseChar2Binary64B_RC(uint64_t* binary_seq, uint64_t size_seq, uint64_t *binary_seq_rc)
{
        int n           = size_seq >> 5;
        int tail        = (size_seq % 32)<<1;
        int tail_rc     = 64 - tail;

        if (tail != 0) {
                uint64_t *p = binary_seq + n;
                uint64_t *q = binary_seq_rc + n;
                while (q != binary_seq_rc) {
                        *q = (*p) >> tail_rc;
                        --p;
                        *q |= (*p) << tail;
                        --q;
                }
                *q = (*p) << tail_rc;
                /* Make sure n is large enough */
                ++n;
        }
        else {
                for (int i = 0; i != n; ++i) {
                        binary_seq_rc[i] = binary_seq[i];
                }
        }

        /* Reverse */
        uint8_t *p_iter         = (uint8_t *)binary_seq;
        uint8_t *p_iter_end     = (uint8_t*)(binary_seq+n);
        uint8_t *q_iter         = (uint8_t *)(binary_seq_rc+n);

        --q_iter;

        while (p_iter != p_iter_end) {
                *p_iter = Array256Swap2bitTable[*q_iter];
                ++p_iter;
                --q_iter;
        }

        for (int i = 0; i <= (size_seq/32); ++i) { PrintBinary(binary_seq[i]); }

        /* Complement */
        for (int i = 0; i != n; ++i) {
                *binary_seq ^= 0xAAAAAAAAAAAAAAAAULL;
                binary_seq++;
        }

        return;
}

void NaiveBaseChar2Binary64B(char *buffer, uint64_t size_buffer, uint64_t* binary_seq)
{
        uint8_t *b = (uint8_t*)binary_seq;
        uint8_t *c = (uint8_t*)buffer;
        for (uint64_t i = 0; i < size_buffer; ++i) {
                if (i%4 == 0 && i != 0) { ++b; }
                *b |= charToDna5[*c]<<((i%4)*2);
                ++c;
        }
}

} /* namespace sbwtio */
