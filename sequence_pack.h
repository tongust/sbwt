#ifndef SBWT_SEQUENCE_PACK_H
#define SBWT_SEQUENCE_PACK_H

#include <vector>
#include <string>
#include <bitset>

using std::vector;
using std::bitset;
using std::string;



namespace sbwtio
{

/// For converting from ASCII to the Dna5 code where A=0, C=1, G=2,
/// T=3, N=4
extern uint8_t charToDna5[256];

/// For converting from ASCII to the reverse-complement Dna5 code where
/// A=3, C=2, G=1, T=0, N=4
extern uint8_t rcCharToDna5[256];

/**
 * The table turns over a bite 2-bit-wise. Such 10-00-11-01 to 01-11-00-10
 * */
extern const uint8_t Array256Swap2bitTable[256];

///
extern const char DnaCharMapReverseComplement[256];
/* Using builtin popcount by default */
#ifndef BUILTIN_POPCOUNT
#define BUILTIN_POPCOUNT
#endif

/**
 * v2, uint64_t, temporary value
 * diff_count, int, Hamming weight
 */
#define HammingWeightDna64(v2) __builtin_popcountll( ((v2<<1) | v2) & 0xAAAAAAAAAAAAAAAAull )

/* Benchmark */
#ifdef BITSET_POPCOUNT
#define HammingWeight32(x) std::bitset<32>(x).count()
#define HammingWeight64(x) std::bitset<64>(x).count()
#endif

/* Speed up 250% in profile */
#ifdef BUILTIN_POPCOUNT
#define HammingWeight32(x) __builtin_popcount(x)
#define HammingWeight64(x) __builtin_popcountll(x)
#endif

/* Speed up 154% in profile */
#ifdef BOWTIE_POPCOUNT
#define HammingWeight64(x) \
        x = x - ((x >> 1) & 0x5555555555555555llu); \
        x = (x & 0x3333333333333333llu) + ((x >> 2) & 0x3333333333333333llu); \
        x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0Fllu; \
        x = x + (x >> 8); \
        x = x + (x >> 16); \
        x = x + (x >> 32); \
        x = (x & 0x3F)
#endif

        template <typename T>
        void PrintBinary(T c, const std::string &end_str = "\n");

/// 0100 0[00]1 - A
/// 0100 0[01]1 - C
/// 0100 0[11]1 - G
/// 0101 0[10]0 - T
/// 0000 0110 ... mask for extracting DNA base pair
#define BASE_EXTRACT_MASK64 0x0606060606060606ULL
#define BASE_EXTRACT_MASK32 0x06060606UL

/// turn 8 chars (1 unsigned 64-bit integer) to 8 base pairs (16-bit)
#define UINT64_TO_8BP(c)\
        c &= BASE_EXTRACT_MASK64;\
        c >>= 1;\
        c |= (c >> 6 ) & 0x0C0C0C0C0C0C0C0CULL;\
        c |= (c >> 12) & 0xF0F0F0F0F0F0F0F0ULL;\
        c &=             0x000000FF000000FFULL;\
        c |= c >> 24


        /* It seems that 'inline' does not work here */
        /* TODO should be included */
        void BaseChar2Binary64B(char*, uint64_t, uint64_t*);
        void BaseChar2Binary8B(char*, uint32_t, uint8_t*);

        /* Reverse Complement version */
        /// 0100 0[00]1 - A
        /// 0100 0[01]1 - C
        /// 0100 0[11]1 - G
        /// 0101 0[10]0 - T
        /// 0100 1[11]0 - N (G)
        void BaseChar2Binary64B_RC(uint64_t* binary_seq, uint64_t size_seq, uint64_t *binary_seq_rc);

        /// Version 1: Streamed from BaseChar2Binary8B (as benchmark)
        void BaseChar2Binary8B_RC(uint8_t *bin/*src*/, uint32_t size, uint8_t *bin_rc/*sink*/);

        /// Version 2: reverse the dnas and then use BaseChar2Binary8B
        void BaseChar2Binary8B_RC(char *buffer, uint32_t size, uint8_t *bin);

        /// Version 3: turned from dnas sequence as BaseChar2Binary8B does.
        void BaseChar2Binary8B_RC_Exter(char *buffer, uint8_t *bin, uint32_t size_char, uint32_t size_8bit, uint32_t mod_4_by2, uint32_t mod_4_r_by2);

        void NaiveBaseChar2Binary64B(char *buffer, uint64_t size_buffer, uint64_t* binary_seq);

} /* namespace sbwtio */
#endif /* SBWT_SEQUENCE_PACK_H */
