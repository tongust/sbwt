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

/*
 * The table turns over a bite 2-bit-wise. Such 10-00-11-01 to 01-11-00-10
 * */
extern const uint8_t Array256Swap2bitTable[256];

/* Using builtin popcount by default */
#ifndef BUILTIN_POPCOUNT
#define BUILTIN_POPCOUNT
#endif

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

/// turn 8 chars (1 unsigned 64-bit integer) to 8 base pairs (16-bit)
#define UINT64_TO_8BP(c)\
        c &= BASE_EXTRACT_MASK64;\
        c >>= 1;\
        c |= (c >> 6 ) & 0x0C0C0C0C0C0C0C0CULL;\
        c |= (c >> 12) & 0xF0F0F0F0F0F0F0F0ULL;\
        c &=             0x000000FF000000FFULL;\
        c |= c >> 24

/* It seems that 'inline' does not work here */
/* TODO should be include */
void BaseChar2Binary64B(char *buffer, uint64_t size_buffer, uint64_t* binary_seq);

/* Reverse Complement version */
/// 0100 0[00]1 - A
/// 0100 0[01]1 - C
/// 0100 0[11]1 - G
/// 0101 0[10]0 - T
/// 0100 1[11]0 - N (G)
inline
void RC_BaseChar2Binary64B(uint64_t* binary_seq, int size_seq, uint64_t *binary_seq_rc);

void NaiveBaseChar2Binary64B(char *buffer, uint64_t size_buffer, uint64_t* binary_seq);

} /* namespace sbwtio */
#endif /* SBWT_SEQUENCE_PACK_H */
