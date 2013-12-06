#include <math.h>
#include <assert.h>
#include "sse_helpers.h"

#define USE_INDEX

#ifdef USE_SSE
#include <emmintrin.h>
#include <x86intrin.h>

const uint32_t cmpestrmFlag = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
const uint32_t cmpestrmMatchFlag = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_POSITIVE_POLARITY; 
#endif

int CountMismatches(const char* __restrict__ read,
        const char* __restrict__ reference,
        int length) {
    int NumMismatches = 0;
    int i = 0;
/*#ifdef USE_SSE
    __m128i dontcarSIMD = _mm_set1_epi8('N');

    int length_a = length - (length % 16);
    for (; i < length_a; i+=16) {
        int toProcess = 16;
        __m128i readSIMD = _mm_lddqu_si128((__m128i* const) &read[i]);
        __m128i inputSIMD = _mm_lddqu_si128((__m128i* const) &reference[i]);
        __m128i readMaskSIMD = _mm_cmpestrm(readSIMD, toProcess, dontcarSIMD, toProcess, cmpestrmFlag);
        __m128i cmpres = _mm_and_si128(readMaskSIMD, _mm_cmpestrm(readSIMD, toProcess, inputSIMD, toProcess, cmpestrmFlag));
        NumMismatches += _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
    }
#endif*/
    for (; i < length; i++) {
        NumMismatches += MismatchPair[(int) read[i]][(int) reference[i]];
    }
    return NumMismatches;
}

int DoInitialSeedAndExtendForward(const Chromosome& chromosome,
                                  int start,
                                  int end,
                                  size_t initExtend,
                                  const std::string& readSeq,
                                  size_t maxMismatches,
                                  std::vector<unsigned int>* PD) {
    if (initExtend > readSeq.size()) {
      return 0;
    }

    char initBase = readSeq[0];
    int initBaseNum = Convert2Num[(int) initBase];

    if (initBase == 'N') {
      return 0;
    }

    const std::string& chromosomeSeq = chromosome.getSeq();

#ifdef USE_SSE
    initExtend = std::min(16, (int) initExtend);
    __m128i readSIMD = _mm_lddqu_si128((__m128i* const) &readSeq[0]);
    __m128i dontCareSIMD = _mm_set1_epi8('N');
    __m128i dontCareMaskSIMD = _mm_cmpestrm(readSIMD, initExtend, dontCareSIMD, initExtend, cmpestrmFlag);
#endif

#ifdef USE_INDEX
    for (const unsigned int *it = chromosome.getPositions(initBaseNum, start); *it < end; it++) {
        unsigned int pos = *it;
        unsigned int nMismatches = 0;
#ifdef USE_SSE
        __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos]);
        __m128i cmpres = _mm_and_si128(dontCareMaskSIMD, _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmFlag));
        nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
#else
        for (unsigned int i = 1; i < initExtend; i++) {
            nMismatches += MismatchPair[(int) readSeq[i]][(int) chromosomeSeq[pos + i]];
        }
#endif
        if (nMismatches <= maxMismatches) {
            PD[nMismatches].push_back(pos + initExtend - 1);
        }
    }
#else
    for (unsigned int pos = start; pos < end; pos++) {
        if (initBase != chromosomeSeq[pos]) continue;
        unsigned int nMismatches = 0;
#ifdef USE_SSE
        __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos]);
        __m128i cmpres = _mm_and_si128(dontCareMaskSIMD, _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmFlag));
        nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
#else
        for (unsigned int i = 1; i < initExtend; i++) {
            nMismatches += MismatchPair[(int) readSeq[i]][(int) chromosomeSeq[pos + i]];
        }
#endif
        if (nMismatches <= maxMismatches) {
            PD[nMismatches].push_back(pos + initExtend - 1);
        }
    }

#endif
    return initExtend;
}

int DoInitialSeedAndExtendReverse(const Chromosome& chromosome,
                                  int start,
                                  int end,
                                  size_t initExtend,
                                  const std::string& readSeq,
                                  size_t maxMismatches,
                                  std::vector<unsigned int>* PD) {
    if (initExtend > readSeq.size()) {
      return 0;
    }
    size_t readLength = readSeq.size();

    char initBase = readSeq[readLength - 1];
    int initBaseNum = Convert2Num[(int) initBase];

    if (initBase == 'N') {
      return 0;
    }

    const std::string& chromosomeSeq = chromosome.getSeq();

#ifdef USE_SSE
    initExtend = std::min(16, (int) initExtend);
    const uint32_t cmpestrmFlag = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;

    __m128i readSIMD = _mm_setzero_si128();
    for (int i = 0; i < initExtend; i++) {
        readSIMD = _mm_slli_si128(readSIMD, 1);
        readSIMD = _mm_insert_epi8(readSIMD, readSeq[readLength - 1 - i], 0);
    }
    __m128i dontCareSIMD = _mm_set1_epi8('N');
    __m128i dontCareMaskSIMD = _mm_cmpestrm(readSIMD, initExtend, dontCareSIMD, initExtend, cmpestrmFlag);
#endif

#ifdef USE_INDEX
    for (const unsigned int *it = chromosome.getPositions(initBaseNum, start); *it < end; it++) {
        unsigned int pos = *it;
        unsigned nMismatches = 0;
#ifdef USE_SSE
        __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos + 1 - initExtend]);
        __m128i cmpres = _mm_and_si128(dontCareMaskSIMD, _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmFlag));
        nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
#else
        for (unsigned int i = 1; i < initExtend; i++) {
            nMismatches += MismatchPair[(int) readSeq[readLength - 1 - i]][(int) chromosomeSeq[pos - i]];
        }
#endif
        if (nMismatches <= maxMismatches) {
            PD[nMismatches].push_back(pos - initExtend + 1);
        }
    }
#else
    for (unsigned int pos = start; pos < end; pos++) {
        if (initBase != chromosomeSeq[pos]) continue;
        unsigned nMismatches = 0;
#ifdef USE_SSE
        __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos + 1 - initExtend]);
        __m128i cmpres = _mm_and_si128(dontCareMaskSIMD, _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmFlag));
        nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
#else
        for (unsigned int i = 1; i < initExtend; i++) {
            nMismatches += MismatchPair[(int) readSeq[readLength - 1 - i]][(int) chromosomeSeq[pos - i]];
        }
#endif
        if (nMismatches <= maxMismatches) {
            PD[nMismatches].push_back(pos - initExtend + 1);
        }
    }
#endif
    return initExtend;
}
