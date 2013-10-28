#include <math.h>
#include <assert.h>
#include "sse_helpers.h"

#ifdef USE_SSE
#include <emmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <x86intrin.h>
#endif

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

    const unsigned int *posS, *posE;
    chromosome.getPositions(initBaseNum, start, end, &posS, &posE);
    const std::string& chromosomeSeq = chromosome.getSeq();

#ifdef USE_SSE
    initExtend = std::min(16, maxExtend);
    const uint32_t cmpestrmFlag = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
    __m128i readSIMD = _mm_lddqu_si128((__m128i* const) &readSeq[0]);
    __m128i dontCareSIMD = _mm_set1_epi8('N');
    __m128i dontCareMaskSIMD = _mm_cmpestrm(readSIMD, initExtend, dontCareSIMD, initExtend, cmpestrmFlag);
    bool hasN = _mm_popcnt_u32(_mm_extract_epi32(dontCareMaskSIMD, 0) < initExtend);

    const unsigned int *posS, *posE;
    chromosome.getPositions(initBaseNum, start, end, &posS, &posE);

    const std::string& chromosomeSeq = chromosome.getSeq();
    if (hasN) {
        for (const unsigned int * it = posS; it != posE; it++) {
            unsigned int pos = *it;
            __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos]);
            __m128i cmpres = _mm_and_si128(dontCareMaskSIMD, _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmflag));
            unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
            if (nMismatches <= maxMismatches) {
                PD[nMismatches].push_back(pos + initExtend - 1);
            }
        }
    } else {
        for (const unsigned int * it = posS; it != posE; it++) {
            unsigned int pos = *it;
            __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos]);
            __m128i cmpres = _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmflag);
            unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
            if (nMismatches <= maxMismatches) {
                PD[nMismatches].push_back(pos + initExtend - 1);
            }
        }
    }
#else
    for (const unsigned int *it = posS; it != posE; it++) {
        unsigned int pos = *it;
        unsigned nMismatches = 0;
        for (unsigned int i = 1; i < initExtend; i++) {
            nMismatches += MismatchPair[(int) readSeq[i]][(int) chromosomeSeq[pos + i]];
        }
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

    const unsigned int *posS, *posE;
    chromosome.getPositions(initBaseNum, start, end, &posS, &posE);
    const std::string& chromosomeSeq = chromosome.getSeq();
#ifdef USE_SSE
    initExtend = std::min(16, maxExtend);
    const uint32_t cmpestrmFlag = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;

    __m128i readSIMD = _mm_setzero_si128();
    for (int i = 0; i < initExtend; i++) {
        readSIMD = _mm_slli_si128(readSIMD, 1);
        readSIMD = _mm_insert_epi8(readSIMD, readSeq[readLength - 1 - i], 0);
    }
    __m128i dontCareSIMD = _mm_set1_epi8('N');
    __m128i dontCareMaskSIMD = _mm_cmpestrm(readSIMD, initExtend, dontCareSIMD, initExtend, cmpestrmFlag);
    bool hasN = _mm_popcnt_u32(_mm_extract_epi32(dontCareMaskSIMD, 0) < initExtend);

    const unsigned int *posS, *posE;
    chromosome.getPositions(initBaseNum, start, end, &posS, &posE);

    const std::string& chromosomeSeq = chromosome.getSeq();
    if (hasN) {
        for (const unsigned int * it = posS; it != posE; it++) {
            unsigned int pos = *it;
            __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos + 1 - initExtend]);
            __m128i cmpres = _mm_and_si128(dontCareMaskSIMD, _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmflag));
            unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
            if (nMismatches <= maxMismatches) {
                PD[nMismatches].push_back(pos - initExtend + 1);
            }
        }
    } else {
        for (const unsigned int * it = posS; it != posE; it++) {
            unsigned int pos = *it;
            __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosomeSeq[pos + 1 - initExtend]);
            __m128i cmpres = _mm_cmpestrm(readSIMD, initExtend, chromosSIMD, initExtend, cmpestrmflag);
            unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
            if (nMismatches <= maxMismatches) {
                PD[nMismatches].push_back(pos - initExtend + 1);
            }
        }
    }
#else
    for (const unsigned int *it = posS; it != posE; it++) {
        unsigned int pos = *it;
        unsigned nMismatches = 0;
        for (unsigned int i = 1; i < initExtend; i++) {
            nMismatches += MismatchPair[(int) readSeq[readLength - 1 - i]][(int) chromosomeSeq[pos - i]];
        }
        if (nMismatches <= maxMismatches) {
            PD[nMismatches].push_back(pos - initExtend + 1);
        }
    }
#endif
    return initExtend;
}
