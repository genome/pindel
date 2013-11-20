#ifndef _SSE_HELPERS_H
#define _SSE_HELPERS_H

#include <string>
#include <vector>
#include "pindel.h"

int CountMismatches(const char* read,
                    const char* reference,
                    int length);
int DoInitialSeedAndExtendForward(const Chromosome& chromosome,
                                  int start,
                                  int end,
                                  size_t initExtend,
                                  const std::string& readSeq,
                                  size_t maxMismatches,
                                  std::vector<unsigned int>* PD);

int DoInitialSeedAndExtendReverse(const Chromosome& chromosome,
                                  int start,
                                  int end,
                                  size_t initExtend,
                                  const std::string& readSeq,
                                  size_t maxMismatches,
                                  std::vector<unsigned int>* PD);
#endif
