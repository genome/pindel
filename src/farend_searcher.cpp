/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by  
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <x86intrin.h>

#include "pindel.h"
#include "searcher.h"
#include "farend_searcher.h"




bool NewUPFarIsBetter(const SortedUniquePoints & UP, const SPLIT_READ& Read) {

    if (UP.MaxLen() < Read.MaxLenFarEnd()) {
       return false;
    }
    if (UP.NumMismatch() < Read.UP_Far.NumMismatch()){
        return true;
    }
    else {
        return false;
    }
}

void FillForward(const unsigned int* positions,
                 int npositions,
                 const __m128i forwardSIMD,
                 const __m128i forwardMaskSIMD,
                 int InitExtend,
                 const std::string& chromosome,
                 std::vector<PosVector>& PD) {
        int npositions_a = npositions - (npositions % 4);
	const uint32_t cmpestrmflag = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
        __m128i InitExtendSIMD = _mm_set1_epi32(InitExtend - 1);
	for (int i = 0; i < npositions_a; i+=4) {
		__m128i posSIMD = _mm_lddqu_si128((__m128i const*) &positions[i]);
		__m128i posExtendedSIMD = _mm_add_epi32(posSIMD, InitExtendSIMD);

		{
			__m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[_mm_extract_epi32(posSIMD, 0)]);
			__m128i cmpres = _mm_and_si128(forwardMaskSIMD, _mm_cmpestrm(forwardSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
			unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
			PD[nMismatches].push_back(_mm_extract_epi32(posExtendedSIMD, 0)); // else
		}
		{
			__m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[_mm_extract_epi32(posSIMD, 1)]);
			__m128i cmpres = _mm_and_si128(forwardMaskSIMD, _mm_cmpestrm(forwardSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
			unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
			PD[nMismatches].push_back(_mm_extract_epi32(posExtendedSIMD, 1)); // else
		}
		{
			__m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[_mm_extract_epi32(posSIMD, 2)]);
			__m128i cmpres = _mm_and_si128(forwardMaskSIMD, _mm_cmpestrm(forwardSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
			unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
			PD[nMismatches].push_back(_mm_extract_epi32(posExtendedSIMD, 2)); // else
		}
		{
			__m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[_mm_extract_epi32(posSIMD, 3)]);
			__m128i cmpres = _mm_and_si128(forwardMaskSIMD, _mm_cmpestrm(forwardSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
			unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
			PD[nMismatches].push_back(_mm_extract_epi32(posExtendedSIMD, 3)); // else
		}
	}
	for (int i = npositions_a; i < npositions; i++) {
                int pos = positions[i];
		__m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[pos]);
		__m128i cmpres = _mm_and_si128(forwardMaskSIMD, _mm_cmpestrm(forwardSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
		unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
		PD[nMismatches].push_back(pos + InitExtend - 1); // else
        }
}

void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions )
{
	const std::string& forwardSeq = Temp_One_Read.getUnmatchedSeq();
        const std::string& reverseSeq = Temp_One_Read.getUnmatchedSeqRev();
        const unsigned readLength = forwardSeq.size();

	// step 1 find out which chromosomes in Regions: set? linear pass of regions
	// step 2 for each identified chromsome, for each regions on the chromosme, do the business.
	const char CurrentBase = forwardSeq[0];
	const char CurrentBaseRC = reverseSeq[readLength - 1];

	if (CurrentBase == 'N' || Temp_One_Read.MaxLenCloseEnd() == 0) return;
	//int CurrentReadLength = Temp_One_Read.getReadLength();
	std::vector <FarEndSearchPerRegion*> WholeGenomeSearchResult;
	unsigned NumberOfHits = 0;
        const int InitExtend = 10;

        if (forwardSeq.size() < InitExtend) {
            return;
        }
        int CurrentBaseNum = Convert2Num[CurrentBase];
        int CurrentBaseRCNum = Convert2Num[CurrentBaseRC];

	const uint32_t cmpestrmflag = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
	__m128i forwardSIMD = _mm_lddqu_si128((__m128i* const) &forwardSeq[0]);

        __m128i reverseSIMD = _mm_setzero_si128();
        for (int i = 0; i < InitExtend; i++) {
            reverseSIMD = _mm_slli_si128(reverseSIMD, 1);
            reverseSIMD = _mm_insert_epi8(reverseSIMD, reverseSeq[readLength - 1 - i], 0);
        }
	__m128i dontcarSIMD = _mm_set1_epi8('N');

	__m128i forwardMaskSIMD = _mm_cmpestrm(forwardSIMD, InitExtend, dontcarSIMD, InitExtend, cmpestrmflag);
	__m128i reverseMaskSIMD = _mm_cmpestrm(reverseSIMD, InitExtend, dontcarSIMD, InitExtend, cmpestrmflag);

        PosVector TmpPos[6];

        unsigned PD_size = std::max(InitExtend, (int)Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED());
	for (unsigned RegionIndex = 0; RegionIndex < Regions.size(); RegionIndex++) {
		FarEndSearchPerRegion* CurrentRegion = new FarEndSearchPerRegion(Regions[RegionIndex].getChromosome(), PD_size, std::max(1, (int)Regions[RegionIndex].getSize()/4));
		const std::string & chromosome = Regions[RegionIndex].getChromosome()->getSeq();

		int Start = Regions[RegionIndex].getStart();
		int End = std::min((unsigned) Regions[RegionIndex].getEnd(), (unsigned) chromosome.size());
		if (Start < 0) Start = End -1;

		for (int i = 0; i < 6; i++) {
			TmpPos[i].clear();
			TmpPos[i].reserve(std::max(End - Start, 1));
		}
		for (int pos = Start; pos < End; pos++) {
			TmpPos[Convert2Num[chromosome[pos]]].push_back(pos);
		}
                //FillForward(&TmpPos[CurrentBaseNum][0], TmpPos[CurrentBaseNum].size(), forwardSIMD, forwardMaskSIMD, InitExtend, chromosome, CurrentRegion->PD_Plus);
                 
		for (PosVector::iterator it = TmpPos[CurrentBaseNum].begin(); it != TmpPos[CurrentBaseNum].end(); it++) {
			int pos = *it;
			__m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[pos]);
			__m128i cmpres = _mm_and_si128(forwardMaskSIMD, _mm_cmpestrm(forwardSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
			unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
			CurrentRegion->PD_Plus[nMismatches].push_back(pos + InitExtend - 1); // else
		}
		for (PosVector::iterator it = TmpPos[CurrentBaseRCNum].begin(); it != TmpPos[CurrentBaseRCNum].end(); it++) {
			int pos = *it;
			__m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[pos + 1 - InitExtend]);
			__m128i cmpres = _mm_and_si128(reverseMaskSIMD, _mm_cmpestrm(reverseSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
			unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
			CurrentRegion->PD_Minus[nMismatches].push_back(pos - InitExtend + 1); // else
		}


		/*for (int pos = Start; pos < End; pos++) {
			if (chromosome[pos] == CurrentBase) {
                                // TODO: Make the SIMD match the MismatchPair based code
                                __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[pos]);
                                __m128i cmpres = _mm_and_si128(forwardMaskSIMD, _mm_cmpestrm(forwardSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
				unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
				CurrentRegion->PD_Plus[nMismatches].push_back(pos + InitExtend - 1); // else
				//CurrentRegion->PD_Plus[0].push_back(pos);
			}
			if (chromosome[pos] == CurrentBaseRC) {
                                __m128i chromosSIMD = _mm_lddqu_si128((__m128i* const) &chromosome[pos + 1 - InitExtend]);
                                __m128i cmpres = _mm_and_si128(reverseMaskSIMD, _mm_cmpestrm(reverseSIMD, InitExtend, chromosSIMD, InitExtend, cmpestrmflag));
				unsigned nMismatches = _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0)); 
				CurrentRegion->PD_Minus[nMismatches].push_back(pos - InitExtend + 1); // else
                                //CurrentRegion->PD_Minus[0].push_back(pos);
			}
		}*/
                for (unsigned i = 0; i < Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(); i++) {
			NumberOfHits += CurrentRegion->PD_Plus[i].size();
                }
                for (unsigned i = 0; i < Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(); i++) {
			NumberOfHits += CurrentRegion->PD_Minus[i].size();
                }
		WholeGenomeSearchResult.push_back(CurrentRegion);
	}

	if (NumberOfHits>0) {
		short BP_Start = 10; // perhaps use global constant like "g_MinimumLengthToReportMatch"
		short BP_End = Temp_One_Read.getReadLengthMinus(); // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
		SortedUniquePoints UP; // temporary container for unique far ends
		CheckBoth(Temp_One_Read, forwardSeq, reverseSeq, WholeGenomeSearchResult, BP_Start, BP_End, InitExtend, UP);

		if ( NewUPFarIsBetter(UP, Temp_One_Read)) { // UP.MaxLen() > Temp_One_Read.MaxLenFarEnd()
			Temp_One_Read.UP_Far.swap(UP);
		}
		UP.clear(); // may not be necessary as this is deleted from the stack anyway
	}
        for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult.size(); RegionIndex++) {
		delete WholeGenomeSearchResult[ RegionIndex ];
	}
}

void SearchFarEndAtPosPerfect( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions )
{
	// step 1 find out which chromosomes in Regions: set? linear pass of regions
	// step 2 for each identified chromsome, for each regions on the chromosme, do the business.
	const char CurrentBase = Temp_One_Read.getUnmatchedSeq()[0];
	const char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];

	if (CurrentBase == 'N' || Temp_One_Read.MaxLenCloseEnd() == 0) return;
	//int CurrentReadLength = Temp_One_Read.getReadLength();

	std::vector <FarEndSearchPerRegion*> WholeGenomeSearchResult;
	unsigned NumberOfHits = 0;


	for (unsigned RegionIndex = 0; RegionIndex < Regions.size(); RegionIndex++) {

		FarEndSearchPerRegion* CurrentRegion = new FarEndSearchPerRegion(Regions[RegionIndex].getChromosome(), Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(), Regions[RegionIndex].getSize());
		const std::string & chromosome = Regions[RegionIndex].getChromosome()->getSeq();

		int Start = Regions[RegionIndex].getStart();
		int End = std::min((unsigned) Regions[RegionIndex].getEnd(), (unsigned) chromosome.size());
		if (Start < 0) Start = End -1;
		for (int pos = Start; pos < End; pos++) {
			if (chromosome[pos] == CurrentBase) {
				CurrentRegion->PD_Plus[0].push_back(pos); // else
			}
			if (chromosome[pos] == CurrentBaseRC) {
				CurrentRegion->PD_Minus[0].push_back(pos);
			}
		}
		NumberOfHits += CurrentRegion->PD_Plus[0].size() + CurrentRegion->PD_Minus[0].size();
		WholeGenomeSearchResult.push_back(CurrentRegion);
	}

	if (NumberOfHits>0) {
		short BP_Start = 10; // perhaps use global constant like "g_MinimumLengthToReportMatch"
		short BP_End = Temp_One_Read.getReadLengthMinus(); // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
		SortedUniquePoints UP; // temporary container for unique far ends
		CheckBothPerfect(Temp_One_Read, Temp_One_Read.getUnmatchedSeq(), Temp_One_Read.getUnmatchedSeqRev(), WholeGenomeSearchResult, BP_Start, BP_End, 1, UP);

		if ( NewUPFarIsBetter(UP, Temp_One_Read)) { // UP.MaxLen() > Temp_One_Read.MaxLenFarEnd()
			Temp_One_Read.UP_Far.swap(UP);
		}
		UP.clear(); // may not be necessary as this is deleted from the stack anyway
	}
        for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult.size(); RegionIndex++) {
		delete WholeGenomeSearchResult[ RegionIndex ];
	}
}


/*void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions )
{
   std::vector<unsigned int> PD_Plus[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
   std::vector<unsigned int> PD_Minus[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
    
   int TotalSize = 0;
   for (unsigned int RegionIndex = 0; RegionIndex < Regions.size(); RegionIndex++) {
      TotalSize += Regions[RegionIndex].getEnd() - Regions[RegionIndex].getStart() + 1;        
   }
    
   for (int CheckIndex = 0; CheckIndex < Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(); CheckIndex++) {
      PD_Plus[CheckIndex].reserve(TotalSize);
      PD_Minus[CheckIndex].reserve(TotalSize);
   }
    
   char CurrentBase = Temp_One_Read.getUnmatchedSeq()[0];
   char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];
	
	if (CurrentBase == 'N' || Temp_One_Read.MaxLenCloseEnd() == 0) return; 
   int CurrentReadLength = Temp_One_Read.getReadLength();
   for (unsigned RegionIndex = 0; RegionIndex < Regions.size(); RegionIndex++) {
      int Start = Regions[RegionIndex].getStart() - CurrentReadLength;
      int End = Regions[RegionIndex].getEnd() + CurrentReadLength;
      for (int pos = Start; pos < End; pos++) {
         if (chromosome.at(pos) == CurrentBase) {
            PD_Plus[0].push_back(pos); // else
         }
         else if (chromosome.at(pos) == CurrentBaseRC) {
            PD_Minus[0].push_back(pos);
         }
      }
   }

	if (PD_Minus[0].size() + PD_Plus[0].size() > 0) {
		short BP_Start = 10; // perhaps use global constant like "g_MinimumLengthToReportMatch"
      short BP_End = Temp_One_Read.getReadLengthMinus(); // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
      SortedUniquePoints UP; // temporary container for unique far ends
      CheckBoth(Temp_One_Read, chromosome, Temp_One_Read.getUnmatchedSeq(), PD_Plus, PD_Minus, BP_Start, BP_End, 1, UP);
        
		if ( UP.MaxLen() > Temp_One_Read.MaxLenFarEnd() ) {
			Temp_One_Read.UP_Far.swap(UP);
		}
      UP.clear(); // may not be necessary as this is deleted from the stack anyway
	}
}*/

