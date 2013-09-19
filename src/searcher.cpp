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

// Pindel header files
#include <emmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <x86intrin.h>
#include <assert.h>
#include "pindel.h"
#include "searcher.h"
#include <cmath>

unsigned int numberOfCompetingPositions( const std::vector < PosVector >& positions, unsigned int startIndex, unsigned int endIndex )
{
	unsigned int sum=0;	
	for (unsigned int j = startIndex; j <= endIndex; j++) {
		sum += positions[j].size();
	}	
	return sum;
}

unsigned int numberOfCompetingPositions( const std::vector < PosVector >& positions, unsigned int maxIndex )
{
	unsigned int sum=0;	
	for (unsigned int j = 0; j <= maxIndex && j < positions.size(); j++) {
   	sum += positions[j].size();
   }	
	return sum;
}

inline bool Matches(const char readBase, const char referenceBase) __attribute__((always_inline));

inline bool Matches( const char readBase, const char referenceBase )
{
        return MatchPair[readBase][referenceBase];
}

/** "CategorizePositions" categorizes the positions in PD_Plus as being extended perfectly or with an (extra) mismatch */  
inline void CategorizePositionsNormal(const char readBase,
		const std::string & chromosomeSeq,
		const PosVector & PD_Input,
		PosVector* PD_Output,
		const int searchDirection) __attribute__((always_inline));

inline void CategorizePositionsNormal(const char readBase,
		const std::string & chromosomeSeq,
		const PosVector & PD_Input,
		PosVector* PD_Output,
		const int searchDirection)
{
	for (int j = 0; j < PD_Input.size(); j++) {
		unsigned int pos = PD_Input[j] + searchDirection;
		PD_Output[MismatchPair[(int)readBase][(int)chromosomeSeq[pos]]].push_back(pos);
	}
}

void CategorizePositionsBoundary(const char readBase,
		const std::string & chromosomeSeq,
		const PosVector & PD_Input,
		PosVector* PD_Output,
		const int searchDirection)
{
	for (int j = 0; j < PD_Input.size(); j++) {
		unsigned int pos = PD_Input[j] + searchDirection;
		if ( !MismatchPair[readBase][chromosomeSeq[pos]] ) {
			PD_Output[0].push_back(pos);
		}
	}
}

inline void CategorizePositions2Normal(char readBase1, char readBase2,
		const std::string & chromosomeSeq,
		const PosVector & PD_Input,
		PosVector* PD_Output,
		const int searchDirection)
{
	for (int j = 0; j < PD_Input.size(); j++) {
		unsigned int pos1 = PD_Input[j] + searchDirection;
                unsigned int pos2 = pos1 + searchDirection;
                int index = MismatchPair[(int)readBase1][(int)chromosomeSeq[pos1]] + MismatchPair[(int)readBase2][(int)chromosomeSeq[pos2]];
		PD_Output[index].push_back(pos2);
	}
}


/** "CategorizePositions" categorizes the positions in PD_Plus as being extended perfectly or with an (extra) mismatch */  
void CategorizePositions(const char readBase, const std::string & chromosomeSeq, const std::vector<PosVector>& PD_Plus, std::vector<PosVector>& PD_Plus_Output, const int numMisMatches, 	
	const int searchDirection,	const int maxNumMismatches )
{
	int SizeOfCurrent = PD_Plus[ numMisMatches ].size();
        const PosVector& PD_Plus_current = PD_Plus[numMisMatches];
        PosVector* PD_Plus_Output_current = &PD_Plus_Output[numMisMatches];

        if (numMisMatches < maxNumMismatches) {
		for (int j = 0; j < SizeOfCurrent; j++) {
			unsigned int pos = PD_Plus_current[j] + searchDirection;
                        (PD_Plus_Output_current + MismatchPair[(int)readBase][(int)chromosomeSeq[pos]])->push_back(pos);
		}
	} else {
		for (int j = 0; j < SizeOfCurrent; j++) {
			unsigned int pos = PD_Plus_current[j] + searchDirection;
			if ( !MismatchPair[readBase][chromosomeSeq[pos]] ) {
				PD_Plus_Output_current->push_back(pos);
			}
		}
        }
}

void ExtendInPlace(char currentChar,
                   const std::string& chromosomeSeq,
                   std::vector<PosVector>& positions,
                   PosVector* TmpPositions,
		   int direction,
                   int minMismatches,
		   int maxMismatches) {
	if (minMismatches > maxMismatches) {
		return;
	}
        TmpPositions[0].clear();
        TmpPositions[1].clear();
	TmpPositions[0].swap(positions[minMismatches]);
	for (int i = minMismatches; i < maxMismatches; i++) {
		TmpPositions[1].swap(positions[i+1]); 
		CategorizePositionsNormal(currentChar, chromosomeSeq, TmpPositions[0], &positions[i], direction);
		TmpPositions[0].swap(TmpPositions[1]);
		TmpPositions[1].clear();
	}
	CategorizePositionsBoundary(currentChar, chromosomeSeq, TmpPositions[0], &positions[maxMismatches], direction);
}

void ExtendInPlace2(char char1, char char2,
                   const std::string& chromosomeSeq,
                   std::vector<PosVector>& positions,
                   PosVector* TmpPositions,
		   int direction,
                   int minMismatches,
		   int maxMismatches) {
	if (minMismatches > maxMismatches) {
		return;
	}
        //assert(positions.size() > maxMismatches + 2);
        TmpPositions[0].clear();
        TmpPositions[1].clear();
        TmpPositions[2].clear();
	TmpPositions[0].swap(positions[minMismatches]);
	TmpPositions[1].swap(positions[minMismatches+1]);
	for (int i = minMismatches; i <= maxMismatches; i++) {
		TmpPositions[2].swap(positions[i+2]); 
		CategorizePositions2Normal(char1, char2, chromosomeSeq, TmpPositions[0], &positions[i], direction);
		TmpPositions[0].swap(TmpPositions[1]);
		TmpPositions[1].swap(TmpPositions[2]);
		TmpPositions[2].clear();
	}
}

void ExtendMatchClose( SPLIT_READ & read, const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector<PosVector>& InputPositions, const short minimumLengthToReportMatch,
		const short BP_End, short CurrentLength,
		SortedUniquePoints &UP, int direction )
{
        if (InputPositions.size() == 0) {
            return;
        }

	{
		const char CurrentChar = ((direction==1) ? readSeq[CurrentLength] : readSeq[read.getReadLengthMinus() - CurrentLength] );
                PosVector TmpPositions[2];
                ExtendInPlace(CurrentChar, chromosomeSeq, InputPositions, TmpPositions, direction, 0, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
	}

	unsigned int Sum = numberOfCompetingPositions( InputPositions, read.getMAX_SNP_ERROR() );

	if (Sum) {
		if ( direction==1 ) {
			CheckLeft_Close(read, chromosomeSeq, readSeq, InputPositions, minimumLengthToReportMatch, BP_End, CurrentLength + 1, UP);
		}
		else {
			CheckRight_Close(read, chromosomeSeq, readSeq, InputPositions, minimumLengthToReportMatch, BP_End, CurrentLength + 1, UP);	
		}
	}
}

void ExtendMatchClosePerfect( SPLIT_READ & read, const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector<PosVector>& InputPositions, const short minimumLengthToReportMatch,
		const short BP_End, const short CurrentLength,
		SortedUniquePoints &UP, int direction )
{   //std::cout << "in ExtendMatchClosePerfect " << std::endl;
	if (InputPositions[0].size() == 0) return;

	const char CurrentChar = ((direction==1) ? readSeq[CurrentLength] : readSeq[read.getReadLengthMinus() - CurrentLength] );
	{
		const std::vector<PosVector> TmpPositions = InputPositions;
		for (int i = 0; i < TmpPositions.size(); i++) { InputPositions[i].clear(); }
		for (int i = 0; i <= userSettings->ADDITIONAL_MISMATCH; i++) {
			CategorizePositions( CurrentChar, chromosomeSeq, TmpPositions, InputPositions, i, direction, userSettings->ADDITIONAL_MISMATCH);
		}
	}
	unsigned int Sum = numberOfCompetingPositions( InputPositions, userSettings->ADDITIONAL_MISMATCH + 1 );
	if (Sum) {
		if ( direction==1 ) {
			CheckLeft_Close_Perfect(read, chromosomeSeq, readSeq, InputPositions, minimumLengthToReportMatch, BP_End, CurrentLength + 1, UP);
		}
		else {
			CheckRight_Close_Perfect(read, chromosomeSeq, readSeq, InputPositions, minimumLengthToReportMatch, BP_End, CurrentLength + 1, UP);
		}
	}
}

unsigned int minimumNumberOfMismatches( const std::vector< PosVector>& mismatches, const unsigned int maxNumberMismatches )
{
	unsigned int numberOfMismatches=0; 
	for (;numberOfMismatches<=maxNumberMismatches; numberOfMismatches++ ) {
		if ( mismatches[ numberOfMismatches ].size() != 0 ) { break; }
	}
	return numberOfMismatches;
}


void CheckLeft_Close (SPLIT_READ & read,
		const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector< PosVector >& Left_PD,
		const short &BP_Left_Start,
		const short &BP_Left_End,
		short CurrentLength, SortedUniquePoints &LeftUP)
{
	int minMismatches = 0;
        int matchRunLength = BP_Left_End;
        PosVector TmpPositions[3];
	while (Left_PD[minMismatches].size() == 0 && minMismatches <= read.getMAX_SNP_ERROR()) {
		minMismatches++;
	}

        for ( ; CurrentLength < BP_Left_Start && minMismatches <= read.getMAX_SNP_ERROR(); CurrentLength++) {
		ExtendInPlace(readSeq[CurrentLength], chromosomeSeq, Left_PD, TmpPositions, 1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());

		minMismatches += (Left_PD[minMismatches].size() == 0)? 1: 0;
        }

	for ( ; CurrentLength <= BP_Left_End && minMismatches <= read.getMAX_SNP_ERROR(); CurrentLength++) {
		if (minMismatches > g_maxMismatch[CurrentLength] ) {
			return; 
		}

		if (minMismatches <= CurrentLength - BP_Left_Start && Left_PD[minMismatches].size() == 1 && matchRunLength >= userSettings->Min_Perfect_Match_Around_BP) {
			const std::string& forwardSeq = read.getUnmatchedSeq();
			const std::string& reverseSeq = read.getUnmatchedSeqRev();
			unsigned int Sum = numberOfCompetingPositions( Left_PD, minMismatches + 1, minMismatches + userSettings->ADDITIONAL_MISMATCH );

			if (Sum == 0) {
				UniquePoint TempOne(g_genome.getChr(read.FragId), CurrentLength, Left_PD[minMismatches][0], FORWARD, ANTISENSE, minMismatches);  
				if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
					LeftUP.push_back (TempOne);
				}
			}
		}

                if (Left_PD[minMismatches].size() > 1 && BP_Left_End - CurrentLength < userSettings->ADDITIONAL_MISMATCH + 1) {
			return;
                }
		if (CurrentLength < BP_Left_End) {
			ExtendInPlace(readSeq[CurrentLength], chromosomeSeq, Left_PD, TmpPositions, 1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
		}
                if (Left_PD[minMismatches].size() == 0) {
                    minMismatches++;
                    matchRunLength = 1;
                } else {
                    matchRunLength++;
                }
		//minMismatches += (Left_PD[minMismatches].size() == 0)? 1: 0;
	}
}

void CheckLeft_Close_Perfect (SPLIT_READ & read,
		const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector< PosVector >& Left_PD,
		const short &BP_Left_Start,
		const short &BP_Left_End,
		short CurrentLength, SortedUniquePoints &LeftUP)
{
        PosVector TmpPositions[2];
        for ( ; CurrentLength < BP_Left_Start && Left_PD[0].size() > 0; CurrentLength++) {
		ExtendInPlace(readSeq[CurrentLength] , chromosomeSeq, Left_PD, TmpPositions, 1, 0, userSettings->ADDITIONAL_MISMATCH);
	}

        for ( ; CurrentLength <= BP_Left_End && Left_PD[0].size() > 0; CurrentLength++) {
		if (Left_PD[0].size() == 1) {
			if (numberOfCompetingPositions(Left_PD, 1, userSettings->ADDITIONAL_MISMATCH) == 0) {
				const std::string& forwardSeq = read.getUnmatchedSeq();
				const std::string& reverseSeq = read.getUnmatchedSeqRev();
				UniquePoint TempOne( g_genome.getChr(read.FragId), CurrentLength, Left_PD[0][0], FORWARD, ANTISENSE, 0);
				if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
					LeftUP.push_back (TempOne);
				} 
			}
		}

		if (CurrentLength < BP_Left_End) {
			ExtendInPlace(readSeq[CurrentLength] , chromosomeSeq, Left_PD, TmpPositions, 1, 0, userSettings->ADDITIONAL_MISMATCH);
		}
	}
}

void CheckRight_Close (SPLIT_READ & read,
		const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector < PosVector >& Right_PD,
		const short &BP_Right_Start,
		const short &BP_Right_End,
		short CurrentLength, SortedUniquePoints &RightUP)
{
        PosVector TmpPositions[2];
	int minMismatches = 0;
	while (Right_PD[minMismatches].size() == 0 && minMismatches <= read.getMAX_SNP_ERROR()) {
		minMismatches++;
	}

        int matchRunLength = BP_Right_Start;
        for ( ; CurrentLength < BP_Right_Start && minMismatches <= read.getMAX_SNP_ERROR(); CurrentLength++) {
			ExtendInPlace(readSeq[read.getReadLengthMinus() - CurrentLength] , chromosomeSeq, Right_PD, TmpPositions, -1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
			
			minMismatches += (Right_PD[minMismatches].size() == 0)? 1: 0;	
        }

	for ( ; CurrentLength <= BP_Right_End && minMismatches <= read.getMAX_SNP_ERROR(); CurrentLength++) {
		if (minMismatches > g_maxMismatch[CurrentLength] ) {
			return; 
		}

		if (minMismatches <= CurrentLength - BP_Right_Start && Right_PD[minMismatches].size() == 1 && matchRunLength >= userSettings->Min_Perfect_Match_Around_BP) {
			const std::string& forwardSeq = read.getUnmatchedSeq();
			const std::string& reverseSeq = read.getUnmatchedSeqRev();
			unsigned int Sum = numberOfCompetingPositions(Right_PD, minMismatches+1, minMismatches+userSettings->ADDITIONAL_MISMATCH );
			if (Sum == 0) {
				UniquePoint TempOne( g_genome.getChr(read.FragId), CurrentLength, Right_PD[minMismatches][0], BACKWARD, SENSE, minMismatches);
				if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
					RightUP.push_back (TempOne);
				} // ###################################
			}
		}

		if (CurrentLength < BP_Right_End) {
			ExtendInPlace(readSeq[read.getReadLengthMinus() - CurrentLength] , chromosomeSeq, Right_PD, TmpPositions, -1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
		}
		//minMismatches += (Right_PD[minMismatches].size() == 0)? 1: 0;
                if (Right_PD[minMismatches].size() == 0) {
                    minMismatches++;
                    matchRunLength = 1;
                } else {
                    matchRunLength++;
                }
	
	}
}

void CheckRight_Close_Perfect (SPLIT_READ & read,
		const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector < PosVector >& Right_PD,
		const short &BP_Right_Start,
		const short &BP_Right_End,
		short CurrentLength, SortedUniquePoints &RightUP)
{
        PosVector TmpPositions[2];
        for ( ; CurrentLength < BP_Right_Start && Right_PD[0].size() > 0; CurrentLength++) {
		ExtendInPlace(readSeq[read.getReadLengthMinus() - CurrentLength] , chromosomeSeq, Right_PD, TmpPositions, -1, 0, userSettings->ADDITIONAL_MISMATCH);
        }

        for ( ; CurrentLength <= BP_Right_End && Right_PD[0].size() > 0; CurrentLength++) {
		if (Right_PD[0].size() == 1) {
			if(numberOfCompetingPositions(Right_PD, 1, userSettings->ADDITIONAL_MISMATCH) == 0) {
				const std::string& forwardSeq = read.getUnmatchedSeq();
				const std::string& reverseSeq = read.getUnmatchedSeqRev();
				UniquePoint TempOne( g_genome.getChr(read.FragId), CurrentLength, Right_PD[0][0], BACKWARD, SENSE, 0);
				if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
					RightUP.push_back (TempOne);
				} // ###################################
			}
		}

		if (CurrentLength < BP_Right_End) {
			ExtendInPlace(readSeq[read.getReadLengthMinus() - CurrentLength] , chromosomeSeq, Right_PD, TmpPositions, -1, 0, userSettings->ADDITIONAL_MISMATCH);
		}
	}
}

bool CheckMismatches (const std::string & TheInput, const std::string & InputReadSeq, const std::string& InputReadSeqRev, const UniquePoint & UP, short & numberOfMismatch)
{
	int Min_Perfect_Match_Around_BP = userSettings->Min_Perfect_Match_Around_BP;
   //std::string CurrentReadSeq;
   const std::string* CurrentReadSeq = (UP.Strand == SENSE)? &InputReadSeq: &InputReadSeqRev;
   short CurrentReadLength = CurrentReadSeq->size ();
   unsigned int Start = 0;
   //std::string BP_On_Read, BP_On_Ref;
   if (UP.Direction == FORWARD) {

      Start = UP.AbsLoc - UP.LengthStr + 1;
      if (UP.LengthStr <= Min_Perfect_Match_Around_BP) {
         return false;
      }
      const char* BP_On_Read_c = &(*CurrentReadSeq)[UP.LengthStr - Min_Perfect_Match_Around_BP];
      const char* BP_On_Ref_c = &TheInput[UP.AbsLoc - Min_Perfect_Match_Around_BP + 1];
      if (strncmp(BP_On_Read_c, BP_On_Ref_c, Min_Perfect_Match_Around_BP) != 0) {
        return false;
      }
   }
   else if (UP.Direction == BACKWARD) {
      Start = UP.AbsLoc + UP.LengthStr - CurrentReadLength;
      if (CurrentReadLength < UP.LengthStr) {
         return false;
      }
      if (CurrentReadLength - UP.LengthStr + Min_Perfect_Match_Around_BP > CurrentReadLength) {
         return false;
      }
      if (UP.AbsLoc + Min_Perfect_Match_Around_BP > TheInput.size()) {
          return false;
      }
      const char* BP_On_Read_c = &(*CurrentReadSeq)[CurrentReadLength - UP.LengthStr];
      const char* BP_On_Ref_c = &TheInput[UP.AbsLoc];
      if (strncmp(BP_On_Read_c, BP_On_Ref_c, Min_Perfect_Match_Around_BP) != 0) {
        return false;
      }
   }
   float MAX_ALLOWED_MISMATCHES = CurrentReadSeq->size () * userSettings->MaximumAllowedMismatchRate;	//

   short NumMismatches = 0;			// Match2N[(short)'A'] = 'N';
   const uint32_t cmpestrmflag       = _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | _SIDD_NEGATIVE_POLARITY;
   __m128i dontcarSIMD = _mm_set1_epi8('N');

   int CurrentReadLength_a = CurrentReadLength - (CurrentReadLength % 16);
   for (short i = 0; i < CurrentReadLength_a; i+=16) {
       int toProcess = 16;
       __m128i readSIMD = _mm_lddqu_si128((__m128i* const) &(*CurrentReadSeq)[i]); // TODO: fix potential seg fault
       __m128i inputSIMD = _mm_lddqu_si128((__m128i* const) &TheInput[Start+i]);
       __m128i readMaskSIMD = _mm_cmpestrm(readSIMD, toProcess, dontcarSIMD, toProcess, cmpestrmflag);
       __m128i cmpres = _mm_and_si128(readMaskSIMD, _mm_cmpestrm(readSIMD, toProcess, inputSIMD, toProcess, cmpestrmflag));
       NumMismatches += _mm_popcnt_u32(_mm_extract_epi32(cmpres, 0));
   }
   for (int i = CurrentReadLength_a; i < CurrentReadLength; i++) {
       NumMismatches += MismatchPair[(int)(*CurrentReadSeq)[i]][(int)TheInput[Start+i]];
   }
   numberOfMismatch = NumMismatches;
   // std::cout << "NumMismatches > MAX_ALLOWED_MISMATCHES " << NumMismatches << " " << MAX_ALLOWED_MISMATCHES << std::endl;
   if ((float)NumMismatches >= MAX_ALLOWED_MISMATCHES) {
      return true;
   }
   else {
      return false;
   }
}
