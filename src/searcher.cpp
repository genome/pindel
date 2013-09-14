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

void ExtendMatchClose( SPLIT_READ & read, const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector<PosVector>& InputPositions, const short minimumLengthToReportMatch,
		const short BP_End, const short CurrentLength,
		SortedUniquePoints &UP, int direction )
{
	const char CurrentChar = ((direction==1) ? readSeq[CurrentLength] : readSeq[read.getReadLengthMinus() - CurrentLength] );
	{
		const std::vector<PosVector> TmpPositions = InputPositions;
		for (unsigned int i = 0; i < TmpPositions.size(); i++) { InputPositions[i].clear(); }
		for (int i = 0; i <= read.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) { 
			CategorizePositions( CurrentChar, chromosomeSeq, TmpPositions, InputPositions, i, direction, read.getTOTAL_SNP_ERROR_CHECKED_Minus() );
		}
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
		const short CurrentLength, SortedUniquePoints &LeftUP)
{
	if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
		if (minimumNumberOfMismatches( Left_PD,read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
			return; 
		}

		const std::string& forwardSeq = read.getUnmatchedSeq();
		const std::string& reverseSeq = read.getUnmatchedSeqRev();
		const int g_maxMismatch_ = std::min((int)std::min(read.getMAX_SNP_ERROR(), (short) g_maxMismatch[CurrentLength]), CurrentLength - BP_Left_Start);
		// put it to LeftUP if unique
		for (short i = 0; i <= g_maxMismatch_; i++) {
			if (Left_PD[i].size() == 1) {
				unsigned int Sum = numberOfCompetingPositions( Left_PD, i + 1, i + userSettings->ADDITIONAL_MISMATCH );

				if (Sum == 0) {
					UniquePoint TempOne(g_genome.getChr(read.FragName), CurrentLength, Left_PD[i][0], FORWARD, ANTISENSE, i );  
					if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
						LeftUP.push_back (TempOne);
					}
				}
			}
                        if (Left_PD[i].size() > 0) {
                            break;
                        }
		}
	}
	if (CurrentLength < BP_Left_End) {
		ExtendMatchClose( read, chromosomeSeq, readSeq, Left_PD, BP_Left_Start, BP_Left_End, CurrentLength, LeftUP, 1 );
	}
}

void CheckLeft_Close_Perfect (SPLIT_READ & read,
		const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector< PosVector >& Left_PD,
		const short &BP_Left_Start,
		const short &BP_Left_End,
		const short CurrentLength, SortedUniquePoints &LeftUP)
{
	if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
		if (minimumNumberOfMismatches( Left_PD,read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
			return;
		}
		const std::string& forwardSeq = read.getUnmatchedSeq();
		const std::string& reverseSeq = read.getUnmatchedSeqRev();
		// put it to LeftUP if unique
		if (Left_PD[0].size() == 1) {
			unsigned int Sum = numberOfCompetingPositions( Left_PD, 1, userSettings->ADDITIONAL_MISMATCH);

			if (Sum == 0) {
				UniquePoint TempOne(g_genome.getChr(read.FragName), CurrentLength, Left_PD[0][0], FORWARD, ANTISENSE, 0 );
				if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
					LeftUP.push_back (TempOne);
				}
			}
		}
	}
	if (CurrentLength < BP_Left_End) {
		ExtendMatchClosePerfect( read, chromosomeSeq, readSeq, Left_PD, BP_Left_Start, BP_Left_End, CurrentLength, LeftUP, 1 );
	}
}

void CheckRight_Close (SPLIT_READ & read,
		const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector < PosVector >& Right_PD,
		const short &BP_Right_Start,
		const short &BP_Right_End,
		const short CurrentLength, SortedUniquePoints &RightUP)
{
	if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End) {
		if (minimumNumberOfMismatches( Right_PD,read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
			return; 
		}
		const std::string& forwardSeq = read.getUnmatchedSeq();
		const std::string& reverseSeq = read.getUnmatchedSeqRev();
		const int g_maxMismatch_ = std::min((int)std::min(read.getMAX_SNP_ERROR(), (short)g_maxMismatch[CurrentLength]), CurrentLength - BP_Right_Start);
		for (short i = 0; i <= g_maxMismatch_; i++) {
			if (Right_PD[i].size () == 1) {
				unsigned int Sum = numberOfCompetingPositions( Right_PD, i+1, i+userSettings->ADDITIONAL_MISMATCH );
				if (Sum == 0) {
					UniquePoint TempOne( g_genome.getChr(read.FragName), CurrentLength, Right_PD[i][0], BACKWARD, SENSE, i);
					if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
						RightUP.push_back (TempOne);
					} // ###################################
				}
			}
                        if (Right_PD[i].size() > 0) {
				break;
			}
		}
	}

	if (CurrentLength < BP_Right_End) {
		ExtendMatchClose( read, chromosomeSeq, readSeq, Right_PD, BP_Right_Start, BP_Right_End, CurrentLength, RightUP, -1 );
	}
}

void CheckRight_Close_Perfect (SPLIT_READ & read,
		const std::string & chromosomeSeq,
		const std::string & readSeq,
		std::vector < PosVector >& Right_PD,
		const short &BP_Right_Start,
		const short &BP_Right_End,
		const short CurrentLength, SortedUniquePoints &RightUP)
{
	if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End) {
		if (minimumNumberOfMismatches( Right_PD,read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
			return;
		}
		const std::string& forwardSeq = read.getUnmatchedSeq();
		const std::string& reverseSeq = read.getUnmatchedSeqRev();
		if (Right_PD[0].size () == 1) {
			unsigned int Sum = numberOfCompetingPositions( Right_PD, 1, userSettings->ADDITIONAL_MISMATCH );
			if (Sum == 0) {
				UniquePoint TempOne( g_genome.getChr(read.FragName), CurrentLength, Right_PD[0][0], BACKWARD, SENSE, 0);
				if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
					RightUP.push_back (TempOne);
				} // ###################################
			}
		}
	}

	if (CurrentLength < BP_Right_End) {
		ExtendMatchClosePerfect( read, chromosomeSeq, readSeq, Right_PD, BP_Right_Start, BP_Right_End, CurrentLength, RightUP, -1 );
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

   for (short i = 0; i < CurrentReadLength; i++) {
      char CurrentReadSeqChar = (*CurrentReadSeq)[i];
      NumMismatches += MismatchPair[CurrentReadSeqChar][TheInput[Start+i]];
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
