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

unsigned int numberOfCompetingPositions( const std::vector < PosVector >& positions, unsigned int maxIndex )
{
	unsigned int sum=0;	
	for (unsigned int j = 0; j <= maxIndex; j++) {
   	sum += positions[j].size();
   }	
	return sum;
}


bool Matches( const char readBase, const char referenceBase )
{
	//std::cout << "read/ref: " << readBase << "," << referenceBase << "\n";
	if (readBase!='N') { return referenceBase == readBase; }
	else { return Match2N[(short) referenceBase] == 'N'; } 
}


/** "CategorizePositions" categorizes the positions in PD_Plus as being extended perfectly or with an (extra) mismatch */  
void CategorizePositions(const char readBase, const std::string & chromosomeSeq, const std::vector<PosVector>& PD_Plus, std::vector<PosVector>& PD_Plus_Output, const int numMisMatches, 	
	const int searchDirection,	const int maxNumMismatches )
{
	int SizeOfCurrent = PD_Plus[ numMisMatches ].size();
 	for (int j = 0; j < SizeOfCurrent; j++) {
      unsigned int pos = PD_Plus[ numMisMatches ][j] + searchDirection;
      if ( Matches( readBase, chromosomeSeq[ pos ] ) ) {
         PD_Plus_Output[ numMisMatches ].push_back(pos);
      }
      else {
         if ( numMisMatches<maxNumMismatches) { PD_Plus_Output[ numMisMatches + 1].push_back(pos); }
      } 
   }
	//std::cout << "From " << PD_Plus[ numMisMatches ].size() << "Good: " << PD_Plus_Output[ numMisMatches ].size() << "Bad: " << PD_Plus_Output[ numMisMatches + 1].size() << std::endl;
}

void ExtendMatchClose( const SPLIT_READ & read, const std::string & chromosomeSeq,
               const std::string & readSeq,
               const std::vector<PosVector> InputPositions, const short minimumLengthToReportMatch,
               const short BP_End, const short CurrentLength,
               SortedUniquePoints &UP, int direction )
{
	std::vector<PosVector> OutputPositions;
	PosVector emptyPosVector;
   OutputPositions.assign( read.getTOTAL_SNP_ERROR_CHECKED(), emptyPosVector);
		
   for (int CheckedIndex = 0; CheckedIndex < read.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
      OutputPositions[CheckedIndex].reserve( InputPositions[CheckedIndex].size()); // this assumes perfect matches and no 'attrition' from higher levels. We may want to test this...
   }
   const char CurrentChar = ((direction==1) ? readSeq[CurrentLength] : readSeq[read.getReadLengthMinus() - CurrentLength] );
		//std::cout << "Matching " << CurrentChar << "\n";
	for (int i = 0; i <= read.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) { 

		CategorizePositions( CurrentChar, chromosomeSeq, InputPositions, OutputPositions, i, direction, read.getTOTAL_SNP_ERROR_CHECKED_Minus() );
   }

   unsigned int Sum = numberOfCompetingPositions( OutputPositions, read.getMAX_SNP_ERROR() );

   if (Sum) {
      const short CurrentLengthOutput = CurrentLength + 1;
		if ( direction==1 ) {
	      CheckLeft_Close(read, chromosomeSeq, readSeq, OutputPositions, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);
		}
		else {
			CheckRight_Close(read, chromosomeSeq, readSeq, OutputPositions, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);	
		}
   }
   else {
     return;
 	} // else-if Sum
}

unsigned int minimumNumberOfMismatches( const std::vector< PosVector>& mismatches, const unsigned int maxNumberMismatches )
{
	unsigned int numberOfMismatches=0; 
	for (;numberOfMismatches<=maxNumberMismatches; numberOfMismatches++ ) {
		if ( mismatches[ numberOfMismatches ].size() != 0 ) { break; }
	}
	return numberOfMismatches;
}


void CheckLeft_Close (const SPLIT_READ & read,
                 const std::string & chromosomeSeq,
                 const std::string & readSeq,
                 const std::vector< PosVector >& Left_PD,
                 const short &BP_Left_Start,
                 const short &BP_Left_End,
                 const short &CurrentLength, SortedUniquePoints &LeftUP)
{
	UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
		if (minimumNumberOfMismatches( Left_PD,read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
			return; 
		}
      // put it to LeftUP if unique
      for (short i = 0; i <= read.getMAX_SNP_ERROR(); i++) {
         if (Left_PD[i].size() == 1 && CurrentLength >= BP_Left_Start + i) {
            unsigned int Sum = numberOfCompetingPositions( Left_PD, i + userSettings->ADDITIONAL_MISMATCH );
				/*std::cout << "In CLC: CurrentLength = " << CurrentLength << ", mismatch count = " << i << ", maxMismatch = " << g_maxMismatch[CurrentLength] << std::endl;
				for (short k=0;k<=read.getMAX_SNP_ERROR(); k++) {
					std::cout << k << "\t" << Left_PD[k].size() << "\n";
				}*/
 
            if (Sum == 1 && (unsigned)i <= g_maxMismatch[CurrentLength] ) {
               UniquePoint TempOne(g_genome.getChr(read.FragName), CurrentLength, Left_PD[i][0], FORWARD, ANTISENSE, i );  
					std::cout << "Saving point\n";            
               if (CheckMismatches(chromosomeSeq, read.getUnmatchedSeq(), TempOne)) {
                  LeftUP.push_back (TempOne);
                  break;
               }
            }
         }
      }
   }
   if (CurrentLength < BP_Left_End) {
		ExtendMatchClose( read, chromosomeSeq, readSeq, Left_PD, BP_Left_Start, BP_Left_End, CurrentLength, LeftUP, 1 );
	}
}


void CheckRight_Close (const SPLIT_READ & read,
                  const std::string & chromosomeSeq,
                  const std::string & readSeq,
                  const std::vector < PosVector >& Right_PD,
                  const short &BP_Right_Start,
                  const short &BP_Right_End,
                  const short &CurrentLength, SortedUniquePoints &RightUP)
{
	UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End) {
		if (minimumNumberOfMismatches( Right_PD,read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
			return; 
		}
      for (short i = 0; i <= read.getMAX_SNP_ERROR(); i++) {
         if (Right_PD[i].size () == 1 && CurrentLength >= BP_Right_Start + i) {
            unsigned int Sum = numberOfCompetingPositions( Right_PD, i+userSettings->ADDITIONAL_MISMATCH );
				/*std::cout << "In CRC: CurrentLength = " << CurrentLength << ", mismatch count = " << i << ", maxMismatch = " << g_maxMismatch[CurrentLength] << std::endl;
				for (short k=0;k<=read.getMAX_SNP_ERROR(); k++) {
					std::cout << k << "\t" << Right_PD[k].size() << "\n";
				}*/
            if (Sum == 1 && (unsigned)i <= g_maxMismatch[CurrentLength] ) {
               UniquePoint TempOne( g_genome.getChr(read.FragName), CurrentLength, Right_PD[i][0], BACKWARD, SENSE, i);
               if (CheckMismatches(chromosomeSeq, read.getUnmatchedSeq(), TempOne)) {
                  RightUP.push_back (TempOne);
                  break;
               } // ###################################
            }
         }
      }
   }

   if (CurrentLength < BP_Right_End) {
		ExtendMatchClose( read, chromosomeSeq, readSeq, Right_PD, BP_Right_Start, BP_Right_End, CurrentLength, RightUP, -1 );
	}
}


bool CheckMismatches (const std::string & TheInput, const std::string & InputReadSeq, const UniquePoint & UP)
{
	int Min_Perfect_Match_Around_BP = UserDefinedSettings::Instance()->Min_Perfect_Match_Around_BP;
   std::string CurrentReadSeq;
   if (UP.Strand == SENSE) {
      CurrentReadSeq = InputReadSeq;
   }
   else {
      CurrentReadSeq = ReverseComplement (InputReadSeq);
   }
   short CurrentReadLength = CurrentReadSeq.size ();
   unsigned int Start = 0;
   std::string BP_On_Read, BP_On_Ref;
   if (UP.Direction == FORWARD) {

      Start = UP.AbsLoc - UP.LengthStr + 1;
      if (UP.LengthStr <= Min_Perfect_Match_Around_BP) {
         return false;
      }
      BP_On_Read = CurrentReadSeq.substr (UP.LengthStr - Min_Perfect_Match_Around_BP, Min_Perfect_Match_Around_BP);
      BP_On_Ref = TheInput.substr (UP.AbsLoc - Min_Perfect_Match_Around_BP + 1, Min_Perfect_Match_Around_BP);

      if (BP_On_Read != BP_On_Ref) {
         return false; //#################################
      }
   }
   else if (UP.Direction == BACKWARD) {
      Start = UP.AbsLoc + UP.LengthStr - CurrentReadLength;
      if (CurrentReadLength < UP.LengthStr) {
         return false;
      }
      BP_On_Read = CurrentReadSeq.substr (CurrentReadLength - UP.LengthStr, Min_Perfect_Match_Around_BP);
      BP_On_Ref = TheInput.substr (UP.AbsLoc, Min_Perfect_Match_Around_BP);
      if (BP_On_Read != BP_On_Ref) {
         return false; //#################################
      }
   }
   float MAX_ALLOWED_MISMATCHES = CurrentReadSeq.size () * UserDefinedSettings::Instance()->MaximumAllowedMismatchRate;	//

   short NumMismatches = 0;			// Match2N[(short)'A'] = 'N';

   for (short i = 0; i < CurrentReadLength; i++) {
      if (CurrentReadSeq[i] == N_char) {
         if (Match2N[(short) TheInput[Start + i]] != N_char) {
            NumMismatches++;
         }
      }
      else {
         if (TheInput[Start + i] != CurrentReadSeq[i]) {
            NumMismatches++;
         }
      }
   }
   if ((float)NumMismatches > MAX_ALLOWED_MISMATCHES) {
      return true;
   }
   else {
      return false;
   }
}
