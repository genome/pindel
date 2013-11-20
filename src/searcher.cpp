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
#include <assert.h>
#include "pindel.h"
#include "searcher.h"
#include "sse_helpers.h"
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

// The following are different flavors of CategorizePositions
inline void CategorizePositionsNormal(int readBase,
        const std::string& chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        int searchDirection) __attribute__((always_inline));

inline void CategorizePositionsNormal(int readBase,
        const std::string & chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        int searchDirection) {
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) {
        unsigned int pos = *it + searchDirection;
        Output[(int) MismatchPair[readBase][(int) chromosomeSeq[pos]]].push_back(pos);
    }
}

void CategorizePositionsBoundary(int readBase,
        const std::string & chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        int searchDirection)
{
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) {
        unsigned int pos = *it + searchDirection;
        if (!MismatchPair[readBase][(int) chromosomeSeq[pos]]) {
            Output[0].push_back(pos);
        }
    }
}

// This one extends two bases at a time, assuming no read bases are N
inline void CategorizePositions2NormalNoN(int readBase1,
        int readBase2,
        const std::string& chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        int searchDirection)
{
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) {
        unsigned int pos1 = *it + searchDirection;
        unsigned int pos2 = pos1 + searchDirection;
        int index = (chromosomeSeq[pos1] == readBase1)? 0: 1;
        index += (chromosomeSeq[pos2] == readBase2)? 0: 1;
        Output[index].push_back(pos2);
    }
}

inline void CategorizePositions2Boundary1NoN(int readBase1,
        int readBase2,
        const std::string& chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        int searchDirection)
{
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) { 
        unsigned int pos1 = *it + searchDirection;
        unsigned int pos2 = pos1 + searchDirection;
        int index = (chromosomeSeq[pos1] == readBase1)? 0: 1;
        index += (chromosomeSeq[pos2] == readBase2)? 0: 1;
        if (index < 2) Output[index].push_back(pos2);
    }
}

inline void CategorizePositions2Boundary2NoN(int readBase1,
        int readBase2,
        const std::string & chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        const int searchDirection)
{
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) {
        unsigned int pos1 = *it + searchDirection;
        unsigned int pos2 = pos1 + searchDirection;
        int index = (chromosomeSeq[pos1] == readBase1)? 0: 1;
        index += (chromosomeSeq[pos2] == readBase2)? 0: 1;
        if (index == 0) Output[index].push_back(pos2);
    }
}

inline void CategorizePositions2Normal(int readBase1,
        int readBase2,
        const std::string& chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        const int searchDirection)
{
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) {
        unsigned int pos1 = *it + searchDirection;
        unsigned int pos2 = pos1 + searchDirection;
        int index = MismatchPair[readBase1][(int) chromosomeSeq[pos1]] + MismatchPair[readBase2][(int) chromosomeSeq[pos2]];
        Output[index].push_back(pos2);
    }
}

inline void CategorizePositions2Boundary1(int readBase1,
        int readBase2,
        const std::string& chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        const int searchDirection)
{
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) {
        unsigned int pos1 = *it + searchDirection;
        unsigned int pos2 = pos1 + searchDirection;
        int index = MismatchPair[readBase1][(int) chromosomeSeq[pos1]] + MismatchPair[readBase2][(int) chromosomeSeq[pos2]];
        if (index < 2) Output[index].push_back(pos2);
    }
}

inline void CategorizePositions2Boundary2(int readBase1,
        int readBase2,
        const std::string& chromosomeSeq,
        const PosVector& Input,
        PosVector* Output,
        const int searchDirection)
{
    for (PosVector::const_iterator it = Input.begin(); it != Input.end(); it++) {
        unsigned int pos1 = *it + searchDirection;
        unsigned int pos2 = pos1 + searchDirection;
        int index = MismatchPair[readBase1][(int) chromosomeSeq[pos1]] + MismatchPair[readBase2][(int) chromosomeSeq[pos2]];
        if (index == 0) Output[index].push_back(pos2);
    }
}

// The original CategorizePositions
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
            if ( !MismatchPair[readBase][(int) chromosomeSeq[pos]] ) {
                PD_Plus_Output_current->push_back(pos);
            }
        }
    }
}

// Extents one base from the positions in the positions array. The idea is to use as less memory as possible and modify the positions array in place.
void ExtendInPlace(int currentChar,
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

// Extends two bases from the positions in the positions array
void ExtendInPlace2(int char1,
        int char2,
        const std::string& chromosomeSeq,
        std::vector<PosVector>& positions,
        PosVector* TmpPositions,
        int direction,
        int minMismatches,
        int maxMismatches) {
    if (minMismatches > maxMismatches - 1) {
        return;
    }
    TmpPositions[0].clear();
    TmpPositions[1].clear();
    TmpPositions[2].clear();
    TmpPositions[0].swap(positions[minMismatches]);
    TmpPositions[1].swap(positions[minMismatches+1]);

    if (char1 == 'N' || char2 == 'N') {
        int i;
        for (i = minMismatches; i < maxMismatches - 1; i++) {
            TmpPositions[2].swap(positions[i+2]); 
            CategorizePositions2Normal(char1, char2, chromosomeSeq, TmpPositions[0], &positions[i], direction);
            TmpPositions[0].swap(TmpPositions[1]);
            TmpPositions[1].swap(TmpPositions[2]);
            TmpPositions[2].clear();
        }
        CategorizePositions2Boundary1(char1, char2, chromosomeSeq, TmpPositions[0], &positions[i], direction);
        TmpPositions[0].swap(TmpPositions[1]);
        i++;

        CategorizePositions2Boundary2(char1, char2, chromosomeSeq, TmpPositions[0], &positions[i], direction);
    } else {
        int i;
        for (i = minMismatches; i < maxMismatches - 1; i++) {
            TmpPositions[2].swap(positions[i+2]); 
            CategorizePositions2NormalNoN(char1, char2, chromosomeSeq, TmpPositions[0], &positions[i], direction);
            TmpPositions[0].swap(TmpPositions[1]);
            TmpPositions[1].swap(TmpPositions[2]);
            TmpPositions[2].clear();
        }
        CategorizePositions2Boundary1NoN(char1, char2, chromosomeSeq, TmpPositions[0], &positions[i], direction);
        TmpPositions[0].swap(TmpPositions[1]);
        i++;

        CategorizePositions2Boundary2NoN(char1, char2, chromosomeSeq, TmpPositions[0], &positions[i], direction);
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
{
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

void CheckLeft_Close (SPLIT_READ & read,
        const std::string & chromosomeSeq,
        const std::string & readSeq,
        std::vector< PosVector >& Left_PD,
        int BP_Left_Start,
        int BP_Left_End,
        int CurrentLength, SortedUniquePoints &LeftUP)
{
    int minMismatches = 0;
    int maxMismatches = read.getMAX_SNP_ERROR();
    PosVector TmpPositions[3];
    while (Left_PD[minMismatches].size() == 0 && minMismatches <= maxMismatches) {
        minMismatches++;
    }

    for ( ; CurrentLength < BP_Left_Start && minMismatches <= maxMismatches; CurrentLength++) {
        ExtendInPlace(readSeq[CurrentLength], chromosomeSeq, Left_PD, TmpPositions, 1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());

        minMismatches += (Left_PD[minMismatches].size() == 0)? 1: 0;
    }

    for ( ; CurrentLength <= BP_Left_End && minMismatches <= maxMismatches; CurrentLength++) {
        if (minMismatches > g_maxMismatch[CurrentLength] ) {
            return; 
        }

        if (minMismatches <= CurrentLength - BP_Left_Start && Left_PD[minMismatches].size() == 1) {
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

        if (CurrentLength < BP_Left_End) {
            ExtendInPlace(readSeq[CurrentLength], chromosomeSeq, Left_PD, TmpPositions, 1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
        }
        if (Left_PD[minMismatches].size() == 0) {
            minMismatches++;
        }
    }
}

void CheckLeft_Close_Perfect (SPLIT_READ & read,
        const std::string & chromosomeSeq,
        const std::string & readSeq,
        std::vector< PosVector >& Left_PD,
        int BP_Left_Start,
        int BP_Left_End,
        int CurrentLength, SortedUniquePoints &LeftUP)
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
        int BP_Right_Start,
        int BP_Right_End,
        int CurrentLength, SortedUniquePoints &RightUP)
{
    PosVector TmpPositions[2];
    int minMismatches = 0;
    while (Right_PD[minMismatches].size() == 0 && minMismatches <= read.getMAX_SNP_ERROR()) {
        minMismatches++;
    }

    for ( ; CurrentLength < BP_Right_Start && minMismatches <= read.getMAX_SNP_ERROR(); CurrentLength++) {
        ExtendInPlace(readSeq[read.getReadLengthMinus() - CurrentLength] , chromosomeSeq, Right_PD, TmpPositions, -1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());

        minMismatches += (Right_PD[minMismatches].size() == 0)? 1: 0;	
    }

    for ( ; CurrentLength <= BP_Right_End && minMismatches <= read.getMAX_SNP_ERROR(); CurrentLength++) {
        if (minMismatches > g_maxMismatch[CurrentLength] ) {
            return; 
        }

        if (minMismatches <= CurrentLength - BP_Right_Start && Right_PD[minMismatches].size() == 1) {
            const std::string& forwardSeq = read.getUnmatchedSeq();
            const std::string& reverseSeq = read.getUnmatchedSeqRev();
            unsigned int Sum = numberOfCompetingPositions(Right_PD, minMismatches+1, minMismatches+userSettings->ADDITIONAL_MISMATCH );
            if (Sum == 0) {
                UniquePoint TempOne( g_genome.getChr(read.FragId), CurrentLength, Right_PD[minMismatches][0], BACKWARD, SENSE, minMismatches);
                if (CheckMismatches(chromosomeSeq, forwardSeq, reverseSeq, TempOne, read.CloseEndMismatch)) {
                    RightUP.push_back (TempOne);
                }
            }
        }

        if (CurrentLength < BP_Right_End) {
            ExtendInPlace(readSeq[read.getReadLengthMinus() - CurrentLength] , chromosomeSeq, Right_PD, TmpPositions, -1, minMismatches, read.getTOTAL_SNP_ERROR_CHECKED_Minus());
        }
        if (Right_PD[minMismatches].size() == 0) {
            minMismatches++;
        }
    }
}

void CheckRight_Close_Perfect (SPLIT_READ & read,
        const std::string & chromosomeSeq,
        const std::string & readSeq,
        std::vector < PosVector >& Right_PD,
        int BP_Right_Start,
        int BP_Right_End,
        int CurrentLength, SortedUniquePoints &RightUP)
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
    const std::string* CurrentReadSeq = (UP.Strand == SENSE)? &InputReadSeq: &InputReadSeqRev;
    int CurrentReadLength = CurrentReadSeq->size ();
    unsigned int Start = 0;
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
        const char* BP_On_Read_c = &(*CurrentReadSeq)[CurrentReadLength - UP.LengthStr];
        const char* BP_On_Ref_c = &TheInput[UP.AbsLoc];
        if (strncmp(BP_On_Read_c, BP_On_Ref_c, Min_Perfect_Match_Around_BP) != 0) {
            return false;
        }
    }
    float MAX_ALLOWED_MISMATCHES = CurrentReadSeq->size() * userSettings->MaximumAllowedMismatchRate;	//
    int NumMismatches = CountMismatches(&(*CurrentReadSeq)[0], &TheInput[Start], CurrentReadLength);
    numberOfMismatch = NumMismatches;
    return NumMismatches >= MAX_ALLOWED_MISMATCHES;
}
