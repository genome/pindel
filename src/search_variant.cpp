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

// System header files
#include <string>
#include <vector>

// Pindel header files
#include "control_state.h"
#include "search_variant.h"
#include "pindel.h"
#include "logdef.h"

SearchVariant::SearchVariant()
{
   Count_Var = 0;
   Count_Var_Plus = 0;
   Count_Var_Minus = 0;
   typeOfVariant = "some type of variant, replace this with the correct name in child class";
}

SearchVariant::~SearchVariant()
{

}

int SearchVariant::Search(ControlState& currentState, const unsigned numBoxes)
{

   std::vector<unsigned> Vars[numBoxes];

// EW150811 DEBUG->
   unsigned int farEndExists = 0;
   unsigned int readsUsed = 0;
   unsigned int bpSum = 0;
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++)  {
      if (currentState.Reads_SR[ReadIndex].Used) {
         readsUsed++;
      }
      if (!currentState.Reads_SR[ReadIndex].UP_Far.empty()) {
         farEndExists++;
         bpSum += currentState.Reads_SR[ReadIndex].UP_Far[ currentState.Reads_SR[ReadIndex].UP_Far.size()-1 ].AbsLoc;
      }
      if (bpSum>1000000000) {
         bpSum -= 1000000000;
      }
   }
   //*logStream << "At start of search deletions:\n\n";
   *logStream << "Reads already used: " << readsUsed << std::endl;
   *logStream << "Far ends already mapped " << farEndExists << std::endl;
   *logStream << "Checksum of far ends: " << bpSum << std::endl;

// <-EW150811 DEBUG

   LOG_INFO(*logStream << "Searching " << typeOfVariant << " ... " << std::endl);
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
      if (currentState.Reads_SR[ReadIndex].Used
            || currentState.Reads_SR[ReadIndex].UP_Far.empty()) {
         continue;
      }
      if (currentState.Reads_SR[ReadIndex].MatchedD == Plus) {
         for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
               <= currentState.Reads_SR[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
            for (unsigned int CloseIndex = 0; CloseIndex
                  < currentState.Reads_SR[ReadIndex].UP_Close.size(); CloseIndex++) {
               if (currentState.Reads_SR[ReadIndex].Used) {
                  break;
               }
               if (currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex]. Mismatches
                     > MAX_SNP_ERROR_index) {
                  continue;
               }
               for (int FarIndex =
                        currentState.Reads_SR[ReadIndex].UP_Far.size() - 1; FarIndex
                     >= 0; FarIndex--) {
                  if (currentState.Reads_SR[ReadIndex].Used) {
                     break;
                  }
                  if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. Mismatches
                        > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. Mismatches
                        + currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex]. Mismatches
                        > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. Direction
                        == Minus) {
                     //							if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc
                     //									== currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                     //											+ 1
                     //									&& currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                     //											+ currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. LengthStr
                     //											< currentState.Reads_SR[ReadIndex].ReadLength)
                     if (decisionBranch1(currentState, ReadIndex, CloseIndex, FarIndex)) {

                        currentState.Reads_SR[ReadIndex].Left
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                             - currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                             + 1;
                        currentState.Reads_SR[ReadIndex].Right
                           = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc
                             + currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. LengthStr
                             - 1;
                        currentState.Reads_SR[ReadIndex].BP
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                             - 1;

                        //								currentState.Reads_SR[ReadIndex].IndelSize
                        //										= currentState.Reads_SR[ReadIndex].ReadLengthMinus
                        //												- (currentState.Reads_SR[ReadIndex].Right
                        //														- currentState.Reads_SR[ReadIndex].Left);
                        currentState.Reads_SR[ReadIndex].IndelSize
                           = calculateIndelSize(currentState,
                                                ReadIndex);

                        //								currentState.Reads_SR[ReadIndex].InsertedStr
                        //										= ReverseComplement(
                        //												currentState.Reads_SR[ReadIndex]. UnmatchedSeq). substr(
                        //												currentState.Reads_SR[ReadIndex].BP
                        //														+ 1,
                        //												currentState.Reads_SR[ReadIndex]. IndelSize);
                        currentState.Reads_SR[ReadIndex].NT_str
                           = getInsertedStr1(currentState, ReadIndex);

                        currentState.Reads_SR[ReadIndex].BPLeft
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                             - g_SpacerBeforeAfter;
                        currentState.Reads_SR[ReadIndex].BPRight
                           = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc
                             - g_SpacerBeforeAfter;
                        if (readTransgressesBinBoundaries(
                                 currentState.Reads_SR[ReadIndex],
                                 currentState.upperBinBorder)) {
                           saveReadForNextCycle(
                              currentState.Reads_SR[ReadIndex],
                              currentState.FutureReads_SR);
                        }
                        else {
                           if (readInSpecifiedRegion(
                                    currentState.Reads_SR[ReadIndex],
                                    currentState.regionStartDefined,
                                    currentState.regionEndDefined,
                                    currentState.startOfRegion,
                                    currentState.endOfRegion)) {
                              Vars[(int) currentState.Reads_SR[ReadIndex]. BPLeft
                                   / BoxSize]. push_back(ReadIndex);
                              currentState.Reads_SR[ReadIndex].Used
                                 = true;
                              Count_Var_Plus++;
                              Count_Var++;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      else if (currentState.Reads_SR[ReadIndex].MatchedD == Minus) {
         for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
               <= currentState.Reads_SR[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
            for (int CloseIndex =
                     currentState.Reads_SR[ReadIndex].UP_Close.size() - 1; CloseIndex
                  >= 0; CloseIndex--) {
               if (currentState.Reads_SR[ReadIndex].Used) {
                  break;
               }
               if (currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex]. Mismatches
                     > MAX_SNP_ERROR_index) {
                  continue;
               }
               for (int FarIndex =
                        currentState.Reads_SR[ReadIndex].UP_Far.size() - 1; FarIndex
                     >= 0; FarIndex--) {
                  if (currentState.Reads_SR[ReadIndex].Used) {
                     break;
                  }
                  if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. Mismatches
                        > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. Mismatches
                        + currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex]. Mismatches
                        > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. Direction
                        == Plus) {
                     //							if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                     //									== currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc
                     //											+ 1
                     //									&& currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].LengthStr
                     //											+ currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                     //											< currentState.Reads_SR[ReadIndex].ReadLength) {
                     if (decisionBranch2(currentState, ReadIndex, CloseIndex, FarIndex)) {

                        currentState.Reads_SR[ReadIndex].Left
                           = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc
                             - currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. LengthStr
                             + 1;
                        currentState.Reads_SR[ReadIndex].Right
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                             + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                             - 1;
                        currentState.Reads_SR[ReadIndex].BP
                           = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. LengthStr
                             - 1;

                        //								currentState.Reads_SR[ReadIndex].IndelSize
                        //										= currentState.Reads_SR[ReadIndex].ReadLengthMinus
                        //												- (currentState.Reads_SR[ReadIndex].Right
                        //														- currentState.Reads_SR[ReadIndex].Left);
                        currentState.Reads_SR[ReadIndex].IndelSize
                           = calculateIndelSize(currentState,
                                                ReadIndex);

                        //								currentState.Reads_SR[ReadIndex].InsertedStr
                        //										= currentState.Reads_SR[ReadIndex].UnmatchedSeq. substr(
                        //												currentState.Reads_SR[ReadIndex].BP
                        //														+ 1,
                        //												currentState.Reads_SR[ReadIndex]. IndelSize);
                        currentState.Reads_SR[ReadIndex].NT_str
                           = getInsertedStr2(currentState,
                                             ReadIndex);

                        currentState.Reads_SR[ReadIndex].BPLeft
                           = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
                             - g_SpacerBeforeAfter;
                        currentState.Reads_SR[ReadIndex].BPRight
                           = currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc
                             - g_SpacerBeforeAfter;

                        if (readTransgressesBinBoundaries(
                                 currentState.Reads_SR[ReadIndex],
                                 currentState.upperBinBorder)) {
                           saveReadForNextCycle(
                              currentState.Reads_SR[ReadIndex],
                              currentState.FutureReads_SR);
                        }
                        else {
                           if (readInSpecifiedRegion(
                                    currentState.Reads_SR[ReadIndex],
                                    currentState.regionStartDefined,
                                    currentState.regionEndDefined,
                                    currentState.startOfRegion,
                                    currentState.endOfRegion)) {
                              Vars[(int) currentState.Reads_SR[ReadIndex]. BPLeft
                                   / BoxSize]. push_back(ReadIndex);
                              currentState.Reads_SR[ReadIndex].Used
                                 = true;
                              Count_Var++;
                              Count_Var_Minus++;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   LOG_INFO(*logStream << "Total: " << Count_Var << "\t+" << Count_Var_Plus << "\t-"
            << Count_Var_Minus << std::endl);

//	std::ofstream SIoutputfile(currentState.SIOutputFilename.c_str(),
//			std::ios::app);
//	SortOutputSI(NumBoxes, currentState.CurrentChr, currentState.Reads_SR, Vars,
//			SIoutputfile);
//	SIoutputfile.close();
   outputResults(currentState, Vars, numBoxes);

   for (unsigned int i = 0; i < numBoxes; i++) {
      Vars[i].clear();
   }

   return EXIT_SUCCESS;
}
