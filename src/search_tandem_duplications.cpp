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
#include "reporter.h"
#include "control_state.h"
#include "logdef.h"

int searchTandemDuplications(ControlState& currentState, unsigned NumBoxes)
{

   static int Count_TD = 0;
   static int Count_TD_Plus = 0;
   static int Count_TD_Minus = 0;

   std::vector<unsigned> TD[NumBoxes];

   LOG_INFO(*logStream << "Searching tandem duplication events ... " << std::endl);
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
      if (currentState.Reads_SR[ReadIndex].Used
            || currentState.Reads_SR[ReadIndex].UP_Far.empty()) {
         continue;
      }
      if (currentState.Reads_SR[ReadIndex].MatchedD == Plus) {
         for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
               <= currentState.Reads_SR[ReadIndex].getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
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

                     if (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].LengthStr
                           + currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].LengthStr
                           == currentState.Reads_SR[ReadIndex].getReadLength()
                           && currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
                           + currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].LengthStr
                           < currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc
                           && currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
                           + currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].LengthStr
                           < currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc) {

                        currentState.Reads_SR[ReadIndex].Right
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                             - currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                             + 1;
                        currentState.Reads_SR[ReadIndex].Left
                           = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                             + currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].LengthStr
                             - 1;
                        currentState.Reads_SR[ReadIndex].BP
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                             - 1;

                        currentState.Reads_SR[ReadIndex].IndelSize
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                             - currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                             + 1;
                        //currentState.Reads_SR[ReadIndex].InsertedStr = "";
                        currentState.Reads_SR[ReadIndex].BPRight
                           = currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc
                             - g_SpacerBeforeAfter;
                        currentState.Reads_SR[ReadIndex].BPLeft
                           = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
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
                              TD[(int) currentState.Reads_SR[ReadIndex]. BPLeft
                                 / BoxSize]. push_back(ReadIndex);
                              currentState.Reads_SR[ReadIndex].Used
                                 = true;
                              Count_TD++;
                              Count_TD_Plus++;
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
               <= currentState.Reads_SR[ReadIndex].getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
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
               for (int FarIndex = 0; FarIndex
                     < (int) currentState.Reads_SR[ReadIndex].UP_Far.size(); FarIndex++) {
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
                     if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                           + currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. LengthStr
                           == currentState.Reads_SR[ReadIndex].getReadLength()
                           && currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                           + currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].LengthStr
                           < currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc
                           && currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                           + currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].LengthStr
                           < currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc) {
                        currentState.Reads_SR[ReadIndex].Right
                           = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                             - currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].LengthStr
                             + 1;
                        currentState.Reads_SR[ReadIndex].Left
                           = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                             + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
                             - 1;
                        currentState.Reads_SR[ReadIndex].BP
                           = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].LengthStr
                             - 1;

                        currentState.Reads_SR[ReadIndex].IndelSize
                           = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                             - currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
                             + 1;
                        //currentState.Reads_SR[ReadIndex].InsertedStr = "";
                        currentState.Reads_SR[ReadIndex].BPRight
                           = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
                             - g_SpacerBeforeAfter;
                        currentState.Reads_SR[ReadIndex].BPLeft
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
                              TD[(int) currentState.Reads_SR[ReadIndex]. BPLeft
                                 / BoxSize]. push_back(ReadIndex);
                              currentState.Reads_SR[ReadIndex].Used
                                 = true;

                              Count_TD++;
                              Count_TD_Minus++;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   LOG_INFO(*logStream << "Total: " << Count_TD << "\t+" << Count_TD_Plus << "\t-"
            << Count_TD_Minus << std::endl);
   std::ofstream TDOutf(currentState.TDOutputFilename.c_str(), std::ios::app);
   SortAndOutputTandemDuplications(NumBoxes, currentState.CurrentChrSeq,
                                   currentState.Reads_SR, TD, TDOutf, false);
   for (unsigned int i = 0; i < NumBoxes; i++) {
      TD[i].clear();
   }

   return EXIT_SUCCESS;
}
