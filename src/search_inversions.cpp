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

int searchInversions(ControlState& currentState, unsigned NumBoxes)
{
    //return 0; 
   static int Count_Inv = 0;
   static int Count_Inv_Plus = 0;
   static int Count_Inv_Minus = 0;

   std::vector<unsigned> Inv[NumBoxes];

   int CloseIndex = 0;
   int FarIndex = 0;
    //*logStream << "Debug Searching inversions ... " << std::endl;
   LOG_INFO(*logStream << "Searching inversions ... " << std::endl);
    //*logStream << "start" << std::endl;
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
      if (currentState.Reads_SR[ReadIndex].Used
            || currentState.Reads_SR[ReadIndex].UP_Far.empty()) {
         continue;
      }
      if (currentState.Reads_SR[ReadIndex].UP_Close[0].Strand
            != currentState.Reads_SR[ReadIndex].UP_Far[0].Strand
            && currentState.Reads_SR[ReadIndex].UP_Close[0].Direction
            == currentState.Reads_SR[ReadIndex].UP_Far[0].Direction) {

         if (currentState.Reads_SR[ReadIndex].MatchedD == Plus) {
            if (currentState.Reads_SR[ReadIndex].UP_Far[0].AbsLoc
                  > currentState.Reads_SR[ReadIndex].UP_Close[currentState.Reads_SR[ReadIndex]. UP_Close.size()
                        - 1].AbsLoc + MIN_IndelSize_Inversion) { // normal situation
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
                     <= currentState.Reads_SR[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
                  for (CloseIndex
                        = (int) currentState.Reads_SR[ReadIndex].UP_Close.size()
                          - 1; CloseIndex >= 0; CloseIndex--) {
                     if (currentState.Reads_SR[ReadIndex].Used) {
                        break;
                     }
                     if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].Mismatches
                           > MAX_SNP_ERROR_index) {
                        continue;
                     }
                     for (FarIndex = 0; FarIndex
                           < (int) currentState.Reads_SR[ReadIndex].UP_Far.size(); FarIndex++) {
                        if (currentState.Reads_SR[ReadIndex].Used) {
                           break;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Direction
                              == Plus) {
                           if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].LengthStr
                                 + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                 == currentState.Reads_SR[ReadIndex].getReadLength()
                                 && currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                 > currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                 + MIN_IndelSize_Inversion) {
                              currentState.Reads_SR[ReadIndex].Left
                                 = (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                    + 1)
                                   - currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr;
                              currentState.Reads_SR[ReadIndex].Right
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                   - currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                   + currentState.Reads_SR[ReadIndex].getReadLength();
                              currentState.Reads_SR[ReadIndex].BP
                                 = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                   - 1;

                              currentState.Reads_SR[ReadIndex].IndelSize
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                   - currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc;
                              currentState.Reads_SR[ReadIndex].NT_str
                                 = "";
                              currentState.Reads_SR[ReadIndex].NT_size
                                 = 0;
                              //currentState.Reads_SR[ReadIndex]. InsertedStr
                              //   = "";
                              currentState.Reads_SR[ReadIndex].BPLeft
                                 = currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc
                                   + 1 - g_SpacerBeforeAfter;
                              currentState.Reads_SR[ReadIndex].BPRight
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
                                    Inv[(int) (currentState.Reads_SR[ReadIndex]. BPLeft + currentState.Reads_SR[ReadIndex]. BPRight)
                                        / BoxSize]. push_back(
                                           ReadIndex);
                                    currentState.Reads_SR[ReadIndex]. Used
                                       = true;
                                    Count_Inv++;
                                    Count_Inv_Plus++;
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
            
            else if (currentState.Reads_SR[ReadIndex]. UP_Far[currentState.Reads_SR[ReadIndex].UP_Far.size()
                     - 1].AbsLoc + MIN_IndelSize_Inversion
                     < currentState.Reads_SR[ReadIndex].UP_Close[0].AbsLoc) {
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
                     <= currentState.Reads_SR[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
                  for (CloseIndex = 0; CloseIndex
                        < (int) currentState.Reads_SR[ReadIndex].UP_Close.size(); CloseIndex++) {
                     if (currentState.Reads_SR[ReadIndex].Used) {
                        break;
                     }
                     if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].Mismatches
                           > MAX_SNP_ERROR_index) {
                        continue;
                     }
                     for (FarIndex
                           = (int) currentState.Reads_SR[ReadIndex].UP_Far.size()
                             - 1; FarIndex >= 0; FarIndex--) {
                        if (currentState.Reads_SR[ReadIndex].Used) {
                           break;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Direction
                              == Plus) {
                            
                           // anchor inside reversed block.
                           if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].LengthStr
                                 + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                 == currentState.Reads_SR[ReadIndex].getReadLength()
                                 && currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                 + MIN_IndelSize_Inversion
                                 < currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc) {
                              currentState.Reads_SR[ReadIndex].Right
                                 = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                   - currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                   + currentState.Reads_SR[ReadIndex].getReadLength();
                              currentState.Reads_SR[ReadIndex].Left
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                   - currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                   + 1;
                              currentState.Reads_SR[ReadIndex].BP
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                   - 1;

                              currentState.Reads_SR[ReadIndex].IndelSize
                                 = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                   - currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc;
                              currentState.Reads_SR[ReadIndex].NT_str
                                 = "";
                              currentState.Reads_SR[ReadIndex].NT_size
                                 = 0;
                              //currentState.Reads_SR[ReadIndex]. InsertedStr
                              //   = "";
                              currentState.Reads_SR[ReadIndex].BPRight
                                 = currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc
                                   - g_SpacerBeforeAfter;
                              currentState.Reads_SR[ReadIndex].BPLeft
                                 = (currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
                                    + 1)
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
                                    Inv[(int) (currentState.Reads_SR[ReadIndex]. BPLeft + currentState.Reads_SR[ReadIndex]. BPRight)
                                        / BoxSize]. push_back(
                                           ReadIndex);
                                    currentState.Reads_SR[ReadIndex]. Used
                                       = true;
                                    Count_Inv++;
                                    Count_Inv_Plus++;
                                 }
                              }
                           }
                             
                        }
                     }
                  }
               }
            }
              
         }
              
         else if (currentState.Reads_SR[ReadIndex].MatchedD == Minus) {
            if (currentState.Reads_SR[ReadIndex]. UP_Close[currentState.Reads_SR[ReadIndex].UP_Close.size()
                  - 1].AbsLoc
                  > currentState.Reads_SR[ReadIndex].UP_Far[0].AbsLoc
                  + MIN_IndelSize_Inversion) {
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
                     <= currentState.Reads_SR[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
                  for (CloseIndex
                        = (int) currentState.Reads_SR[ReadIndex].UP_Close.size()
                          - 1; CloseIndex >= 0; CloseIndex--) {
                     if (currentState.Reads_SR[ReadIndex].Used) {
                        break;
                     }
                     if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].Mismatches
                           > MAX_SNP_ERROR_index) {
                        continue;
                     }
                     for (FarIndex = 0; FarIndex
                           < (int) currentState.Reads_SR[ReadIndex].UP_Far.size(); FarIndex++) {
                        if (currentState.Reads_SR[ReadIndex].Used) {
                           break;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Direction
                              == Minus) {
                           if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                 + currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                 == currentState.Reads_SR[ReadIndex].getReadLength()
                                 && currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                 > currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                 + MIN_IndelSize_Inversion) {
                              currentState.Reads_SR[ReadIndex].Left
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                   + currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                   - currentState.Reads_SR[ReadIndex].getReadLength();
                              currentState.Reads_SR[ReadIndex].Right
                                 = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                   + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                   - 1;
                              currentState.Reads_SR[ReadIndex].BP
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                   - 1;

                              currentState.Reads_SR[ReadIndex].IndelSize
                                 = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                   - currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc;
                              currentState.Reads_SR[ReadIndex].NT_str
                                 = "";
                              currentState.Reads_SR[ReadIndex].NT_size
                                 = 0;
                              //currentState.Reads_SR[ReadIndex]. InsertedStr
                              //   = "";
                              currentState.Reads_SR[ReadIndex].BPLeft
                                 = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
                                   - g_SpacerBeforeAfter;
                              currentState.Reads_SR[ReadIndex].BPRight
                                 = currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc
                                   - 1 - g_SpacerBeforeAfter;
                              if (readInSpecifiedRegion(
                                       currentState.Reads_SR[ReadIndex],
                                       currentState.regionStartDefined,
                                       currentState.regionEndDefined,
                                       currentState.startOfRegion,
                                       currentState.endOfRegion)) {
                                 Inv[(int) (currentState.Reads_SR[ReadIndex]. BPLeft + currentState.Reads_SR[ReadIndex]. BPRight)
                                     / BoxSize]. push_back(
                                        ReadIndex);
                                 currentState.Reads_SR[ReadIndex]. Used
                                    = true;

                                 Count_Inv++;
                                 Count_Inv_Minus++;
                              }
                           }
                        }
                     }
                  }
               }
            }
             
            else if (currentState.Reads_SR[ReadIndex].UP_Close[0].AbsLoc
                     + MIN_IndelSize_Inversion
                     < currentState.Reads_SR[ReadIndex].UP_Far[currentState.Reads_SR[ReadIndex]. UP_Far.size()
                           - 1].AbsLoc) {
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
                     <= currentState.Reads_SR[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
                  for (unsigned int CloseIndex = 0; CloseIndex
                        < currentState.Reads_SR[ReadIndex].UP_Close.size(); CloseIndex++) {
                     if (currentState.Reads_SR[ReadIndex].Used) {
                        break;
                     }
                     if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].Mismatches
                           > MAX_SNP_ERROR_index) {
                        continue;
                     }
                     for (int FarIndex =
                              currentState.Reads_SR[ReadIndex].UP_Far.size()
                              - 1; FarIndex >= 0; FarIndex--) {
                        if (currentState.Reads_SR[ReadIndex].Used) {
                           break;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Mismatches
                              + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. Mismatches
                              > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].Direction
                              == Minus) {
                            
                           // anchor inside reversed block.
                           if (currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                 + currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                 == currentState.Reads_SR[ReadIndex].getReadLength()
                                 && currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                 + MIN_IndelSize_Inversion
                                 < currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc) {
                              currentState.Reads_SR[ReadIndex].Right
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                   + currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex]. LengthStr
                                   - 1;
                              currentState.Reads_SR[ReadIndex].Left
                                 = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
                                   + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                   - currentState.Reads_SR[ReadIndex].getReadLength();
                              currentState.Reads_SR[ReadIndex].BP
                                 = currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. LengthStr
                                   - 1;

                              currentState.Reads_SR[ReadIndex].IndelSize
                                 = currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].AbsLoc
                                   - currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex]. AbsLoc;
                              currentState.Reads_SR[ReadIndex].NT_str
                                 = "";
                              currentState.Reads_SR[ReadIndex].NT_size
                                 = 0;
                              //currentState.Reads_SR[ReadIndex]. InsertedStr
                              //   = "";
                              currentState.Reads_SR[ReadIndex].BPLeft
                                 = currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc
                                   - g_SpacerBeforeAfter;
                              currentState.Reads_SR[ReadIndex].BPRight
                                 = currentState.Reads_SR[ReadIndex].UP_Far[FarIndex].AbsLoc
                                   - 1 - g_SpacerBeforeAfter;
                              if (readInSpecifiedRegion(
                                       currentState.Reads_SR[ReadIndex],
                                       currentState.regionStartDefined,
                                       currentState.regionEndDefined,
                                       currentState.startOfRegion,
                                       currentState.endOfRegion)) {
                                 Inv[(int) (currentState.Reads_SR[ReadIndex]. BPLeft + currentState.Reads_SR[ReadIndex]. BPRight)
                                     / BoxSize]. push_back(
                                        ReadIndex);
                                 currentState.Reads_SR[ReadIndex]. Used
                                    = true;

                                 Count_Inv++;
                                 Count_Inv_Minus++;
                              }
                           }
                             
                        }
                     }
                  }
               }
            }
            
         }
      }
   }
    //*logStream << "end" << std::endl;
   LOG_INFO(*logStream << "Total: " << Count_Inv << "\t+" << Count_Inv_Plus << "\t-"
            << Count_Inv_Minus << std::endl);
   std::ofstream InversionOutf(currentState.InversionOutputFilename.c_str(),
                               std::ios::app);
   SortOutputInv(NumBoxes, currentState.CurrentChrSeq, currentState.Reads_SR, Inv,
                 InversionOutf);
   for (unsigned int i = 0; i < NumBoxes; i++) {
      Inv[i].clear();
   }

   return EXIT_SUCCESS;
}
