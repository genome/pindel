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
#include "logstream.h"
#include "reporter.h"
#include "control_state.h"
#include "logdef.h"
#include "search_tandem_duplications.h"

int searchTandemDuplications(ControlState& currentState, unsigned NumBoxes, const SearchWindow& window)
{

   static int Count_TD = 0;
   static int Count_TD_Plus = 0;
   static int Count_TD_Minus = 0;

   std::vector<unsigned> TD[NumBoxes];
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();


   LOG_INFO(*logStream << "Searching tandem duplication events ... " << std::endl);
    //std::cout << "std1" << std::endl;
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
		SPLIT_READ& currentRead = currentState.Reads_SR[ReadIndex];
       //std::cout << "std1a" << std::endl;
      if (currentRead.Used || currentRead.UP_Far.empty() || currentRead.FragName != currentRead.FarFragName) {
         continue;
      }
       //std::cout << "std1b" << std::endl;
      if (currentRead.MatchedD == Plus) {
          //std::cout << "in plus" << std::endl;
         for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
             //std::cout << "std1c" << std::endl;
            for (unsigned int CloseIndex = 0; CloseIndex < currentRead.UP_Close.size(); CloseIndex++) {
                //std::cout << "std1d" << std::endl;
               if (currentRead.Used /* || currentRead. BPLeft == 0*/) {
                  break;
               }
               if (currentRead.UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                  continue;
               }
                //std::cout << "std1e " << currentRead.UP_Far.size() << std::endl;
               for (int FarIndex = currentRead.UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
                  //std::cout << "std1f" << std::endl;
                  if (currentRead.Used) {
                     break;
                  }
                  if (currentRead.UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentRead.UP_Far[FarIndex].Mismatches + currentRead.UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                   //std::cout << "std1g" << std::endl;
                  if (currentRead.UP_Far[FarIndex].Direction  == Minus) {

                     if (currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr == currentRead.getReadLength() && 
                         currentRead.UP_Far[FarIndex].AbsLoc + currentRead.UP_Far[FarIndex].LengthStr < currentRead.UP_Close[CloseIndex].AbsLoc && 
								 currentRead.UP_Far[FarIndex].AbsLoc + currentRead.UP_Close[CloseIndex].LengthStr < currentRead.UP_Close[CloseIndex].AbsLoc) {

                        currentRead.Right = currentRead.UP_Close[CloseIndex].AbsLoc - currentRead.UP_Close[CloseIndex].LengthStr + 1;
                        currentRead.Left = currentRead. UP_Far[FarIndex].AbsLoc + currentRead. UP_Far[FarIndex].LengthStr - 1;
                        currentRead.BP = currentRead.UP_Close[CloseIndex].LengthStr - 1;

                        currentRead.IndelSize = currentRead.UP_Close[CloseIndex].AbsLoc - currentRead. UP_Far[FarIndex].AbsLoc + 1;
                        currentRead.BPRight = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                        currentRead.BPLeft = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
                         if (currentRead. BPLeft == 0) continue;
                         LeftMostTD(currentState, currentRead, window);
                        if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                           saveReadForNextCycle( currentRead, currentState.FutureReads_SR);
                        }
                        else {
                           if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                              TD[(int) currentRead. BPLeft / BoxSize]. push_back(ReadIndex);
                              currentRead.Used = true;
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
      else if (currentRead.MatchedD == Minus) {
         //std::cout << "in Minus" << std::endl;
         for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
            //std::cout << "in Minus 1" << std::endl;
            for (int CloseIndex = currentRead.UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
                //std::cout << "in Minus 2" << std::endl;
               if (currentRead.Used /*|| currentRead. BPLeft == 0*/) {
                  break;
               }
               if (currentRead.UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                  continue;
               }
                //std::cout << "in Minus 3" << std::endl;
               for (int FarIndex = 0; FarIndex < (int) currentRead.UP_Far.size(); FarIndex++) {
                  if (currentRead.Used) {
                     break;
                  }
                  if (currentRead.UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentRead.UP_Far[FarIndex].Mismatches + currentRead.UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentRead.UP_Far[FarIndex]. Direction == Plus) {
                      //std::cout << "BPLeft" << currentRead. BPLeft << std::endl;
                     if (currentRead.UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr == currentRead.getReadLength() && 
                         currentRead.UP_Close[CloseIndex].AbsLoc + currentRead.UP_Close[CloseIndex].LengthStr < currentRead.UP_Far[FarIndex]. AbsLoc && 
                         currentRead.UP_Close[CloseIndex].AbsLoc + currentRead.UP_Far[FarIndex].LengthStr < currentRead.UP_Far[FarIndex]. AbsLoc) {
                        currentRead.Right = currentRead. UP_Far[FarIndex].AbsLoc - currentRead. UP_Far[FarIndex].LengthStr + 1;
                        currentRead.Left = currentRead.UP_Close[CloseIndex].AbsLoc + currentRead.UP_Close[CloseIndex].LengthStr - 1;
                        currentRead.BP = currentRead. UP_Far[FarIndex].LengthStr - 1;

                        currentRead.IndelSize = currentRead. UP_Far[FarIndex].AbsLoc - currentRead.UP_Close[CloseIndex].AbsLoc + 1;
                        currentRead.BPRight = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
                        currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                         if (currentRead. BPLeft == 0) continue;
                         LeftMostTD(currentState, currentRead, window);
                        if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                           saveReadForNextCycle( currentRead, currentState.FutureReads_SR);
                        }
                        else {
                           if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                              TD[(int) currentRead. BPLeft / BoxSize]. push_back(ReadIndex);
                              currentRead.Used = true;

                              Count_TD++;
                              Count_TD_Minus++;
                           }
                        }
                     }
                  }
               }
                //std::cout << "in Minus 4" << std::endl;
            }
         }
      }
   }
    //std::cout << "std2" << std::endl;
   LOG_INFO(*logStream << "Total: " << Count_TD << "\t+" << Count_TD_Plus << "\t-"  << Count_TD_Minus << std::endl);
   std::ofstream TDOutf(userSettings->getTDOutputFilename().c_str(), std::ios::app);
   SortAndOutputTandemDuplications(currentState, NumBoxes, window.getChromosome()->getSeq(), currentState.Reads_SR, TD, TDOutf, false);
    //std::cout << "std3" << std::endl;
   for (unsigned int i = 0; i < NumBoxes; i++) {
      TD[i].clear();
   }
   // std::cout << "std4" << std::endl;
   return EXIT_SUCCESS;
}

void LeftMostTD(ControlState& currentState, SPLIT_READ & currentRead, const SearchWindow& window) {
    const std::string & TheInput = window.getChromosome()->getSeq();
    
    unsigned int PosIndex = currentRead.BPLeft + g_SpacerBeforeAfter;
    unsigned int original_PosIndex = PosIndex;
    //unsigned int Start = PosIndex + 1;
    unsigned int End = currentRead.BPRight + g_SpacerBeforeAfter - 1;
    while (TheInput[PosIndex] == TheInput[End]) {
        --PosIndex;
        --End;
    }
    short DIFF = original_PosIndex - PosIndex;
    if (DIFF > 0) {
        if (DIFF >= currentRead.BP) DIFF = currentRead.BP - 1;
        currentRead.BPLeft -= DIFF;
        currentRead.BPRight -= DIFF;
        currentRead.BP -= DIFF;
    }
}

