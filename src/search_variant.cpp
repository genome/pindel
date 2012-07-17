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
#include "logstream.h"
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

int SearchVariant::Search(ControlState& currentState, const unsigned numBoxes, const SearchWindow& window)
{
   std::vector<unsigned> Vars[numBoxes];

   unsigned int farEndExists = 0;
   unsigned int readsUsed = 0;
   unsigned int bpSum = 0;
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++)  {
		SPLIT_READ& currentRead = currentState.Reads_SR[ReadIndex];
      if (currentRead.Used) {
         readsUsed++;
      }
      if (!currentRead.UP_Far.empty()) {
         farEndExists++;
         bpSum += currentRead.UP_Far[ currentRead.UP_Far.size()-1 ].AbsLoc;
      }
      if (bpSum>1000000000) {
         bpSum -= 1000000000;
      }
   }
   *logStream << "Reads already used: " << readsUsed << std::endl;
   *logStream << "Far ends already mapped " << farEndExists << std::endl;
   *logStream << "Checksum of far ends: " << bpSum << std::endl;

	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
   LOG_INFO(*logStream << "Searching " << typeOfVariant << " ... " << std::endl);
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
		SPLIT_READ& currentRead = currentState.Reads_SR[ReadIndex];
      if (currentRead.Used || currentRead.UP_Far.empty()) {
         continue;
      }
      if (currentRead.MatchedD == Plus) {
         for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
            for (unsigned int CloseIndex = 0; CloseIndex
                  < currentRead.UP_Close.size(); CloseIndex++) {
               if (currentRead.Used) {
                  break;
               }
               if (currentRead.UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                  continue;
               }
               for (int FarIndex = currentRead.UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
                  if (currentRead.Used) {
                     break;
                  }
                  if (currentRead.UP_Far[FarIndex]. Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentRead.UP_Far[FarIndex]. Mismatches + currentRead.UP_Close[CloseIndex]. Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentRead.UP_Far[FarIndex]. Direction == Minus) {
                     if (decisionBranch1(currentRead, CloseIndex, FarIndex)) {

                        currentRead.Left = currentRead. UP_Close[CloseIndex].AbsLoc - currentRead. UP_Close[CloseIndex].LengthStr + 1;
                        currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc + currentRead.UP_Far[FarIndex]. LengthStr - 1;
                        currentRead.BP = currentRead. UP_Close[CloseIndex].LengthStr - 1;

                        currentRead.IndelSize = calculateIndelSize(currentRead);

                        currentRead.NT_str = getInsertedStr1(currentRead);

                        currentRead.BPLeft = currentRead. UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                        currentRead.BPRight = currentRead.UP_Far[FarIndex]. AbsLoc - g_SpacerBeforeAfter;
                        if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                           saveReadForNextCycle( currentRead, currentState.FutureReads_SR);
                        }
                        else {
                           if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                              Vars[(int) currentRead. BPLeft / BoxSize]. push_back(ReadIndex);
                              currentRead.Used = true;
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
      else if (currentRead.MatchedD == Minus) {
         for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
            for (int CloseIndex = currentRead.UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
               if (currentRead.Used) {
                  break;
               }
               if (currentRead.UP_Close[CloseIndex]. Mismatches > MAX_SNP_ERROR_index) {
                  continue;
               }
               for (int FarIndex = currentRead.UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
                  if (currentRead.Used) {
                     break;
                  }
                  if (currentRead.UP_Far[FarIndex]. Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentRead.UP_Far[FarIndex]. Mismatches + currentRead.UP_Close[CloseIndex]. Mismatches > MAX_SNP_ERROR_index) {
                     continue;
                  }
                  if (currentRead.UP_Far[FarIndex]. Direction == Plus) {
                     if (decisionBranch2(currentRead, CloseIndex, FarIndex)) {

                        currentRead.Left = currentRead.UP_Far[FarIndex]. AbsLoc - currentRead.UP_Far[FarIndex]. LengthStr + 1;
                        currentRead.Right = currentRead. UP_Close[CloseIndex].AbsLoc + currentRead. UP_Close[CloseIndex].LengthStr - 1;
                        currentRead.BP = currentRead.UP_Far[FarIndex]. LengthStr - 1;


                        currentRead.IndelSize = calculateIndelSize(currentRead);

                        currentRead.NT_str = getInsertedStr2(currentRead);

                        currentRead.BPLeft = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
                        currentRead.BPRight = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;

                        if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                           saveReadForNextCycle( currentRead, currentState.FutureReads_SR);
                        }
                        else {
                           if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                              Vars[(int) currentRead. BPLeft / BoxSize]. push_back(ReadIndex);
                              currentRead.Used = true;
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
   LOG_INFO(*logStream << "Total: " << Count_Var << "\t+" << Count_Var_Plus << "\t-" << Count_Var_Minus << std::endl);


   outputResults(currentState, Vars, numBoxes, window);

   for (unsigned int i = 0; i < numBoxes; i++) {
      Vars[i].clear();
   }

   return EXIT_SUCCESS;
}
