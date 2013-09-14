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

int searchTandemDuplicationsNT(ControlState& currentState, unsigned NumBoxes, const SearchWindow& window)
{

   static int Count_TD_NT = 0;
   static int Count_TD_NT_Plus = 0;
   static int Count_TD_NT_Minus = 0;

   std::vector<unsigned> TD_NT[NumBoxes];

   int CloseIndex = 0;
   int FarIndex = 0;

	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   LOG_INFO(*logStream << "Searching tandem duplication events with non-template sequence ... " << std::endl);
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
		SPLIT_READ& currentRead = currentState.Reads_SR[ReadIndex];
      if (currentRead.Used || currentRead.UP_Far.empty() || currentRead.FragName != currentRead.FarFragName) {
         continue;
      }
      CloseIndex = currentRead.UP_Close.size() - 1;
      FarIndex = currentRead.UP_Far.size() - 1;
      if (currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr >= currentRead.getReadLength()) {
         continue;
      }
      if (currentRead.UP_Far[FarIndex].Mismatches + currentRead.UP_Close[CloseIndex].Mismatches > 
          (short) (1 + userSettings->Seq_Error_Rate * (currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr))) {
         continue;
      }
      if (currentRead.MatchedD == Plus) {
         if (currentRead.UP_Far[FarIndex].Direction == Minus) {
            if (currentRead.UP_Far[FarIndex].AbsLoc + currentRead.UP_Far[FarIndex].LengthStr < currentRead.UP_Close[CloseIndex].AbsLoc && 
                currentRead.UP_Far[FarIndex].AbsLoc + currentRead.UP_Close[CloseIndex].LengthStr < currentRead.UP_Close[CloseIndex].AbsLoc && 
                currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr > userSettings->Min_Num_Matched_Bases) {

               currentRead.Right = currentRead.UP_Close[CloseIndex].AbsLoc - currentRead.UP_Close[CloseIndex].LengthStr + 1;
               currentRead.Left = currentRead.UP_Far[FarIndex].AbsLoc + currentRead.UP_Far[FarIndex].LengthStr - 1;
               currentRead.BP = currentRead.UP_Close[CloseIndex].LengthStr - 1;

               currentRead.IndelSize = currentRead.UP_Close[CloseIndex].AbsLoc - currentRead.UP_Far[FarIndex].AbsLoc + 1;
               currentRead.NT_size = currentRead.getReadLength() - currentRead.UP_Close[CloseIndex].LengthStr - currentRead.UP_Far[FarIndex].LengthStr;
               currentRead.NT_str = currentRead.getUnmatchedSeqRev().substr( currentRead.BP + 1, currentRead.NT_size);
               currentRead.BPRight = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
               currentRead.BPLeft = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
               if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                  saveReadForNextCycle(currentRead, currentState.FutureReads_SR);
               }
               else {
                  if ( readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                     TD_NT[(int) currentRead. BPLeft / BoxSize]. push_back(ReadIndex);
                     currentRead.Used = true;
                     Count_TD_NT++;
                     Count_TD_NT_Plus++;
                  }
               }
            }
         }
      }
      else if (currentRead.MatchedD == Minus) {
         if (currentRead.UP_Far[FarIndex].Direction == Plus) {
            if (currentRead.UP_Close[CloseIndex].AbsLoc + currentRead.UP_Close[CloseIndex].LengthStr < currentRead.UP_Far[FarIndex].AbsLoc && 
                currentRead.UP_Close[CloseIndex].AbsLoc + currentRead.UP_Far[FarIndex].LengthStr < currentRead.UP_Far[FarIndex].AbsLoc && 
                currentRead.UP_Far[FarIndex].LengthStr + currentRead.UP_Close[CloseIndex].LengthStr > userSettings->Min_Num_Matched_Bases) {

               currentRead.Right = currentRead.UP_Far[FarIndex].AbsLoc - currentRead.UP_Far[FarIndex].LengthStr + 1;
               currentRead.Left = currentRead.UP_Close[CloseIndex].AbsLoc + currentRead.UP_Close[CloseIndex].LengthStr - 1;
               currentRead.BP = currentRead.UP_Far[FarIndex].LengthStr - 1;

               currentRead.IndelSize = currentRead.UP_Far[FarIndex].AbsLoc - currentRead.UP_Close[CloseIndex].AbsLoc + 1;
               currentRead.NT_size = currentRead.getReadLength() - currentRead.UP_Close[CloseIndex].LengthStr - currentRead.UP_Far[FarIndex].LengthStr;
               currentRead.NT_str = currentRead.getUnmatchedSeq(). substr( currentRead.BP + 1, currentRead.NT_size);

               currentRead.BPRight = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
               currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
               if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                  saveReadForNextCycle(currentRead, currentState.FutureReads_SR);
               }
               else {
                  if ( readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                     TD_NT[(int) currentRead. BPLeft / BoxSize]. push_back(ReadIndex);
                     currentRead.Used = true;

                     Count_TD_NT++;
                     Count_TD_NT_Minus++;
                  }
               }
            }
         }
      }
   }
   LOG_INFO(*logStream << "Total: " << Count_TD_NT << "\t+" << Count_TD_NT_Plus << "\t-" << Count_TD_NT_Minus << std::endl);
   std::ofstream TDOutf(userSettings->getTDOutputFilename().c_str(), std::ios::app);
   SortAndOutputTandemDuplications(currentState, NumBoxes, window.getChromosome()->getSeq(), currentState.Reads_SR, TD_NT, TDOutf, true);
   TDOutf.close();
   for (unsigned int i = 0; i < NumBoxes; i++) {
      TD_NT[i].clear();
   }

   return EXIT_SUCCESS;
}
