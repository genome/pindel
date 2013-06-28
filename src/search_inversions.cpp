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
#include "search_inversions.h"


// EW: McCabe complexity through the roof here; extract some blocks
int searchInversions(ControlState& currentState, unsigned NumBoxes, const SearchWindow& window)
{
   static int Count_Inv = 0;
   static int Count_Inv_Plus = 0;
   static int Count_Inv_Minus = 0;

   std::vector<unsigned> Inv[NumBoxes];

	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
	int MIN_IndelSize_Inversion = userSettings->MIN_IndelSize_Inversion;

   int CloseIndex = 0;
   int FarIndex = 0;
   LOG_INFO(*logStream << "Searching inversions ... " << std::endl);
   for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads_SR.size(); ReadIndex++) {
		SPLIT_READ& currentRead = currentState.Reads_SR[ReadIndex];
      if (currentRead.Used || currentRead.UP_Far.empty() || currentRead.FragName != currentRead.FarFragName) {
         continue;
      }
      if (currentRead.UP_Close[0].Strand != currentRead.UP_Far[0].Strand
            && currentRead.UP_Close[0].Direction == currentRead.UP_Far[0].Direction) {

         if (currentRead.MatchedD == Plus) {
            if (currentRead.UP_Far[0].AbsLoc > currentRead.getLastAbsLocCloseEnd() + MIN_IndelSize_Inversion) { // normal situation
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
                  for (CloseIndex = (int) currentRead.UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
                     if (currentRead.Used) {
                        break;
                     }
                     if (currentRead. UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                        continue;
                     }
                     for (FarIndex = 0; FarIndex < (int) currentRead.UP_Far.size(); FarIndex++) {
                        if (currentRead.Used) {
                           break;
                        }
                        if (currentRead. UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Mismatches + currentRead. UP_Close[CloseIndex]. Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Direction == Plus) {
                           if (currentRead. UP_Far[FarIndex].LengthStr + currentRead. UP_Close[CloseIndex]. LengthStr == currentRead.getReadLength()
                                 && currentRead. UP_Far[FarIndex].AbsLoc > currentRead. UP_Close[CloseIndex]. AbsLoc + MIN_IndelSize_Inversion) {
                              currentRead.Left = (currentRead. UP_Close[CloseIndex]. AbsLoc + 1) - currentRead. UP_Close[CloseIndex]. LengthStr;
                              currentRead.Right = currentRead. UP_Far[FarIndex].AbsLoc - currentRead. UP_Far[FarIndex]. LengthStr + currentRead.getReadLength();
                              currentRead.BP = currentRead. UP_Close[CloseIndex]. LengthStr - 1;

                              currentRead.IndelSize = currentRead. UP_Far[FarIndex].AbsLoc - currentRead. UP_Close[CloseIndex]. AbsLoc;
                              currentRead.NT_str = "";
                              currentRead.NT_size = 0;

                              currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc + 1 - g_SpacerBeforeAfter;
                              currentRead.BPRight = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
                               LeftMostINV(currentState, currentRead, window);
                              if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                                 saveReadForNextCycle( currentRead, currentState.FutureReads_SR);
                              }
                              else {
                                 if (readInSpecifiedRegion( currentRead,userSettings->getRegion())) {
                                    Inv[(int) (currentRead. BPLeft + currentRead. BPRight) / BoxSize]. push_back( ReadIndex);
                                    currentRead.Used = true;
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
            
            else if (currentRead. UP_Far[currentRead.UP_Far.size() - 1].AbsLoc + MIN_IndelSize_Inversion  < currentRead.UP_Close[0].AbsLoc) {
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
                  for (CloseIndex = 0; CloseIndex < (int) currentRead.UP_Close.size(); CloseIndex++) {
                     if (currentRead.Used) {
                        break;
                     }
                     if (currentRead. UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                        continue;
                     }
                     for (FarIndex = (int) currentRead.UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
                        if (currentRead.Used) {
                           break;
                        }
                        if (currentRead. UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Mismatches + currentRead. UP_Close[CloseIndex]. Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Direction == Plus) {
                            
                           // anchor inside reversed block.
                           if (currentRead. UP_Far[FarIndex].LengthStr + currentRead. UP_Close[CloseIndex]. LengthStr == currentRead.getReadLength()
                                 && currentRead. UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion < currentRead. UP_Close[CloseIndex].AbsLoc) {
                              currentRead.Right = currentRead. UP_Close[CloseIndex]. AbsLoc - currentRead. UP_Close[CloseIndex]. LengthStr + currentRead.getReadLength();
                              currentRead.Left = currentRead.UP_Far[FarIndex].AbsLoc - currentRead.UP_Far[FarIndex].LengthStr + 1;
                              currentRead.BP = currentRead. UP_Far[FarIndex]. LengthStr - 1;

                              currentRead.IndelSize = currentRead. UP_Close[CloseIndex]. AbsLoc - currentRead. UP_Far[FarIndex].AbsLoc;
                              currentRead.NT_str = "";
                              currentRead.NT_size = 0;

                              currentRead.BPRight = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                              currentRead.BPLeft = (currentRead.UP_Far[FarIndex].AbsLoc + 1) - g_SpacerBeforeAfter;
                               LeftMostINV(currentState, currentRead, window);
                              if (readTransgressesBinBoundaries( currentRead, window.getEnd())) {
                                 saveReadForNextCycle( currentRead, currentState.FutureReads_SR);
                              }
                              else {
                                 if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                                    Inv[(int) (currentRead. BPLeft + currentRead. BPRight) / BoxSize]. push_back( ReadIndex);
                                    currentRead.Used = true;
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
              
         else if (currentRead.MatchedD == Minus) {
            if (currentRead. UP_Close[currentRead.UP_Close.size() - 1].AbsLoc > currentRead.UP_Far[0].AbsLoc + MIN_IndelSize_Inversion) {
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
                  for (CloseIndex = (int) currentRead.UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
                     if (currentRead.Used) {
                        break;
                     }
                     if (currentRead. UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) {
                        continue;
                     }
                     for (FarIndex = 0; FarIndex < (int) currentRead.UP_Far.size(); FarIndex++) {
                        if (currentRead.Used) {
                           break;
                        }
                        if (currentRead. UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Mismatches + currentRead. UP_Close[CloseIndex]. Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Direction == Minus) {
                           if (currentRead. UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr == currentRead.getReadLength()
                                 && currentRead.UP_Close[CloseIndex]. AbsLoc > currentRead.UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion) {
                              currentRead.Left = currentRead. UP_Far[FarIndex].AbsLoc + currentRead. UP_Far[FarIndex]. LengthStr - currentRead.getReadLength();
                              currentRead.Right = currentRead. UP_Close[CloseIndex]. AbsLoc + currentRead. UP_Close[CloseIndex]. LengthStr - 1;
                              currentRead.BP = currentRead. UP_Far[FarIndex]. LengthStr - 1;

                              currentRead.IndelSize = currentRead. UP_Close[CloseIndex]. AbsLoc - currentRead. UP_Far[FarIndex].AbsLoc;
                              currentRead.NT_str = "";
                              currentRead.NT_size = 0;

                              currentRead.BPLeft = currentRead.UP_Far[FarIndex].AbsLoc - g_SpacerBeforeAfter;
                              currentRead.BPRight = currentRead.UP_Close[CloseIndex].AbsLoc - 1 - g_SpacerBeforeAfter;
                               LeftMostINV(currentState, currentRead, window);
                              if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                                 Inv[(int) (currentRead. BPLeft + currentRead. BPRight) / BoxSize]. push_back( ReadIndex);
                                 currentRead. Used = true;
                                 Count_Inv++;
                                 Count_Inv_Minus++;
                              }
                           }
                        }
                     }
                  }
               }
            }
             
            else if (currentRead.UP_Close[0].AbsLoc + MIN_IndelSize_Inversion < currentRead.UP_Far[currentRead. UP_Far.size() - 1].AbsLoc) {
               for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= currentRead.getMAX_SNP_ERROR(); MAX_SNP_ERROR_index++) {
                  for (unsigned int CloseIndex = 0; CloseIndex < currentRead.UP_Close.size(); CloseIndex++) {
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
                        if (currentRead. UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Mismatches + currentRead. UP_Close[CloseIndex]. Mismatches > MAX_SNP_ERROR_index) {
                           continue;
                        }
                        if (currentRead. UP_Far[FarIndex].Direction == Minus) {
                            
                           // anchor inside reversed block.
                           if (currentRead. UP_Close[CloseIndex].LengthStr + currentRead. UP_Far[FarIndex]. LengthStr == currentRead.getReadLength()
                                 && currentRead. UP_Close[CloseIndex]. AbsLoc + MIN_IndelSize_Inversion < currentRead. UP_Far[FarIndex].AbsLoc) {
                              currentRead.Right = currentRead. UP_Far[FarIndex].AbsLoc + currentRead. UP_Far[FarIndex]. LengthStr - 1;
                              currentRead.Left = currentRead. UP_Close[CloseIndex]. AbsLoc + currentRead. UP_Close[CloseIndex]. LengthStr - currentRead.getReadLength();
                              currentRead.BP = currentRead. UP_Close[CloseIndex]. LengthStr - 1;

                              currentRead.IndelSize = currentRead. UP_Far[FarIndex].AbsLoc - currentRead. UP_Close[CloseIndex]. AbsLoc;
                              currentRead.NT_str = "";
                              currentRead.NT_size = 0;

                              currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                              currentRead.BPRight = currentRead.UP_Far[FarIndex].AbsLoc - 1 - g_SpacerBeforeAfter;
                               LeftMostINV(currentState, currentRead, window);
                              if (readInSpecifiedRegion( currentRead, userSettings->getRegion())) {
                                 Inv[(int) (currentRead. BPLeft + currentRead. BPRight) / BoxSize]. push_back( ReadIndex);
                                 currentRead. Used = true;

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

   LOG_INFO(*logStream << "Total: " << Count_Inv << "\t+" << Count_Inv_Plus << "\t-" << Count_Inv_Minus << std::endl);
   std::ofstream InversionOutf(userSettings->getINVOutputFilename().c_str(), std::ios::app);
   SortOutputInv(currentState, NumBoxes, window.getChromosome()->getSeq(), currentState.Reads_SR, Inv, InversionOutf);
   for (unsigned int i = 0; i < NumBoxes; i++) {
      Inv[i].clear();
   }

   return EXIT_SUCCESS;
}

void LeftMostINV(ControlState& currentState, SPLIT_READ & currentRead, const SearchWindow& window) {
    const std::string & TheInput = window.getChromosome()->getSeq();
    
    unsigned int PosIndex = currentRead.BPLeft + g_SpacerBeforeAfter + 1;
    unsigned int original_PosIndex = PosIndex;
    //unsigned int Start = PosIndex + 1;
    unsigned int End = currentRead.BPRight + g_SpacerBeforeAfter - 1;
    while (TheInput[PosIndex] == Convert2RC4N[(short) TheInput[End]]) {
        ++PosIndex;
        --End;
    }
    short DIFF = PosIndex - original_PosIndex;
    if (DIFF > 0) {
        if (currentRead.MatchedD == Plus) { // BP++
            if (DIFF >= currentRead.BP) DIFF = currentRead.BP - 1;
        }
        else { //BP--
            if (DIFF + currentRead.BP >= currentRead.getReadLength()) DIFF = currentRead.getReadLength() - currentRead.BP - 1;
            currentRead.BPLeft += DIFF;
            currentRead.BPRight -= DIFF;
            currentRead.BP += DIFF;
        }
    } // Convert2RC4N[(short) CurrentBase];
}

