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

int searchTandemDuplicationsNT(ControlState& currentState, unsigned NumBoxes) {

	static int Count_TD_NT = 0;
	static int Count_TD_NT_Plus = 0;
	static int Count_TD_NT_Minus = 0;

	std::vector<unsigned> TD_NT[NumBoxes];

	int CloseIndex = 0;
	int FarIndex = 0;

	LOG_INFO(std::cout
			<< "Searching tandem duplication events with non-template sequence ... "
			<< std::endl);
	for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads.size(); ReadIndex++) {
		if (currentState.Reads[ReadIndex].Used
				|| currentState.Reads[ReadIndex].UP_Far.empty())
			continue;
		CloseIndex = currentState.Reads[ReadIndex].UP_Close.size() - 1;
		FarIndex = currentState.Reads[ReadIndex].UP_Far.size() - 1;
		if (currentState.Reads[ReadIndex].UP_Far[FarIndex].LengthStr
				+ currentState.Reads[ReadIndex].UP_Close[CloseIndex].LengthStr
				>= currentState.Reads[ReadIndex].ReadLength)
			continue;
		if (currentState.Reads[ReadIndex].UP_Far[FarIndex].Mismatches
				+ currentState.Reads[ReadIndex].UP_Close[CloseIndex].Mismatches
				> (short) (1
						+ Seq_Error_Rate
								* (currentState.Reads[ReadIndex].UP_Far[FarIndex].LengthStr
										+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr)))
			continue;
		if (currentState.Reads[ReadIndex].MatchedD == Plus) {
			if (currentState.Reads[ReadIndex].UP_Far[FarIndex].Direction
					== Minus) {
				if (currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
						+ currentState.Reads[ReadIndex].ReadLength
						< currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
						&& currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
								+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
								> Min_Num_Matched_Bases) {
					currentState.Reads[ReadIndex].Right
							= currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
									- currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
									+ 1;
					currentState.Reads[ReadIndex].Left
							= currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
									+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
									- 1;
					currentState.Reads[ReadIndex].BP
							= currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
									- 1;

					currentState.Reads[ReadIndex].IndelSize
							= currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
									- currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
									+ 1;
					currentState.Reads[ReadIndex].NT_size
							= currentState.Reads[ReadIndex].ReadLength
									- currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
									- currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr;
					currentState.Reads[ReadIndex].NT_str
							= ReverseComplement(
									currentState.Reads[ReadIndex]. UnmatchedSeq). substr(
									currentState.Reads[ReadIndex].BP + 1,
									currentState.Reads[ReadIndex].NT_size);
					currentState.Reads[ReadIndex].InsertedStr = "";
					currentState.Reads[ReadIndex].BPRight
							= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									- SpacerBeforeAfter;
					currentState.Reads[ReadIndex].BPLeft
							= currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
									- SpacerBeforeAfter;
					if (readTransgressesBinBoundaries(
							currentState.Reads[ReadIndex],
							currentState.upperBinBorder)) {
						saveReadForNextCycle(currentState.Reads[ReadIndex],
								currentState.FutureReads);
					} else {
						if (currentState.Reads[ReadIndex].NT_size
								<= Max_Length_NT && readInSpecifiedRegion(
								currentState.Reads[ReadIndex],
								currentState.startOfRegion,
								currentState.endOfRegion)) {
							TD_NT[(int) currentState.Reads[ReadIndex]. BPLeft
									/ BoxSize]. push_back(ReadIndex);
							currentState.Reads[ReadIndex].Used = true;
							Count_TD_NT++;
							Count_TD_NT_Plus++;
						}
					}
				}
			}
		} else if (currentState.Reads[ReadIndex].MatchedD == Minus) {
			if (currentState.Reads[ReadIndex].UP_Far[FarIndex].Direction
					== Plus) {
				if (currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
						+ currentState.Reads[ReadIndex].ReadLength
						< currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
						&& currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
								+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
								> Min_Num_Matched_Bases) {

					currentState.Reads[ReadIndex].Right
							= currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
									- currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
									+ 1;
					currentState.Reads[ReadIndex].Left
							= currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
									+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
									- 1;
					currentState.Reads[ReadIndex].BP
							= currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
									- 1;

					currentState.Reads[ReadIndex].IndelSize
							= currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
									- currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
									+ 1;
					currentState.Reads[ReadIndex].NT_size
							= currentState.Reads[ReadIndex].ReadLength
									- currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
									- currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr;
					currentState.Reads[ReadIndex].NT_str
							= currentState.Reads[ReadIndex].UnmatchedSeq. substr(
									currentState.Reads[ReadIndex].BP + 1,
									currentState.Reads[ReadIndex].NT_size);
					currentState.Reads[ReadIndex].InsertedStr = "";
					currentState.Reads[ReadIndex].BPRight
							= currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
									- SpacerBeforeAfter;
					currentState.Reads[ReadIndex].BPLeft
							= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									- SpacerBeforeAfter;
					if (readTransgressesBinBoundaries(
							currentState.Reads[ReadIndex],
							currentState.upperBinBorder)) {
						saveReadForNextCycle(currentState.Reads[ReadIndex],
								currentState.FutureReads);
					} else {
						if (currentState.Reads[ReadIndex].NT_size
								<= Max_Length_NT && readInSpecifiedRegion(
								currentState.Reads[ReadIndex],
								currentState.startOfRegion,
								currentState.endOfRegion)) {
							TD_NT[(int) currentState.Reads[ReadIndex]. BPLeft
									/ BoxSize]. push_back(ReadIndex);
							currentState.Reads[ReadIndex].Used = true;

							Count_TD_NT++;
							Count_TD_NT_Minus++;
						}
					}
				}
			}
		}
	}
	LOG_INFO(std::cout << "Total: " << Count_TD_NT << "\t+" << Count_TD_NT_Plus << "\t-"
			<< Count_TD_NT_Minus << std::endl);
	std::ofstream TDOutf(currentState.TDOutputFilename.c_str(), std::ios::app);
	SortAndOutputTandemDuplications(NumBoxes, currentState.CurrentChr,
			currentState.Reads, TD_NT, TDOutf, true);
	TDOutf.close();
	for (unsigned int i = 0; i < NumBoxes; i++)
		TD_NT[i].clear();

	return EXIT_SUCCESS;
}
