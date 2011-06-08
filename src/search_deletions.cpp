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

int searchDeletions(ControlState& currentState, unsigned NumBoxes) {

	static int Count_D = 0;
	static int Count_D_Plus = 0;
	static int Count_D_Minus = 0;

	std::vector<unsigned> Deletions[NumBoxes];

	LOG_INFO(std::cout << "Searching deletion events ... " << std::endl);
	for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads.size(); ReadIndex++) {
		if (currentState.Reads[ReadIndex].Used
				|| currentState.Reads[ReadIndex].UP_Far.empty())
			continue;

		if (currentState.Reads[ReadIndex].MatchedD == Plus) {
			for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
					<= currentState.Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
				for (unsigned int CloseIndex = 0; CloseIndex
						< currentState.Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
					if (currentState.Reads[ReadIndex].Used)
						break;
					if (currentState.Reads[ReadIndex].UP_Close[CloseIndex]. Mismatches
							> MAX_SNP_ERROR_index)
						continue;
					for (int FarIndex =
							currentState.Reads[ReadIndex].UP_Far.size() - 1; FarIndex
							>= 0; FarIndex--) {
						if (currentState.Reads[ReadIndex].Used)
							break;
						if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. Mismatches
								> MAX_SNP_ERROR_index)
							continue;
						if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. Mismatches
								+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. Mismatches
								> MAX_SNP_ERROR_index)
							continue;
						if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. Direction
								== Minus) {
							if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
									+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
									== currentState.Reads[ReadIndex].ReadLength
									&& currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
											> currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
													+ 1) {
								currentState.Reads[ReadIndex].Left
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												- currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												+ 1;
								currentState.Reads[ReadIndex].Right
										= currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
												+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
												- 1;
								currentState.Reads[ReadIndex].BP
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												- 1;

								currentState.Reads[ReadIndex].IndelSize
										= (currentState.Reads[ReadIndex].Right
												- currentState.Reads[ReadIndex].Left)
												- currentState.Reads[ReadIndex].ReadLengthMinus;
								currentState.Reads[ReadIndex].BPLeft
										= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
												- SpacerBeforeAfter;
								currentState.Reads[ReadIndex].BPRight
										= currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
												- SpacerBeforeAfter;
								if (readTransgressesBinBoundaries(
										currentState.Reads[ReadIndex],
										currentState.upperBinBorder)) {
									saveReadForNextCycle(
											currentState.Reads[ReadIndex],
											currentState.FutureReads);
								} else {
									if (readInSpecifiedRegion(
											currentState.Reads[ReadIndex],
											currentState.startOfRegion,
											currentState.endOfRegion)) {
										Deletions[(int) currentState.Reads[ReadIndex]. BPLeft
												/ BoxSize]. push_back(ReadIndex);
										currentState.Reads[ReadIndex].Used
												= true;
										Count_D++;
										Count_D_Plus++;
									}
								}
							}
						}
					}
				}
			}
		} else if (currentState.Reads[ReadIndex].MatchedD == Minus) {
			for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
					<= currentState.Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
				for (int CloseIndex =
						currentState.Reads[ReadIndex].UP_Close.size() - 1; CloseIndex
						>= 0; CloseIndex--) {
					if (currentState.Reads[ReadIndex].Used)
						break;
					if (currentState.Reads[ReadIndex].UP_Close[CloseIndex]. Mismatches
							> MAX_SNP_ERROR_index)
						continue;
					for (int FarIndex =
							currentState.Reads[ReadIndex].UP_Far.size() - 1; FarIndex
							>= 0; FarIndex--) {
						if (currentState.Reads[ReadIndex].Used)
							break;
						if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. Mismatches
								> MAX_SNP_ERROR_index)
							continue;
						if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. Mismatches
								+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. Mismatches
								> MAX_SNP_ERROR_index)
							continue;
						if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. Direction
								== Plus) {
							if (currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
									+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
									== currentState.Reads[ReadIndex].ReadLength
									&& currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
											> currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
													+ 1) {

								currentState.Reads[ReadIndex].Left
										= currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
												- currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
												+ 1;
								currentState.Reads[ReadIndex].Right
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												- 1;
								currentState.Reads[ReadIndex].BP
										= currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
												- 1;

								currentState.Reads[ReadIndex].IndelSize
										= (currentState.Reads[ReadIndex].Right
												- currentState.Reads[ReadIndex].Left)
												- currentState.Reads[ReadIndex].ReadLengthMinus;
								currentState.Reads[ReadIndex].BPLeft
										= currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
												- SpacerBeforeAfter;
								currentState.Reads[ReadIndex].BPRight
										= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
												- SpacerBeforeAfter;
								if (readTransgressesBinBoundaries(
										currentState.Reads[ReadIndex],
										currentState.upperBinBorder)) {
									saveReadForNextCycle(
											currentState.Reads[ReadIndex],
											currentState.FutureReads);
								} else {
									if (readInSpecifiedRegion(
											currentState.Reads[ReadIndex],
											currentState.startOfRegion,
											currentState.endOfRegion)) {
										Deletions[(int) currentState.Reads[ReadIndex]. BPLeft
												/ BoxSize]. push_back(ReadIndex);
										currentState.Reads[ReadIndex].Used
												= true;
										Count_D++;
										Count_D_Minus++;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	LOG_INFO(std::cout << "Total: " << Count_D << "\t+" << Count_D_Plus << "\t-"
			<< Count_D_Minus << std::endl);
	std::ofstream DeletionOutf(currentState.DeletionOutputFilename.c_str(),
			std::ios::app);
	SortOutputD(NumBoxes, currentState.CurrentChr, currentState.Reads,
			Deletions, DeletionOutf);
	DeletionOutf.close();

	for (unsigned int i = 0; i < NumBoxes; i++)
		Deletions[i].clear();

	return EXIT_SUCCESS;
}
