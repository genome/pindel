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

SearchVariant::SearchVariant() {
	Count_Var = 0;
	Count_Var_Plus = 0;
	Count_Var_Minus = 0;
	typeOfVariant = "some type of variant, replace this with the correct name in child class";
}

SearchVariant::~SearchVariant() {

}

int SearchVariant::Search(ControlState& currentState, const unsigned numBoxes) {

	std::vector<unsigned> Vars[numBoxes];

	LOG_INFO(std::cout << "Searching " << typeOfVariant << " ... " << std::endl);
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
							//							if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
							//									== currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
							//											+ 1
							//									&& currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
							//											+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
							//											< currentState.Reads[ReadIndex].ReadLength)
							if (decisionBranch1(currentState, ReadIndex, CloseIndex, FarIndex)) {

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

								//								currentState.Reads[ReadIndex].IndelSize
								//										= currentState.Reads[ReadIndex].ReadLengthMinus
								//												- (currentState.Reads[ReadIndex].Right
								//														- currentState.Reads[ReadIndex].Left);
								currentState.Reads[ReadIndex].IndelSize
										= calculateIndelSize(currentState,
												ReadIndex);

								//								currentState.Reads[ReadIndex].InsertedStr
								//										= ReverseComplement(
								//												currentState.Reads[ReadIndex]. UnmatchedSeq). substr(
								//												currentState.Reads[ReadIndex].BP
								//														+ 1,
								//												currentState.Reads[ReadIndex]. IndelSize);
								currentState.Reads[ReadIndex].InsertedStr
										= getInsertedStr1(currentState, ReadIndex);

								currentState.Reads[ReadIndex].BPLeft
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												- SpacerBeforeAfter;
								currentState.Reads[ReadIndex].BPRight
										= currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
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
											currentState.regionStartDefined,
											currentState.regionEndDefined,
											currentState.startOfRegion,
											currentState.endOfRegion)) {
										Vars[(int) currentState.Reads[ReadIndex]. BPLeft
												/ BoxSize]. push_back(ReadIndex);
										currentState.Reads[ReadIndex].Used
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
							//							if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
							//									== currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
							//											+ 1
							//									&& currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
							//											+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
							//											< currentState.Reads[ReadIndex].ReadLength) {
							if (decisionBranch2(currentState, ReadIndex, CloseIndex, FarIndex)) {

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

								//								currentState.Reads[ReadIndex].IndelSize
								//										= currentState.Reads[ReadIndex].ReadLengthMinus
								//												- (currentState.Reads[ReadIndex].Right
								//														- currentState.Reads[ReadIndex].Left);
								currentState.Reads[ReadIndex].IndelSize
										= calculateIndelSize(currentState,
												ReadIndex);

								//								currentState.Reads[ReadIndex].InsertedStr
								//										= currentState.Reads[ReadIndex].UnmatchedSeq. substr(
								//												currentState.Reads[ReadIndex].BP
								//														+ 1,
								//												currentState.Reads[ReadIndex]. IndelSize);
								currentState.Reads[ReadIndex].InsertedStr
										= getInsertedStr2(currentState,
												ReadIndex);

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
											currentState.regionStartDefined,
											currentState.regionEndDefined,
											currentState.startOfRegion,
											currentState.endOfRegion)) {
										Vars[(int) currentState.Reads[ReadIndex]. BPLeft
												/ BoxSize]. push_back(ReadIndex);
										currentState.Reads[ReadIndex].Used
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
	LOG_INFO(std::cout << "Total: " << Count_Var << "\t+" << Count_Var_Plus << "\t-"
			<< Count_Var_Minus << std::endl);

//	std::ofstream SIoutputfile(currentState.SIOutputFilename.c_str(),
//			std::ios::app);
//	SortOutputSI(NumBoxes, currentState.CurrentChr, currentState.Reads, Vars,
//			SIoutputfile);
//	SIoutputfile.close();
	outputResults(currentState, Vars, numBoxes);

	for (unsigned int i = 0; i < numBoxes; i++)
		Vars[i].clear();

	return EXIT_SUCCESS;
}
