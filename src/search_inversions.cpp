// Pindel header files
#include "reporter.h"
#include "ControlState.h"

int searchInversions(ControlState& currentState, unsigned NumBoxes) {

	static int Count_Inv = 0;
	static int Count_Inv_Plus = 0;
	static int Count_Inv_Minus = 0;

	std::vector<unsigned> Inv[NumBoxes];

	int CloseIndex = 0;
	int FarIndex = 0;

	std::cout << "Searching inversions ... " << std::endl;
	for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads.size(); ReadIndex++) {
		if (currentState.Reads[ReadIndex].Used
				|| currentState.Reads[ReadIndex].UP_Far.empty())
			continue;
		if (currentState.Reads[ReadIndex].UP_Close[0].Strand
				!= currentState.Reads[ReadIndex].UP_Far[0].Strand
				&& currentState.Reads[ReadIndex].UP_Close[0].Direction
						== currentState.Reads[ReadIndex].UP_Far[0].Direction) {

			if (currentState.Reads[ReadIndex].MatchedD == Plus) {
				if (currentState.Reads[ReadIndex].UP_Far[0].AbsLoc
						> currentState.Reads[ReadIndex].UP_Close[currentState.Reads[ReadIndex]. UP_Close.size()
								- 1].AbsLoc + MIN_IndelSize_Inversion) { // normal situation
					for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
							<= currentState.Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
						for (CloseIndex
								= (int) currentState.Reads[ReadIndex].UP_Close.size()
										- 1; CloseIndex >= 0; CloseIndex--) {
							if (currentState.Reads[ReadIndex].Used)
								break;
							if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex].Mismatches
									> MAX_SNP_ERROR_index)
								continue;
							for (FarIndex = 0; FarIndex
									< (int) currentState.Reads[ReadIndex].UP_Far.size(); FarIndex++) {
								if (currentState.Reads[ReadIndex].Used)
									break;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Direction
										== Plus) {
									if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
											+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
											== currentState.Reads[ReadIndex].ReadLength
											&& currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
													> currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
															+ MIN_IndelSize_Inversion) {
										currentState.Reads[ReadIndex].Left
												= (currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
														+ 1)
														- currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr;
										currentState.Reads[ReadIndex].Right
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
														- currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
														+ currentState.Reads[ReadIndex]. ReadLength;
										currentState.Reads[ReadIndex].BP
												= currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
														- 1;

										currentState.Reads[ReadIndex].IndelSize
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
														- currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc;
										currentState.Reads[ReadIndex].NT_str
												= "";
										currentState.Reads[ReadIndex].NT_size
												= 0;
										currentState.Reads[ReadIndex]. InsertedStr
												= "";
										currentState.Reads[ReadIndex].BPLeft
												= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
														+ 1 - SpacerBeforeAfter;
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
												Inv[(int) currentState.Reads[ReadIndex]. BPLeft
														/ BoxSize]. push_back(
														ReadIndex);
												currentState.Reads[ReadIndex]. Used
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
				} else if (currentState.Reads[ReadIndex]. UP_Far[currentState.Reads[ReadIndex].UP_Far.size()
						- 1].AbsLoc + MIN_IndelSize_Inversion
						< currentState.Reads[ReadIndex].UP_Close[0].AbsLoc) {
					for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
							<= currentState.Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
						for (CloseIndex = 0; CloseIndex
								< (int) currentState.Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
							if (currentState.Reads[ReadIndex].Used)
								break;
							if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex].Mismatches
									> MAX_SNP_ERROR_index)
								continue;
							for (FarIndex
									= (int) currentState.Reads[ReadIndex].UP_Far.size()
											- 1; FarIndex >= 0; FarIndex--) {
								if (currentState.Reads[ReadIndex].Used)
									break;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Direction
										== Plus) {
									// anchor inside reversed block.
									if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
											+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
											== currentState.Reads[ReadIndex].ReadLength
											&& currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
													+ MIN_IndelSize_Inversion
													< currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc) {
										currentState.Reads[ReadIndex].Right
												= currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
														- currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
														+ currentState.Reads[ReadIndex]. ReadLength;
										currentState.Reads[ReadIndex].Left
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
														- currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
														+ 1;
										currentState.Reads[ReadIndex].BP
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
														- 1;

										currentState.Reads[ReadIndex].IndelSize
												= currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
														- currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc;
										currentState.Reads[ReadIndex].NT_str
												= "";
										currentState.Reads[ReadIndex].NT_size
												= 0;
										currentState.Reads[ReadIndex]. InsertedStr
												= "";
										currentState.Reads[ReadIndex].BPRight
												= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
														- SpacerBeforeAfter;
										currentState.Reads[ReadIndex].BPLeft
												= (currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
														+ 1)
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
												Inv[(int) currentState.Reads[ReadIndex]. BPLeft
														/ BoxSize]. push_back(
														ReadIndex);
												currentState.Reads[ReadIndex]. Used
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
			} else if (currentState.Reads[ReadIndex].MatchedD == Minus) {
				if (currentState.Reads[ReadIndex]. UP_Close[currentState.Reads[ReadIndex].UP_Close.size()
						- 1].AbsLoc
						> currentState.Reads[ReadIndex].UP_Far[0].AbsLoc
								+ MIN_IndelSize_Inversion) {
					for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
							<= currentState.Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
						for (CloseIndex
								= (int) currentState.Reads[ReadIndex].UP_Close.size()
										- 1; CloseIndex >= 0; CloseIndex--) {
							if (currentState.Reads[ReadIndex].Used)
								break;
							if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex].Mismatches
									> MAX_SNP_ERROR_index)
								continue;
							for (FarIndex = 0; FarIndex
									< (int) currentState.Reads[ReadIndex].UP_Far.size(); FarIndex++) {
								if (currentState.Reads[ReadIndex].Used)
									break;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Direction
										== Minus) {
									if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
											+ currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
											== currentState.Reads[ReadIndex].ReadLength
											&& currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
													> currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
															+ MIN_IndelSize_Inversion) {
										currentState.Reads[ReadIndex].Left
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
														+ currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
														- currentState.Reads[ReadIndex]. ReadLength;
										currentState.Reads[ReadIndex].Right
												= currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
														+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
														- 1;
										currentState.Reads[ReadIndex].BP
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
														- 1;

										currentState.Reads[ReadIndex].IndelSize
												= currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
														- currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc;
										currentState.Reads[ReadIndex].NT_str
												= "";
										currentState.Reads[ReadIndex].NT_size
												= 0;
										currentState.Reads[ReadIndex]. InsertedStr
												= "";
										currentState.Reads[ReadIndex].BPLeft
												= currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
														- SpacerBeforeAfter;
										currentState.Reads[ReadIndex].BPRight
												= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
														- 1 - SpacerBeforeAfter;
										if (readInSpecifiedRegion(
												currentState.Reads[ReadIndex],
												currentState.startOfRegion,
												currentState.endOfRegion)) {
											Inv[(int) currentState.Reads[ReadIndex]. BPLeft
													/ BoxSize]. push_back(
													ReadIndex);
											currentState.Reads[ReadIndex]. Used
													= true;

											Count_Inv++;
											Count_Inv_Minus++;
										}
									}
								}
							}
						}
					}
				} else if (currentState.Reads[ReadIndex].UP_Close[0].AbsLoc
						+ MIN_IndelSize_Inversion
						< currentState.Reads[ReadIndex].UP_Far[currentState.Reads[ReadIndex]. UP_Far.size()
								- 1].AbsLoc) {
					for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index
							<= currentState.Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
						for (unsigned int CloseIndex = 0; CloseIndex
								< currentState.Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
							if (currentState.Reads[ReadIndex].Used)
								break;
							if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex].Mismatches
									> MAX_SNP_ERROR_index)
								continue;
							for (int FarIndex =
									currentState.Reads[ReadIndex].UP_Far.size()
											- 1; FarIndex >= 0; FarIndex--) {
								if (currentState.Reads[ReadIndex].Used)
									break;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Mismatches
										+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. Mismatches
										> MAX_SNP_ERROR_index)
									continue;
								if (currentState.Reads[ReadIndex]. UP_Far[FarIndex].Direction
										== Minus) {
									// anchor inside reversed block.
									if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
											+ currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
											== currentState.Reads[ReadIndex].ReadLength
											&& currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
													+ MIN_IndelSize_Inversion
													< currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc) {
										currentState.Reads[ReadIndex].Right
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
														+ currentState.Reads[ReadIndex]. UP_Far[FarIndex]. LengthStr
														- 1;
										currentState.Reads[ReadIndex].Left
												= currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc
														+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
														- currentState.Reads[ReadIndex]. ReadLength;
										currentState.Reads[ReadIndex].BP
												= currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. LengthStr
														- 1;

										currentState.Reads[ReadIndex].IndelSize
												= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
														- currentState.Reads[ReadIndex]. UP_Close[CloseIndex]. AbsLoc;
										currentState.Reads[ReadIndex].NT_str
												= "";
										currentState.Reads[ReadIndex].NT_size
												= 0;
										currentState.Reads[ReadIndex]. InsertedStr
												= "";
										currentState.Reads[ReadIndex].BPLeft
												= currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
														- SpacerBeforeAfter;
										currentState.Reads[ReadIndex].BPRight
												= currentState.Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
														- 1 - SpacerBeforeAfter;
										if (readInSpecifiedRegion(
												currentState.Reads[ReadIndex],
												currentState.startOfRegion,
												currentState.endOfRegion)) {
											Inv[(int) currentState.Reads[ReadIndex]. BPLeft
													/ BoxSize]. push_back(
													ReadIndex);
											currentState.Reads[ReadIndex]. Used
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
	std::cout << "Total: " << Count_Inv << "\t+" << Count_Inv_Plus << "\t-"
			<< Count_Inv_Minus << std::endl;
	std::ofstream InversionOutf(currentState.InversionOutputFilename.c_str(),
			std::ios::app);
	SortOutputInv(NumBoxes, currentState.CurrentChr, currentState.Reads, Inv,
			InversionOutf);
	for (unsigned int i = 0; i < NumBoxes; i++)
		Inv[i].clear();

	return EXIT_SUCCESS;
}
