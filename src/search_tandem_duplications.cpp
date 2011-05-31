// Pindel header files
#include "reporter.h"
#include "ControlState.h"

int searchTandemDuplications(ControlState& currentState, unsigned NumBoxes) {

	static int Count_TD = 0;
	static int Count_TD_Plus = 0;
	static int Count_TD_Minus = 0;

	std::vector<unsigned> TD[NumBoxes];

	std::cout << "Searching tandem duplication events ... " << std::endl;
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
									+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
									== currentState.Reads[ReadIndex].ReadLength
									&& currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
											+ currentState.Reads[ReadIndex].ReadLength
											< currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc) {

								currentState.Reads[ReadIndex].Right
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												- currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												+ 1;
								currentState.Reads[ReadIndex].Left
										= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
												+ currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
												- 1;
								currentState.Reads[ReadIndex].BP
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												- 1;

								currentState.Reads[ReadIndex].IndelSize
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												- currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
												+ 1;
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
									saveReadForNextCycle(
											currentState.Reads[ReadIndex],
											currentState.FutureReads);
								} else {
									if (readInSpecifiedRegion(
											currentState.Reads[ReadIndex],
											currentState.startOfRegion,
											currentState.endOfRegion)) {
										TD[(int) currentState.Reads[ReadIndex]. BPLeft
												/ BoxSize]. push_back(ReadIndex);
										currentState.Reads[ReadIndex].Used
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
					for (int FarIndex = 0; FarIndex
							< (int) currentState.Reads[ReadIndex].UP_Far.size(); FarIndex++) {
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
							if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
									+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
									== currentState.Reads[ReadIndex].ReadLength
									&& currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
											+ currentState.Reads[ReadIndex].ReadLength
											< currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc) {
								currentState.Reads[ReadIndex].Right
										= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
												- currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
												+ 1;
								currentState.Reads[ReadIndex].Left
										= currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												- 1;
								currentState.Reads[ReadIndex].BP
										= currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
												- 1;

								currentState.Reads[ReadIndex].IndelSize
										= currentState.Reads[ReadIndex]. UP_Far[FarIndex].AbsLoc
												- currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												+ 1;
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
									saveReadForNextCycle(
											currentState.Reads[ReadIndex],
											currentState.FutureReads);
								} else {
									if (readInSpecifiedRegion(
											currentState.Reads[ReadIndex],
											currentState.startOfRegion,
											currentState.endOfRegion)) {
										TD[(int) currentState.Reads[ReadIndex]. BPLeft
												/ BoxSize]. push_back(ReadIndex);
										currentState.Reads[ReadIndex].Used
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
	std::cout << "Total: " << Count_TD << "\t+" << Count_TD_Plus << "\t-"
			<< Count_TD_Minus << std::endl;
	std::ofstream TDOutf(currentState.TDOutputFilename.c_str(), std::ios::app);
	SortAndOutputTandemDuplications(NumBoxes, currentState.CurrentChr,
			currentState.Reads, TD, TDOutf, false);
	for (unsigned int i = 0; i < NumBoxes; i++)
		TD[i].clear();

	return EXIT_SUCCESS;
}
