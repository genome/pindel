// Pindel header files
#include "reporter.h"
#include "ControlState.h"

int searchShortInsertions(ControlState& currentState, unsigned NumBoxes) {

	static int Count_SI = 0;
	static int Count_SI_Plus = 0;
	static int Count_SI_Minus = 0;

	std::vector<unsigned> SIs[NumBoxes];

	std::cout << "Searching short insertions ... " << std::endl;
	for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads.size(); ReadIndex++) {
		if (currentState.Reads[ReadIndex].Used
				|| currentState.Reads[ReadIndex].UP_Far.empty())
			continue;
		if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
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
								if (currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
										== currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
												+ 1
										&& currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
												< currentState.Reads[ReadIndex].ReadLength) {

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
											= currentState.Reads[ReadIndex].ReadLengthMinus
													- (currentState.Reads[ReadIndex].Right
															- currentState.Reads[ReadIndex].Left);
									currentState.Reads[ReadIndex].InsertedStr
											= ReverseComplement(
													currentState.Reads[ReadIndex]. UnmatchedSeq). substr(
													currentState.Reads[ReadIndex].BP
															+ 1,
													currentState.Reads[ReadIndex]. IndelSize);
									currentState.Reads[ReadIndex].NT_str = "";
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
												currentState.startOfRegion,
												currentState.endOfRegion)) {
											SIs[(int) currentState.Reads[ReadIndex]. BPLeft
													/ BoxSize]. push_back(
													ReadIndex);
											currentState.Reads[ReadIndex].Used
													= true;
											Count_SI_Plus++;
											Count_SI++;
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
								if (currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
										== currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
												+ 1
										&& currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
												+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
												< currentState.Reads[ReadIndex].ReadLength) {

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
											= currentState.Reads[ReadIndex].ReadLengthMinus
													- (currentState.Reads[ReadIndex].Right
															- currentState.Reads[ReadIndex].Left);
									currentState.Reads[ReadIndex].InsertedStr
											= currentState.Reads[ReadIndex].UnmatchedSeq. substr(
													currentState.Reads[ReadIndex].BP
															+ 1,
													currentState.Reads[ReadIndex]. IndelSize);
									currentState.Reads[ReadIndex].NT_str = "";
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
											SIs[(int) currentState.Reads[ReadIndex]. BPLeft
													/ BoxSize]. push_back(
													ReadIndex);
											currentState.Reads[ReadIndex].Used
													= true;
											Count_SI++;
											Count_SI_Minus++;
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
	std::cout << "Total: " << Count_SI << "\t+" << Count_SI_Plus << "\t-"
			<< Count_SI_Minus << std::endl;
	std::ofstream SIoutputfile(currentState.SIOutputFilename.c_str(),
			std::ios::app);
	SortOutputSI(NumBoxes, currentState.CurrentChr, currentState.Reads, SIs,
			SIoutputfile);
	SIoutputfile.close();
	for (unsigned int i = 0; i < NumBoxes; i++)
		SIs[i].clear();

	return EXIT_SUCCESS;
}
