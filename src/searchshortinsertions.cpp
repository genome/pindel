#include "searchshortinsertions.h"
#include "pindel.h"
#include "reporter.h"

SearchShortInsertion::SearchShortInsertion() {
	typeOfVariant = "short insertions";
}

SearchShortInsertion::~SearchShortInsertion() {

}

bool SearchShortInsertion::decisionBranch1(ControlState& currentState, unsigned ReadIndex, unsigned int CloseIndex, int FarIndex) {
	return
		currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
			== currentState.Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1
		&& currentState.Reads[ReadIndex].UP_Close[CloseIndex].LengthStr
			+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
			< currentState.Reads[ReadIndex].ReadLength;
}

bool SearchShortInsertion::decisionBranch2(ControlState& currentState, unsigned ReadIndex, unsigned int CloseIndex, int FarIndex) {
	return
		currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
			== currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc + 1
		&& currentState.Reads[ReadIndex]. UP_Far[FarIndex].LengthStr
			+ currentState.Reads[ReadIndex]. UP_Close[CloseIndex].LengthStr
			< currentState.Reads[ReadIndex].ReadLength;
}

unsigned int SearchShortInsertion::calculateIndelSize(ControlState& currentState, unsigned ReadIndex) {
	return currentState.Reads[ReadIndex].ReadLengthMinus
			- (currentState.Reads[ReadIndex].Right - currentState.Reads[ReadIndex].Left);
}

std::string SearchShortInsertion::getInsertedStr1(ControlState& currentState, unsigned ReadIndex) {
	return ReverseComplement(
		currentState.Reads[ReadIndex]. UnmatchedSeq). substr(
		currentState.Reads[ReadIndex].BP
				+ 1,
		currentState.Reads[ReadIndex]. IndelSize);
}

std::string SearchShortInsertion::getInsertedStr2(ControlState& currentState, unsigned ReadIndex) {
	return currentState.Reads[ReadIndex].UnmatchedSeq. substr(
		currentState.Reads[ReadIndex].BP + 1,
		currentState.Reads[ReadIndex]. IndelSize);
}

void SearchShortInsertion::outputResults(ControlState& currentState, std::vector<unsigned> Vars[], const unsigned NumBoxes) {
	std::ofstream SIoutputfile(currentState.SIOutputFilename.c_str(), std::ios::app);
	SortOutputSI(NumBoxes, currentState.CurrentChr, currentState.Reads, Vars, SIoutputfile);
	SIoutputfile.close();
}
