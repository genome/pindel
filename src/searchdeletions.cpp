#include "searchdeletions.h"
#include "pindel.h"
#include "reporter.h"

SearchDeletions::SearchDeletions() {
	typeOfVariant = "deletions";
}

SearchDeletions::~SearchDeletions() {

}

bool SearchDeletions::decisionBranch1(ControlState& currentState,
		unsigned ReadIndex, unsigned int CloseIndex, int FarIndex) {
	return currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
			+ currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
			== currentState.Reads[ReadIndex].ReadLength
			&& currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
					> currentState.Reads[ReadIndex].UP_Close[CloseIndex]. AbsLoc
							+ 1;
}

bool SearchDeletions::decisionBranch2(ControlState& currentState,
		unsigned ReadIndex, unsigned int CloseIndex, int FarIndex) {
	return currentState.Reads[ReadIndex].UP_Close[CloseIndex]. LengthStr
			+ currentState.Reads[ReadIndex].UP_Far[FarIndex]. LengthStr
			== currentState.Reads[ReadIndex].ReadLength
			&& currentState.Reads[ReadIndex]. UP_Close[CloseIndex].AbsLoc
					> currentState.Reads[ReadIndex].UP_Far[FarIndex]. AbsLoc
							+ 1;
}

unsigned int SearchDeletions::calculateIndelSize(ControlState& currentState,
		unsigned ReadIndex) {
	return (currentState.Reads[ReadIndex].Right
			- currentState.Reads[ReadIndex].Left)
			- currentState.Reads[ReadIndex].ReadLengthMinus;
}

std::string SearchDeletions::getInsertedStr1(ControlState& currentState,
		unsigned ReadIndex) {
	return "";
}

std::string SearchDeletions::getInsertedStr2(ControlState& currentState,
		unsigned ReadIndex) {
	return "";
}

void SearchDeletions::outputResults(ControlState& currentState,
		std::vector<unsigned> Vars[], const unsigned NumBoxes) {
	std::ofstream DeletionOutf(currentState.DeletionOutputFilename.c_str(),
			std::ios::app);
	SortOutputD(NumBoxes, currentState.CurrentChr, currentState.Reads, Vars,
			DeletionOutf);
	DeletionOutf.close();
}
