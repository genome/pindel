#ifndef SEARCHVARIANT_H
#define SEARCHVARIANT_H

#include <string>
#include <vector>
#include "ControlState.h"

class SearchVariant {
public:
	SearchVariant();
	~SearchVariant();

	int Search(ControlState& currentState, const unsigned numBoxes);

protected:
	int Count_Var;
	int Count_Var_Plus;
	int Count_Var_Minus;

	std::string typeOfVariant;

	virtual bool decisionBranch1(ControlState& currentState, unsigned ReadIndex, unsigned int CloseIndex, int FarIndex) = 0;

	virtual bool decisionBranch2(ControlState& currentState, unsigned ReadIndex, unsigned int CloseIndex, int FarIndex) = 0;

	virtual unsigned int calculateIndelSize(ControlState& currentState, unsigned ReadIndex) = 0;

	virtual std::string getInsertedStr1(ControlState& currentState, unsigned ReadIndex) = 0;
	virtual std::string getInsertedStr2(ControlState& currentState, unsigned ReadIndex) = 0;

	virtual void outputResults(ControlState& currentState, std::vector<unsigned> Vars[], const unsigned NumBoxes) = 0;

private:
	SearchVariant(const SearchVariant&);

};

#endif // SEARCHVARIANT_H
