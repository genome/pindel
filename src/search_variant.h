#ifndef SEARCHVARIANT_H
#define SEARCHVARIANT_H

#include "bddata.h"


class SearchVariant {
public:
	SearchVariant();
	virtual ~SearchVariant();

	int Search(BDData & g_bdData, ControlState& currentState, const unsigned numBoxes, const SearchWindow& window);

protected:
	int Count_Var;
	int Count_Var_Plus;
	int Count_Var_Minus;

	std::string typeOfVariant;

	virtual bool decisionBranch1(const SPLIT_READ& read, unsigned int CloseIndex, int FarIndex) = 0;

	virtual bool decisionBranch2(const SPLIT_READ& read, unsigned int CloseIndex, int FarIndex) = 0;

	virtual unsigned int calculateIndelSize(const SPLIT_READ& read ) = 0;

	virtual std::string getInsertedStr1(const SPLIT_READ& read ) = 0;
	virtual std::string getInsertedStr2(const SPLIT_READ& read ) = 0;

	virtual void outputResults(BDData & g_bdData, ControlState& currentState, std::vector<unsigned> Vars[], const unsigned NumBoxes, const SearchWindow& currentWindow) = 0;

private:
	SearchVariant(const SearchVariant&);

};

#endif // SEARCHVARIANT_H
