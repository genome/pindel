#ifndef SEARCHDELETIONS_H
#define SEARCHDELETIONS_H

class SearchDeletions : public SearchVariant {
public:
	SearchDeletions();
	virtual ~SearchDeletions();

protected:
	int Count_Var;
	int Count_Var_Plus;
	int Count_Var_Minus;

	std::string typeOfVariant;

	bool decisionBranch1(ControlState& currentState, unsigned ReadIndex, unsigned int CloseIndex, int FarIndex);

	bool decisionBranch2(ControlState& currentState, unsigned ReadIndex, unsigned int CloseIndex, int FarIndex);

	unsigned int calculateIndelSize(ControlState& currentState, unsigned ReadIndex);

	std::string getInsertedStr1(ControlState& currentState, unsigned ReadIndex);

	std::string getInsertedStr2(ControlState& currentState, unsigned ReadIndex);

	void outputResults(ControlState& currentState, std::vector<unsigned> Vars[], const unsigned NumBoxes);

private:
	SearchDeletions(const SearchDeletions&);

};

#endif // SEARCHDELETIONS_H
