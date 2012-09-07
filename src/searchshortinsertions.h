#ifndef SEARCHSHORTINSERTIONS_H
#define SEARCHSHORTINSERTIONS_H

class SearchShortInsertions : public SearchVariant {
public:
	SearchShortInsertions();
	virtual ~SearchShortInsertions();

protected:
	int Count_Var;
	int Count_Var_Plus;
	int Count_Var_Minus;

	std::string typeOfVariant;

	bool decisionBranch1(const SPLIT_READ& read, unsigned int CloseIndex, int FarIndex);

	bool decisionBranch2(const SPLIT_READ& read, unsigned int CloseIndex, int FarIndex);

	unsigned int calculateIndelSize(const SPLIT_READ& read);

	std::string getInsertedStr1(const SPLIT_READ& read);

	std::string getInsertedStr2(const SPLIT_READ& read);

	void outputResults(BDData & g_bdData, ControlState& currentState, std::vector<unsigned> Vars[], const unsigned NumBoxes, const SearchWindow& currentWindow);

private:
	SearchShortInsertions(const SearchShortInsertions&);

};

#endif // SEARCHSHORTINSERTIONS_H
