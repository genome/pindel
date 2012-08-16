#ifndef SEARCH_TANDEM_DUPLICATIONS_H
#define	SEARCH_TANDEM_DUPLICATIONS_H

int searchTandemDuplications(ControlState& currentState, unsigned NumBoxes, const SearchWindow& currentWindow);
void LeftMostTD(ControlState& currentState, SPLIT_READ & currentRead, const SearchWindow& window);
#endif /* SEARCH_TANDEM_DUPLICATIONS_H */
