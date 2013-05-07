/* BDData.h

	Contains all data necessary to store the information from breakdancer events, and check of a read whether it corresponds to one or more events

	Eric-Wubbo Lameijer, section of Molecular Epidemiology, LUMC, 10-08-2011, e.m.w.lameijer@lumc.nl
*/

#ifndef BDDATA_H_
#define BDDATA_H_

#include <string>
#include <vector>

#include "control_state.h"

typedef std::vector< SearchWindow > SearchWindowCluster;
typedef std::vector< BreakDancerEvent>::iterator BDIterator;

class BDData {

public:
	BDData();
	void createRegionCluster(const BDIterator& startOfEventList, const BDIterator& endOfEventList, SearchWindowCluster& newCluster);
	void loadBDFile(const std::string& filename);
	void loadRegion( const SearchWindow& searchWindow );
	void UpdateBD(ControlState & currentState);
	bool isBreakDancerEvent( const unsigned int leftPosition, const unsigned int rightPosition ) const;
	// returns positions belonging to the the complementary breakdancer calls
	const SearchWindowCluster& getCorrespondingSearchWindowCluster( const SPLIT_READ& read ) const;

	unsigned GetBDSize_external();
    
	unsigned GetBDSize_total();

private: 
	unsigned int *m_breakDancerMask;
	std::vector<BreakDancerEvent> m_bdEvents, m_bdEvents_external; 

	std::vector<SearchWindowCluster> m_regionsToScanCollection;
	//std::string m_currentChrName;
	SearchWindow m_currentWindow;
};

bool Compare2Str (const std::string & first, const std::string & second);

#endif /* BDDATA_H_ */
