/* BDData.h

	Contains all data necessary to store the information from breakdancer events, and check of a read whether it corresponds to one or more events

	Eric-Wubbo Lameijer, section of Molecular Epidemiology, LUMC, 10-08-2011, e.m.w.lameijer@lumc.nl
*/

#ifndef BDDATA_H_
#define BDDATA_H_

#include <string>
#include <vector>

#include "control_state.h"

class BDData {

public:
	BDData();
	void loadBDFile(const std::string& filename);
	void loadChromosome( const std::string& chromosomeName, const unsigned int chromosomeSize );

	// returns positions belonging to the the complementary breakdancer calls
	void getCorrespondingEvents( const SPLIT_READ& read, std::vector<unsigned int>& complementaryPositions ) const; 

private: 
   int *m_breakDancerMask;
	std::vector<BreakDancer> m_bdEvents; 
	std::string m_currentChrName;
};

#endif /* BDDATA_H_ */
