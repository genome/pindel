#include <fstream>
#include "bddata.h"	



BDData::BDData()
{
	m_breakDancerMask = NULL;
	m_currentChrName = "";
}


void BDData::loadBDFile(const std::string& filename)
{
	std::ifstream bdFile( filename.c_str() );
	if (!bdFile.good()) {
		std::cout << "Error: cannot load breakdancer file '" << filename << std::endl;
		exit( EXIT_FAILURE );
	}
	char firstChar;

	std::string tempLine;
	std::string tempStringItem;
	BreakDancer temp_BD_event;

	m_bdEvents.push_back(temp_BD_event); // empty event since we work with +EventIndex and -EventIndex to identifty starts and ends of bdEvent with 
												// index EventIndex; and that wouldn't work for 0 (+0==-0)

	while (!bdFile.eof()) {
		firstChar= bdFile.peek();
		if (firstChar == '#') {
			std::getline(bdFile, tempLine);
		} else {
         bdFile >> temp_BD_event.ChrName_A >> temp_BD_event.POS1 >> tempStringItem 
                >> temp_BD_event.ChrName_B >> temp_BD_event.POS2 >> tempStringItem;
         std::getline(bdFile, tempLine);  // get rest of line

         temp_BD_event.POS1 += g_SpacerBeforeAfter;
         temp_BD_event.POS2 += g_SpacerBeforeAfter;

			m_bdEvents.push_back(temp_BD_event);
		}
	}
	std::cout << "BreakDancer events: " << m_bdEvents.size() - 1 << std::endl;
}


void BDData::loadChromosome( const std::string& chromosomeName, const unsigned int chromosomeSize )
{
	if (m_currentChrName != chromosomeName ) { // so don't overwrite current chromosome data if not necessary
	
		m_currentChrName = chromosomeName;

	  	// create a new array to house the breakdancer data pertaining to the current chromosome
		if (m_breakDancerMask != NULL) {
			delete[] m_breakDancerMask;
		}

   	m_breakDancerMask = new int[ chromosomeSize ];
		
   
   	// set all positions to 0 (=not covered by any breakdancer reads) by default
   	for (unsigned i = 0; i < chromosomeSize; i++) {
			m_breakDancerMask[i] = 0;
   	}
	
		// give positions which are covered by a breakdancer read the index of the breakdancer read.

	   for (unsigned i = 1; i < m_bdEvents.size(); i++) {
			if ( m_bdEvents[i].ChrName_A == chromosomeName && m_bdEvents[i].ChrName_B == chromosomeName ) {

				// sets reads around POS1 of the breakdancer item to +i		
				unsigned int Start = m_bdEvents[i].POS1 - 200;
				unsigned int End = m_bdEvents[i].POS1 + 200;
	   	   for (unsigned j = Start; j < End; j++) {
	  				m_breakDancerMask[j] = i;
				}

				// sets reads around POS2 of the breakdancer item to -i	
				Start = m_bdEvents[i].POS2 - 200;
				End = m_bdEvents[i].POS2 + 200;
				for (unsigned j = Start; j < End; j++) {
					m_breakDancerMask[j] = i * (-1);
				} 
			}
		}
	}
}


void BDData::getCorrespondingEvents( const SPLIT_READ& read, std::vector<unsigned int>& complementaryPositions ) const
{
	if (m_breakDancerMask == NULL || m_currentChrName != read.FragName ) {
		std::cout << "Error: chromosome data has not been set with BDData::setChromosome\n";
		exit(EXIT_FAILURE);
	}
	unsigned int firstMappingPoint = read.getLastAbsLocCloseEnd();
	int bdIndex = m_breakDancerMask[ firstMappingPoint ];
	if ( bdIndex != 0 ) { // there is an overlapping BD event
		unsigned int secondMappingPoint = 0;
		if (bdIndex < 0 ) {
			secondMappingPoint = m_bdEvents[ -bdIndex ].POS1; 
		} else {
			secondMappingPoint = m_bdEvents[ bdIndex ].POS2;
		}		
		complementaryPositions.push_back( secondMappingPoint );
	}

}

/* 'isBreakdancerEvent' returns whether an event between leftPosition and rightPosition in the current chromosome 
	is confirmed by a breakdancer result. */
bool BDData::isBreakDancerEvent( const unsigned int leftPosition, const unsigned int rightPosition ) const
{
   unsigned int rawLeftPosition = leftPosition + g_SpacerBeforeAfter;
   unsigned int rawRightPosition = rightPosition + g_SpacerBeforeAfter;

	if ( m_breakDancerMask[ rawLeftPosition ]!=0  && 
		 m_breakDancerMask[ rawLeftPosition ] == -m_breakDancerMask[ rawRightPosition ] ) {
		return true;
	} else {
		return false;		
	}
}
