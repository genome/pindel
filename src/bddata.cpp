#include <algorithm>
#include <fstream>
#include "bddata.h"

BDData::BDData()
{
   m_breakDancerMask = NULL;
   //m_currentChrName = "";
}

bool isNumber(const std::string & input) {
    for (unsigned int index = 0; index < input.size(); index++) {
        switch (input[index]) {
            case '0':
                break;
            case '1':
                break;
            case '2':
                break;
            case '3':
                break;
            case '4':
                break;
            case '5':
                break;
            case '6':
                break;
            case '7':
                break;
            case '8':
                break;
            case '9':
                break;
            default:
                return false;
                //break;
        }
    }
    return true;
}


bool AtLeast6Fields(const std::string & input) {
    unsigned NumSpaceFields = 0;
    bool InSpace = false;
    unsigned InputSize = input.size();
    if (input[0] == '\t' || input[0] == ' ') return false; // must start with non-space char.
    for (unsigned index = 1; index < InputSize; index++) {
        if (input[index] == '\t' || input[index] == ' ') { // if current char is space, set InSpace to true
            InSpace = true;
        }
        else if (InSpace == true) { // not a space, NumSpaceFields++ if previous char is a space
            NumSpaceFields++;
            InSpace = false;
        }
    }
    if (NumSpaceFields >= 5) return true;
    else return false;
}


short CheckBreakDancerFileFormat(const std::string& filename) {
    std::ifstream CheckbdFileFirst( filename.c_str() );
    std::ifstream CheckbdFileSecond( filename.c_str() );
    char firstChar = 'k';
    std::string tempLine, errorLine;
    std::string Pos1, Pos2;
    
    while (!CheckbdFileFirst.eof()) {
        CheckbdFileFirst >> firstChar;
        //CheckbdFileFirst >> tempLine;
        //std::cout << "firstChar:" << firstChar << ":" << tempLine << std::endl;
        if (firstChar == '#') {
            std::getline(CheckbdFileFirst, tempLine);
            std::getline(CheckbdFileSecond, tempLine);
            //std::cout << "#" << tempLine << std::endl;
        }
        else {
            std::getline(CheckbdFileFirst, errorLine);
            //std::cout << "else " << errorLine << std::endl;
            if (AtLeast6Fields(errorLine)) {
                CheckbdFileSecond >> tempLine >> Pos1 >> tempLine
                >> tempLine >> Pos2 >> tempLine;
                std::getline(CheckbdFileSecond, tempLine);  // get rest of line
                if (isNumber(Pos1) && isNumber(Pos2)) {
                    
                }
                else {
                    std::cout << "something is wrong with this line \"" << firstChar << errorLine << "\"" << std::endl;
                    return 1;
                }
            }
            else {
                std::cout << "something is wrong with this line \"" << firstChar << errorLine << "\"" << std::endl;
                return 1;
            }
        }
    }
    return 0;
}





void BDData::loadBDFile(const std::string& filename)
{
	if (CheckBreakDancerFileFormat(filename) == 1) {
        std::cout << "\nIgnore breakdancer file due to an error in the BreakDancer file format.\n" << std::endl;
        return; // if lines not start with #, there must be at least 6 fields and NO. 2 and No. 5 must be int.
    }
   std::ifstream bdFile( filename.c_str() );
   if (!bdFile.good()) {
      std::cout << "Error: cannot load breakdancer file '" << filename << std::endl;
      exit( EXIT_FAILURE );
   }

   char firstChar;
   std::string tempLine;

   while (!bdFile.eof()) {
      firstChar= bdFile.peek();
      if (firstChar == '#') {
         std::getline(bdFile, tempLine);
      }
      else {
			std::string firstChrName, secondChrName;
		   std::string tempStringItem;
			unsigned int firstPos, secondPos;

         bdFile >> firstChrName >> firstPos >> tempStringItem
                >> secondChrName >> secondPos >> tempStringItem;
         std::getline(bdFile, tempLine);  // get rest of line

         firstPos += g_SpacerBeforeAfter; // ??? ask Kai
         secondPos += g_SpacerBeforeAfter;
			BreakDancerCoordinate firstBDCoordinate( firstChrName, firstPos );
			BreakDancerCoordinate secondBDCoordinate( secondChrName, secondPos ); 
			
         m_bdEvents.push_back(BreakDancerEvent( firstBDCoordinate, secondBDCoordinate ));
         m_bdEvents.push_back(BreakDancerEvent( secondBDCoordinate, firstBDCoordinate ));			
      }
   }
	sort( m_bdEvents.begin(), m_bdEvents.end(), sortOnFirstBDCoordinate ); 
   std::cout << "BreakDancer events: " << m_bdEvents.size()/2 << std::endl;
}


bool regionsOverlap( const BDIterator& firstRegionIt, const BDIterator& secondRegionIt )
{
	return ( ( secondRegionIt->second.chromosomeName == firstRegionIt->second.chromosomeName )
		&& ( secondRegionIt->second.startOfWindow() <= firstRegionIt->second.endOfWindow() + 1 ) );
}

void BDData::createRegionCluster(const BDIterator& startOfEventList, const BDIterator& endOfEventList, SearchWindowCluster& newCluster)
{
	std::vector<BreakDancerEvent> relevantSubcluster( startOfEventList, endOfEventList );
	sort( relevantSubcluster.begin(), relevantSubcluster.end(), sortOnSecondBDCoordinate );

	newCluster.clear();
	for (BDIterator eventIter=relevantSubcluster.begin(); eventIter!=relevantSubcluster.end(); eventIter++ ) {
		// NOTE: below code will be removed once we start working on interchromosomal translocations
		if (eventIter->second.chromosomeName != m_currentWindow.getChromosomeName() ) {
			std::cout << "Possible translocation from chromosome " << m_currentWindow.getChromosomeName() << " to chromosome " << eventIter->second.chromosomeName << "\n";
		}
		else {
			SearchWindow currentEventWindow( eventIter->second.chromosomeName, eventIter->second.startOfWindow(), eventIter->second.endOfWindow() );
			while ( eventIter+1!=relevantSubcluster.end() && regionsOverlap( eventIter, eventIter+1 ) ) {
				eventIter++;
				currentEventWindow.setEnd( eventIter->second.endOfWindow() );
			}
			newCluster.push_back( currentEventWindow );
		}
	}
}


void BDData::loadRegion( const SearchWindow& searchWindow  )
{
	// check how we handle borders/reads near borders at the moment
	
	m_currentWindow = searchWindow;
	const int INSERT_SIZE = 1000; // will need to adjust this!
	if (m_currentWindow.getStart() >= 3*INSERT_SIZE) {
		m_currentWindow.setStart( m_currentWindow.getStart() - 3*INSERT_SIZE );
	}
	else {
		m_currentWindow.setStart( 0 ); 
	}
	m_currentWindow.setEnd( m_currentWindow.getEnd() + 3*INSERT_SIZE );
	//std::cout << "Starting LoadRegion\n";
	BreakDancerCoordinate emptyBDCoord("",0);
	BreakDancerCoordinate startOfWindowBDCoord( m_currentWindow.getChromosomeName(), m_currentWindow.getStart() );
	BDIterator startRegionInBDEvents = lower_bound( m_bdEvents.begin(), m_bdEvents.end(), BreakDancerEvent(startOfWindowBDCoord, emptyBDCoord), sortOnFirstBDCoordinate );
	BreakDancerCoordinate endOfWindowBDCoord( m_currentWindow.getChromosomeName(), m_currentWindow.getEnd() );
	BDIterator endRegionInBDEvents = upper_bound( m_bdEvents.begin(), m_bdEvents.end(), BreakDancerEvent(endOfWindowBDCoord, emptyBDCoord), sortOnFirstBDCoordinate ); 
   
	delete[] m_breakDancerMask; // removing NULL is also safe
   m_breakDancerMask = new unsigned int[ m_currentWindow.getSize() ];

	SearchWindowCluster emptyCluster;	
	m_regionsToScanCollection.clear();
	m_regionsToScanCollection.push_back( emptyCluster );
	//std::cout << "Continuing LoadRegion\n";		
	BDIterator startOfEventList = startRegionInBDEvents;
	BDIterator endOfEventList = startRegionInBDEvents;
	int index = 0;
	
	for ( unsigned int position=m_currentWindow.getStart(); position< m_currentWindow.getEnd(); position++ ) {
		bool changed = false;
	//std::cout << "At position " << position << "\n";	

		// remove events that have been passed already
		for ( BDIterator eventIter=startOfEventList; eventIter<endOfEventList; eventIter++ ) {
			if ( position > eventIter->first.endOfWindow() ) {
				startOfEventList++;
				changed = true;
			}
			else { // <= to end position of event
				break;
			}
		}
	
		// add new events
		for (BDIterator eventIter=endOfEventList; eventIter<endRegionInBDEvents; eventIter++ ) {
			if ( position < eventIter->first.startOfWindow() ) {
				break;
			}
			else if ( position == eventIter->first.startOfWindow() ) {
				endOfEventList++;
				changed = true;
			}
		}
		// so if there are no events suitable, save 0 as the value
		if ( startOfEventList == endOfEventList ) {
			m_breakDancerMask[ position-m_currentWindow.getStart() ] = 0;
		}
		else {
			if (changed) {
				index++;
				SearchWindowCluster newCluster;			
				createRegionCluster( startOfEventList, endOfEventList, newCluster);
				m_regionsToScanCollection.push_back( newCluster );
			}
			m_breakDancerMask[ position-m_currentWindow.getStart() ] = index;
			//if (index%50==0) { std::cout << "Mask making: position " << position << " has index " << index << "\n"; }
		}	
	}
	
  /* if (m_currentChrName != chromosomeName ) { // so don't overwrite current chromosome data if not necessary
      m_currentChrName = chromosomeName;

      // create a new array to house the breakdancer data pertaining to the current chromosome
      delete[] m_breakDancerMask; // removing NULL is also safe

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
   }*/
}


const SearchWindowCluster& BDData::getCorrespondingSearchWindowCluster( const SPLIT_READ& read ) const
{
	if ( read.getLastAbsLocCloseEnd() - m_currentWindow.getStart() < 0 || read.getLastAbsLocCloseEnd() - m_currentWindow.getStart() > m_currentWindow.getSize() ) {
		std::cout << "Coordinate of the close end " << read.getLastAbsLocCloseEnd() << " is bad.\n";
		exit( EXIT_FAILURE );		
	}
	unsigned int clusterIndex = m_breakDancerMask[ read.getLastAbsLocCloseEnd() - m_currentWindow.getStart() ];
	return m_regionsToScanCollection[ clusterIndex ];
   /*if (m_breakDancerMask == NULL || m_currentChrName != read.FragName ) {
       //std::cout <<  m_currentChrName << " " << read.FragName  << std::endl;
      std::cout << "Error: chromosome data has not been set with BDData::setChromosome\n";
      exit(EXIT_FAILURE);
   }
   unsigned int firstMappingPoint = read.getLastAbsLocCloseEnd();
   int bdIndex = m_breakDancerMask[ firstMappingPoint ];
   if ( bdIndex != 0 ) { // there is an overlapping BD event
      unsigned int secondMappingPoint = 0;
      if (bdIndex < 0 ) {
         secondMappingPoint = m_bdEvents[ -bdIndex ].POS1;
      }
      else {
         secondMappingPoint = m_bdEvents[ bdIndex ].POS2;
      }
      complementaryPositions.push_back( secondMappingPoint );
   }*/

}

bool haveCommonBDEvent( const SearchWindowCluster& event, const unsigned int position, const std::string& chromosomeName )
{
	for (unsigned int index=0; index<event.size(); index++) {
		if (event[ index ].encompasses( chromosomeName, position )) {
			return true;	
		}
	}
	return false;
} 

/* 'isBreakdancerEvent' returns whether an event between leftPosition and rightPosition in the current chromosome
	is confirmed by a breakdancer result. */
bool BDData::isBreakDancerEvent( const unsigned int leftPosition, const unsigned int rightPosition ) const
{
   unsigned int rawLeftPosition = leftPosition + g_SpacerBeforeAfter;
   unsigned int rawRightPosition = rightPosition + g_SpacerBeforeAfter;

   if ( m_breakDancerMask[ rawLeftPosition - m_currentWindow.getStart() ]!=0  &&
		  m_breakDancerMask[ rawRightPosition - m_currentWindow.getStart() ]!=0 ) {
        
		return haveCommonBDEvent( m_regionsToScanCollection[ m_breakDancerMask[ rawLeftPosition - m_currentWindow.getStart() ] ],
				rawRightPosition - m_currentWindow.getStart(), m_currentWindow.getChromosomeName() );
   }
   else {
      return false;
   }
}
