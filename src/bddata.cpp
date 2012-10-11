#include <algorithm>
#include <fstream>
#include "bddata.h"

BDData::BDData() : 	m_currentWindow( &g_dummyChromosome, 0, 0 )
{
   m_breakDancerMask = NULL;

   //m_currentChrName = "";
}

bool isNumber(const std::string & input) 
{
   for (unsigned int index = 0; index < input.size(); index++) {
		if (input[ index ]<'0' || input[ index ]>'9') return false;
	}	
 /*       switch (input[index]) {
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
    }*/
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
        firstChar = CheckbdFileFirst.peek();// >> firstChar;
        if (firstChar == '#') {
            safeGetline(CheckbdFileFirst, tempLine);
            
            safeGetline(CheckbdFileSecond, tempLine);
        }
        else {
            safeGetline(CheckbdFileFirst, errorLine);
            //safeGetline(CheckbdFileSecond, errorLine);
            //std::cout << "errorLine/" << errorLine << "/" << std::endl;
            
            
            if (AtLeast6Fields(errorLine)) {
                //std::cout << "at least 6 fields" << std::endl;
                CheckbdFileSecond >> tempLine >> Pos1 >> tempLine
                >> tempLine >> Pos2 >> tempLine;
                safeGetline(CheckbdFileSecond, tempLine);  // get rest of line
                //std::cout << errorLine << " " << Pos1 << " " << Pos2 << std::endl;
                if (isNumber(Pos1) && isNumber(Pos2)) {
                    
                }
                else {
                    std::cout << "1 something is wrong with this line \"" << firstChar << errorLine << "\"" << std::endl;
                    return 1;
                }
            }
            else if (errorLine != "") {
                std::cout << "2 something is wrong with this line \"" << firstChar << errorLine << "\"" << std::endl;
                return 1;
            }
        }
    }
    return 0;
}


void BDData::loadBDFile(const std::string& filename)
{
    //std::cout << " BDData::loadBDFile " << std::endl;
	if (CheckBreakDancerFileFormat(filename) == 1) {
        //std::cout << "\nIgnore breakdancer file due to an error in the BreakDancer file format.\n" << std::endl;
        return; // if lines not start with #, there must be at least 6 fields and NO. 2 and No. 5 must be int.
    }
   std::ifstream bdFile( filename.c_str() );
   if (!bdFile.good()) {
      //std::cout << "Error: cannot load breakdancer file '" << filename << std::endl;
      exit( EXIT_FAILURE );
   }

   char firstChar;
   std::string tempLine;

   while (!bdFile.eof()) {
      firstChar= bdFile.peek();
      if (firstChar == '#') {
         safeGetline(bdFile, tempLine);
      }
      else {
			std::string firstChrName, secondChrName;
		   std::string tempStringItem;
			unsigned int firstPos, secondPos;

         bdFile >> firstChrName >> firstPos >> tempStringItem
                >> secondChrName >> secondPos >> tempStringItem;
          //std::cout << firstChrName << " " << firstPos << " " << tempStringItem << " " << secondChrName << " " << secondPos << std::endl;
         safeGetline(bdFile, tempLine);  // get rest of line

         firstPos += g_SpacerBeforeAfter; // ??? ask Kai
         secondPos += g_SpacerBeforeAfter;
			if ( firstChrName!="" && secondChrName!="" ) {
				BreakDancerCoordinate firstBDCoordinate( firstChrName, firstPos );
				BreakDancerCoordinate secondBDCoordinate( secondChrName, secondPos ); 
			
   	      m_bdEvents.push_back(BreakDancerEvent( firstBDCoordinate, secondBDCoordinate ));
   	      m_bdEvents.push_back(BreakDancerEvent( secondBDCoordinate, firstBDCoordinate ));			
			}
      }
   }
	sort( m_bdEvents.begin(), m_bdEvents.end(), sortOnFirstBDCoordinate ); 
   std::cout << "BreakDancer events: " << m_bdEvents.size()/2 << std::endl;
}

bool Compare2Str (const std::string & first, const std::string & second) {
    if (first == second) return false;
    unsigned Length_first = first.size();
    unsigned Length_second = second.size();
    unsigned Length = std::min(Length_first, Length_second);
    for (unsigned index = 0; index < Length; index++) {
        if (first[index] > second[index]) return true;
        if (first[index] < second[index]) return false;
    }
    if (Length_first < Length_second) return true;
    else if (Length_first >= Length_second) return false;
    return false;
}

void SortRPByChrPos(const std::vector <RP_READ> & Reads_RP, std::vector <unsigned> & RP_index) {
    unsigned DistanceCutoff = 1000;
    unsigned Num_RP_Reads = Reads_RP.size();
    for (unsigned index = 0; index < Num_RP_Reads; index++) {
        RP_index.push_back(index);
    }
    unsigned TempIndex;
    bool Exchange;
    for (unsigned first = 0; first < Num_RP_Reads - 1; first++) {
        for (unsigned second = first + 1; second < Num_RP_Reads; second++) {
            if (Reads_RP[RP_index[second]].PosA > Reads_RP[RP_index[first]].PosA + DistanceCutoff
                && Reads_RP[RP_index[second]].PosA > Reads_RP[RP_index[first]].PosB + DistanceCutoff
                && Reads_RP[RP_index[second]].PosB > Reads_RP[RP_index[first]].PosA + DistanceCutoff
                && Reads_RP[RP_index[second]].PosB > Reads_RP[RP_index[first]].PosB + DistanceCutoff) break;
            Exchange = false;
            if (Compare2Str(Reads_RP[RP_index[first]].ChrNameA, Reads_RP[RP_index[second]].ChrNameA)) Exchange = true;
            else if (Reads_RP[RP_index[first]].ChrNameA == Reads_RP[RP_index[second]].ChrNameA
                     && Reads_RP[RP_index[first]].PosA > Reads_RP[RP_index[second]].PosA) Exchange = true;
            else if (Reads_RP[RP_index[first]].ChrNameA == Reads_RP[RP_index[second]].ChrNameA
                     && Reads_RP[RP_index[first]].PosA == Reads_RP[RP_index[second]].PosA
                     && Compare2Str(Reads_RP[RP_index[first]].ChrNameB, Reads_RP[RP_index[second]].ChrNameB)) Exchange = true;
            else if (Reads_RP[RP_index[first]].ChrNameA == Reads_RP[RP_index[second]].ChrNameA
                     && Reads_RP[RP_index[first]].PosA == Reads_RP[RP_index[second]].PosA
                     && Reads_RP[RP_index[first]].ChrNameB == Reads_RP[RP_index[second]].ChrNameB
                     && Reads_RP[RP_index[first]].PosB > Reads_RP[RP_index[second]].PosB) Exchange = true;
            if (Exchange) {
                TempIndex = RP_index[first];
                RP_index[first] = RP_index[second];
                RP_index[second] = TempIndex;
            }
        }
    }
}

void ModifyRP(std::vector <RP_READ> & Reads_RP, std::vector <unsigned> & RP_index) {
    unsigned Range, PosA_start, PosA_end, PosB_start, PosB_end, Stop_Pos_A, Stop_Pos_B;
    int Start_Index_Second;
    for (unsigned first = 0; first < Reads_RP.size(); first++) {
        RP_READ & Current_first = Reads_RP[RP_index[first]];
        Range = (unsigned)(Current_first.InsertSize * 1.1);
        if (Current_first.PosA > Range)
        PosA_start = Current_first.PosA - Range;
        else PosA_start = 1;
        PosA_end = Current_first.PosA + Range;
        if (Current_first.PosB > Range)
        PosB_start = Current_first.PosB - Range;
        else PosB_start = 1;
        PosB_end = Current_first.PosB + Range;
        Stop_Pos_A = PosA_end + Range;
        Stop_Pos_B = PosB_end + Range;
        Start_Index_Second = (int)first - 200;
        if (Start_Index_Second < 0) Start_Index_Second = 0;
        //std::cout << "Before: " << Current_first.DA << " " << Current_first.PosA << " " << Current_first.DB << " " << Current_first.PosB << std::endl;
        for (unsigned second = Start_Index_Second; second < Reads_RP.size(); second++) {
            RP_READ & Current_second = Reads_RP[RP_index[second]];
            if (Current_second.PosA > Stop_Pos_A && Current_second.PosA > Stop_Pos_B && Current_second.PosB > Stop_Pos_A && Current_second.PosB > Stop_Pos_B ) break;
            if (Current_first.ChrNameA == Current_second.ChrNameA
                && Current_first.ChrNameB == Current_second.ChrNameB
                && Current_first.DA == Current_second.DA
                && Current_first.DB == Current_second.DB) {
                if (Current_second.PosA > PosA_start
                    && Current_second.PosA < PosA_end
                    && Current_second.PosB > PosB_start
                    && Current_second.PosB < PosB_end) { // same cluster
                    if ((Current_first.DA == '+' && Current_first.PosA < Current_second.PosA)
                        || (Current_first.DA == '-' && Current_first.PosA > Current_second.PosA)) {
                        Current_first.PosA = Current_second.PosA;
                    }
                    if ((Current_first.DB == '+' && Current_first.PosB < Current_second.PosB)
                        || (Current_first.DB == '-' && Current_first.PosB > Current_second.PosB)) {
                        Current_first.PosB = Current_second.PosB;
                    }
                }
            }
        }
        //std::cout << "After: " << Current_first.DA << " " << Current_first.PosA << " " << Current_first.DB << " " << Current_first.PosB << std::endl;
    }
    for (unsigned first = 0; first < Reads_RP.size(); first++) {
        if (Reads_RP[first].DA == '+') Reads_RP[first].PosA = Reads_RP[first].PosA + Reads_RP[first].ReadLength;
        if (Reads_RP[first].DB == '+') Reads_RP[first].PosB = Reads_RP[first].PosB + Reads_RP[first].ReadLength;
        //std::cout << "Final: " << Reads_RP[first].ChrNameA << " " << Reads_RP[first].DA << " " << Reads_RP[first].PosA << "\t" << Reads_RP[first].ChrNameB << " " << Reads_RP[first].DB << " " << Reads_RP[first].PosB << std::endl;
    }
}

void Summarize(std::vector <RP_READ> & Reads_RP) {
    unsigned Cutoff = 3;
	if (Reads_RP.size()==0) return;
    for (unsigned int first = 0; first < Reads_RP.size() - 1; first++) {
        if (Reads_RP[first].Visited == true) continue;
        for (unsigned second = first + 1; second < Reads_RP.size(); second++) {
            if (Reads_RP[second].Visited == true) continue;
            if (Reads_RP[first].ChrNameA == Reads_RP[second].ChrNameA
                && Reads_RP[first].ChrNameB == Reads_RP[second].ChrNameB
                && Reads_RP[first].PosA == Reads_RP[second].PosA
                && Reads_RP[first].PosB == Reads_RP[second].PosB
                && Reads_RP[first].DA == Reads_RP[second].DA
                && Reads_RP[first].DB == Reads_RP[second].DB) {
                Reads_RP[first].NumberOfIdentical++;
                Reads_RP[second].Visited = true;
            }
        }
        if (Reads_RP[first].NumberOfIdentical >= Cutoff) {
            Reads_RP[first].Report = true;
            //std::cout << Reads_RP[first].ChrNameA << " " << Reads_RP[first].PosA << " " << Reads_RP[first].ChrNameB << " " << Reads_RP[first].PosB << std::endl;
        }
    }
}

void BDData::UpdateBD(ControlState & currentState) {
    std::cout << "Discovery RP: " << currentState.Reads_RP_Discovery.size() << std::endl;
    std::vector <unsigned> RP_index;
    SortRPByChrPos(currentState.Reads_RP_Discovery, RP_index);
    std::cout << "sorting RP complete." << std::endl;
    ModifyRP(currentState.Reads_RP_Discovery, RP_index);
    std::cout << "Modify RP complete." << std::endl;
    Summarize(currentState.Reads_RP_Discovery);
    #pragma omp parallel default(shared)
    {
    #pragma omp for
        for (int read_index = 0; read_index < (int)currentState.Reads_RP_Discovery.size(); read_index++) {
            if (currentState.Reads_RP_Discovery[read_index].Report == false) continue;
    
            std::string firstChrName = currentState.Reads_RP_Discovery[read_index].ChrNameA;
            std::string secondChrName = currentState.Reads_RP_Discovery[read_index].ChrNameB;
            unsigned int firstPos = currentState.Reads_RP_Discovery[read_index].PosA + g_SpacerBeforeAfter;
            unsigned int secondPos  = currentState.Reads_RP_Discovery[read_index].PosB + g_SpacerBeforeAfter;
            if ( firstChrName!="" && secondChrName!="" ) {
                BreakDancerCoordinate firstBDCoordinate( firstChrName, firstPos );
                BreakDancerCoordinate secondBDCoordinate( secondChrName, secondPos );
                #pragma omp critical
                {
                    m_bdEvents.push_back(BreakDancerEvent( firstBDCoordinate, secondBDCoordinate ));
                    m_bdEvents.push_back(BreakDancerEvent( secondBDCoordinate, firstBDCoordinate ));
                }
                //std::cout << "adding " << firstChrName << " " << firstPos - g_SpacerBeforeAfter << "\t" << secondChrName << " " << secondPos - g_SpacerBeforeAfter <<  " to breakdancer events. " << abs((int)secondPos - (int)firstPos) << std::endl;
            }
        
        }
    }
    //currentState.Reads_RP_Discovery.clear();
    std::cout << "sorting BD..." << std::endl;
	sort( m_bdEvents.begin(), m_bdEvents.end(), sortOnFirstBDCoordinate );
    std::cout << "sorting BD... done." << std::endl;
}


bool regionsOverlap( const BDIterator& firstRegionIt, const BDIterator& secondRegionIt )
{
	return ( ( secondRegionIt->second.getChromosomeName() == firstRegionIt->second.getChromosomeName() )
		&& ( secondRegionIt->second.startOfWindow() <= firstRegionIt->second.endOfWindow() + 1 ) );
}

void BDData::createRegionCluster(const BDIterator& startOfEventList, const BDIterator& endOfEventList, SearchWindowCluster& newCluster)
{
	std::vector<BreakDancerEvent> relevantSubcluster( startOfEventList, endOfEventList );
	sort( relevantSubcluster.begin(), relevantSubcluster.end(), sortOnSecondBDCoordinate );
		//std::cout << "Making cluster" << "\n";
	newCluster.clear();
	for (BDIterator eventIter=relevantSubcluster.begin(); eventIter!=relevantSubcluster.end(); eventIter++ ) {
		//std::cout << "FC: " << eventIter->first.getChromosomeName() << "FS: " << eventIter->first.startOfWindow() << "FE: " << eventIter->first.endOfWindow() << 
		//	"SC: " << eventIter->second.getChromosomeName() << "SS: " << eventIter->second.startOfWindow()<< "SE: " << eventIter->second.endOfWindow() << "\n";
		// NOTE: below code will be removed once we start working on interchromosomal translocations
		/*if (eventIter->second.getChromosomeName() != m_currentWindow.getChromosomeName() ) {
			std::cout << "Possible translocation from chromosome " << m_currentWindow.getChromosomeName() << " to chromosome " << eventIter->second.getChromosomeName() << "\n";
		}
		else {*/

			SearchWindow currentEventWindow( eventIter->second.getChromosome(), eventIter->second.startOfWindow(), eventIter->second.endOfWindow() );
			while ( eventIter+1!=relevantSubcluster.end() && regionsOverlap( eventIter, eventIter+1 ) ) {
				eventIter++;
				currentEventWindow.setEnd( eventIter->second.endOfWindow() );
			}
			//std::cout << "Pushing back\n";
			newCluster.push_back( currentEventWindow );
		//}
	}
}


void BDData::loadRegion( const SearchWindow& searchWindow  )
{
	// check how we handle borders/reads near borders at the moment
	//std::cout << "Starting loadRegion...\n";
	m_currentWindow = searchWindow;
	//std::cout << "m_currentWindow start " << m_currentWindow.getStart() << " end " << m_currentWindow.getEnd() << "size" << m_currentWindow.getSize()<< "\n";
	const int INSERT_SIZE = 1000; // will need to adjust this!
	if (m_currentWindow.getStart() >= 3*INSERT_SIZE) {
		m_currentWindow.setStart( m_currentWindow.getStart() - 3*INSERT_SIZE );
	}
	else {
		m_currentWindow.setStart( 0 ); 
	}
	//std::cout << "chromosomePtr=" << m_currentWindow.getChromosome() << "\n";	
	//std::cout << "m_currentWindow start " << m_currentWindow.getChromosome()->getName() << " z " << m_currentWindow.getStart() << " end " << m_currentWindow.getEnd() << "size" << m_currentWindow.getSize()<< "\n";
	/*for (unsigned int ki=0; ki<m_bdEvents.size(); ki++ ) {
		std::cout << m_bdEvents[ki].first.getChromosomeName() <<" " << m_bdEvents[ki].first.position << " - " << m_bdEvents[ki].second.position << "-";
		std::cout << m_bdEvents[ki].second.getChromosomeName() <<" " << m_bdEvents[ki].second.position << " - " << m_bdEvents[ki].second.position << "\n";
	}*/
//std::cout << "X\n";
//std::cout << "m_currentWindow startA " << m_currentWindow.getChromosome()->getName() << "\n";
	m_currentWindow.setEnd( m_currentWindow.getEnd() + 3*INSERT_SIZE );
//std::cout << "m_currentWindow startB " << m_currentWindow.getChromosome()->getName() << "\n";
	BreakDancerCoordinate emptyBDCoord("",0);
//std::cout << "m_currentWindow startC " << m_currentWindow.getChromosome()->getName() << "\n";
//std::cout << "Y\n";
	BreakDancerCoordinate startOfWindowBDCoord( m_currentWindow.getChromosome(), m_currentWindow.getStart() );
//std::cout << "Z\n";
	BDIterator startRegionInBDEvents = lower_bound( m_bdEvents.begin(), m_bdEvents.end(), BreakDancerEvent(startOfWindowBDCoord, emptyBDCoord), sortOnFirstBDCoordinate );
//std::cout << "Z1\n";
	BreakDancerCoordinate endOfWindowBDCoord( m_currentWindow.getChromosome(), m_currentWindow.getEnd() );
//std::cout << "Last coordinate: " << endOfWindowBDCoord.position << std::endl;
	BreakDancerEvent endEvent(endOfWindowBDCoord, emptyBDCoord);
//std::cout << "Last coordinate(2): " << endEvent.first.position << std::endl;	
	BDIterator endRegionInBDEvents = upper_bound( m_bdEvents.begin(), m_bdEvents.end(), endEvent, sortOnFirstBDCoordinate ); 
//std::cout << "The last is at..." <<   endRegionInBDEvents-m_bdEvents.begin() << "\n";;
//std::cout<< "Investigating: " << endRegionInBDEvents->first.position << " - " << startRegionInBDEvents->first.position << "of" << m_bdEvents.end()-m_bdEvents.begin() << "\n"; 
	/*std::cout << "After selecting:\n";

	for (BDIterator it=startRegionInBDEvents; it!=endRegionInBDEvents; it++ ) {
		std::cout << it->first.getChromosomeName() <<" " << it->first.position << " - " << it->second.position << "\n";
	}*/


	delete[] m_breakDancerMask; // removing NULL is also safe
   m_breakDancerMask = new unsigned int[ m_currentWindow.getSize() ];

	SearchWindowCluster emptyCluster;	
	m_regionsToScanCollection.clear();
	m_regionsToScanCollection.push_back( emptyCluster );
	////std::cout << "Continuing LoadRegion\n";		
	BDIterator startOfEventList = startRegionInBDEvents;
	BDIterator endOfEventList = startRegionInBDEvents;
	int index = 0;
	
	for ( unsigned int position=m_currentWindow.getStart(); position< m_currentWindow.getEnd(); position++ ) {
		bool changed = false;
	////std::cout << "At position " << position << "\n";	

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
                //std::cout << "I should be making an event now!\n";
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
                                //std::cout << "I should be making a cluster now!\n";
				createRegionCluster( startOfEventList, endOfEventList, newCluster);
				m_regionsToScanCollection.push_back( newCluster );
                                //std::cout << "With " << newCluster.size() << "events\n";
			}
			m_breakDancerMask[ position-m_currentWindow.getStart() ] = index;
			//if (index%50==0) { //std::cout << "Mask making: position " << position << " has index " << index << "\n"; }
		}	
	}
}


const SearchWindowCluster& BDData::getCorrespondingSearchWindowCluster( const SPLIT_READ& read ) const
{
	if ( read.getLastAbsLocCloseEnd() - m_currentWindow.getStart() < 0 || read.getLastAbsLocCloseEnd() - m_currentWindow.getStart() > m_currentWindow.getSize() ) {
		//std::cout << "Coordinate of the close end " << read.getLastAbsLocCloseEnd() << " is bad.\n";
		//std::cout << "Start: " << m_currentWindow.getStart() << " end " << m_currentWindow.getEnd() << " size " << m_currentWindow.getSize() << "\n";
		exit( EXIT_FAILURE );		
	}
	unsigned int clusterIndex = m_breakDancerMask[ read.getLastAbsLocCloseEnd() - m_currentWindow.getStart() ];
	return m_regionsToScanCollection[ clusterIndex ];
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
		//std::cout << "Found breakdancerEvent near " << leftPosition << "-" << rightPosition << "\n";
        
		return haveCommonBDEvent( m_regionsToScanCollection[ m_breakDancerMask[ rawLeftPosition - m_currentWindow.getStart() ] ],
				rawRightPosition - m_currentWindow.getStart(), m_currentWindow.getChromosomeName() );
   }
   else {
      return false;
   }
}

unsigned BDData::GetBDSize(){return m_bdEvents.size();}
