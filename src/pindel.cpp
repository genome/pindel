/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


// System header files
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <omp.h>
#include <set>

// Pindel header files
#include "logstream.h"
#include "pindel.h"
#include "fn_parameters.h"
#include "bddata.h"
#include "reader.h"
#include "searcher.h"
#include "reporter.h"
#include "parameter.h"
#include "control_state.h"
#include "search_deletions_nt.h"
#include "search_inversions.h"
#include "search_inversions_nt.h"
#include "search_MEI.h"
#include "search_tandem_duplications.h"
#include "search_tandem_duplications_nt.h"
#include "read_buffer.h"
#include "farend_searcher.h"
#include "search_variant.h"
#include "searchshortinsertions.h"
#include "searchdeletions.h"
#include "user_defined_settings.h"
#include "logdef.h"
#include "assembly.h"
#include "genotyping.h"
#include "ifstream_line_reader.h"
#include "gz_line_reader.h"

/*v Kai Ye update 0.2.4h, Oct 31 2011, update for MOSAIK */
/*v EW update 0.2.4j, Pindel will now abort when insert size is set too small. */
/*v Kai/EW update 0.2.4k, pindel will now give the consensus inserted sequence instead of always the first one. */
/*v EW update 0.2.4l, -Werror removed for user convenience, also made some typecasts in reader.cpp explicit; also reset -l and -k options to true by default. */
/*v Kai update 0.2.4m, removed an error in invertion calling. */
/* EW update 0.2.4n, improved VCF creator: less memory, more speed, removed bug. */
/* EW/Kai/Matthijs: update 0.2.4o: does not report short inversions as deletions anymore, also displays the reference for short inversions correctly, slightly changed sensitivity, small memory leak fixed. */
/* Kai/EW: update 0.2.4p: sorts reads in output, number of unique calls should be correct now, pindel now gives error if using config file that does not exist or has other problems */
/* EW: update 0.2.4q: also works with -c all instead of -c ALL */
/* Kai: update 0.2.4q: use ChrName and ChrSeq for clarity; start to include assembly */
/* Kai: min_perfect_match_around_BP to 5 */
/* EW: update 0.2.4s: bugfix for -p option of Pindel0.2.4r */
/* EW: update 0.2.4t, updates now shown in RELEASE document in trunk directory */

const std::string Pindel_Version_str = "Pindel version 0.2.4t, August 13 2012.";

std::ofstream g_logFile;

int g_binIndex = -1; // global variable for the bin index, as I cannot easily pass an extra parameter to the diverse functions
unsigned int g_maxPos = 0; // to calculate which is the last position in the chromosome, and hence to calculate the number of bins
short g_MinClose = 8;

std::set<std::string> g_sampleNames;

short Before, After;

BDData g_bdData;

unsigned int CountIndels = 0;
const int alphs = 4;
const char alphabet[alphs] = { 'A', 'C', 'G', 'T' };

unsigned long long int TheMax = 0;
const short MAX_MISMATCHES = 4;
float ExtraDistanceRate = 0.1;

double Const_S = 0.0;
double LOG14 = log10(0.25);
unsigned int BoxSize = 10000; // 10k is fine
const double Min_Filter_Ratio = 0.5;
unsigned int SPACERSIZE = 1;
unsigned int OriginalNumRead = 0;
const std::string NonACGT = "$";
short MIN_Len_Match = 4;
unsigned int NumberOfSIsInstances = 0;
unsigned int NumberOfDeletionsInstances = 0;
unsigned int NumberOfDIInstances = 0;
unsigned int g_numberOfInvInstances = 0;
unsigned int NumberOfTDInstances = 0;
short g_reportLength = 1;
char Match[256];
char Match2N[256];
char Convert2RC[256];
char Convert2RC4N[256];
char Cap2LowArray[256];
bool FirstChr = true;
const double InsertSizeExtra = 2;
unsigned int CONS_Chr_Size;
unsigned int DSizeArray[15];
int g_maxInsertSize=0;

std::string BreakDancerMask;
std::string CurrentChrMask;
std::vector<Parameter *> parameters;

// #########################################################

int NumRead2ReportCutOff_BP = 2;
int MaxRangeIndex = 9; // 5 or 6 or 7 or maximum 8      //# // user
double MaximumAllowedMismatchRate = 0.1; //#  // user
int Max_Length_NT = 30; // user
const bool ReportSVReads = false;
const bool ReportLargeInterChrSVReads = false;
const unsigned short Indel_SV_cutoff = 50;
unsigned int WINDOW_SIZE = 10000000;
const int AROUND_REGION_BUFFER = 10000; // how much earlier reads should be selected if only a region of the chromosome needs be specified.
// #########################################################

unsigned int Distance = 300;
short MinFar_D = 8; //atoi(argv[3]);
const short MaxDI = 30;

UniquePoint::UniquePoint( const short lengthStr, const unsigned int absLoc, const char direction, const char strand, const short mismatches ) :
	LengthStr( lengthStr ), AbsLoc( absLoc ), Direction( direction ), Strand( strand ), Mismatches( mismatches )
{
}

bool SearchWindow::encompasses( const std::string& chromosomeName, const unsigned int position ) const
{
	return ( ( m_chromosomeName == chromosomeName ) && ( position >= m_currentStart ) && ( position <= m_currentEnd ) );
}

SearchWindow::SearchWindow(const std::string& chromosomeName, const int regionStart, const int regionEnd ) : m_chromosomeName( chromosomeName )
{
	m_currentStart = regionStart;
	m_currentEnd = regionEnd;
}

LoopingSearchWindow::LoopingSearchWindow(const SearchRegion* region, const int chromosomeSize, const int binSize, const std::string& chromosomeName ) : 
	SearchWindow(chromosomeName,0,chromosomeSize), m_BIN_SIZE( binSize )
{
	if (region->isStartDefined()) {
		m_officialStart = region->getStart();
		// if the user defines a region, you need to start with reads before that, but not before the start of the chromosome
		m_globalStart = std::max( 0 , region->getStart() - AROUND_REGION_BUFFER ); 
	}
	else {
		m_officialStart = m_globalStart = 0;
	}

	if (region->isEndDefined()) {
		m_officialEnd = region->getEnd();
		// if the user defines a region, you need to end with reads after that, but not after the end of the chromosome
		m_globalEnd = std::min( chromosomeSize , region->getEnd() + AROUND_REGION_BUFFER ); 
	}
	else {
		m_officialEnd = m_globalEnd = chromosomeSize;
	}

	m_currentStart = m_globalStart;
	m_displayedStart = m_officialStart;
	updateEndPositions();
}


void LoopingSearchWindow::updateEndPositions()
{
	m_currentEnd = m_currentStart + m_BIN_SIZE;
	if (m_currentEnd > m_globalEnd ) { 
		m_currentEnd = m_globalEnd; 
	}
	m_displayedEnd = m_displayedStart + m_BIN_SIZE;
	if (m_displayedEnd > m_officialEnd ) {
		m_displayedEnd = m_officialEnd;
	}
}


void LoopingSearchWindow::next()
{
	m_currentStart += m_BIN_SIZE;
	m_displayedStart += m_BIN_SIZE;
	updateEndPositions();
}

std::string LoopingSearchWindow::display() const
{
	std::stringstream ss;
	if (m_displayedStart < m_displayedEnd) {
      ss << "\nLooking at chromosome " << m_chromosomeName << " bases " << m_displayedStart << " to " << m_displayedEnd << ".\n";
   }
   else {
      ss << "Checking out reads near the borders of the specified regions for extra evidence.\n";
   }
	return ss.str();
}

bool LoopingSearchWindow::finished() const
{
	// ugly hack for speed purposes when using Pindel-formatted input
	if (UserDefinedSettings::Instance()->pindelFilesAsInput() &&  m_currentStart >= g_maxPos ) { return true; }
	return ( m_currentStart > m_globalEnd );
}

unsigned int SPLIT_READ::getLastAbsLocCloseEnd() const
{
    return UP_Close[ UP_Close.size()-1 ].AbsLoc;
}

bool SPLIT_READ::goodFarEndFound() const
{
    return UP_Far.MaxLen() + UP_Close.MaxLen() >= UnmatchedSeq.size();
}

bool SPLIT_READ::hasCloseEnd() const
{
    return !UP_Close.empty();
}

unsigned int SortedUniquePoints::MaxLen() const
{
    if (m_positions.size() == 0 ) {
        return 0;
    }
    else {
        int lastElementIndex = m_positions.size()-1;
        return m_positions[ lastElementIndex ].LengthStr;
    }
}

unsigned int SPLIT_READ::MaxLenFarEnd() const
{
    return UP_Far.MaxLen();
}

unsigned int SPLIT_READ::MaxLenCloseEnd() const
{
    return UP_Close.MaxLen();
}


bool doOutputBreakdancerEvents()
{
    return ( UserDefinedSettings::Instance()->breakdancerOutputFilename != "" && parameters[findParameter("-b",parameters)]->isSet());
}

void outputBreakDancerEvent( const std::string& chromosomeName, const int leftPosition, const int rightPosition,
                             const int svSize, const std::string& svType, const int svCounter)
{
    std::ofstream breakDancerOutputFile(UserDefinedSettings::Instance()->breakdancerOutputFilename.c_str(), std::ios::app);
    breakDancerOutputFile << chromosomeName << "\t" << leftPosition << "\t" << rightPosition << "\t" <<
                          svSize << "\t" << svType << "\t" << svCounter << "\n";
}

void reportBreakDancerEvent( const std::string& chromosomeName, const int leftPosition, const int rightPosition,
                             const int svSize, const std::string& svType, const int svCounter)
{
    if (doOutputBreakdancerEvents() && g_bdData.isBreakDancerEvent( leftPosition, rightPosition) ) {
        outputBreakDancerEvent(chromosomeName, leftPosition,rightPosition, svSize, svType, svCounter);
    }
}


struct Region {
    Region() {
        start = 0;
        end = 0;
    }
    unsigned start;
    unsigned end;
};


short WhetherRemoveDuplicates;

std::string TempLine_DB_Unique;

std::vector<Region>
Merge(const std::vector<Region> &AllRegions);

bool readTransgressesBinBoundaries(SPLIT_READ & read, const unsigned int &upperBinBorder)
{
    return (read.BPRight > upperBinBorder - 2 * read.InsertSize);
}

/** 'readInSpecifiedRegion' if a region is specified, check if the read is in it. */
bool readInSpecifiedRegion(const SPLIT_READ & read, // in: the read
                           const SearchRegion* region
                          )
{
    bool passesFilter = true;

    // if a start position has been defined, and the breakpoint is before it
    if (region->isStartDefined() && (read.BPLeft + 1 < (unsigned int) region->getStart())) {
        passesFilter = false;
    }

    // if an end of the region has been specified
    if (region->isEndDefined() && (read.BPLeft + 1 > (unsigned int) region->getEnd())) {
        passesFilter = false;

    }
    // no region specified, so all positions are okay
    return passesFilter;
}

void saveReadForNextCycle(SPLIT_READ & read,
                          std::vector<SPLIT_READ> &futureReads)
{
    futureReads.push_back(read);
    read.Used = true; // as it cannot be used for this round of analyses anymore
}

std::ostream& operator<<(std::ostream& os, const UniquePoint& up )
{
    os << "LengthStr: " << up.LengthStr << " * Absloc: " << up.AbsLoc << " * Direction: " << up.Direction << " * Strand: " << up.Strand;
    os << " * Mismatches: " << up.Mismatches << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const SPLIT_READ& splitRead)
{
    os << "Fragname: " << splitRead.FragName << std::endl;
    os << "Name: " << splitRead.Name << std::endl;
    os << "UnmatchedSeq: " << splitRead.getUnmatchedSeq() << std::endl;
    os << "MatchedD: " << splitRead.MatchedD << " * MatchedRelPos: " << splitRead.MatchedRelPos << " * MS: " << splitRead.MS << " * ";
    os << "InsertSize: " << splitRead.InsertSize << std::endl;
    os << "Tag: " << splitRead.Tag << std::endl;
    os << "UP_Close: " ;
    for (unsigned int i=0; i<splitRead.UP_Close.size(); i++) {
        os << "[" << i << "]=" << splitRead.UP_Close[i] << " ";
    }
    os << std::endl;
    os << "UP_Far: " ;
    for (unsigned int i=0; i<splitRead.UP_Far.size(); i++) {
        os << "[" << i << "]=" << splitRead.UP_Far[i] << " ";
    }
    os << std::endl;
    os << std::endl;
    os << "ReadLength: " << splitRead.getReadLength() << std::endl;
    os << "ReadLengthMinus: " << splitRead.getReadLengthMinus() << std::endl;
    os << "MAX_SNP_ERROR:" << splitRead.getMAX_SNP_ERROR() << std::endl;
    os << "TOTAL_SNP_ERROR_CHECKED:" << splitRead.getTOTAL_SNP_ERROR_CHECKED() << std::endl;
    os << "getTOTAL_SNP_ERROR_CHECKED_Minus():" << splitRead.getTOTAL_SNP_ERROR_CHECKED_Minus() << std::endl;
    os << "BP:" << splitRead.BP << std::endl;
    os << "Left:" << splitRead.Left << std::endl;
    os << "Right:" << splitRead.Right << std::endl;
    os << "BPLeft:" << splitRead.BPLeft << std::endl;
    os << "BPRight:" << splitRead.BPRight << std::endl;
    os << "IndelSize:" << splitRead.IndelSize << std::endl;
    os << "UniqueRead:" << splitRead.UniqueRead << std::endl;
    os << "NT_str:" << splitRead.NT_str  << std::endl;
    os << "NT_size:" << splitRead.NT_size << std::endl;
    return os;
}

bool fileExists(const std::string& filename )
{
    std::ifstream test(filename.c_str());
    bool exists = (bool)test;
    test.close();
    return exists;
}

void readBamConfigFile(std::string& bamConfigFilename, ControlState& currentState )
{
    int sampleCounter=0;
    std::ifstream BamConfigFile( bamConfigFilename.c_str() );
    if (BamConfigFile) {
        while (BamConfigFile.good()) {
				bam_info tempBamInfo;
            BamConfigFile >> tempBamInfo.BamFile >> tempBamInfo.InsertSize;
            if (!BamConfigFile.good()) break;
            tempBamInfo.Tag = "";
            BamConfigFile >> tempBamInfo.Tag;
            if (tempBamInfo.Tag=="") {
                *logStream << "Missing tag in line '" << tempBamInfo.BamFile << "\t" << tempBamInfo.InsertSize << "' in configuration file " << bamConfigFilename << "\n";
                exit(EXIT_FAILURE);
            }
            g_sampleNames.insert( tempBamInfo.Tag );
            if (! fileExists( tempBamInfo.BamFile )) {
                *logStream << "I cannot find the file '"<< tempBamInfo.BamFile << "'. referred to in configuration file '" << bamConfigFilename << "'. Please change the BAM configuration file.\n\n";
                exit(EXIT_FAILURE);
            }
            if (! fileExists( tempBamInfo.BamFile+".bai" )) {
                *logStream << "I cannot find the bam index-file '"<< tempBamInfo.BamFile << ".bai' that should accompany the file " << tempBamInfo.BamFile << " mentioned in the configuration file " << bamConfigFilename << ". Please run samtools index on " <<
                          tempBamInfo.BamFile << ".\n\n";
                exit(EXIT_FAILURE);
            }

            //copy kai and throw crap into useless variable
				std::string restOfLine;
            std::getline(BamConfigFile, restOfLine);
            currentState.bams_to_parse.push_back(tempBamInfo);
            sampleCounter++;
        } // while
        if (sampleCounter==0) {
            *logStream << "Could not find any samples in the sample file '" << bamConfigFilename
                      << "'. Please run Pindel again with a config-file of the specified type (format 'A.bam	<insert-size>	sample_label)\n\n";
            exit( EXIT_FAILURE );
        }
    }
    else {
        // no config-file defined
        *logStream << "BAM configuration file '" << bamConfigFilename << "' does not exist. Please run Pindel again with an existing config-file (format 'A.bam	insert-size	sample_label')\n\n";
        exit( EXIT_FAILURE );
    }
}

void readPindelConfigFile(std::string& pindelConfigFilename, std::vector<std::string>& pindelfilesToParse )
{
	int sampleCounter=0;
   std::ifstream pindelConfigFile( pindelConfigFilename.c_str() );
	if (pindelConfigFile) {
		while (pindelConfigFile.good()) {
			std::string pindelFilename;
			pindelConfigFile >> pindelFilename;
			if (!pindelConfigFile.good()) break;

			if (! fileExists( pindelFilename )) {
				*logStream << "I cannot find the file '"<< pindelFilename << "'. referred to in configuration file '" << pindelConfigFilename << "'. Please change the Pindel-configurationfile.\n\n";
				exit(EXIT_FAILURE);
			}

			std::string restOfLine;
			std::getline( pindelConfigFile, restOfLine);
			pindelfilesToParse.push_back( pindelFilename );
			sampleCounter++;
		} // while
		if (sampleCounter==0) {
			*logStream << "Could not find any samples in the sample file '" << pindelConfigFilename
                      << "'. Please run Pindel again with a config-file of the specified type (format 'coloA.pindel_input<newline>colo_tumor.pindel_input<newline>...')\n\n";
			exit( EXIT_FAILURE );
		}
	}
	else {
		// no config-file defined
		*logStream << "Pindel configuration file '" << pindelConfigFilename << "' does not exist. Please run Pindel again with an existing config-file (format 'coloA.pindel_input<newline>colo_tumor.pindel_input<newline>...')\n\n";
		exit( EXIT_FAILURE );
	}
}


LineReader *getLineReaderByFilename(const char *filename)
{
	const std::string strFilename = filename;
	const size_t len = strFilename.length();
	const std::string suffix = strFilename.substr(len - 3, 3);
	
	if (len > 3 && suffix == ".gz")
		return new GZLineReader(filename);
	else
		return new IfstreamLineReader(filename);
}

void TestFileForOutput( const char* filename )
{
   std::ofstream outputfileTest( filename );
   if (!outputfileTest) {
      LOG_ERROR(*logStream << "Sorry, cannot write to the file: " << filename << std::endl);
      exit( EXIT_FAILURE );
   }
   outputfileTest.close();
}

void TestFileForOutput( const std::string& filename )
{
	TestFileForOutput( filename.c_str() );
}

void CheckWhetherFasta( const std::string& filename )
{
	std::ifstream FastaFile( filename.c_str() );
   char FirstCharOfFasta;
   FastaFile >> FirstCharOfFasta;
   if (FirstCharOfFasta != '>') {
      *logStream << "The reference genome must be in fasta format!" << std::endl;
      exit( EXIT_FAILURE );
   }
	FastaFile.close();
}



void init(int argc, char *argv[], ControlState& currentState )
{
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
	logStream=&std::cout;

    if (userSettings->NumRead2ReportCutOff == 1) {
        userSettings->BalanceCutoff = 300000000;
    }

    // define all the parameters you have
    defineParameters( parameters );

    // now read the parameters from the command line
    readParameters(argc, argv, parameters);

	if (userSettings->logFilename != "" ) {
		g_logFile.open( userSettings->logFilename.c_str() );
		logStream = &g_logFile;
	}

    *logStream << Pindel_Version_str << std::endl;

    if (argc <= 1) { // the user has not given any parameters
        printHelp( parameters );
        exit ( EXIT_FAILURE);
    }

    // check parameters
    if (!checkParameters( parameters )) {
        exit ( EXIT_FAILURE);
    }
    if (userSettings->FLOAT_WINDOW_SIZE > 5000.0) {
        LOG_ERROR(*logStream << "Window size of " << userSettings->FLOAT_WINDOW_SIZE << " million bases is too large" << std::endl);
        exit ( EXIT_FAILURE);
    }
    else if (userSettings->FLOAT_WINDOW_SIZE > 100.0) {
        LOG_ERROR(*logStream << "Window size of " << userSettings->FLOAT_WINDOW_SIZE
                  << " million bases is rather large; this may produce bad::allocs or segmentation faults. If that happens, either try to reduce the window size or deactivate the searching for breakpoints and long insertions by adding the command-line options \"-l false -k false\"." << std::endl);
    }
    WINDOW_SIZE = (unsigned int)(1000000 * userSettings->FLOAT_WINDOW_SIZE);

    // if all parameters are okay, open the files

    if (userSettings->singlePindelFileAsInput()) {
		currentState.lineReader= getLineReaderByFilename(userSettings->pindelFilename.c_str());
        currentState.inf_Pindel_Reads = new PindelReadReader(*currentState.lineReader);
    }

   if (userSettings->pindelConfigFileAsInput()) {
		readPindelConfigFile( userSettings->pindelConfigFilename, currentState.pindelfilesToParse );
	}

    if (userSettings->bamFilesAsInput()) {
        readBamConfigFile( userSettings->bamConfigFilename, currentState );
    }

    bool BreakDancerDefined = parameters[findParameter("-b",parameters)]->isSet();
    if (BreakDancerDefined) {
        g_bdData.loadBDFile(userSettings->breakdancerFilename);
    }

    bool AssemblyInputDefined = parameters[findParameter("-z",parameters)]->isSet();
    if (AssemblyInputDefined) {
        currentState.inf_AssemblyInput.open(userSettings->inf_AssemblyInputFilename.c_str());
    }
    
    bool GenotypingInputDefined = parameters[findParameter("-g",parameters)]->isSet();
    if (GenotypingInputDefined) {
        currentState.inf_GenotypingInput.open(userSettings->inf_GenotypingInputFilename.c_str());
    }

    omp_set_num_threads(userSettings->numThreads);

    if (userSettings->MaxRangeIndex > 9) {
       LOG_ERROR(*logStream
                  << "Maximal range index (-x) exceeds the allowed value (9) - resetting to 9."
                  << std::endl);
        userSettings->MaxRangeIndex = 9;
    }

	TestFileForOutput( userSettings->getSIOutputFilename() );
	TestFileForOutput( userSettings->getDOutputFilename() );
	TestFileForOutput( userSettings->getTDOutputFilename() );
	TestFileForOutput( userSettings->getINVOutputFilename() );
	TestFileForOutput( userSettings->getLIOutputFilename() );
	TestFileForOutput( userSettings->getBPOutputFilename() );
	TestFileForOutput( userSettings->getCloseEndOutputFilename() );
	CheckWhetherFasta( userSettings->referenceFilename );

   if ( userSettings->breakdancerOutputFilename != "" ) {
		TestFileForOutput( userSettings->breakdancerOutputFilename );
   }

    Match[(short) 'A'] = 'A';
    Match[(short) 'C'] = 'C';
    Match[(short) 'G'] = 'G';
    Match[(short) 'T'] = 'T';
    Match[(short) 'N'] = 'X';
    Match[(short) '$'] = '$';
    Match2N[(short) 'A'] = 'N';
    Match2N[(short) 'C'] = 'N';
    Match2N[(short) 'G'] = 'N';
    Match2N[(short) 'T'] = 'N';
    Match2N[(short) 'N'] = 'X';
    Match2N[(short) '$'] = '$';
    Convert2RC[(short) 'A'] = 'T';
    Convert2RC[(short) 'C'] = 'G';
    Convert2RC[(short) 'G'] = 'C';
    Convert2RC[(short) 'T'] = 'A';
    Convert2RC[(short) 'N'] = 'X';
    Convert2RC[(short) '$'] = '$';
    Convert2RC4N[(short) 'A'] = 'T';
    Convert2RC4N[(short) 'C'] = 'G';
    Convert2RC4N[(short) 'G'] = 'C';
    Convert2RC4N[(short) 'T'] = 'A';
    Convert2RC4N[(short) 'N'] = 'N';
    Cap2LowArray[(short) 'A'] = 'a';
    Cap2LowArray[(short) 'C'] = 'c';
    Cap2LowArray[(short) 'G'] = 'g';
    Cap2LowArray[(short) 'T'] = 't';
    Cap2LowArray[(short) 'N'] = 'n';
    Cap2LowArray[(short) '$'] = 'n';

    std::string Spacer = "";
    for (unsigned int i = 0; i < g_SpacerBeforeAfter; i++) {
        Spacer += "N";
    }

    Distance = 300;

    DSizeArray[0] = 0;
    DSizeArray[1] = 128;
	for (int dIndex=2; dIndex<15; dIndex++ ) {
		DSizeArray[ dIndex ] = DSizeArray[ dIndex-1 ] * 4;
	}

    if (!userSettings->getRegion()->isTargetChromosomeDefined() && AssemblyInputDefined == false && GenotypingInputDefined == false) {
        *logStream << "Looping over all chromosomes." << std::endl;
    }
}


void SearchFarEnd( const std::string& chromosome, SPLIT_READ& read)
{
   const int START_SEARCH_SPAN = 128;

   // when using bins, some reads may already have been assigned far ends already if they were members of the previous bins; they
   // can be skipped here
   if (read.Investigated) {
       return;
   }

	const std::vector< SearchWindow>& searchCluster =  g_bdData.getCorrespondingSearchWindowCluster( read );
	if (searchCluster.size()!=0) {
		//std::cout << "Starting read searching with cluster size " << searchCluster.size() << " at read " << read.getLastAbsLocCloseEnd() << std::endl;
      SearchFarEndAtPos( chromosome, read, searchCluster);
		//std::cout << "Ending read searching with cluster size " << searchCluster.size() << std::endl;
		if (read.goodFarEndFound()) {
			read.Investigated = true;
	      return;
		}
   }

   UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

    // if breakdancer does not find the event, or not find an event we trust, we turn to regular pattern matching
    int searchSpan=START_SEARCH_SPAN;
    int centerOfSearch = read.getLastAbsLocCloseEnd();
	//std::cout << "Starting regular searching"  << " at read " << read.getLastAbsLocCloseEnd() << "\n";
    for (int rangeIndex=1; rangeIndex<=userSettings->MaxRangeIndex; rangeIndex++ ) {
        SearchFarEndAtPos( chromosome, read, centerOfSearch, searchSpan );
        searchSpan *= 4;
        if (read.goodFarEndFound()) {
				//std::cout << "Ending regular searching\n";
				read.Investigated = true;
            return;
        }
    }
	read.Investigated = true;
}

void ReportCloseMappedReads( const std::vector<SPLIT_READ>& reads )
{
	std::ofstream CloseEndMappedOutput( UserDefinedSettings::Instance()->getCloseEndOutputFilename().c_str(), std::ios::app);
	int TotalNumReads = reads.size();
   for (int Index = 0; Index < TotalNumReads; Index++) {
		const SPLIT_READ& currentRead = reads[ Index ];
   	CloseEndMappedOutput 
			<< currentRead.Name
         << "\n" << currentRead.getUnmatchedSeq()
         << "\n" << currentRead.MatchedD
         << "\t" << currentRead.FragName
         << "\t" << currentRead.MatchedRelPos
         << "\t" << currentRead.MS 
         << "\t" << currentRead.InsertSize 
         << "\t" << currentRead.Tag << "\n";
   }
}

void ReportCloseAndFarEndCounts( const std::vector<SPLIT_READ>& reads )
{
	unsigned Count_Far = 0;
   unsigned Count_Used = 0;
   unsigned Count_Unused = 0;

   for (int Index = reads.size() - 1; Index >= 0; Index--) {
		const SPLIT_READ& currentRead = reads[ Index ];
      if (!currentRead.UP_Far.empty()) {
         Count_Far++;
      }
      if (currentRead.Used) {
         Count_Used++;
      }
   }
	Count_Unused = reads.size() - Count_Far;
   *logStream << "Total: " << reads.size() << ";\tClose_end_found " << reads.size() << ";\tFar_end_found " << Count_Far << ";\tUsed\t" << 
 		Count_Used << ".\n\n";
   *logStream << "For LI and BP: " << Count_Unused << "\n\n";
}

void SearchFarEnds( const std::string chromosomeSeq, std::vector<SPLIT_READ>& reads)
{
	#pragma omp parallel default(shared)
   {
	   #pragma omp for
      for (int readIndex= 0; readIndex < (int)reads.size(); readIndex++ ) {
         SearchFarEnd( chromosomeSeq, reads[readIndex] );
      }
   }
   *logStream << "Far end searching completed for this window." << std::endl;
}


void SearchSVs( ControlState& currentState, const int NumBoxes, const SearchWindow& currentWindow )
{					
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   SearchDeletions searchD;
   searchD.Search(currentState, NumBoxes, currentWindow);

   searchIndels(currentState, NumBoxes, currentWindow);

   if (userSettings->Analyze_TD) {
   	searchTandemDuplications(currentState, NumBoxes, currentWindow);
   	searchTandemDuplicationsNT(currentState, NumBoxes, currentWindow);
   }

   if (userSettings->Analyze_INV) {
   	searchInversions(currentState, NumBoxes, currentWindow);
   	searchInversionsNT(currentState, NumBoxes, currentWindow);
   }

   SearchShortInsertions searchSI;
   searchSI.Search(currentState, NumBoxes, currentWindow);

	ReportCloseAndFarEndCounts( currentState.Reads_SR );                 

   if (userSettings->Analyze_LI) {
   	SortOutputLI(currentState.CurrentChrSeq, currentState.Reads_SR, currentWindow, userSettings->getLIOutputFilename());
   }

   if (userSettings->Analyze_BP) {
   	SortOutputRest(currentState.CurrentChrSeq, currentState.Reads_SR, currentWindow, userSettings->getBPOutputFilename());
   }
}


int main(int argc, char *argv[])
{
	//TODO: These are counters that are only used in individual steps. They should be moved to separate functions later.
   //Below are variables used for cpu time measurement
   time_t Time_Load_S, Time_Load_E, Time_Mine_E, Time_Sort_E;
   Time_Load_S = time(NULL);
   unsigned int AllLoadings = 0;
   unsigned int AllSortReport = 0;
   ControlState currentState;
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   init(argc, argv, currentState );

	std::ifstream FastaFile( userSettings->referenceFilename.c_str() );
   char FirstCharOfFasta;
   FastaFile >> FirstCharOfFasta;
    /* Start of shortcut to genotyping */ // currentState.inf_AssemblyInput.open(par.inf_AssemblyInputFilename.c_str());
	bool GenotypingInputDefined = parameters[findParameter("-g",parameters)]->isSet();
   if (GenotypingInputDefined) {
      doGenotyping(currentState, FastaFile );
      exit(EXIT_SUCCESS);
   }
    
   bool AssemblyInputDefined = parameters[findParameter("-z",parameters)]->isSet();
   if (AssemblyInputDefined) {
      doAssembly(currentState, FastaFile );
      exit(EXIT_SUCCESS);
   }

    // If -q parameter given, search for mobile element insertions and quit.
    if (parameters[findParameter("-q", parameters)]->isSet()) {
        exit(searchMEImain(currentState, FastaFile, userSettings));
    }

   /* Normal pindel functioning: search SVs*/

   // Get a new chromosome again and again until you have visited the specified chromosome or the file ends
   // CurrentChrName stores the name of the chromosome.

   bool SpecifiedChrVisited = false;

   while (SpecifiedChrVisited == false && FastaFile >> currentState.CurrentChrName && !FastaFile.eof()) {

      *logStream << "Processing chromosome: " << currentState.CurrentChrName << std::endl;
        //TODO: check with Kai what's the use of this line.
        // dangerous, there may be no other elements on the fasta header line
	   std::string emptystr;
      std::getline(FastaFile, emptystr);
      if (userSettings->loopOverAllChromosomes()) {
      	GetOneChrSeq(FastaFile, currentState.CurrentChrSeq, true);
      }
      else if (currentState.CurrentChrName == userSettings->getRegion()->getTargetChromosomeName()) {   // just one chr and this is the correct one
         GetOneChrSeq(FastaFile, currentState.CurrentChrSeq, true);
         SpecifiedChrVisited = true;
      }
      else {   // not build up sequence
         GetOneChrSeq(FastaFile, currentState.CurrentChrSeq, false);
         *logStream << "Skipping chromosome: " << currentState.CurrentChrName << std::endl;
         continue;
      }

      CONS_Chr_Size = currentState.CurrentChrSeq.size() - 2 * g_SpacerBeforeAfter; // #################
      g_maxPos = 0; // #################
      *logStream << "Chromosome Size: " << CONS_Chr_Size << std::endl;
      CurrentChrMask.resize(currentState.CurrentChrSeq.size());

      for (unsigned int i = 0; i < currentState.CurrentChrSeq.size(); i++) {
         CurrentChrMask[i] = 'N';
      }
      BoxSize = currentState.CurrentChrSeq.size() / 30000;
        if (BoxSize == 0) BoxSize = 1;
        unsigned NumBoxes = (unsigned) (currentState.CurrentChrSeq.size() * 2 / BoxSize) + 1; // box size
        (*logStream << "NumBoxes: " << NumBoxes << "\tBoxSize: " << BoxSize << std::endl);

        /* 3.2 apply sliding windows to input datasets starts. This is the 2nd level while loop */
        g_binIndex = 0; // to start with 0... 
    
        LoopingSearchWindow currentWindow( userSettings->getRegion(), CONS_Chr_Size, WINDOW_SIZE, currentState.CurrentChrName ); 
			

        // loop over one chromosome
        do {
            /* 3.2.1 preparation starts */

				*logStream << currentWindow.display();
			SearchWindow currentWindow_cs = currentWindow.makePindelCoordinateCopy(); // _cs means computer science coordinates
	     g_bdData.loadRegion( currentWindow_cs );

            if (Time_Load_S == 0) {
                Time_Load_S = time(NULL);
            }
            get_SR_Reads(currentState, currentWindow ); 
            Time_Mine_E = time(NULL);

            if (currentState.Reads_SR.size() ) {
                *logStream << "There are " << currentState.Reads_SR.size() << " reads for this chromosome region." << std::endl; // what region?

                if (userSettings->reportCloseMappedReads() ) {
						 ReportCloseMappedReads( currentState.Reads_SR );       
                }
                Time_Load_E = time(NULL);
                if (!userSettings->reportOnlyCloseMappedReads) {
						SearchFarEnds( currentState.CurrentChrSeq, currentState.Reads_SR );
						SearchSVs( currentState, NumBoxes, currentWindow );
                }
                Time_Sort_E = time(NULL);

                AllLoadings += (unsigned int) difftime(Time_Load_E, Time_Load_S);
                AllSortReport += (unsigned int) difftime(Time_Sort_E, Time_Load_E);
                currentState.Reads_SR.clear();
                *logStream << "There are " << currentState.FutureReads_SR. size()  << " reads saved for the next cycle.\n" << std::endl;
                currentState.Reads_SR.swap(currentState.FutureReads_SR);
            }
            else {
                *logStream << "There are no reads for this bin.\n";
            }
            Time_Load_S = 0;
				currentWindow.next();
            g_binIndex++;

        } // do {
        while (!currentWindow.finished());

    } // while ( loopOverAllChromosomes && chromosomeIndex < chromosomes.size() );

    *logStream << "Loading genome sequences and reads: " << AllLoadings << " seconds." << std::endl;
    *logStream << "Mining, Sorting and output results: " << AllSortReport << " seconds." << std::endl;
    exit( EXIT_SUCCESS) ;
} //main

std::vector<std::string> ReverseComplement(
    const std::vector<std::string> &InputPatterns)
{
    std::vector < std::string > OutputPatterns; // = InputPatterns;
    unsigned int NumPattern = InputPatterns.size();
    OutputPatterns.resize(NumPattern);
    for (unsigned int i = 0; i < NumPattern; i++) {
        OutputPatterns[i] = ReverseComplement(InputPatterns[i]);
    }
    return OutputPatterns;
}

std::string Reverse(const std::string & InputPattern)
{
    std::string OutputPattern = InputPattern;
    unsigned int LenPattern = InputPattern.size();
    for (unsigned int j = 0; j < LenPattern; j++) {
        OutputPattern[j] = InputPattern[j];
    }
    return OutputPattern;
}

std::string ReverseComplement(const std::string & InputPattern)
{
    std::string OutputPattern = InputPattern;

    unsigned int LenPattern = InputPattern.size();

    for (unsigned int j = 0; j < LenPattern; j++)
        OutputPattern[j] = Convert2RC4N[(unsigned int) InputPattern[LenPattern - j - 1]];

    return OutputPattern;
}

std::string Cap2Low(const std::string & input)
{
    std::string output = input;
    for (unsigned int i = 0; i < output.size(); i++) {
        output[i] = Cap2LowArray[(short) input[i]];
    }
    return output;
}

bool ReportEvent(const std::vector<SPLIT_READ> &Deletions,
                 const unsigned int &S, const unsigned int &E)
{
    short ReadLength = Deletions[S].getReadLength() - Deletions[S].NT_size;
    short Min_Length = (short) ((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
    short Max_Length = (short) (ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
    bool LeftMin = false;
    bool LeftMax = false;
    bool RightMin = false;
    bool RightMax = false;
    for (unsigned i = S; i <= E; i++) {
        ReadLength = Deletions[i].getReadLength() - Deletions[i].NT_size;
        Min_Length = (short) ((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
        Max_Length = (short) (ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
        if (Deletions[i].BP <= Min_Length) {
            LeftMin = true;
        }
        if (Deletions[i].getReadLength() - Deletions[i].BP - Deletions[i].NT_size <= Min_Length) {
            RightMin = true;
        }
        if (Deletions[i].BP >= Max_Length) {
            LeftMax = true;
        }
        if (Deletions[i].getReadLength() - Deletions[i].BP - Deletions[i].NT_size >= Max_Length) {
            RightMax = true;
        }
    }

    if (LeftMin && LeftMax && RightMin && RightMax) {
        return true;
    }
    else {
        return false;
    }
}

void GetRealStart4Deletion(const std::string & TheInput,
                           unsigned int &RealStart, unsigned int &RealEnd)
{
    unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;
    unsigned int Start = PosIndex + 1;
    unsigned int End = RealEnd + g_SpacerBeforeAfter - 1;
    while (TheInput[PosIndex] == TheInput[End]) {
        --PosIndex;
        --End;
    }
    RealStart = PosIndex - g_SpacerBeforeAfter;
    PosIndex = RealEnd + g_SpacerBeforeAfter;
    while (TheInput[PosIndex] == TheInput[Start]) {
        ++PosIndex;
        ++Start;
    }
    RealEnd = PosIndex - g_SpacerBeforeAfter;
}

void GetRealStart4Insertion(const std::string & TheInput,
                            const std::string & InsertedStr, unsigned int &RealStart,
                            unsigned int &RealEnd)
{
    unsigned int IndelSize = InsertedStr.size();
    unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;

    for (int i = IndelSize - 1; i >= 0; i--) {
        if (TheInput[PosIndex] == InsertedStr[i]) {
            PosIndex--;
        }
        else {
            break;
        }
    }
    if (PosIndex == RealStart + g_SpacerBeforeAfter - IndelSize) {
        while (TheInput[PosIndex] == TheInput[PosIndex + IndelSize]) {
            PosIndex--;
        }
    }
    RealStart = PosIndex - g_SpacerBeforeAfter;
    PosIndex = RealEnd + g_SpacerBeforeAfter;
    for (unsigned int i = 0; i < IndelSize; i++) {
        if (TheInput[PosIndex] == InsertedStr[i]) {
            PosIndex++;
        }
        else {
            break;
        }
    }
    if (PosIndex == RealEnd + g_SpacerBeforeAfter + IndelSize) {
        while (TheInput[PosIndex] == TheInput[PosIndex - IndelSize]) {
            PosIndex++;
        }
    }
    RealEnd = PosIndex - g_SpacerBeforeAfter;
}

std::vector<Region> Merge(const std::vector<Region> &AllRegions)
{
    return AllRegions;
}

void GetCloseEndInner(const std::string & CurrentChrSeq, SPLIT_READ & Temp_One_Read)
{
    std::string CurrentReadSeq;
    std::vector<unsigned int> PD[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
    if (Temp_One_Read.InsertSize > g_maxInsertSize) {
        g_maxInsertSize = Temp_One_Read.InsertSize;
    }
    for (int CheckIndex = 0; CheckIndex < Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(); CheckIndex++) {
        PD[CheckIndex].reserve(3 * Temp_One_Read.InsertSize);
    }
    SortedUniquePoints UP;
    int Start, End;
    short BP_Start; // = MinClose;
    short BP_End; // = ReadLength - MinClose;

    Temp_One_Read.UP_Close.clear();
    BP_Start = g_MinClose;
    BP_End = Temp_One_Read.getReadLengthMinus();
    if (Temp_One_Read.MatchedD == Plus) {
        CurrentReadSeq = ReverseComplement(Temp_One_Read.getUnmatchedSeq());
        Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter;
        End = Start + 3 * Temp_One_Read.InsertSize;
        char LeftChar;
        LeftChar = CurrentReadSeq[0];
        if (LeftChar != 'N') {
            {
                for (int pos = Start; pos < End; pos++) {
                    if (CurrentChrSeq[pos] == LeftChar) {
                        PD[0].push_back(pos);
                    }
                }
            }
        }

        if (PD[0].size())
            CheckLeft_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
        if (UP.empty()) {}
        else {
            Temp_One_Read.Used = false;
            Temp_One_Read.UP_Close.swap(UP);
            UP.clear();
        }
    }
    else if (Temp_One_Read.MatchedD == Minus) {
        CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
        End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter;
        Start = End - 3 * Temp_One_Read.InsertSize;
        char RightChar;
        RightChar = CurrentReadSeq[Temp_One_Read.getReadLengthMinus()];
        if (RightChar != 'N') {
            for (int pos = Start; pos < End; pos++) {
                if (CurrentChrSeq[pos] == RightChar) {
                    PD[0].push_back(pos);
                }
            }
        }

        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
        CheckRight_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD,
                         BP_Start, BP_End, 1, UP);
        LOG_DEBUG(*logStream << UP.size() << std::endl);
        if (UP.empty()) {}
        else {
            Temp_One_Read.Used = false;
            Temp_One_Read.UP_Close.swap(UP);
            UP.clear();
        }
    }
    return;
}

void GetCloseEnd(const std::string & CurrentChrSeq, SPLIT_READ & Temp_One_Read)
{

    GetCloseEndInner( CurrentChrSeq, Temp_One_Read );

    if (Temp_One_Read.UP_Close.size()==0) { // no good close ends found
        Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );

        GetCloseEndInner( CurrentChrSeq, Temp_One_Read );

    }
}

void CheckBoth(const SPLIT_READ & OneRead, const std::string & TheInput,
               const std::string & CurrentReadSeq,
               const std::vector<unsigned int> PD_Plus[],
               const std::vector<unsigned int> PD_Minus[], const short minimumLengthToReportMatch,
               const short BP_End, const short CurrentLength,
               SortedUniquePoints &UP)
{   
	int Sum;

	UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   if (CurrentLength >= minimumLengthToReportMatch && CurrentLength <= BP_End) {
      for (short numberOfMismatches = 0; numberOfMismatches <= OneRead.getMAX_SNP_ERROR(); numberOfMismatches++) {
         if (PD_Plus[numberOfMismatches].size() + PD_Minus[numberOfMismatches].size() == 1 && CurrentLength >= minimumLengthToReportMatch + numberOfMismatches) {
            Sum = 0;
            if (userSettings->ADDITIONAL_MISMATCH > 0) { // what if this is ADDITIONAL_MISMATCH is 0? Do you save anything then?
					// what if j +ADD_MISMATCH exceeds the max allowed number of mismatches (the PD_element does not exist?)
				
					// feeling the need to comment here - so factor out?
					// only report reads if there are no reads with fewer mismatches or only one or two more mismatches
               for (short mismatchCount = 0; mismatchCount <= numberOfMismatches + userSettings->ADDITIONAL_MISMATCH; mismatchCount++) {
                  Sum += PD_Plus[mismatchCount].size() + PD_Minus[mismatchCount].size();
               }
               if (Sum == 1 && numberOfMismatches <= (short) (userSettings->Seq_Error_Rate * CurrentLength + 1)) {
							// why I love constructors
						UniquePoint MatchPosition;
                  if (PD_Plus[numberOfMismatches].size() == 1) {
							UniquePoint PlusMatch( CurrentLength, PD_Plus[numberOfMismatches][0], FORWARD, SENSE, numberOfMismatches );
							MatchPosition = PlusMatch; 
                  }
                  else {
							UniquePoint MinMatch( CurrentLength, PD_Minus[numberOfMismatches][0], BACKWARD, ANTISENSE, numberOfMismatches );
							MatchPosition = MinMatch; 
					   }
                  if (CheckMismatches(TheInput, OneRead.getUnmatchedSeq(), MatchPosition)) {
                     UP.push_back (MatchPosition);
                     break;
                  } // if CheckMismatches
               } // if Sum==1
            } // if AdditionalMismatches
        	} // if sumsize ==1
      } // for-loop
	} // if length of match is sufficient to be reportable

   if (CurrentLength < BP_End) {
		// the below lines give me an odd deja-vu-feeling
      std::vector<unsigned int> PD_Plus_Output[OneRead.getTOTAL_SNP_ERROR_CHECKED()];
      std::vector<unsigned int> PD_Minus_Output[OneRead.getTOTAL_SNP_ERROR_CHECKED()];
		
      for (int CheckedIndex = 0; CheckedIndex < OneRead.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
         PD_Plus_Output[CheckedIndex].reserve( PD_Plus[CheckedIndex].size());
         PD_Minus_Output[CheckedIndex].reserve( PD_Minus[CheckedIndex].size());
      }
      const char CurrentChar = CurrentReadSeq[CurrentLength];
      const char CurrentCharRC = Convert2RC4N[(short) CurrentChar];
      unsigned int pos=0;
      int SizeOfCurrent=0;
		for (int i = 0; i < OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) {
			
         SizeOfCurrent = PD_Plus[i].size();
         if (CurrentChar == 'N') {
            for (int j = 0; j < SizeOfCurrent; j++) {
               pos = PD_Plus[i][j] + 1;
               if (Match2N[(short) TheInput[pos]] == 'N') {
                  PD_Plus_Output[i].push_back(pos);
               }
               else {
                  PD_Plus_Output[i + 1].push_back(pos);
               }
            }
		   }
         else { // currentChar != 'N'; note that these two loops could be fused if we had a Match function or array that would function properly for normal bases and Ns
            for (int j = 0; j < SizeOfCurrent; j++) {
               pos = PD_Plus[i][j] + 1;
               if (TheInput[pos] == CurrentChar) {
                  PD_Plus_Output[i].push_back(pos);
               }
               else {
                  PD_Plus_Output[i + 1].push_back(pos);
               }
            }
         }
			
			// code repetition; I'd experiment with factoring these two out in a separate function	and check performance
         SizeOfCurrent = PD_Minus[i].size();
         if (CurrentCharRC == 'N') {
            for (int j = 0; j < SizeOfCurrent; j++) {
               pos = PD_Minus[i][j] - 1;
               if (Match2N[(short) TheInput[pos]] == 'N') {
                  PD_Minus_Output[i].push_back(pos);
               }
               else {
                  PD_Minus_Output[i + 1].push_back(pos);
               }
            }
         }
         else {
            for (int j = 0; j < SizeOfCurrent; j++) {
               pos = PD_Minus[i][j] - 1;
               if (TheInput[pos] == CurrentCharRC) {
                  PD_Minus_Output[i].push_back(pos);
               }
               else {
                  PD_Minus_Output[i + 1].push_back(pos);
               }
            }
         }
      }

		// really strong DEJA-VU-feeling; prime candidate for deduplicating code and testing performance afterwards.
      SizeOfCurrent = PD_Plus[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()].size();
      if (CurrentChar == 'N') {
         for (int j = 0; j < SizeOfCurrent; j++) {
            pos = PD_Plus[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()][j] + 1;
            if (Match2N[(short) TheInput[pos]] == 'N') {
               PD_Plus_Output[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()]. push_back( pos);
				}
         }
      }
      else {
         for (int j = 0; j < SizeOfCurrent; j++) {
            pos = PD_Plus[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()][j] + 1;
            if (TheInput[pos] == CurrentChar) {
               PD_Plus_Output[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()]. push_back( pos);
				}
         }
      }
      SizeOfCurrent = PD_Minus[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()].size();
      if (CurrentCharRC == 'N') {
         for (int j = 0; j < SizeOfCurrent; j++) {
            pos = PD_Minus[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()][j] - 1;
            if (Match2N[(short) TheInput[pos]] == 'N') {
               PD_Minus_Output[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()]. push_back( pos);
				}
         }
      }
      else {
         for (int j = 0; j < SizeOfCurrent; j++) {
            pos = PD_Minus[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()][j] - 1;
            if (TheInput[pos] == CurrentCharRC) {
               PD_Minus_Output[OneRead.getTOTAL_SNP_ERROR_CHECKED_Minus()]. push_back( pos);
				}
         }
      }

		// this loop looks familiar; candidate for factoring out mini-function?
      Sum = 0;
      for (int i = 0; i <= OneRead.getMAX_SNP_ERROR(); i++) {
         Sum += PD_Plus_Output[i].size() + PD_Minus_Output[i].size();
      }
      if (Sum) {
         const short CurrentLengthOutput = CurrentLength + 1;
         CheckBoth(OneRead, TheInput, CurrentReadSeq, PD_Plus_Output, PD_Minus_Output, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);
      }
      else {
         return;
   	} // else-if Sum
 	} // if length < BP_End
 	else {
	   return;
   } // if length>=BP_End
}

void CleanUniquePoints(SortedUniquePoints &Input_UP)
{
    SortedUniquePoints TempUP; //vector <UniquePoint> UP_Close; UP_Far
    UniquePoint LastUP = Input_UP[Input_UP.size() - 1];
    char LastDirection = LastUP.Direction;
    char LastStrand = LastUP.Strand;
    unsigned int Terminal;

    if (LastDirection == FORWARD) {
        Terminal = LastUP.AbsLoc - LastUP.LengthStr;
        for (unsigned i = 0; i < Input_UP.size(); i++) {
           if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
              if (Terminal == Input_UP[i].AbsLoc - Input_UP[i].LengthStr) {
                 TempUP.push_back(Input_UP[i]);
              }
           }
        }
    }
    else if (LastDirection == BACKWARD) {
       Terminal = LastUP.AbsLoc + LastUP.LengthStr;
       for (unsigned i = 0; i < Input_UP.size(); i++) {
          if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
             if (Terminal == Input_UP[i].AbsLoc + Input_UP[i].LengthStr) {
                TempUP.push_back(Input_UP[i]);
             }
          }
       }
    }
   Input_UP.clear();
   Input_UP = TempUP;
}


