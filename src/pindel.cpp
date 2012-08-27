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
#include <sstream>

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
#include "searcher.h"
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

const Chromosome g_dummyChromosome("","");
Genome g_genome;
std::ofstream g_logFile;

int g_binIndex = -1; // global variable for the bin index, as I cannot easily pass an extra parameter to the diverse functions
unsigned int g_maxPos = 0; // to calculate which is the last position in the chromosome, and hence to calculate the number of bins
short g_MinClose = 8;
std::set<std::string> g_sampleNames;
short Before, After;
BDData g_bdData;
const int alphs = 4;
const char alphabet[alphs] = { 'A', 'C', 'G', 'T' };
unsigned int BoxSize = 10000; // 10k is fine
const double Min_Filter_Ratio = 0.5;
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
unsigned int DSizeArray[15];
int g_maxInsertSize=0;
std::string CurrentChrMask;
std::vector<Parameter *> parameters;



// #########################################################

int NumRead2ReportCutOff_BP = 2;
const int g_MAX_RANGE_INDEX = 9; // 5 or 6 or 7 or maximum 8      //# // user
unsigned int WINDOW_SIZE = 10000000;
const unsigned int AROUND_REGION_BUFFER = 10000; // how much earlier reads should be selected if only a region of the chromosome needs be specified.
// #########################################################

const short MaxDI = 30;

// from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
void safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    /*std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\r':
            c = sb->sgetc();
            if(c == '\n')
                sb->sbumpc();
            return is;
        case '\n':
        case EOF:
            return is;
        default:
            t += (char)c;
        }
    }*/
	getline( is, t );
	/*unsigned int lastIndex = t.size()-1;
	while (lastIndex>=0 && t[ lastIndex ]=='\r' ) {
		t.resize( lastIndex );
		lastIndex--;
	} */
	//return is;
}


void SPLIT_READ::setUnmatchedSeq( const std::string unmatchedSeq ) 
{ 
	//const double EPSILON = 0.00001; // to compensate for downrounding off errors (0.04 = 0.03999999 on some computers)
	
	UnmatchedSeq = unmatchedSeq; 
	if (UnmatchedSeq.size()>0) {
		/*for (unsigned int x=0; x<UnmatchedSeq.size();x++ ) {
			std::cout << "seq[" << x << "]="<< int(UnmatchedSeq[ x ]) << "('" << UnmatchedSeq[ x ] << "')\n";
		}*/
		unsigned int lastCharIndex = UnmatchedSeq.size()-1;
		while (!isalnum( UnmatchedSeq[ lastCharIndex ] )) {
			UnmatchedSeq.resize( lastCharIndex );
			lastCharIndex--;
		} 
		/*for (unsigned int x=0; x<UnmatchedSeq.size();x++ ) {
			std::cout << "resseq[" << x << "]="<< int(UnmatchedSeq[ x ]) << "('" << UnmatchedSeq[ x ] << "')\n";
		}*/
	}

	ReadLength = UnmatchedSeq.size();
	ReadLengthMinus = ReadLength - 1;

	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

	MAX_SNP_ERROR = g_maxMismatch[ ReadLength ];
	TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + userSettings->ADDITIONAL_MISMATCH;
	TOTAL_SNP_ERROR_CHECKED = TOTAL_SNP_ERROR_CHECKED_Minus + 1;
}

const Chromosome* Genome::addChromosome( Chromosome* newChromosome ) 
{ 
	for (unsigned int i=0; i<m_chromosomes.size(); i++ ) {
		if ( m_chromosomes[ i ]->getName() == newChromosome->getName() ) {
			delete m_chromosomes[ i ];
			m_chromosomes[ i ] = newChromosome;
			return newChromosome;
		}
	}
	m_chromosomes.push_back( newChromosome );
	return newChromosome;
}

const Chromosome* Genome::getChr( unsigned int index ) const 
{ 
	if (index<m_chromosomes.size()) return m_chromosomes[ index ]; 
	else return NULL; 
}

const Chromosome* Genome::getChr( const std::string& chromosomeName ) const 
{ 
	for (unsigned int i=0; i<m_chromosomes.size(); i++ ) {
		if ( m_chromosomes[ i ]->getName() == chromosomeName ) {
			return m_chromosomes[ i ];
		}
	}
	return NULL;
}

void Genome::clear()
{
	for (unsigned int i=0; i<m_chromosomes.size(); i++ ) {
		delete m_chromosomes[ i ];
	}
	m_fullMode = false; 
	m_currentChromosomeIndex = -1;
	m_referenceFile.close();
}

void Genome::reset() {
    m_referenceFile.clear();
    m_referenceFile.seekg(0, std::ios::beg);
    m_currentChromosomeIndex = -1;
}

void Genome::load( const std::string& referenceFileName )
{
	clear();
	m_referenceFile.open( referenceFileName.c_str() );
}

void Genome::loadAll(const std::string& referenceFileName)
{
	load( referenceFileName );
	while (!m_referenceFile.eof() ) {
		loadChromosome();
	}
	m_fullMode = true;
}

const Chromosome* Genome::getNextChromosome()
{
	if (m_fullMode) {
		// all chromosomes loaded 
		if ( (unsigned int)(m_currentChromosomeIndex+1)<m_chromosomes.size() ) {
			m_currentChromosomeIndex++;
			return m_chromosomes[ m_currentChromosomeIndex ];
		}
		else {
			return NULL;
		}
	}
	else {
		// need to load next chromosome
		if (m_currentChromosomeIndex>=0) {
			// there is a previous chromosome
			const std::string oldName = m_chromosomes[ m_currentChromosomeIndex ]->getName();
			delete m_chromosomes[ m_currentChromosomeIndex ];
			m_chromosomes[ m_currentChromosomeIndex ] = new Chromosome( oldName, "");
		}
		m_currentChromosomeIndex++;
		return loadChromosome();
	}
}



const Chromosome* Genome::loadChromosome()
{	
	std::string chromosomeName = "";
	std::string chromosomeSeq = "";
	if (m_referenceFile.eof()) {
		return NULL;
	}
	char chrIndicatorChar; // '>' preceeds all chromosome names in a FASTA file
	m_referenceFile >> chrIndicatorChar;
	if (chrIndicatorChar != '>' ) {
		std::cout << "Error: fasta line starts with " << chrIndicatorChar << " instead of '>'. Aborting.\n";
		exit( EXIT_FAILURE );
		return NULL;
	}
	m_referenceFile >> chromosomeName;
	std::string restOfChromosomeHeadline;
	safeGetline( m_referenceFile, restOfChromosomeHeadline );
	char ch;
	do {
		m_referenceFile >> ch;
		if (ch=='>') {
			// next chromosome
			m_referenceFile.putback( ch );
			break;
		} 
		else {
			ch = toupper( ch );
			if (ch!='A' && ch!='C' && ch!='G' && ch!='T') { 
				ch='N';
			}
			chromosomeSeq += ch;
		}		
	} while (!m_referenceFile.eof() );

   std::string Spacer = "";
   for (unsigned i = 0; i < g_SpacerBeforeAfter; i++) {
      Spacer += "N";
   }
   chromosomeSeq = Spacer + chromosomeSeq + Spacer;
	Chromosome* newChromosome= new Chromosome( chromosomeName, chromosomeSeq );
	return addChromosome( newChromosome );
}

UniquePoint::UniquePoint( const Chromosome* chromosome_ptr, const short lengthStr, const unsigned int absLoc, const char direction, const char strand, const short mismatches ) :
	chromosome_p( chromosome_ptr), LengthStr( lengthStr ), AbsLoc( absLoc ), Direction( direction ), Strand( strand ), Mismatches( mismatches )
{
}

SearchWindow& SearchWindow::operator=(const SearchWindow& other )
{
   if (this != &other) {
      this->SearchWindow::~SearchWindow(); // explicit non-virtual destructor
      new (this) SearchWindow( other ); // placement new
   }
   return *this;
}

bool SearchWindow::encompasses( const std::string& chromosomeName, const unsigned int position ) const
{
	return ( ( m_chromosome->getName() == chromosomeName ) && ( position >= m_currentStart ) && ( position <= m_currentEnd ) );
}

SearchWindow::SearchWindow( const Chromosome* chromosome, const unsigned int regionStart, const unsigned int regionEnd ) : m_chromosome( chromosome )
{
	m_currentStart = regionStart;
	m_currentEnd = regionEnd;
}

SearchWindow::SearchWindow(const SearchWindow& other) : 
	m_chromosome(other.m_chromosome), m_currentStart( other.m_currentStart ), m_currentEnd( other.m_currentEnd) 
{}

LoopingSearchWindow::LoopingSearchWindow(const SearchRegion* region, const Chromosome* chromosome, const int binSize ) : 
	SearchWindow(chromosome,0,chromosome->getBiolSize()), m_BIN_SIZE( binSize )
{
	if (region->isStartDefined()) {
		m_officialStart = region->getStart();
		// if the user defines a region, you need to start with reads before that, but not before the start of the chromosome
		m_globalStart = ( region->getStart() >= AROUND_REGION_BUFFER ? region->getStart() - AROUND_REGION_BUFFER : 0); 
	}
	else {
		m_officialStart = m_globalStart = 0;
	}

	if (region->isEndDefined()) {
		m_officialEnd = region->getEnd();
		// if the user defines a region, you need to end with reads after that, but not after the end of the chromosome
		m_globalEnd = std::min( chromosome->getBiolSize() , region->getEnd() + AROUND_REGION_BUFFER ); 
	}
	else {
		m_officialEnd = m_globalEnd = chromosome->getBiolSize();
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
      ss << "\nLooking at chromosome " << m_chromosome->getName() << " bases " << m_displayedStart << " to " << m_displayedEnd << ".\n";
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
    os << "\nName: " << splitRead.Name << std::endl;
    os << "Fragname: " << splitRead.FragName << std::endl;
    os << "FarFragname: " << splitRead.FarFragName << std::endl;
    
    os << "UnmatchedSeq: " << splitRead.getUnmatchedSeq() << std::endl;
    os << "MatchedD: " << splitRead.MatchedD << " * MatchedRelPos: " << splitRead.MatchedRelPos << " * MS: " << splitRead.MS << " * ";
    os << "InsertSize: " << splitRead.InsertSize << std::endl;
    os << "Tag: " << splitRead.Tag << std::endl;
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
    //os << std::endl;

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
            safeGetline(BamConfigFile, restOfLine);
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
			safeGetline( pindelConfigFile, restOfLine);
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

/** 'probOfReadWithTheseErrors' gives the probability that a read of length 'length' will have number of errors 'numberOfErrors' if the sequencing error rate is 'seqErrorRate'. */
double probOfReadWithTheseErrors(const unsigned int length, const unsigned int numberOfErrors, const double seqErrorRate )
{
	double chanceOfCorrect = 1.0-seqErrorRate;
	unsigned int numberOfCorrectBases = length - numberOfErrors;
	double matchedPart = pow( chanceOfCorrect, numberOfCorrectBases );

	double mismatchedPart = 1.0;
	for (unsigned int i=0; i<numberOfErrors; i++ ) {
		mismatchedPart *= (((length-i) * seqErrorRate) / (numberOfErrors-i) );
	}
	return matchedPart * mismatchedPart;
}

std::vector<unsigned int> g_maxMismatch;

/**  'createProbTable' fills the g_probMismatch[ length ] table with the maximum amount of differences between read and reference are acceptable.
	For example: if sequencing error rate is 1%, and SNP rate is 0.1%, total error rate is 1.1%. Throwing away all 100-base reads with 1 error would throw away over 1.1% of reads, 
	which is acceptable if the sensitivity is set to 95% (0.95), but not if it is set to 1%. Uses binomial formula. */
void createProbTable(const double seqErrorRate, const double sensitivity) 
{
	const unsigned int MAX_LENGTH = 500;
	g_maxMismatch.assign( MAX_LENGTH, 0);
	for (unsigned int length=0; length<MAX_LENGTH; length++) {
		double totalErrorProb = 0.0;
		for (unsigned int numberOfErrors=0; numberOfErrors<=length; numberOfErrors++ ) {
			totalErrorProb += probOfReadWithTheseErrors( length, numberOfErrors, seqErrorRate );
			if (totalErrorProb > sensitivity ) {
				g_maxMismatch[ length ] = numberOfErrors;
				//g_maxMisMatch.push_back( numberOfErrors );
				//std::cout << length << " bases has max errors" << numberOfErrors << "\n";
				break; // break out of this length, up to the next
			}
		}
	}
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
	createProbTable(0.001+userSettings->Seq_Error_Rate, userSettings->sensitivity); 
	std::string fastaFilename( userSettings->referenceFilename.c_str() );
	if (userSettings->reportInterchromosomalEvents) {
		g_genome.loadAll( fastaFilename );
	}
	else {
		g_genome.load( fastaFilename );
	}

   bool BreakDancerDefined = parameters[findParameter("-b",parameters)]->isSet();
   if (BreakDancerDefined) {
      g_bdData.loadBDFile(userSettings->breakdancerFilename);
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



    bool AssemblyInputDefined = parameters[findParameter("-z",parameters)]->isSet();
    if (AssemblyInputDefined) {
        currentState.inf_AssemblyInput.open(userSettings->inf_AssemblyInputFilename.c_str());
    }
    
    bool GenotypingInputDefined = parameters[findParameter("-g",parameters)]->isSet();
    if (GenotypingInputDefined) {
        currentState.inf_GenotypingInput.open(userSettings->inf_GenotypingInputFilename.c_str());
    }

    omp_set_num_threads(userSettings->numThreads);

    if (userSettings->MaxRangeIndex > g_MAX_RANGE_INDEX) {
       LOG_ERROR(*logStream
                  << "Maximal range index (-x) exceeds the allowed value (" << g_MAX_RANGE_INDEX << ") - resetting to " << g_MAX_RANGE_INDEX << ".\n" );
        userSettings->MaxRangeIndex = g_MAX_RANGE_INDEX;
    }

	if (userSettings->ADDITIONAL_MISMATCH<1) {
       LOG_ERROR(*logStream << "Number of additional mismatches (-a) is less than the allowed value (1) - resetting to 1.\n" );
        userSettings->ADDITIONAL_MISMATCH = 1;
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

    //Distance = 300;

    DSizeArray[0] = 0;
    DSizeArray[1] = 128;
	for (int dIndex=2; dIndex<15; dIndex++ ) {
		DSizeArray[ dIndex ] = DSizeArray[ dIndex-1 ] * 4;
	}

    if (!userSettings->getRegion()->isTargetChromosomeDefined() && AssemblyInputDefined == false && GenotypingInputDefined == false) {
        *logStream << "Looping over all chromosomes." << std::endl;
    }
}


void SearchFarEnd( const std::string& chromosome, SPLIT_READ& read, const Chromosome& currentChromosome)
{
	const int START_SEARCH_SPAN = 128;

   // when using bins, some reads may already have been assigned far ends already if they were members of the previous bins; they
   // can be skipped here
   /*if (read.Investigated) {
       return;
   }*/

//std::cout << "CurrentChromosome name /seq =" << currentChromosome.getName() << " = " << currentChromosome.getSeq() << "\n";
	const std::vector< SearchWindow>& searchCluster =  g_bdData.getCorrespondingSearchWindowCluster( read );
  	 //std::cout << read.Name << " - " << read.getLastAbsLocCloseEnd() << "\n"; 
	if (searchCluster.size()!=0) {
      //  std::cout << "Breakdancer input is not empty " << searchCluster.size() << std::endl;
      SearchFarEndAtPos( chromosome, read, searchCluster);
		if (read.goodFarEndFound()) {
			//read.Investigated = true;
	      return;
		}
   }

   UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   // if breakdancer does not find the event, or not find an event we trust, we turn to regular pattern matching
   unsigned int searchSpan=START_SEARCH_SPAN;
   unsigned int centerOfSearch = read.getLastAbsLocCloseEnd();
	
   for (int rangeIndex=1; rangeIndex<=userSettings->MaxRangeIndex; rangeIndex++ ) {
		// note: searching the central range again and again may seem overkill, but since Pindel at the moment wants an unique best point, you can't skip the middle part
		// may be stuff for future changes/refactorings though
		std::vector< SearchWindow > aroundCESearchCluster;
      unsigned int Start;
      if (centerOfSearch > searchSpan+g_SpacerBeforeAfter) {
			Start = centerOfSearch - searchSpan;
		}
      else {
			Start = g_SpacerBeforeAfter;
		}
		//std::cout << "Start is (abs) " << Start << " (rel): " << Start - g_SpacerBeforeAfter << "\n";
		//std::cout << "End is (abs) " << centerOfSearch+searchSpan << " (rel): " << centerOfSearch+searchSpan-10000000 << "Size chrom=" << chromosome.size() << "\n";		
      //std::cout << Start << "FirstStart " << centerOfSearch+searchSpan-10000000 << "<COS" << centerOfSearch-10000000 << " span " << searchSpan<< std::endl;
		SearchWindow regularWindow( &currentChromosome, Start, centerOfSearch+searchSpan );
      //std::cout << Start << " " << centerOfSearch+searchSpan-10000000 << std::endl;
		aroundCESearchCluster.push_back( regularWindow );
      SearchFarEndAtPos( chromosome, read, aroundCESearchCluster );
      searchSpan *= 4;
      if (read.goodFarEndFound()) {
			//read.Investigated = true;
         return;
      }
   }
	//read.Investigated = true;
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

void SearchFarEnds( const std::string chromosomeSeq, std::vector<SPLIT_READ>& reads, const Chromosome& currentChromosome)
{
	#pragma omp parallel default(shared)
   {
	   #pragma omp for
      for (int readIndex= 0; readIndex < (int)reads.size(); readIndex++ ) {
         SearchFarEnd( chromosomeSeq, reads[readIndex], currentChromosome );
          
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
   	SortOutputLI(currentWindow.getChromosome()->getSeq(), currentState.Reads_SR, currentWindow, userSettings->getLIOutputFilename());
   }

   if (userSettings->Analyze_BP) {
   	SortOutputRest(currentWindow.getChromosome()->getSeq(), currentState.Reads_SR, currentWindow, userSettings->getBPOutputFilename());
   }
}

class TimerItem {

public:
	TimerItem( const std::string& id );
	void stop();
	void restart();
	const std::string getReport() const;
	const std::string& getId() const { return m_id; }	

private:		
	std::string m_id;
	time_t m_countSoFar;
	time_t m_lastStart;
};

TimerItem::TimerItem( const std::string& id ) 
{
	m_id = id;
	m_countSoFar = 0;
	m_lastStart = time( NULL );
}

void TimerItem::stop()
{
	m_countSoFar += (time( NULL ) - m_lastStart);
}

void TimerItem::restart()
{
	m_lastStart = time( NULL );	
}

const std::string TimerItem::getReport() const
{
	std::stringstream ss;
	ss <<  m_id << " " << m_countSoFar << " seconds.\n";
	return ss.str();
}


class Timer {

public:
	void switchTo( const std::string& itemName );
	void reportAll(std::ostream& os);
	Timer() : m_currentItemIndex( -1 ) {}

private:
	std::vector< TimerItem > m_timerItems;
	int m_currentItemIndex;
};

void Timer::switchTo( const std::string& itemName )
{
	if (m_currentItemIndex!=-1) {
		m_timerItems[ m_currentItemIndex ].stop();
	}
	for (unsigned int itemIter=0; itemIter<m_timerItems.size(); itemIter++ ) {
		if ( m_timerItems[ itemIter].getId() == itemName ) {
			m_timerItems[ itemIter].restart();
			m_currentItemIndex = (int)itemIter;
			return;
		}
	}
	// no existing element found? new element needs to be constructed
	TimerItem newItem( itemName );
	m_timerItems.push_back( newItem );
	m_currentItemIndex = m_timerItems.size() - 1;
}

void Timer::reportAll(std::ostream& os)
{
	if (m_currentItemIndex!=-1) {
		m_timerItems[ m_currentItemIndex ].stop();
	}
	for (std::vector< TimerItem>::iterator itemIter=m_timerItems.begin(); itemIter!=m_timerItems.end(); itemIter++ ) {
		os << itemIter->getReport();
	}	
}

void UpdateFarFragName(std::vector <SPLIT_READ> & input) {
    for (unsigned index = 0; index < input.size(); index++) {
        if (input[index].UP_Far.size()) {
            input[index].FarFragName = input[index].UP_Far[0].chromosome_p->getName();
            input[index].MatchedFarD = input[index].UP_Far[0].Strand;
        }
    }
}



int main(int argc, char *argv[])
{
	//TODO: These are counters that are only used in individual steps. They should be moved to separate functions later.
   //Below are variables used for cpu time measurement
   time_t Time_Load_S, Time_Load_E, Time_Mine_E, Time_Sort_E;
	Timer timer;
	timer.switchTo("Initializing pindel");
	
   Time_Load_S = time(NULL);
   unsigned int AllLoadings = 0;
   unsigned int AllSortReport = 0;
   ControlState currentState;
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   init(argc, argv, currentState );

    /* Start of shortcut to genotyping */ // currentState.inf_AssemblyInput.open(par.inf_AssemblyInputFilename.c_str());
	/*bool GenotypingInputDefined = parameters[findParameter("-g",parameters)]->isSet();
   if (GenotypingInputDefined) {
      doGenotyping(currentState, FastaFile );
      exit(EXIT_SUCCESS);
   }
    
   bool AssemblyInputDefined = parameters[findParameter("-z",parameters)]->isSet();
   if (AssemblyInputDefined) {
      doAssembly(currentState, FastaFile );
      exit(EXIT_SUCCESS);
   }*/

    // If -q parameter given, search for mobile element insertions and quit.
    if (parameters[findParameter("-q", parameters)]->isSet()) {
        exit(searchMEImain(currentState, g_genome, userSettings));
    }

   /* Normal pindel functioning: search SVs*/

   // Get a new chromosome again and again until you have visited the specified chromosome or the file ends
   // CurrentChrName stores the name of the chromosome.

   bool SpecifiedChrVisited = false;

   do {
		timer.switchTo("Loading chromosomes");
		const Chromosome* currentChromosome = g_genome.getNextChromosome();
		if ( (currentChromosome == NULL ) || ( SpecifiedChrVisited==true )) {
			break;
		}
		if (!userSettings->loopOverAllChromosomes()) {
			if (currentChromosome->getName() == userSettings->getRegion()->getTargetChromosomeName() ) {
				SpecifiedChrVisited = true;
			}
			else {
         	*logStream << "Skipping chromosome: " << currentChromosome->getName() << std::endl;
				continue;
			}
		}

      *logStream << "Processing chromosome: " << currentChromosome->getName() << std::endl;

      g_maxPos = 0; // #################
      *logStream << "Chromosome Size: " << currentChromosome->getBiolSize() << std::endl;
      CurrentChrMask.resize(currentChromosome->getCompSize());

      for (unsigned int i = 0; i < currentChromosome->getCompSize(); i++) {
         CurrentChrMask[i] = 'N';
      }
      BoxSize = currentChromosome->getCompSize()/ 30000;
      if (BoxSize == 0) BoxSize = 1;
      unsigned NumBoxes = (unsigned) (currentChromosome->getCompSize() * 2 / BoxSize) + 1; // box size
      (*logStream << "NumBoxes: " << NumBoxes << "\tBoxSize: " << BoxSize << std::endl);

      g_binIndex = 0; // to start with 0... 
    
      LoopingSearchWindow currentWindow( userSettings->getRegion(), currentChromosome, WINDOW_SIZE ); 

      // loop over one chromosome
      do {
      	/* 3.2.1 preparation starts */
			*logStream << currentWindow.display();
			SearchWindow currentWindow_cs = currentWindow.makePindelCoordinateCopy(); // _cs means computer science coordinates
	     	//g_bdData.loadRegion( currentWindow_cs );
         if (Time_Load_S == 0) {
            Time_Load_S = time(NULL);
         }
			timer.switchTo("Reading in reads + matching close ends");//
         if (userSettings->bamFilesAsInput()) {
             get_RP_Reads_Discovery(currentState, currentWindow );
             g_bdData.UpdateBD(currentState);
				
         }
          g_bdData.loadRegion( currentWindow_cs );
          //std::cout << "g_bdData.size() " << g_bdData.GetBDSize() << std::endl;
         get_SR_Reads(currentState, currentWindow );
         Time_Mine_E = time(NULL);

         if (currentState.Reads_SR.size() ) {
            *logStream << "There are " << currentState.Reads_SR.size() << " reads for this chromosome region." << std::endl; // what region?

            if (userSettings->reportCloseMappedReads() ) {
				 	ReportCloseMappedReads( currentState.Reads_SR );       
            }
            Time_Load_E = time(NULL);
            if (!userSettings->reportOnlyCloseMappedReads) {
					timer.switchTo("Searching far ends");
					SearchFarEnds( currentChromosome->getSeq(), currentState.Reads_SR, *currentChromosome );
                    UpdateFarFragName(currentState.Reads_SR);
                    //std::cout << "before currentState.InterChromosome_SR size " << currentState.InterChromosome_SR.size() << std::endl;
                    if (userSettings->reportInterchromosomalEvents) {
                        for (unsigned index = 0; index < currentState.Reads_SR.size(); index++) {
                            if (currentState.Reads_SR[index].UP_Far.size()) {
                                //std::cout << currentState.Reads_SR[index].FragName << " " << currentState.Reads_SR[index].FarFragName << std::endl;
                                if (currentState.Reads_SR[index].FragName != currentState.Reads_SR[index].FarFragName)
                                    currentState.InterChromosome_SR.push_back(currentState.Reads_SR[index]);
                            }
                        }
                        //std::cout << "currentState.InterChromosome_SR.size() = " << currentState.InterChromosome_SR.size() << std::endl;
                    }

                    
					currentState.Reads_SR.insert( currentState.Reads_SR.end(), currentState.FutureReads_SR.begin(), currentState.FutureReads_SR.end() );
					currentState.FutureReads_SR.clear();
					timer.switchTo("Searching and reporting variations");
					SearchSVs( currentState, NumBoxes, currentWindow );
            }
            Time_Sort_E = time(NULL);

            AllLoadings += (unsigned int) difftime(Time_Load_E, Time_Load_S);
            AllSortReport += (unsigned int) difftime(Time_Sort_E, Time_Load_E);
            currentState.Reads_SR.clear();
            *logStream << "There are " << currentState.FutureReads_SR. size()  << " reads saved for the next cycle.\n" << std::endl;
            //currentState.Reads_SR.swap(currentState.FutureReads_SR);
         }
         else {
             *logStream << "There are no reads for this bin.\n";
         }
         Time_Load_S = 0;
			currentWindow.next();
         g_binIndex++;
      } while (!currentWindow.finished());

   } while (true);

    if (userSettings->reportInterchromosomalEvents)
       SortAndReportInterChromosomalEvents(currentState, g_genome, userSettings);
    
	timer.reportAll( *logStream );
   //*logStream << "Loading genome sequences and reads: " << AllLoadings << " seconds." << std::endl;
   //*logStream << "Mining, Sorting and output results: " << AllSortReport << " seconds." << std::endl;
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

// rotates a string back: KAI -> IKA 
void rotateBack( std::string& str )
{
	char lastChar = str[str.size()-1];
	str = lastChar + str.substr(0,str.size()-1);
}

void rotateForward( std::string& str )
{
	char firstChar = str[0];
	str = str.substr(1) + firstChar;
}

// the true principle here would be that you can move the insertion backwards, 
void GetRealStart4Insertion(const std::string & chromosomeSeq,
                            std::string & InsertedStr, unsigned int &RealStart,
                            unsigned int &RealEnd)
{

	unsigned int lastPosAfterInsertion_comp = RealEnd + g_SpacerBeforeAfter;
	//std::cout << "AInsertedStr: " << InsertedStr << ", start= " << RealStart  << chromosomeSeq[ g_SpacerBeforeAfter + RealStart] << ", end= " << RealEnd << chromosomeSeq[ g_SpacerBeforeAfter + RealEnd] << std::endl;
	//for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
	//std::cout << "\n";
	while ( chromosomeSeq[ lastPosAfterInsertion_comp ] == InsertedStr[ 0 ] ) {
		rotateForward( InsertedStr );
		lastPosAfterInsertion_comp++;
	}
	RealEnd = lastPosAfterInsertion_comp - g_SpacerBeforeAfter;
	//std::cout << "BInsertedStr: " << InsertedStr << ", start= " << RealStart << ", end= " << RealEnd << std::endl;
	//for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
	//std::cout << "\n";
	unsigned int lastPosBeforeInsertion_comp = lastPosAfterInsertion_comp-1;
	while ( chromosomeSeq[ lastPosBeforeInsertion_comp ] == InsertedStr[ InsertedStr.size()-1 ] ) {
		rotateBack( InsertedStr );
		lastPosBeforeInsertion_comp--;
	}
	RealStart = lastPosBeforeInsertion_comp - g_SpacerBeforeAfter;
	//std::cout << "CInsertedStr: " << InsertedStr << ", start= " << RealStart << ", end= " << RealEnd << std::endl;
	//for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
	//std::cout << "\n";

}
	//std::cout << "RS: " << RealStart << " RE: " << RealEnd << " IS: " <<  InsertedStr << "\n";
	//std::cout << "----\n";

   /*unsigned int IndelSize = InsertedStr.size();
   unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;
   unsigned int original_RealStart = RealStart;

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
    unsigned DIFF = RealStart - original_RealStart;
    InsertedStr = InsertedStr.substr(0, IndelSize - DIFF) + InsertedStr.substr(IndelSize, DIFF);*/
//}


/*void GetRealStart4Insertion(const std::string & TheInput,
                            std::string & InsertedStr, unsigned int &RealStart,
                            unsigned int &RealEnd)
{
    unsigned int IndelSize = InsertedStr.size();
    unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;
    unsigned int original_RealStart = RealStart;

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
    unsigned DIFF = RealStart - original_RealStart;
    InsertedStr = InsertedStr.substr(0, IndelSize - DIFF) + InsertedStr.substr(IndelSize, DIFF);
}*/

std::vector<Region> Merge(const std::vector<Region> &AllRegions)
{
    return AllRegions;
}

void GetCloseEndInner(const std::string & CurrentChrSeq, SPLIT_READ & Temp_One_Read)
{
    std::string CurrentReadSeq;
    //std::vector<unsigned int> PD[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
	std::vector<PosVector> PD;
	PosVector emptyPosVector;
   PD.assign( Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(), emptyPosVector);
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

        if (PD[0].size()) {
            CheckLeft_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
			}
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
		//std::cout << "Starting to fit the close end with character" << RightChar << "\n";
        if (RightChar != 'N') {
            for (int pos = Start; pos < End; pos++) {
                if (CurrentChrSeq[pos] == RightChar) {
                    PD[0].push_back(pos);
                }
            }
        }
        //std::cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl;
        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
        CheckRight_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
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


unsigned CountElements(const FarEndSearchPerRegion * OneRegionSearchResult_input, int levels) {
    unsigned Sum = 0;
    
    for (int LevelIndex = 0; LevelIndex < levels; LevelIndex++) {
        Sum += OneRegionSearchResult_input->PD_Plus[LevelIndex].size() + OneRegionSearchResult_input->PD_Minus[LevelIndex].size();
    }
    return Sum;
}

void ExtendMatch( const SPLIT_READ & read,
                 const std::string & readSeq,
                 const std::vector <FarEndSearchPerRegion*> & WholeGenomeSearchResult_input,
                 const short minimumLengthToReportMatch,
                 const short BP_End, const short CurrentLength,
                 SortedUniquePoints &UP )
{
	const char CurrentChar = readSeq[CurrentLength];
   const char CurrentCharRC = Convert2RC4N[(short) CurrentChar];
   bool AllEmpty = true;
   std::vector <FarEndSearchPerRegion*> WholeGenomeSearchResult_output;
   for (unsigned IndexOfRegion = 0; IndexOfRegion < WholeGenomeSearchResult_input.size(); IndexOfRegion++) {
      const FarEndSearchPerRegion* CurrentRegion_input = WholeGenomeSearchResult_input[IndexOfRegion];
      unsigned int Max_size = 0;
      for (int CheckedIndex = 0; CheckedIndex < read.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
         if (Max_size < CurrentRegion_input->PD_Plus[CheckedIndex].size()) Max_size = CurrentRegion_input->PD_Plus[CheckedIndex].size();
         if (Max_size < CurrentRegion_input->PD_Minus[CheckedIndex].size()) Max_size = CurrentRegion_input->PD_Minus[CheckedIndex].size();
      }
      const std::string & chromosomeSeq = CurrentRegion_input->CurrentChromosome->getSeq();
      FarEndSearchPerRegion* CurrentRegion_output=new FarEndSearchPerRegion (CurrentRegion_input->CurrentChromosome, read.getTOTAL_SNP_ERROR_CHECKED(), Max_size);
      for (int i = 0; i <= read.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) {
         CategorizePositions( CurrentChar, chromosomeSeq, CurrentRegion_input->PD_Plus, CurrentRegion_output->PD_Plus, i, 1, read.getTOTAL_SNP_ERROR_CHECKED_Minus() );
         CategorizePositions( CurrentCharRC, chromosomeSeq, CurrentRegion_input->PD_Minus, CurrentRegion_output->PD_Minus, i, -1, read.getTOTAL_SNP_ERROR_CHECKED_Minus() );
      }
      if (CountElements(CurrentRegion_output, read.getTOTAL_SNP_ERROR_CHECKED())) {
         AllEmpty = false;
         WholeGenomeSearchResult_output.push_back(CurrentRegion_output);
      }
   }
	/*std::cout << "Matching " << CurrentLength << " length and char " << CurrentChar << ", " << CurrentCharRC << std::endl;
	for (unsigned IndexOfRegion = 0; IndexOfRegion < WholeGenomeSearchResult_input.size(); IndexOfRegion++) {
		std::cout << "Region index: " << IndexOfRegion << "\n";
		for (int CheckedIndex = 0; CheckedIndex < read.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
			std::cout << "R["<< CheckedIndex << "]=" << WholeGenomeSearchResult_input[IndexOfRegion].PD_Plus[CheckedIndex].size() << "/" << WholeGenomeSearchResult_input[IndexOfRegion].PD_Minus[CheckedIndex].size() << std::endl;
		}
	}
	std::cout << "UP-size:" << UP.size() << "\n";*/
	// this loop looks familiar; candidate for factoring out mini-function?
   if (AllEmpty == false ) {
      const short CurrentLengthOutput = CurrentLength + 1;
      CheckBoth(read, readSeq, WholeGenomeSearchResult_output, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);
   }
   else {
 	} // else-if Sum
 	for (unsigned int i=0; i<WholeGenomeSearchResult_output.size();i++) {
		delete WholeGenomeSearchResult_output[ i ];
   }
}

unsigned int minimumNumberOfMismatches( const std::vector <FarEndSearchPerRegion*> & WholeGenomeSearchResult_input, const unsigned int maxNumberMismatches )
{
	unsigned int Sum=0; 
	unsigned int numberOfMismatches=0;
	for (;numberOfMismatches<=maxNumberMismatches; numberOfMismatches++ ) {
		for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
         Sum += WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[numberOfMismatches].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[numberOfMismatches].size();
      }
		if ( Sum != 0 ) { break; }
	}
	return numberOfMismatches;
}

void CheckBoth(const SPLIT_READ & read,
               const std::string & readSeq,
               const std::vector <FarEndSearchPerRegion*> & WholeGenomeSearchResult_input,
               const short minimumLengthToReportMatch,
               const short BP_End,
               const short CurrentLength,
               SortedUniquePoints &UP)
{
	UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
   unsigned NumberOfMatchPositionsWithLessMismatches = 0;
   int Sum = 0;
	
	if (CurrentLength >= minimumLengthToReportMatch && CurrentLength <= BP_End) {
		if (minimumNumberOfMismatches( WholeGenomeSearchResult_input, read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
			return; 
		}
      for (short numberOfMismatches = 0; numberOfMismatches <= read.getMAX_SNP_ERROR(); numberOfMismatches++) {
         if (NumberOfMatchPositionsWithLessMismatches) break;
			
         Sum = 0;
         for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
            Sum += WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[numberOfMismatches].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[numberOfMismatches].size();
         }
         NumberOfMatchPositionsWithLessMismatches = Sum;
         if (Sum == 1 && CurrentLength >= minimumLengthToReportMatch + numberOfMismatches) {
            Sum = 0;
            if (userSettings->ADDITIONAL_MISMATCH > 0) { // what if this is ADDITIONAL_MISMATCH is 0? Do you save anything then?
					
					// what if j +ADD_MISMATCH exceeds the max allowed number of mismatches (the PD_element does not exist?)
                    
					// feeling the need to comment here - so factor out?
					// only report reads if there are no reads with fewer mismatches or only one or two more mismatches
               unsigned int regionWithMatch = 0;
               for (short mismatchCount = 0; mismatchCount <= numberOfMismatches + userSettings->ADDITIONAL_MISMATCH; mismatchCount++) {
                  for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
                     unsigned int hitsInRegion= WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[mismatchCount].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[mismatchCount].size();
						 	Sum += hitsInRegion;
							if (hitsInRegion>0) {
								regionWithMatch = RegionIndex;
							}
                  }
               }
				/*if (read.Name=="@read_6990/2" ) {
				std::cout << "In CFE: CurrentLength = " << CurrentLength << ", mismatch count = " << numberOfMismatches << ", maxMismatch = " << g_maxMismatch[CurrentLength] << std::endl;
				for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) { 
					std::cout << "Region " << RegionIndex << " is " <<  WholeGenomeSearchResult_input[ RegionIndex ].CurrentChromosome->getName() << ":" << WholeGenomeSearchResult_input[ RegionIndex ].PD_Plus[0].size() << "-" <<
								WholeGenomeSearchResult_input[ RegionIndex ].PD_Minus[0].size()<< "\n";
				}
				for (short k=0;k<=read.getMAX_SNP_ERROR(); k++) {
					std::cout << k << "\t" << WholeGenomeSearchResult_input[0].PD_Plus[k].size() + WholeGenomeSearchResult_input[0].PD_Minus[k].size() << "\n";
				}}*/
             if (Sum == 1 && (unsigned)numberOfMismatches <= g_maxMismatch[CurrentLength] ) {
                        // why I love constructors
								UniquePoint MatchPosition;
								const FarEndSearchPerRegion* hitRegion = WholeGenomeSearchResult_input[ regionWithMatch ];
                        if (WholeGenomeSearchResult_input[ regionWithMatch ]->PD_Plus[numberOfMismatches].size() == 1) {
									UniquePoint PlusMatch( hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Plus[numberOfMismatches][0], FORWARD, SENSE, numberOfMismatches );
									MatchPosition = PlusMatch;
                        }
                        else {
									UniquePoint MinMatch(  hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Minus[numberOfMismatches][0], BACKWARD, ANTISENSE, numberOfMismatches );
									MatchPosition = MinMatch;
                        }

                        if (CheckMismatches(WholeGenomeSearchResult_input[regionWithMatch]->CurrentChromosome->getSeq(), readSeq, MatchPosition)) {
                            UP.push_back (MatchPosition);
                            break;
                        } // if CheckMismatches
                    } // if Sum==1
                } // if AdditionalMismatches
        		} // if sumsize ==1
        } // for-loop
	} // if length of match is sufficient to be reportable
    
   if (CurrentLength < BP_End) {
		ExtendMatch( read, readSeq, WholeGenomeSearchResult_input, minimumLengthToReportMatch, BP_End, CurrentLength, UP );
	}
}

/*void ExtendMatch( const SPLIT_READ & read, const std::string & chromosomeSeq,
               const std::string & readSeq,
               const std::vector<unsigned int> PD_Plus[],
               const std::vector<unsigned int> PD_Minus[], const short minimumLengthToReportMatch,
               const short BP_End, const short CurrentLength,
               SortedUniquePoints &UP )
{
	std::vector<unsigned int> PD_Plus_Output[read.getTOTAL_SNP_ERROR_CHECKED()];
   std::vector<unsigned int> PD_Minus_Output[read.getTOTAL_SNP_ERROR_CHECKED()];
		
   for (int CheckedIndex = 0; CheckedIndex < read.getTOTAL_SNP_ERROR_CHECKED(); CheckedIndex++) {
      PD_Plus_Output[CheckedIndex].reserve( PD_Plus[CheckedIndex].size()); // this assumes perfect matches and no 'attrition' from higher levels. We may want to test this...
      PD_Minus_Output[CheckedIndex].reserve( PD_Minus[CheckedIndex].size());
   }
   const char CurrentChar = readSeq[CurrentLength];
   const char CurrentCharRC = Convert2RC4N[(short) CurrentChar];

	for (int i = 0; i <= read.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) {
		CategorizePositions( CurrentChar, chromosomeSeq, PD_Plus, PD_Plus_Output, i, 1, read.getTOTAL_SNP_ERROR_CHECKED_Minus() );
		CategorizePositions( CurrentCharRC, chromosomeSeq, PD_Minus, PD_Minus_Output, i, -1, read.getTOTAL_SNP_ERROR_CHECKED_Minus() );
   }

	// this loop looks familiar; candidate for factoring out mini-function?
   unsigned int Sum = 0;
   for (int i = 0; i <= read.getMAX_SNP_ERROR(); i++) {
      Sum += PD_Plus_Output[i].size() + PD_Minus_Output[i].size();
   }
   if (Sum) {
      const short CurrentLengthOutput = CurrentLength + 1;
      CheckBoth(read, chromosomeSeq, readSeq, PD_Plus_Output, PD_Minus_Output, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);
   }
   else {
     return;
 	} // else-if Sum
}

void CheckBoth(const SPLIT_READ & read, const std::string & chromosomeSeq,
               const std::string & readSeq,
               const std::vector<unsigned int> PD_Plus[],
               const std::vector<unsigned int> PD_Minus[], const short minimumLengthToReportMatch,
               const short BP_End, const short CurrentLength,
               SortedUniquePoints &UP)
{   
	UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   if (CurrentLength >= minimumLengthToReportMatch && CurrentLength <= BP_End) {
      for (short numberOfMismatches = 0; numberOfMismatches <= read.getMAX_SNP_ERROR(); numberOfMismatches++) {
         if (PD_Plus[numberOfMismatches].size() + PD_Minus[numberOfMismatches].size() == 1 && CurrentLength >= minimumLengthToReportMatch + numberOfMismatches) {
            int Sum = 0;
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
                  if (CheckMismatches(chromosomeSeq,readSeq, MatchPosition)) {
                     UP.push_back (MatchPosition);
                     break;
                  } // if CheckMismatches
               } // if Sum==1
            } // if AdditionalMismatches
        	} // if sumsize ==1
      } // for-loop
	} // if length of match is sufficient to be reportable

   if (CurrentLength < BP_End) {
		ExtendMatch( read, chromosomeSeq, readSeq, PD_Plus, PD_Minus, minimumLengthToReportMatch, BP_End, CurrentLength, UP );
	}
}*/


void CleanUniquePoints(SortedUniquePoints &Input_UP)
{
    SortedUniquePoints TempUP; //vector <UniquePoint> UP_Close; UP_Far
    UniquePoint LastUP = Input_UP[Input_UP.size() - 1];
    char LastDirection = LastUP.Direction;
    char LastStrand = LastUP.Strand;
    const Chromosome * LastChr = LastUP.chromosome_p;
    //LastChr == .chromosome_p
    unsigned int Terminal;

    if (LastDirection == FORWARD) {
        Terminal = LastUP.AbsLoc - LastUP.LengthStr;
        for (unsigned i = 0; i < Input_UP.size(); i++) {
           if (LastChr != Input_UP[i].chromosome_p) continue;
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
          if (LastChr != Input_UP[i].chromosome_p) continue; 
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


