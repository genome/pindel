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
#include <map>
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
//#include "genotyping.h"
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

const std::string Pindel_Version_str = "Pindel version 0.2.5b9, 20160729.";

const Chromosome g_dummyChromosome("","");
Genome g_genome;
std::ofstream g_logFile;


std::vector <ChrNameAndSizeAndIndex> g_ChrNameAndSizeAndIndex;

int g_binIndex = -1; // global variable for the bin index, as I cannot easily pass an extra parameter to the diverse functions
unsigned int g_maxPos = 0; // to calculate which is the last position in the chromosome, and hence to calculate the number of bins
short g_MinClose = 8;
std::set<std::string> g_sampleNames;
std::map<std::string,unsigned> g_SampleName2Index;
std::map<std::string,unsigned> g_ChrName2Ploidy;
std::vector <RefCoveragePerPosition> g_RefCoverageRegion;
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
unsigned g_RegionStart, g_RegionEnd;
char Match[256];
char Match2N[256];
char Convert2RC[256];
char Convert2RC4N[256];
char Cap2LowArray[256];
bool FirstChr = true;
unsigned int DSizeArray[15];
int g_maxInsertSize=0;
unsigned g_NumberOfGapAlignedReads = 0;
std::string CurrentChrMask;
std::vector<Parameter *> parameters;

UserDefinedSettings* userSettings;


// #########################################################

int NumRead2ReportCutOff_BP = 2;
const int g_MAX_RANGE_INDEX = 9; // 5 or 6 or 7 or maximum 8      //# // user
unsigned int WINDOW_SIZE = 10000000;
const unsigned int AROUND_REGION_BUFFER = 10000; // how much earlier reads should be selected if only a region of the chromosome needs be specified.
// #########################################################

const short MaxDI = 30;


// Note: in case one needs to handle differnt line delimuters (like crlf)
// from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
void safeGetline(std::istream& is, std::string& t)
{
   t.clear();
   getline( is, t );
}


void SPLIT_READ::setUnmatchedSeq( const std::string & unmatchedSeq )
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

   //UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   MAX_SNP_ERROR = g_maxMismatch[ ReadLength ];
   TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + userSettings->ADDITIONAL_MISMATCH;
   TOTAL_SNP_ERROR_CHECKED = TOTAL_SNP_ERROR_CHECKED_Minus + 1;
}

const Chromosome* Genome::addChromosome( Chromosome* newChromosome )
{
   for (unsigned int i=0; i<m_chromosomes.size(); i++ ) {
      if ( m_chromosomes[ i ]->getName() == newChromosome->getName() ) {
         delete m_chromosomes[ i ];
         std::cout << "delete one chromsome." << std::endl;
         m_chromosomes[ i ] = newChromosome;
         return newChromosome;
      }
   }
   m_chromosomes.push_back( newChromosome );
   return newChromosome;
}

const Chromosome* Genome::getChr( unsigned int index ) const
{
   if (index<m_chromosomes.size()) {
      return m_chromosomes[ index ];
   } else {
      return NULL;
   }
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

short Genome::getChrID( const std::string& chromosomeName )
{
   for (unsigned int i=0; i<m_chromosomes.size(); i++ ) {
      if ( m_chromosomes[ i ]->getName() == chromosomeName ) {
         return i;
      }
   }
   return -1;
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

void Genome::reset()
{
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
      //std::cout << "m_chromosomes.size() " << m_chromosomes.size() << " " << m_chromosomes[m_chromosomes.size() - 1]->getName() << std::endl;
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
      } else {
         return NULL;
      }
   } else {
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
      } else {
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

LoopingSearchWindow::LoopingSearchWindow(const SearchRegion* region, const Chromosome* chromosome, const int binSize) :
   SearchWindow(chromosome,0,chromosome->getBiolSize()), m_BIN_SIZE( binSize )
{

   bool noOverlapWithChromosome = false;
   if (region->isStartDefined()) {
      m_officialStart = region->getStart();
      if (m_officialStart > chromosome->getBiolSize()) {
         noOverlapWithChromosome = true;
      }
      // if the user defines a region, you need to start with reads before that, but not before the start of the chromosome
      m_globalStart = ( region->getStart() >= AROUND_REGION_BUFFER ? region->getStart() - AROUND_REGION_BUFFER : 0);
   } else {
      m_officialStart = m_globalStart = 0;
   }

   if (region->isEndDefined()) {
      m_officialEnd = region->getEnd();
      if (m_officialEnd < 0 ) {
         noOverlapWithChromosome = true;
      }
      // if the user defines a region, you need to end with reads after that, but not after the end of the chromosome
      m_globalEnd = std::min( chromosome->getBiolSize() , region->getEnd() + AROUND_REGION_BUFFER );
   } else {
      m_officialEnd = m_globalEnd = chromosome->getBiolSize();
   }

   if (noOverlapWithChromosome) {
      std::cout << "Error: the region to scan (" << m_officialStart << ", " << m_officialEnd << ") does not overlap with the "
                << "chromosome (positions 0 to " << chromosome->getBiolSize() << std::endl;
      exit( EXIT_FAILURE );
   }

   m_currentStart = m_globalStart;
   m_displayedStart = m_officialStart;
   updateEndPositions();

}

LoopingSearchWindow::LoopingSearchWindow(const SearchRegion* region, const Chromosome* chromosome, const int binSize, const unsigned Bed_start, const unsigned Bed_end ) :
   SearchWindow(chromosome,0,chromosome->getBiolSize()), m_BIN_SIZE( binSize )
{
   m_officialStart = Bed_start;
   m_globalStart = ( Bed_start >= AROUND_REGION_BUFFER ? Bed_start - AROUND_REGION_BUFFER : 0);

   m_officialEnd = Bed_end;
   m_globalEnd = std::min( chromosome->getBiolSize() , Bed_end + AROUND_REGION_BUFFER );

   m_currentStart = m_globalStart;
   m_displayedStart = m_officialStart;
   updateEndPositions();
   /*
       bool noOverlapWithChromosome = false;
       if (region->isStartDefined()) {
           m_officialStart = region->getStart();
           if (m_officialStart > chromosome->getBiolSize()) {
               noOverlapWithChromosome = true;
           }
           // if the user defines a region, you need to start with reads before that, but not before the start of the chromosome
           m_globalStart = ( region->getStart() >= AROUND_REGION_BUFFER ? region->getStart() - AROUND_REGION_BUFFER : 0);
       }
       else {
           m_officialStart = m_globalStart = 0;
       }

       if (region->isEndDefined()) {
           m_officialEnd = region->getEnd();
           if (m_officialEnd < 0 ) {
               noOverlapWithChromosome = true;
           }
           // if the user defines a region, you need to end with reads after that, but not after the end of the chromosome
           m_globalEnd = std::min( chromosome->getBiolSize() , region->getEnd() + AROUND_REGION_BUFFER );
       }
       else {
           m_officialEnd = m_globalEnd = chromosome->getBiolSize();
       }

       if (noOverlapWithChromosome) {
           std::cout << "Error: the region to scan (" << m_officialStart << ", " << m_officialEnd << ") does not overlap with the "
           << "chromosome (positions 0 to " << chromosome->getBiolSize() << std::endl;
           exit( EXIT_FAILURE );
       }

       m_currentStart = m_globalStart;
       m_displayedStart = m_officialStart;
       updateEndPositions();
   */
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
      ss << "\nLooking at chromosome " << m_chromosome->getName() << " bases " << m_displayedStart << " to " << m_displayedEnd << " of the bed region: chromosome " << m_chromosome->getName() << ":" << m_officialStart << "-" << m_officialEnd << " \n";
   } else {
      ss << "Checking out reads near the borders of the specified regions for extra evidence.\n";
   }
   return ss.str();
}

bool LoopingSearchWindow::finished() const
{
   // ugly hack for speed purposes when using Pindel-formatted input
   if (userSettings->pindelFilesAsInput() &&  m_currentStart >= g_maxPos ) {
      return true;
   }
   return ( m_currentStart > m_globalEnd );
}

unsigned int SPLIT_READ::getLastAbsLocCloseEnd() const
{
   return UP_Close[ UP_Close.size()-1 ].AbsLoc;
}

bool SPLIT_READ::goodFarEndFound() const
{
   return ((UP_Far.MaxLen() + UP_Close.MaxLen() >= UnmatchedSeq.size()) );
}

bool SPLIT_READ::hasCloseEnd() const
{
   return !UP_Close.empty();
}

unsigned int SortedUniquePoints::MaxLen() const
{
   if (m_positions.size() == 0 ) {
      return 0;
   } else {
      int lastElementIndex = m_positions.size()-1;
      return m_positions[ lastElementIndex ].LengthStr;
   }
}

unsigned int SortedUniquePoints::NumMismatch() const
{
   if (m_positions.size() == 0 ) {
      return 1000;
   } else {
      int lastElementIndex = m_positions.size()-1;
      return m_positions[ lastElementIndex ].Mismatches;
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
   return ( userSettings->breakdancerOutputFilename != "" && parameters[findParameter("-b",parameters)]->isSet());
}

void outputBreakDancerEvent( const std::string& chromosomeName, const int leftPosition, const int rightPosition,
                             const int svSize, const std::string& svType, const int svCounter)
{
   std::ofstream breakDancerOutputFile(userSettings->breakdancerOutputFilename.c_str(), std::ios::app);
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
   Region()
   {
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

/** 'convertToXBaiFileName' converts the name of the BAM file "bamFileName", for example x.bam, into x.bai
    (or at least returns "x.bai" as a result" **/
//convertToXBaiFilename
std::string convertToXBaiFilename( const std::string& bamFileName )
{
   std::string outputName = bamFileName;
   int nameLength = outputName.length();
   outputName[ nameLength - 1 ] = 'i'; // "x.bam" -> "x.bai"
   return outputName;
}

void readBamConfigFile(std::string& bamConfigFilename, ControlState& currentState )
{
   int sampleCounter=0;
   std::ifstream BamConfigFile( bamConfigFilename.c_str() );
   if (BamConfigFile) {
      while (BamConfigFile.good()) {
         bam_info tempBamInfo;
         BamConfigFile >> tempBamInfo.BamFile >> tempBamInfo.InsertSize;
         if (!BamConfigFile.good()) {
            break;
         }
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
         if (! fileExists( tempBamInfo.BamFile+".bai" ) && ! fileExists( convertToXBaiFilename( tempBamInfo.BamFile ))) {
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
   } else {
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
         if (!pindelConfigFile.good()) {
            break;
         }

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
   } else {
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

   if (len > 3 && suffix == ".gz") {
      return new GZLineReader(filename);
   } else {
      return new IfstreamLineReader(filename);
   }
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
            g_maxMismatch[ length ] = numberOfErrors + 1;
            //g_maxMisMatch.push_back( numberOfErrors );
            //std::cout << length << " bases has max errors \t" << g_maxMismatch[length] << "\n";
            break; // break out of this length, up to the next
         }
      }
   }
   g_maxMismatch[ 0 ] = 0;
   g_maxMismatch[ 1 ] = 0;
   g_maxMismatch[ 2 ] = 0;
   g_maxMismatch[ 3 ] = 0;
}



void init(int argc, char *argv[], ControlState& currentState )
{
   //std::cout << "1" << std::endl;
   //UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
   logStream=&std::cout;
//std::cout << "2" << std::endl;
   //if (userSettings->NumRead2ReportCutOff == 1) {
   //    userSettings->BalanceCutoff = 300000000;
   //}
//std::cout << "2" << std::endl;
   // define all the parameters you have
   defineParameters( parameters );
//std::cout << "3" << std::endl;
   // now read the parameters from the command line
   readParameters(argc, argv, parameters);
//std::cout << "4" << std::endl;
   if (userSettings->logFilename != "" ) {
      g_logFile.open( userSettings->logFilename.c_str() );
      logStream = &g_logFile;
   }
//std::cout << "5" << std::endl;
   *logStream << Pindel_Version_str << std::endl;

   if (argc <= 1) { // the user has not given any parameters
      printHelp( parameters );
      exit ( EXIT_FAILURE);
   }
//std::cout << "6" << std::endl;
   // check parameters
   if (!checkParameters( parameters )) {
      exit ( EXIT_FAILURE);
   }
//std::cout << "7" << std::endl;
   createProbTable(0.001+userSettings->Seq_Error_Rate, userSettings->sensitivity);
//std::cout << "7a" << std::endl;
   std::string fastaFilename( userSettings->referenceFilename.c_str() );
   std::cout << "Loading reference genome ..." << std::endl;
   //if (userSettings->reportInterchromosomalEvents) {
   //std::cout << "7b" << std::endl;
   g_genome.loadAll( fastaFilename );
   //std::cout << "7c" << std::endl;
   //}
   //else {
   //std::cout << "7d" << std::endl;
   //	g_genome.load( fastaFilename );
   //std::cout << "7e" << std::endl;
   //}
   std::cout << "Loading reference genome done." << std::endl;
//std::cout << "8" << std::endl;
   bool BreakDancerDefined = parameters[findParameter("-b",parameters)]->isSet();
   if (BreakDancerDefined) {
      g_bdData.loadBDFile(userSettings->breakdancerFilename);
   }

//std::cout << "9" << std::endl;
   if (userSettings->FLOAT_WINDOW_SIZE > 5000.0) {
      LOG_ERROR(*logStream << "Window size of " << userSettings->FLOAT_WINDOW_SIZE << " million bases is too large" << std::endl);
      exit ( EXIT_FAILURE);
   } else if (userSettings->FLOAT_WINDOW_SIZE > 100.0) {
      LOG_ERROR(*logStream << "Window size of " << userSettings->FLOAT_WINDOW_SIZE
                << " million bases is rather large; this may produce bad::allocs or segmentation faults. If that happens, either try to reduce the window size or deactivate the searching for breakpoints and long insertions by adding the command-line options \"-l false -k false\"." << std::endl);
   }
   WINDOW_SIZE = (unsigned int)(1000000 * userSettings->FLOAT_WINDOW_SIZE);

//std::cout << "10" << std::endl;

   // if all parameters are okay, open the files

   if (userSettings->singlePindelFileAsInput()) {
      currentState.lineReader= getLineReaderByFilename(userSettings->pindelFilename.c_str());
      currentState.inf_Pindel_Reads = new PindelReadReader(*currentState.lineReader);
   }
//std::cout << "11" << std::endl;
   if (userSettings->pindelConfigFileAsInput()) {
      readPindelConfigFile( userSettings->pindelConfigFilename, currentState.pindelfilesToParse );
   }

   if (userSettings->bamFilesAsInput()) {
      readBamConfigFile( userSettings->bamConfigFilename, currentState );
   }

//std::cout << "12" << std::endl;

   bool AssemblyInputDefined = parameters[findParameter("-z",parameters)]->isSet();
   if (AssemblyInputDefined) {
      currentState.inf_AssemblyInput.open(userSettings->inf_AssemblyInputFilename.c_str());
   }

   bool GenotypingInputDefined = parameters[findParameter("-g",parameters)]->isSet();
   if (GenotypingInputDefined) {
      currentState.inf_GenotypingInput.open(userSettings->inf_GenotypingInputFilename.c_str());
   }
//std::cout << "13" << std::endl;

   omp_set_num_threads(userSettings->numThreads);
   g_MinClose = userSettings->minClose;
   //std::cout << "minClose = " << g_MinClose << std::endl;
//std::cout << "14" << std::endl;
   if (userSettings->MaxRangeIndex > g_MAX_RANGE_INDEX) {
      LOG_ERROR(*logStream
                << "Maximal range index (-x) exceeds the allowed value (" << g_MAX_RANGE_INDEX << ") - resetting to " << g_MAX_RANGE_INDEX << ".\n" );
      userSettings->MaxRangeIndex = g_MAX_RANGE_INDEX;
   }
//std::cout << "15" << std::endl;
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

//std::cout << "16" << std::endl;

   if ( userSettings->breakdancerOutputFilename != "" ) {
      TestFileForOutput( userSettings->breakdancerOutputFilename );
   }
//std::cout << "17" << std::endl;

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

//std::cout << "18" << std::endl;

   std::string Spacer = "";
   for (unsigned int i = 0; i < g_SpacerBeforeAfter; i++) {
      Spacer += "N";
   }
//std::cout << "19" << std::endl;
   //Distance = 300;


   DSizeArray[0] = 0;
   DSizeArray[1] = 128;
   for (int dIndex=2; dIndex<15; dIndex++ ) {
      DSizeArray[ dIndex ] = DSizeArray[ dIndex-1 ] * 4;
   }

   //if (!userSettings->getRegion()->isTargetChromosomeDefined() && AssemblyInputDefined == false && GenotypingInputDefined == false) {
   //    *logStream << "Looping over all chromosomes." << std::endl;
   //}
//std::cout << "20" << std::endl;
}


void SearchFarEnd( const std::string& chromosome, SPLIT_READ& read, const Chromosome& currentChromosome)
{
   //std::cout << "entering SearchFarEnd" << std::endl;
   const int START_SEARCH_SPAN = 64;
   //std::cout << "getting searchCluster" << std::endl;
   const std::vector< SearchWindow>& searchCluster =  g_bdData.getCorrespondingSearchWindowCluster( read );
   //std::cout << "searchCluster size: " << searchCluster.size() << std::endl;
   if (searchCluster.size()!=0) {
      //std::cout << "Breakdancer input is not empty " << searchCluster.size() << std::endl;
      //for (unsigned index = 0; index < searchCluster.size(); index++)
      //	searchCluster[index].display();
      SearchFarEndAtPos( chromosome, read, searchCluster); // SearchFarEndAtPos
      //std::cout << "finished" << std::endl;
      if (read.goodFarEndFound()) {
         //read.Investigated = true;
         //std::cout << "return" << std::endl;
         return;
      }
      //else SearchFarEndAtPosPerfect( chromosome, read, searchCluster);
   }
   //std::cout << "SearchFarEnd	2" << std::endl;
   //UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   // if breakdancer does not find the event, or not find an event we trust, we turn to regular pattern matching
   unsigned int searchSpan=START_SEARCH_SPAN;
   unsigned int centerOfSearch = read.getLastAbsLocCloseEnd();
   //std::cout << "SearchFarEnd	3" << std::endl;
   for (int rangeIndex=1; rangeIndex<=userSettings->MaxRangeIndex + 1; rangeIndex++ ) {
      //std::cout << "rangeIndex " << rangeIndex << "\t";
      // note: searching the central range again and again may seem overkill, but since Pindel at the moment wants an unique best point, you can't skip the middle part
      // may be stuff for future changes/refactorings though
      std::vector< SearchWindow > aroundCESearchCluster;
      unsigned int Start, End;
      if (centerOfSearch > searchSpan+g_SpacerBeforeAfter) {
         Start = centerOfSearch - searchSpan;
      } else {
         Start = g_SpacerBeforeAfter;
      }
      if (centerOfSearch + searchSpan + g_SpacerBeforeAfter < chromosome.size()) {
         End = centerOfSearch + searchSpan;
      } else {
         End = chromosome.size() - g_SpacerBeforeAfter;
      }
      //std::cout << "Start is (abs) " << Start << " (rel): " << Start - g_SpacerBeforeAfter << "\n";
      //std::cout << "End is (abs) " << End << " (rel): " << End - g_SpacerBeforeAfter << " Size chrom= " << chromosome.size() << "\n";
      //std::cout << Start << "FirstStart " << centerOfSearch+searchSpan-10000000 << "<COS" << centerOfSearch-10000000 << " span " << searchSpan<< std::endl;
      SearchWindow regularWindow( &currentChromosome, Start, End );
      //if (rangeIndex == 1 && searchCluster.size()) {
      //	std::cout << rangeIndex << "\t";
      //regularWindow.display();
      //}

      //std::cout << Start << " " << centerOfSearch+searchSpan-10000000 << std::endl;
      aroundCESearchCluster.clear();
      aroundCESearchCluster.push_back( regularWindow );
      //std::cout << rangeIndex << "\tSearchFarEndAtPos" << std::endl;
      SearchFarEndAtPos( chromosome, read, aroundCESearchCluster ); // SearchFarEndAtPosPerfect
      //std::cout << "end\tSearchFarEndAtPos" << std::endl;
      if (read.goodFarEndFound()) {
         //read.Investigated = true;
         return;
      }
      //else SearchFarEndAtPosPerfect( chromosome, read, searchCluster);
      if (read.goodFarEndFound()) {
         //read.Investigated = true;
         //std::cout << "SearchFarEnd	found ###################3" << std::endl;
         return;
      }
      searchSpan *= 4;
   }
   //std::cout << std::endl;
   //read.Investigated = true;
   //std::cout << "leaving SearchFarEnd" << std::endl;
}

void ReportCloseMappedReads( const std::vector<SPLIT_READ>& reads )
{
   std::ofstream CloseEndMappedOutput( userSettings->getCloseEndOutputFilename().c_str(), std::ios::app);
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

void SearchFarEnds( const std::string & chromosomeSeq, std::vector<SPLIT_READ>& reads, const Chromosome& currentChromosome)
{
   //std::cout << "report per 1k reads, not sequential due to openmp" << std::endl;
   #pragma omp parallel default(shared)
   {
      #pragma omp for
      for (int readIndex= 0; readIndex < (int)reads.size(); readIndex++ ) {
         //std::cout << "readIndex: " << readIndex << std::endl;
         //if (readIndex % 1000 == 0)
         //	std::cout << "readIndex: " << readIndex << std::endl;
         //	std::cout << "readIndex: " << readIndex << "\t" << reads[readIndex].Name << "\t"
         //		<< reads[readIndex].FragName << "\t" <<  reads[readIndex].MatchedD << "\t"
         //		<< reads[readIndex].MatchedRelPos << "\t" << reads[readIndex].MS << "\t"
         //		<< reads[readIndex].UnmatchedSeq << std::endl;
         if (reads[readIndex].MapperSplit == false) {
            SearchFarEnd( chromosomeSeq, reads[readIndex], currentChromosome );
         }
         //else std::cout << "skip far end search" << std::endl;
         //std::cout << reads[readIndex] << std::endl;
      }
   }

   *logStream << "Far end searching completed for this window." << std::endl;
}


void SearchSVs(ControlState& currentState, const int NumBoxes, const SearchWindow& currentWindow )
{
   //std::cout << "1" << std::endl;
   //UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   SearchDeletions searchD;
   searchD.Search(g_bdData, currentState, NumBoxes, currentWindow);
//std::cout << "2" << std::endl;
   searchIndels(currentState, NumBoxes, currentWindow);
//std::cout << "3" << std::endl;
   if (userSettings->Analyze_TD) {
      searchTandemDuplications(currentState, NumBoxes, currentWindow);
      searchTandemDuplicationsNT(currentState, NumBoxes, currentWindow);
   }
//std::cout << "4" << std::endl;
   if (userSettings->Analyze_INV) {
      searchInversions(currentState, NumBoxes, currentWindow);
      searchInversionsNT(currentState, NumBoxes, currentWindow);
   }
//std::cout << "5" << std::endl;

   SearchShortInsertions searchSI;
   searchSI.Search(g_bdData, currentState, NumBoxes, currentWindow);
//std::cout << "6" << std::endl;
   ReportCloseAndFarEndCounts( currentState.Reads_SR );
//std::cout << "7" << std::endl;
   if (userSettings->Analyze_LI) {
      SortOutputLI(currentState, currentWindow.getChromosome()->getSeq(), currentState.Reads_SR, currentWindow, userSettings->getLIOutputFilename());
   }
//std::cout << "8" << std::endl;
   //if (userSettings->Analyze_BP) {
   //SortOutputRest(currentState, currentWindow.getChromosome()->getSeq(), currentState.Reads_SR, currentWindow, userSettings->getBPOutputFilename());
   //}
//std::cout << "9" << std::endl;
}

class TimerItem
{

public:
   TimerItem( const std::string& id );
   void stop();
   void restart();
   const std::string getReport() const;
   const std::string& getId() const
   {
      return m_id;
   }

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


class Timer
{

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

void UpdateFarFragName(std::vector <SPLIT_READ> & input)
{
   for (unsigned index = 0; index < input.size(); index++) {
      if (input[index].UP_Far.size()) {
         input[index].FarFragName = input[index].UP_Far[0].chromosome_p->getName();
         input[index].MatchedFarD = input[index].UP_Far[0].Strand;
      }
   }
}

short UpdateRefReadCoverage(ControlState& currentState, const SearchWindow& currentWindow )
{
   //std::set<std::string> g_sampleNames;
   //std::map<std::string,unsigned> g_SampleName2Index;
   std::cout << "There are " << g_sampleNames.size() << " samples." << std::endl;
   g_SampleName2Index.clear();
   std::set<std::string>::iterator it;
   unsigned index = 0;
   for (it=g_sampleNames.begin(); it != g_sampleNames.end(); it++ ) {
      g_SampleName2Index.insert ( std::pair<std::string, unsigned>(*it, index) );
      index++;
   }
   std::cout << "SampleName2Index done" << std::endl;
   // std::vector <RefCoveragePerPosition> g_RefCoverageRegion;
   unsigned Start = currentWindow.getStart();
   unsigned End = currentWindow.getEnd();
   unsigned Length = End - Start + 1;
   unsigned NumberOfSamples = g_sampleNames.size();
   std::cout << "declaring g_RefCoverageRegion for " << g_sampleNames.size() << " samples and " << Length << " positions." << std::endl;
   g_RefCoverageRegion.clear();
   RefCoveragePerPosition OnePos;
   for (unsigned SampleIndex = 0; SampleIndex < NumberOfSamples; SampleIndex++) {
      OnePos.RefCoveragePerSample.push_back(0);
   }
   for (unsigned index = 0; index < Length; index++) {
      g_RefCoverageRegion.push_back(OnePos);
   }
   //std::cout << "declare g_RefCoverageRegion done" << std::endl;

   #pragma omp parallel default(shared)
   {
      #pragma omp for
      for (int readIndex = 0; readIndex < (int)currentState.RefSupportingReads.size(); readIndex++) {
         if (currentState.RefSupportingReads[readIndex].Pos < Start || currentState.RefSupportingReads[readIndex].Pos + currentState.RefSupportingReads[readIndex].ReadLength > End) {
            continue;
         }
         unsigned SampleID = g_SampleName2Index.find(currentState.RefSupportingReads[readIndex].Tag)->second;
         unsigned PosStart = currentState.RefSupportingReads[readIndex].Pos - Start;
         #pragma omp critical
         {
            for (unsigned posIndex = 1; posIndex < (unsigned)currentState.RefSupportingReads[readIndex].ReadLength - 1; posIndex++) {
               g_RefCoverageRegion[PosStart + posIndex].RefCoveragePerSample[SampleID]++;
            }
         }
      }//RefSupportingReads
   }
   //std::cout << "update g_RefCoverageRegion done" << std::endl;
   /*
   for (unsigned PosIndex = 0; PosIndex < Length; PosIndex++) {
       std::cout << Start + PosIndex;
       it=g_sampleNames.begin();
       for (unsigned SampleIndex = 0; SampleIndex < NumberOfSamples; SampleIndex++) {
           std::cout << "\t" << *it << ":" << g_RefCoverageRegion[PosIndex].RefCoveragePerSample[SampleIndex];
           it++;
       }
       std::cout << "\n";
   }*/
   return 0;
}

short init_g_ChrNameAndSizeAndIndex(std::string RefIndexFileName)
{
   ChrNameAndSizeAndIndex OneChr;
   std::string TempStr;
   std::ifstream RefIndexInput(RefIndexFileName.c_str());
   if (!RefIndexInput) {
      return 1;
   }
   short ChrCount = 0;
   while (RefIndexInput >> OneChr.ChrName >> OneChr.ChrSize) {
      getline(RefIndexInput, TempStr);
      OneChr.ChrIndex = ChrCount;
      ChrCount++;
      g_ChrNameAndSizeAndIndex.push_back(OneChr);
   }
   return 0;
}

unsigned getChrIndex(std::string & ChrName)
{
   for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
      if (g_ChrNameAndSizeAndIndex[index].ChrName == ChrName) {
         return index;
      }
   }
   return g_ChrNameAndSizeAndIndex.size();
}

unsigned getChrSize(std::string & ChrName)
{
   for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
      if (g_ChrNameAndSizeAndIndex[index].ChrName == ChrName) {
         return g_ChrNameAndSizeAndIndex[index].ChrSize;
      }
   }
   return 0;
}

short CheckChrName(std::string ChrName)
{
   for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
      if (g_ChrNameAndSizeAndIndex[index].ChrName == ChrName) {
         return 1;
      }
   }
   return 0;
}

void CleanUpBedRecord(std::vector <BED> & include, std::vector <BED> & exclude)
{
   if (exclude.size() == 0) {
      return;
   }
   //std::cout << "before include.size() " << include.size() << std::endl;
   for (unsigned include_index = 0; include_index < include.size(); include_index++) {
      //std::cout << "include " << include[include_index].ChrName << "\t" << include[include_index].Start << "\t" << include[include_index].End << std::endl;
      for (unsigned exclude_index = 0; exclude_index < exclude.size(); exclude_index++) {
         //std::cout << "exclude " << include[include_index].ChrName << "\t" << include[include_index].Start << "\t" << include[include_index].End << std::endl;
         if (include[include_index].Start == include[include_index].End) {
            break;
         }
         if (include[include_index].ChrName != exclude[exclude_index].ChrName) {
            continue;
         } else { // same chromosome
            if (include[include_index].Start > exclude[exclude_index].End || exclude[exclude_index].Start > include[include_index].End) {
               continue;
            }
            // "exclude" contains current "include", set start = end
            if (exclude[exclude_index].Start <= include[include_index].Start && include[include_index].End <= exclude[exclude_index].End) { // "exclude" contains current "include", set start = end
               include[include_index].End = include[include_index].Start;
            }
            // "include" contains current "exclude", break into two regions
            else if (include[include_index].Start < exclude[exclude_index].Start && exclude[exclude_index].End < include[include_index].End) {

               BED NewOne;
               NewOne.ChrName = include[include_index].ChrName;
               NewOne.Start  = exclude[exclude_index].End;
               NewOne.End = include[include_index].End;
               include.push_back(NewOne);
               include[include_index].End = exclude[exclude_index].Start;
               //std::cout << include[include_index].ChrName << " " << include[include_index].Start << " " << include[include_index].End << "\t" << NewOne.ChrName << " " << NewOne.Start << " " << NewOne.End << std::endl;
            }
            // intersect I
            else if (exclude[exclude_index].Start <= include[include_index].Start && include[include_index].Start < exclude[exclude_index].End && exclude[exclude_index].End < include[include_index].End) {
               include[include_index].Start = exclude[exclude_index].End;
            }
            // intersect II
            else if (include[include_index].Start < exclude[exclude_index].Start && exclude[exclude_index].Start < include[include_index].End && include[include_index].End < exclude[exclude_index].End) {
               include[include_index].End = exclude[exclude_index].Start;
            }
         }
      }
   }
   //std::cout << "after include.size() " << include.size() << std::endl;
   //for (unsigned index = 0; index < include.size(); index++) {
   //	std::cout << include[index].ChrName << " " << include[index].Start << " " << include[index].End << std::endl;
   //}
   std::vector <BED> result;
   //std::cout << "a list of regions after excluding regions: " << std::endl;
   for (unsigned include_index = 0; include_index < include.size(); include_index++) {
      if (include[include_index].Start != include[include_index].End) {
         result.push_back(include[include_index]);
         //std::cout << include[include_index].ChrName << "\t" << include[include_index].Start << "\t" << include[include_index].End << std::endl;
      }
   }
   //std::cout << "merging the list if there is overlap: " << std::endl;
   for (unsigned first = 0; first < result.size() - 1; first++) {
      for (unsigned second = first + 1; second < result.size(); second++) {
         if (result[first].ChrName != result[second].ChrName) {
            continue;
         } else { // same chromosome

            // non overlap
            if (result[first].Start > result[second].End || result[second].Start > result[first].End) {
               continue;
            }
            // "second" contains current "first", set start = end
            if (result[second].Start <= result[first].Start && result[first].End <=result[second].End) { // "second" contains current "first", set start = end
               result[first].End = result[first].Start;
               break;
            }
            // "first" contains current "second", discard second
            else if (result[first].Start <= result[second].Start && result[second].End <= result[first].End) {
               result[second].Start = result[second].End;
               break;
            }
            // intersect I
            else if (result[second].Start <= result[first].Start && result[first].Start <= result[second].End && result[second].End <= result[first].End) {
               result[first].Start = result[second].Start;
               result[second].Start = result[second].End;
            }
            // intersect II
            else if (result[first].Start <= result[second].Start && result[second].Start <= result[first].End && result[first].End <= result[second].End) {
               result[first].End = result[second].End;
               result[second].Start = result[second].End;
            }
         }
      }
   }
   std::vector <BED> final;
   for (unsigned index = 0; index < result.size(); index++)
      if (result[index].Start != result[index].End) {
         final.push_back(result[index]);
      }
   // sorting
   //std::cout << "final regions" << std::endl;
   bool exchange;
   for (unsigned first = 0; first < final.size() - 1; first++) {
      for (unsigned second = first + 1; second < final.size(); second++) {
         exchange = false;
         unsigned first_ChrIndex = getChrIndex(final[first].ChrName);
         unsigned second_ChrIndex = getChrIndex(final[second].ChrName);
         if (first_ChrIndex < second_ChrIndex) {
            continue;
         } else if (first_ChrIndex > second_ChrIndex) {
            exchange = true;
         } else { // same chr
            if (final[first].Start > final[second].Start) {
               exchange = true;
            }
         }
         if (exchange) {
            BED temp = final[first];
            final[first] = final[second];
            final[second] = temp;
         }
      }
      unsigned CurrentSize = getChrSize(final[first].ChrName);
      if (CurrentSize < final[first].End) {
         final[first].End  = CurrentSize;
      }
      //std::cout << final[first].ChrName << "\t" << final[first].Start << "\t" << final[first].End << std::endl;
   }
   exclude.clear();
   std::cout << "\nfinal list of regions:" << std::endl;
   for (unsigned first = 0; first < final.size(); first++) {
      std::cout << "\t" << final[first].ChrName << "\t" << final[first].Start << "\t" << final[first].End << std::endl;
   }
   std::cout << std::endl;
   include = final;
}

struct InterChrCall {
   char AnchorD;
   char FirstD;
   std::string FirstChrName;
   unsigned FirstPos;
   char SecondD;
   std::string SecondChrName;
   unsigned SecondPos;
   unsigned NumSupport;
   std::string InsertedSequence;
};

void MergeInterChr(ControlState& currentState, UserDefinedSettings *usersettings)
{
   unsigned cutoff = 2;
   std::ifstream INT_input(usersettings->getINTOutputFilename().c_str());
   InterChrCall one;
   std::vector <InterChrCall> All;
   std::string tempstr;
   while (INT_input >> tempstr >> one.AnchorD >> one.FirstChrName >> one.FirstPos >> one.FirstD >> one.SecondChrName >> one.SecondPos >> one.SecondD >> one.InsertedSequence >> tempstr >> one.NumSupport) {
      //if (one.FirstPos != 0 && one.SecondPos != 0)
      All.push_back(one);
      //std::cout << "getting " << one.FirstChrName << "\t" << one.FirstPos << "\t" << one.SecondChrName << "\t" << one.SecondPos << "\t" << one.NumSupport << std::endl;
   }
   std::ofstream INToutputfile((usersettings->getINTOutputFilename() + "_final").c_str());
   if (All.size() == 0) {
      return;
   } else if (All.size() < 2) {
      if (All[0].NumSupport >= cutoff * 2) {
         INToutputfile << All[0].FirstChrName << "\t" << All[0].FirstPos << "\t" << All[0].SecondChrName << "\t"
                       << All[0].SecondPos << "\t" << All[0].InsertedSequence << "\t" << All[0].NumSupport << "\t"
                       << All[0].AnchorD << "\t" << All[0].FirstChrName << "\t" << All[0].FirstPos << "\t" << All[0].FirstD << "\t" << All[0].SecondChrName << "\t"
                       << All[0].SecondPos << "\t" << All[0].SecondD << "\t" << All[0].InsertedSequence << "\t" << All[0].NumSupport << std::endl;
      }
   }
   bool reported;
   for (unsigned index_a = 0; index_a < All.size(); index_a++) {
      reported = false;
      for (unsigned index_b = index_a; index_b < All.size(); index_b++) {
         if (index_a == index_b) {
            continue;
         }
         if (All[index_a].FirstChrName == All[index_b].FirstChrName && All[index_a].SecondChrName == All[index_b].SecondChrName) {
	   if (abs((long long int) All[index_a].FirstPos - (long long int) All[index_b].FirstPos) < 10 && abs((long long int) (All[index_a].SecondPos - All[index_b].SecondPos)) < 10 && All[index_a].NumSupport + All[index_b].NumSupport >= cutoff) {

               INToutputfile << "chr\t" << All[index_a].FirstChrName << "\tpos\t" << unsigned((All[index_a].FirstPos + All[index_b].FirstPos) / 2) << "\tchr\t" << All[index_a].SecondChrName << "\tpos\t"
                             << unsigned((All[index_a].SecondPos + All[index_b].SecondPos) / 2) << "\tseq\t" << All[index_a].InsertedSequence << "\tsupport\t" << All[index_a].NumSupport + All[index_b].NumSupport << "\tINFOR\t"
                             << All[index_a].AnchorD << "\t" << All[index_a].FirstChrName << "\t" << All[index_a].FirstPos << "\t" << All[index_a].FirstD << "\t" << All[index_a].SecondChrName << "\t"
                             << All[index_a].SecondPos << "\t" << All[index_a].SecondD << "\t" << All[index_a].InsertedSequence << "\t" << All[index_a].NumSupport << "\t" << All[index_b].AnchorD << "\t"
                             << All[index_b].FirstChrName << "\t" << All[index_b].FirstPos << "\t" << All[index_b].FirstD << "\t" << All[index_b].SecondChrName << "\t" << All[index_b].SecondPos << "\t"
                             << All[index_b].SecondD << "\t" << All[index_b].InsertedSequence << "\t" << All[index_b].NumSupport << std::endl;
               reported = true;
               break;
            }
         }
      }
      if (reported == false) {
         if (All[index_a].NumSupport >= cutoff * 2) {
            INToutputfile << "chr\t" << All[index_a].FirstChrName << "\tpos\t" << All[index_a].FirstPos << "\tchr\t" << All[index_a].SecondChrName << "\tpos\t"
                          << All[index_a].SecondPos << "\tseq\t" << All[index_a].InsertedSequence << "\tsupport\t" << All[index_a].NumSupport << "\tINFOR\t"
                          << All[index_a].AnchorD << "\t" << All[index_a].FirstChrName << "\t" << All[index_a].FirstPos << "\t" << All[index_a].FirstD << "\t" << All[index_a].SecondChrName << "\t"
                          << All[index_a].SecondPos << "\t" << All[index_a].SecondD << "\t" << All[index_a].InsertedSequence << "\t" << All[index_a].NumSupport << std::endl;
         }
      }
   }
}

int main(int argc, char *argv[])
{
   //TODO: These are counters that are only used in individual steps. They should be moved to separate functions later.
   //Below are variables used for cpu time measurement
   time_t Time_Load_S, Time_Load_E, Time_Sort_E;
   Timer timer;
   timer.switchTo("Initializing pindel");

   Time_Load_S = time(NULL);
   unsigned int AllLoadings = 0;
   unsigned int AllSortReport = 0;

   ControlState currentState;
   userSettings = UserDefinedSettings::Instance();
   std::cout << "Initializing parameters..." << std::endl;
   init(argc, argv, currentState );
   std::cout << "Initializing parameters done." << std::endl;



   if (init_g_ChrNameAndSizeAndIndex(userSettings->getRefFilename() + ".fai") == 1) {
      std::cout << "Please use samtools to index your reference file.\n .fai is missing.\n" << std::endl;
      return 1;
   }

   if (userSettings->loopOverAllChromosomes()) { // WGS
      if (parameters[findParameter("-j",parameters)]->isSet()) { // if a bed file is provided to processing
         if (userSettings->inf_InclusiveBedFileName != "") {
            BED OneBedRecord;
            std::ifstream inf_IncludeBed( userSettings->inf_InclusiveBedFileName.c_str() );
            std::string restofbedline;
            while (inf_IncludeBed >> OneBedRecord.ChrName >> OneBedRecord.Start >> OneBedRecord.End) {
               getline(inf_IncludeBed, restofbedline);
               if (OneBedRecord.Start > OneBedRecord.End) {
                  unsigned tempint = OneBedRecord.Start;
                  OneBedRecord.Start = OneBedRecord.End;
                  OneBedRecord.End = tempint;
               }
               unsigned ChrSize = g_genome.getChr(OneBedRecord.ChrName)->getBiolSize();
               if (OneBedRecord.End > ChrSize) {
                  OneBedRecord.End = ChrSize;
               }
               currentState.IncludeBed.push_back(OneBedRecord);
               //std::cout << currentState.IncludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
            }
         }
      } else {
         for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
            BED OneBedRecord;
            OneBedRecord.ChrName = g_ChrNameAndSizeAndIndex[index].ChrName;
            OneBedRecord.Start = 1;
            OneBedRecord.End = g_ChrNameAndSizeAndIndex[index].ChrSize;
            currentState.IncludeBed.push_back(OneBedRecord);
            //std::cout << currentState.IncludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
         }
      }
   } else { // one region
      if (parameters[findParameter("-j",parameters)]->isSet()) { // if a bed file is provided to processing
         if (userSettings->inf_InclusiveBedFileName != "") {

            std::string ChrName = userSettings->getRegion()->getTargetChromosomeName();
            unsigned Start = userSettings->getRegion()->getStart();
            unsigned End = userSettings->getRegion()->getEnd();
            unsigned ChrSize = g_genome.getChr(ChrName)->getBiolSize();
            if (End > ChrSize) {
               End = ChrSize;
            }
            BED OneBedRecord;
            std::ifstream inf_IncludeBed( userSettings->inf_InclusiveBedFileName.c_str() );
            std::string restofbedline;
            while (inf_IncludeBed >> OneBedRecord.ChrName >> OneBedRecord.Start >> OneBedRecord.End) {
               getline(inf_IncludeBed, restofbedline);
               if (OneBedRecord.Start > OneBedRecord.End) {
                  unsigned tempint = OneBedRecord.Start;
                  OneBedRecord.Start = OneBedRecord.End;
                  OneBedRecord.End = tempint;
               }
               if (OneBedRecord.ChrName != ChrName) {
                  continue;   // not in the same chromosome
               }
               if (OneBedRecord.Start > End) {
                  continue;   // no overlap
               }
               if (Start > OneBedRecord.End) {
                  continue;   // no overlap
               }
               if (OneBedRecord.Start < Start) {
                  OneBedRecord.Start = Start;
               }
               if (OneBedRecord.End > End) {
                  OneBedRecord.End = End;
               }
               currentState.IncludeBed.push_back(OneBedRecord);
               //std::cout << currentState.IncludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
            }
         }
      } else {
         BED OneBedRecord;
         OneBedRecord.ChrName = userSettings->getRegion()->getTargetChromosomeName();
         OneBedRecord.Start = userSettings->getRegion()->getStart();
         OneBedRecord.End = userSettings->getRegion()->getEnd();

         unsigned ChrSize = g_genome.getChr(OneBedRecord.ChrName)->getBiolSize();
         if (OneBedRecord.End > ChrSize) {
            OneBedRecord.End = ChrSize;
         }
         currentState.IncludeBed.push_back(OneBedRecord);
      }
   }


   //std::cout << "currentState.IncludeBed.size() " << currentState.IncludeBed.size() << std::endl;
   //for (unsigned index = 0; index < currentState.IncludeBed.size(); index++) {
   //	std::cout << currentState.IncludeBed.size() << "\t" << currentState.IncludeBed[index].ChrName << "\t" << currentState.IncludeBed[index].Start << "\t" << currentState.IncludeBed[index].End << std::endl;
   //}

   //std::cout << "inf_ExclusiveBedFileName" << userSettings->inf_ExclusiveBedFileName << std::endl;
   if (parameters[findParameter("-J",parameters)]->isSet()) { // if a bed file is provided to exclude regions for processing
      if (userSettings->inf_ExclusiveBedFileName != "") {
         BED OneBedRecord;
         std::ifstream inf_ExcludeBed( userSettings->inf_ExclusiveBedFileName.c_str() );
         std::string restofbedline;
         while (inf_ExcludeBed >> OneBedRecord.ChrName >> OneBedRecord.Start >> OneBedRecord.End) {
            getline(inf_ExcludeBed, restofbedline);
            if (OneBedRecord.Start > OneBedRecord.End) {
               unsigned tempint = OneBedRecord.Start;
               OneBedRecord.Start = OneBedRecord.End;
               OneBedRecord.End = tempint;
            }
            currentState.ExcludeBed.push_back(OneBedRecord);
            //std::cout << currentState.ExcludeBed.size() << "\t" << OneBedRecord.ChrName << "\t" << OneBedRecord.Start << "\t" << OneBedRecord.End << std::endl;
         }
      }
   }

   CleanUpBedRecord(currentState.IncludeBed, currentState.ExcludeBed);

   //return 0;

   if (!userSettings->loopOverAllChromosomes()) {
      if (CheckChrName(userSettings->getRegion()->getTargetChromosomeName()) == 0) {
         std::cout << "Please check chromosome name in the reference file, BAM files and command line. \n Make sure that they are consistent.\n" << std::endl;
         return 1;
      }
   }


   /* Start of shortcut to genotyping *///
   currentState.inf_AssemblyInput.open(userSettings->inf_AssemblyInputFilename.c_str());
   bool GenotypingInputDefined = parameters[findParameter("-g",parameters)]->isSet();

   //std::ifstream FastaFile(userSettings->getRefFilename().c_str());
   if (GenotypingInputDefined) {
      //doGenotyping(currentState, userSettings );
      exit(EXIT_SUCCESS);
   }

   bool AssemblyInputDefined = parameters[findParameter("-z",parameters)]->isSet();
   if (AssemblyInputDefined) {
      //doAssembly(currentState, userSettings );
      exit(EXIT_SUCCESS);
   }

   // If -q parameter given, search for mobile element insertions and quit.
   if (parameters[findParameter("-q", parameters)]->isSet()) {
      exit(searchMEImain(currentState, g_genome, userSettings));
   }

   if (parameters[findParameter("-Y", parameters)]->isSet()) {
      std::string ChrName, TempStr;
      unsigned Ploidy;
      std::ifstream PloidyFile( userSettings->PloidyFileName.c_str() );
      std::map<std::string,unsigned>::iterator it;
      while (PloidyFile >> ChrName >> Ploidy) {
         std::cout << ChrName << "\t" << Ploidy << std::endl;
         getline(PloidyFile, TempStr);
         it = g_ChrName2Ploidy.find(ChrName);
         if (it != g_ChrName2Ploidy.end()) {
            std::cout << "This chromosome " << ChrName << " already seen. Please check " << userSettings->PloidyFileName << std::endl;
            return 1;
         } else {
            g_ChrName2Ploidy.insert(std::pair<std::string, unsigned>(ChrName, Ploidy) );
         }
      }
   }

   std::ofstream RPoutputfile(userSettings->getRPOutputFilename().c_str());

   /* Normal pindel functioning: search SVs*/

   // Get a new chromosome again and again until you have visited the specified chromosome or the file ends
   // CurrentChrName stores the name of the chromosome.
   //std::cout << "loading reference genome" << std::endl;

   std::string CurrentChrName;
   std::string PreviousChrName = "";

   for (unsigned bed_index = 0; bed_index < currentState.IncludeBed.size(); bed_index++) {
      if (parameters[findParameter("-j",parameters)]->isSet()) {
         currentState.CleanUPReads();
      }
      //std::cout << "here" << std::endl;
      timer.switchTo("Loading chromosomes");
      //const Chromosome* currentChromosome = g_genome.getNextChromosome();
      std::string Bed_ChrName = currentState.IncludeBed[bed_index].ChrName;
      unsigned Bed_start = currentState.IncludeBed[bed_index].Start;
      unsigned Bed_end = currentState.IncludeBed[bed_index].End;

      const Chromosome* currentChromosome = g_genome.getChr(Bed_ChrName);

      if (currentChromosome == NULL) {
         std::cout << "There is no " << CurrentChrName << " in the reference file." << std::endl;
         return 1;
      }

      *logStream << "Processing region: " << Bed_ChrName << "\t" << Bed_start << "\t" << Bed_end << std::endl;
      //std::cout << "CurrentChrName: " << currentChromosome->getName() << std::endl;

      g_maxPos = 0; // #################
      *logStream << "Chromosome Size: " << currentChromosome->getBiolSize() << std::endl;
      CurrentChrMask.resize(currentChromosome->getCompSize());
      //std::cout << "currentChromosome->getCompSize()" << currentChromosome->getCompSize() << std::endl;
      for (unsigned int i = 0; i < currentChromosome->getCompSize(); i++) {
         CurrentChrMask[i] = 'N';
      }
      BoxSize = currentChromosome->getCompSize()/ 30000;
      if (BoxSize == 0) {
         BoxSize = 1;
      }
      unsigned NumBoxes = (unsigned) (currentChromosome->getCompSize() * 2 / BoxSize) + 1; // box size
      (*logStream << "NumBoxes: " << NumBoxes << "\tBoxSize: " << BoxSize << std::endl);
      //std::cout << "NumBoxes: " << NumBoxes << "\tBoxSize: " << BoxSize << std::endl;

      g_binIndex = 0; // to start with 0...
      userSettings->getRegion()->SetRegion(Bed_ChrName, Bed_start, Bed_end);
      LoopingSearchWindow currentWindow( userSettings->getRegion(), currentChromosome, WINDOW_SIZE, Bed_start, Bed_end );

      // loop over one bed region
      do {
         //std::cout << "test 1" << std::endl;
         /* 3.2.1 preparation starts */
         *logStream << currentWindow.display();
         //currentState.CURRENT_WINDOW = currentWindow;
         g_RegionStart = currentWindow.getStart();
         g_RegionEnd = currentWindow.getEnd();
         //std::cout << "test 2" << std::endl;

         SearchWindow currentWindow_cs = currentWindow.makePindelCoordinateCopy(); // _cs means computer science coordinates
         //std::cout << "test 3" << std::endl;
         //g_bdData.loadRegion( currentWindow_cs );
         if (Time_Load_S == 0) {
            Time_Load_S = time(NULL);
         }

         timer.switchTo("Reading in reads + matching close ends");//

         //std::cout << "test 4" << std::endl;
         if (userSettings->bamFilesAsInput() && userSettings->SearchDiscordantReadPair) {
            //std::cout << "test 4a" << std::endl;
            //std::cout << "Before" << std::endl;
            get_RP_Reads_Discovery(currentState, currentWindow );
            //std::cout << "test 4b" << std::endl;

            //std::cout << "After" << std::endl;
            g_bdData.UpdateBD(currentState);
            std::cout << "external BD events: " << g_bdData.GetBDSize_external() << " Added BreakDancer-like events: " << (g_bdData.GetBDSize_total() - g_bdData.GetBDSize_external()) / 2 << std::endl;

         }

         //std::cout << "test 5" << std::endl;
         //
         // std::cout << " 1 " << std::endl;
         g_bdData.loadRegion( currentWindow_cs );
         // std::cout << " 2 " << std::endl;

         //std::cout << "test 6" << std::endl;
         //std::cout << "g_bdData.size() " << g_bdData.GetBDSize() << std::endl;
         //std::cout << "Before" << std::endl;
         //g_ReadSeq2Index.clear();
         get_SR_Reads(currentState, currentWindow );


         //std::cout << "g_ReadSeq2Index.size(): " << g_ReadSeq2Index.size() << std::endl;
         //g_ReadSeq2Index.clear();
         std::cout << "There are " << currentState.RefSupportingReads.size() << " reads supporting the reference allele." << std::endl;
         //if (userSettings->bamFilesAsInput())
         UpdateRefReadCoverage(currentState, currentWindow);

         //std::cout << "test 7" << std::endl;

         //for (unsigned index = 0; index < currentState.Reads_SR.size(); index++) {
         //if (currentState.Reads_SR[index].Name == "@M01144:44:000000000-A6N99:1:1101:18676:4723/1")
         //  			std::cout << currentState.Reads_SR[index].Name << std::endl;
         //}
         //Time_Mine_E = time(NULL);
         if (currentState.Reads_SR.size() ) {
            *logStream << "There are " << currentState.Reads_SR.size() << " split-reads for this chromosome region.\n" << std::endl; // what region?
            std::cout << "There are " << g_NumberOfGapAlignedReads << " split-reads mapped by aligner." << std::endl;
            g_NumberOfGapAlignedReads = 0;
            if (userSettings->reportCloseMappedReads() ) {
               *logStream << "report closeMappedReads" << std::endl;
               ReportCloseMappedReads( currentState.Reads_SR );
            }
            Time_Load_E = time(NULL);
            if (!userSettings->reportOnlyCloseMappedReads) {
               timer.switchTo("Searching far ends");
               *logStream << "search far ends" << std::endl;
               SearchFarEnds( currentChromosome->getSeq(), currentState.Reads_SR, *currentChromosome );
               *logStream << "update FarFragName" << std::endl;
               UpdateFarFragName(currentState.Reads_SR);
               *logStream << "update FarFragName done" << std::endl;
               /*
                               			for (unsigned index = 0; index < currentState.Reads_SR.size(); index++) {
               						if (currentState.Reads_SR[index].Name.substr(1, currentState.Reads_SR[index].Name.length() - 1) == "M02294:134:000000000-AEBFT:1:2107:22870:24838" || currentState.Reads_SR[index].Name.substr(1, currentState.Reads_SR[index].Name.length() - 1) == "M02294:134:000000000-AEBFT:1:2107:26950:22204" || currentState.Reads_SR[index].Name.substr(1, currentState.Reads_SR[index].Name.length() - 1) == "M02294:134:000000000-AEBFT:1:2107:22870:24838" || currentState.Reads_SR[index].Name.substr(1, currentState.Reads_SR[index].Name.length() - 1) == "M02294:134:000000000-AEBFT:1:2107:12036:24928" || currentState.Reads_SR[index].Name.substr(1, currentState.Reads_SR[index].Name.length() - 1) == "M02294:134:000000000-AEBFT:1:2110:3389:16550")
                               			        	std::cout << currentState.Reads_SR[index];
                               			}
               */
               /*
                               			for (unsigned index = 0; index < currentState.Reads_SR.size(); index++) {
               						if (currentState.Reads_SR[index].Name.substr(1, currentState.Reads_SR[index].Name.length() - 1) == "M02294:134:000000000-AEBFT:1:2107:22870:24838")
                               			        	std::cout << currentState.Reads_SR[index];
                               			}
               */
               //std::cout << "before currentState.InterChromosome_SR size " << currentState.InterChromosome_SR.size() << std::endl;
               if (userSettings->reportInterchromosomalEvents) {
                  *logStream << "save interchromsome SR" << std::endl;
                  for (unsigned index = 0; index < currentState.Reads_SR.size(); index++) {
                     if (currentState.Reads_SR[index].UP_Far.size()) {
                        //std::cout << currentState.Reads_SR[index].FragName << " " << currentState.Reads_SR[index].FarFragName << std::endl;
                        if (currentState.Reads_SR[index].FragName != currentState.Reads_SR[index].FarFragName) {
                           currentState.InterChromosome_SR.push_back(currentState.Reads_SR[index]);
                        }
                     }
                  }
                  //std::cout << "currentState.InterChromosome_SR.size() = " << currentState.InterChromosome_SR.size() << std::endl;
                  //std::cout << "currentState.Reads_RP_Discovery_InterChr.size() = " << currentState.Reads_RP_Discovery_InterChr.size() << std::endl;
               }

               if (PreviousChrName == CurrentChrName && currentState.FutureReads_SR.size()) {
                  currentState.Reads_SR.insert( currentState.Reads_SR.end(), currentState.FutureReads_SR.begin(), currentState.FutureReads_SR.end() );
               }
               //currentState.FutureReads_SR.clear();
               timer.switchTo("Searching and reporting variations");
               *logStream << "Searching and reporting variations" << std::endl;
               SearchSVs(currentState, NumBoxes, currentWindow );
            }
            Time_Sort_E = time(NULL);

            AllLoadings += (unsigned int) difftime(Time_Load_E, Time_Load_S);
            AllSortReport += (unsigned int) difftime(Time_Sort_E, Time_Load_E);
            //currentState.Reads_SR.clear();
            //currentState.OneEndMappedReads.clear();
            *logStream << "There are " << currentState.FutureReads_SR. size()  << " split-reads saved for the next cycle.\n" << std::endl;
            //currentState.Reads_SR.swap(currentState.FutureReads_SR);
         } else {
            *logStream <<  "no reads " << std::endl;
            //std::cout << " no reads here" << std::endl;
         }
         std::cout << "InterChromosome_SR.size(): " << currentState.InterChromosome_SR.size() << std::endl;
         if (userSettings->reportInterchromosomalEvents && currentState.InterChromosome_SR.size()) {
            SortAndReportInterChromosomalEvents(currentState, g_genome, userSettings);
         }
         {

            /*			std::vector <SPLIT_READ> InputReads_SR, Reads_SR, FutureReads_SR, InterChromosome_SR, OneEndMappedReads;
            			std::vector <RPVector> Reads_RP;
            			std::vector <RP_READ> Reads_RP_Discovery, Reads_RP_Discovery_InterChr;
            			std::vector <REF_READ> RefSupportingReads;

            			currentState.InterChromosome_SR.swap(InterChromosome_SR);
            			currentState.InputReads_SR.swap(InputReads_SR);
            			currentState.Reads_SR.swap(Reads_SR);
            			currentState.FutureReads_SR.swap(FutureReads_SR);
            			currentState.OneEndMappedReads.swap(OneEndMappedReads);
            			currentState.Reads_RP.swap(Reads_RP);
            			currentState.Reads_RP_Discovery.swap(Reads_RP_Discovery);
            			currentState.RefSupportingReads.swap(RefSupportingReads);
            			currentState.Reads_RP_Discovery_InterChr.swap(Reads_RP_Discovery_InterChr);
            */

            currentState.InterChromosome_SR.clear();
            currentState.InputReads_SR.clear();
            currentState.Reads_SR.clear();
            currentState.FutureReads_SR.clear();
            currentState.OneEndMappedReads.clear();
            currentState.Reads_RP.clear();
            currentState.Reads_RP_Discovery.clear();
            currentState.RefSupportingReads.clear();
            currentState.Reads_RP_Discovery_InterChr.clear();

         }

         //std::cout << "after 1" << std::endl;
         /*
         			std::vector <SPLIT_READ> InputReads_SR, Reads_SR, FutureReads_SR, InterChromosome_SR, OneEndMappedReads;
             			std::vector <RPVector> Reads_RP;
             			std::vector <RP_READ> Reads_RP_Discovery, Reads_RP_Discovery_InterChr;
             			std::vector <REF_READ> RefSupportingReads;
         */

         Time_Load_S = 0;
         currentWindow.next();
         g_binIndex++;
         //std::cout << "after 2" << std::endl;
      } while (!currentWindow.finished());
      // name of previous chromosome
      PreviousChrName = CurrentChrName;
      std::cout << "PreviousChrName: " << CurrentChrName << std::endl;
   }

   // If -q parameter given, search for dispersed duplications.
   if (parameters[findParameter("-q", parameters)]->isSet()) {
      int DDresult = searchMEImain(currentState, g_genome, userSettings);
      if (DDresult != 0) {
         exit(DDresult);
      }
   }

   MergeInterChr(currentState, userSettings);

   //std::cout << "before report int " << std::endl;



   //std::cout << "before based on ploidy " << std::endl;
   //if (g_sampleNames.size() == 1 && userSettings->IndelCorrection && currentState.Reads_SR.size()) {
   //	GetConsensusBasedOnPloidy(currentState, g_genome, userSettings);
   //}

   timer.reportAll( *logStream );
   //*logStream << "Loading genome sequences and reads: " << AllLoadings << " seconds." << std::endl;
   //*logStream << "Mining, Sorting and output results: " << AllSortReport << " seconds." << std::endl;
   exit( EXIT_SUCCESS) ;
} //main

std::vector<std::string> ReverseComplement(const std::vector<std::string> &InputPatterns)
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

   for (unsigned int j = 0; j < LenPattern; j++) {
      OutputPattern[j] = Convert2RC4N[(unsigned int) InputPattern[LenPattern - j - 1]];
   }

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
   //return true;
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
   } else {
      return false;
   }
}

void GetRealStart4Deletion(const std::string & TheInput,
                           unsigned int &RealStart, unsigned int &RealEnd)
{

   if (TheInput.size() < RealStart || TheInput.size() < RealEnd) {
      return;
   }
   unsigned int PosIndex = RealStart + g_SpacerBeforeAfter;
   unsigned int Start = PosIndex + 1;
   unsigned int End = RealEnd + g_SpacerBeforeAfter - 1;
   while (TheInput[PosIndex] == TheInput[End] && TheInput[PosIndex] != 'N') {
      --PosIndex;
      --End;
   }
   RealStart = PosIndex - g_SpacerBeforeAfter;
   PosIndex = RealEnd + g_SpacerBeforeAfter;
   while (TheInput[PosIndex] == TheInput[Start] && TheInput[PosIndex] != 'N') {
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

   if (chromosomeSeq.size() < RealStart || chromosomeSeq.size() < RealEnd) {
      return;
   }
   unsigned int lastPosAfterInsertion_comp = RealEnd + g_SpacerBeforeAfter;
   //std::cout << "AInsertedStr: " << InsertedStr << ", start= " << RealStart  << chromosomeSeq[ g_SpacerBeforeAfter + RealStart] << ", end= " << RealEnd << chromosomeSeq[ g_SpacerBeforeAfter + RealEnd] << std::endl;
   //for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
   //std::cout << "\n";
   while ( chromosomeSeq[ lastPosAfterInsertion_comp ] == InsertedStr[ 0 ] && chromosomeSeq[ lastPosAfterInsertion_comp ] != 'N') {
      rotateForward( InsertedStr );
      lastPosAfterInsertion_comp++;
   }
   RealEnd = lastPosAfterInsertion_comp - g_SpacerBeforeAfter;
   //std::cout << "BInsertedStr: " << InsertedStr << ", start= " << RealStart << ", end= " << RealEnd << std::endl;
   //for (int x=-5; x<5; x++ ) { std::cout << chromosomeSeq[ g_SpacerBeforeAfter + RealStart + x] ; }
   //std::cout << "\n";
   unsigned int lastPosBeforeInsertion_comp = lastPosAfterInsertion_comp-1;
   while ( chromosomeSeq[ lastPosBeforeInsertion_comp ] == InsertedStr[ InsertedStr.size()-1 ] && chromosomeSeq[ lastPosBeforeInsertion_comp ] != 'N') {
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

void GetCloseEndInner(const std::string & CurrentChrSeq, SPLIT_READ & Temp_One_Read, int RangeIndex)
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
      PD[CheckIndex].reserve(4 * Temp_One_Read.InsertSize);
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
      Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - RangeIndex * Temp_One_Read.InsertSize; /////////////
      End = Start + (2 * RangeIndex + 1) * Temp_One_Read.InsertSize;
      //std::cout << "1+" << Start << " " << End << std::endl;
      //Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - 2 * RangeIndex * Temp_One_Read.getReadLength(); /////////////
      //End = Start + 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();

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

      CheckLeft_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
      if (UP.empty()) {}
      else {
         Temp_One_Read.Used = false;
         Temp_One_Read.UP_Close.swap(UP);
         UP.clear();
      }
   } else if (Temp_One_Read.MatchedD == Minus) {

      CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
      End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex) * Temp_One_Read.InsertSize; /////////
      Start = End - (2 * RangeIndex + 1) * Temp_One_Read.InsertSize;
      //std::cout << "1-" << Start << " " << End << std::endl;
      //End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + 2 * (RangeIndex) * Temp_One_Read.getReadLength(); /////////
      //Start = End - 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();

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
//        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
      CheckRight_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
//        LOG_DEBUG(*logStream << UP.size() << std::endl);
      if (UP.empty()) {}
      else {
         Temp_One_Read.Used = false;
         Temp_One_Read.UP_Close.swap(UP);
         UP.clear();
      }
   }
   /*
   if (Temp_One_Read.UP_Close.size() == 0) {
       if (Temp_One_Read.MatchedD == Plus) {
           CurrentReadSeq = ReverseComplement(Temp_One_Read.getUnmatchedSeq());
           Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////////
           End = Start + (RangeIndex + 1) * Temp_One_Read.InsertSize;
           //std::cout << "2+" << Start << " " << End << std::endl;
           //Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - 2 * (RangeIndex - 1) * Temp_One_Read.getReadLength(); /////////////
           //End = Start + 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();
           if (Start + 10000 < End) {
               std::cout << "warning: in GetCloseEndInner Start + 10000 < End, slow" << std::endl;
           }
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
           End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////
           Start = End - (RangeIndex + 1) * Temp_One_Read.InsertSize;
           //std::cout << "2-" << Start << " " << End << std::endl;
           //End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + 2 * (RangeIndex - 1) * Temp_One_Read.getReadLength(); /////////
           //Start = End - 2 * (RangeIndex + 1) * Temp_One_Read.getReadLength();
           if (Start + 10000 < End) {
               std::cout << "warning: in GetCloseEndInner Start + 10000 < End, slow" << std::endl;
           }
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
           //        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
           CheckRight_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
           //        LOG_DEBUG(*logStream << UP.size() << std::endl);
           if (UP.empty()) {}
           else {
               Temp_One_Read.Used = false;
               Temp_One_Read.UP_Close.swap(UP);
               UP.clear();
           }
       }
   }
    */

   return;
}

void GetCloseEndInnerPerfectMatch(const std::string & CurrentChrSeq, SPLIT_READ & Temp_One_Read, unsigned RangeIndex)
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
      Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - RangeIndex * Temp_One_Read.InsertSize; /////////////
      End = Start + (RangeIndex + 1) * Temp_One_Read.InsertSize;
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
         CheckLeft_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
      }
      if (UP.empty()) {}
      else {
         Temp_One_Read.Used = false;
         Temp_One_Read.UP_Close.swap(UP);
         UP.clear();
      }
   } else if (Temp_One_Read.MatchedD == Minus) {

      CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
      End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex) * Temp_One_Read.InsertSize; /////////
      Start = End - (RangeIndex + 1) * Temp_One_Read.InsertSize;
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
//        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
      if (PD[0].size()) {
         CheckRight_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
      }
//        LOG_DEBUG(*logStream << UP.size() << std::endl);
      if (UP.empty()) {}
      else {
         Temp_One_Read.Used = false;
         Temp_One_Read.UP_Close.swap(UP);
         UP.clear();
      }
   }

   if (Temp_One_Read.UP_Close.size()) {
      if (Temp_One_Read.MatchedD == Plus) {
         CurrentReadSeq = ReverseComplement(Temp_One_Read.getUnmatchedSeq());
         Start = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter - (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////////
         End = Start + (RangeIndex + 1) * Temp_One_Read.InsertSize;
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
            CheckLeft_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP); // LengthStr
         }
         if (UP.empty()) {}
         else {
            Temp_One_Read.Used = false;
            Temp_One_Read.UP_Close.swap(UP);
            UP.clear();
         }
      } else if (Temp_One_Read.MatchedD == Minus) {

         CurrentReadSeq = Temp_One_Read.getUnmatchedSeq();
         End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter + (RangeIndex - 1) * Temp_One_Read.InsertSize; /////////
         Start = End - (RangeIndex + 1) * Temp_One_Read.InsertSize;
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
         //        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
         if (PD[0].size()) {
            CheckRight_Close_Perfect(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD, BP_Start, BP_End, 1, UP);
         }
         //        LOG_DEBUG(*logStream << UP.size() << std::endl);
         if (UP.empty()) {}
         else {
            Temp_One_Read.Used = false;
            Temp_One_Read.UP_Close.swap(UP);
            UP.clear();
         }
      }
   }
   return;
}

void GetCloseEnd(const std::string & CurrentChrSeq, SPLIT_READ & Temp_One_Read)
{
   const int MaxRange = 2; // 3 insert size away
   //if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
   //    std::cout << "found @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
   //}
   for (int RangeIndex = 0; RangeIndex < MaxRange; RangeIndex++) {
      GetCloseEndInner( CurrentChrSeq, Temp_One_Read, RangeIndex);
      //std::cout << "\nfirst: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
      if (Temp_One_Read.UP_Close.size()==0) { // no good close ends found
         //if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
         //    std::cout << "found step 1 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
         //}
         //std::cout << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
         Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );
         GetCloseEndInner( CurrentChrSeq, Temp_One_Read, RangeIndex );
         //std::cout << "second: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
      }
      //if (Temp_One_Read.UP_Close.size()==0) {
      /*if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
          std::cout << "found step 2 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
      }
       */
      //GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read, RangeIndex );
      //std::cout << "third: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << std::endl;
      //}

      //if (Temp_One_Read.UP_Close.size()==0) { // no good close ends found
      /*if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
          std::cout << "found step 3 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
      }
       */
      //Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );

      //GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read, RangeIndex );
      //std::cout << "fourth: " << Temp_One_Read.Name << " " << Temp_One_Read.UP_Close.size() << "\n" <<  std::endl;
      //}
      /*
      if (Temp_One_Read.Name == "@M01144:44:000000000-A6N99:1:1104:14364:9012/2") {
          std::cout << "found step 4 @M01144:44:000000000-A6N99:1:1104:14364:9012/2" << std::endl;
      }
       */
      if (Temp_One_Read.hasCloseEnd()) {
         break;
      }
   }



   /*
   	GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read );

   	SortedUniquePoints First_UP = Temp_One_Read.UP_Close; // backup

   	Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );

   	GetCloseEndInnerPerfectMatch( CurrentChrSeq, Temp_One_Read );



   	if (First_UP.size() + Temp_One_Read.UP_Close.size()) {
   		std::cout << First_UP.size() << " " << Temp_One_Read.UP_Close.size() << std::endl;
   		if (First_UP.size() > Temp_One_Read.UP_Close.size())
   			Temp_One_Read.UP_Close = First_UP;
   	}
   	else {
   		GetCloseEndInner( CurrentChrSeq, Temp_One_Read );
   		First_UP = Temp_One_Read.UP_Close;
   		Temp_One_Read.setUnmatchedSeq( ReverseComplement( Temp_One_Read.getUnmatchedSeq() ) );
   		GetCloseEndInner( CurrentChrSeq, Temp_One_Read );
   		if (First_UP.size() > Temp_One_Read.UP_Close.size())
   			Temp_One_Read.UP_Close = First_UP;
   	}
   */
}


unsigned CountElements(const FarEndSearchPerRegion * OneRegionSearchResult_input, int levels)
{
   unsigned Sum = 0;

   for (int LevelIndex = 0; LevelIndex < levels; LevelIndex++) {
      Sum += OneRegionSearchResult_input->PD_Plus[LevelIndex].size() + OneRegionSearchResult_input->PD_Minus[LevelIndex].size();
   }
   return Sum;
}

void ExtendMatchPerfect(SPLIT_READ & read,
                        const std::string & readSeq,
                        const std::vector <FarEndSearchPerRegion*> & WholeGenomeSearchResult_input,
                        const short minimumLengthToReportMatch,
                        const short BP_End, const short CurrentLength,
                        SortedUniquePoints &UP )
{
//UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
   const char CurrentChar = readSeq[CurrentLength];
   const char CurrentCharRC = Convert2RC4N[(short) CurrentChar];
   bool AllEmpty = true;
   std::vector <FarEndSearchPerRegion*> WholeGenomeSearchResult_output;
   for (unsigned IndexOfRegion = 0; IndexOfRegion < WholeGenomeSearchResult_input.size(); IndexOfRegion++) {
      const FarEndSearchPerRegion* CurrentRegion_input = WholeGenomeSearchResult_input[IndexOfRegion];
      unsigned int Max_size = 0;
      for (int CheckedIndex = 0; CheckedIndex <= userSettings -> ADDITIONAL_MISMATCH; CheckedIndex++) {
         if (Max_size < CurrentRegion_input->PD_Plus[CheckedIndex].size()) {
            Max_size = CurrentRegion_input->PD_Plus[CheckedIndex].size();
         }
         if (Max_size < CurrentRegion_input->PD_Minus[CheckedIndex].size()) {
            Max_size = CurrentRegion_input->PD_Minus[CheckedIndex].size();
         }
      }
      const std::string & chromosomeSeq = CurrentRegion_input->CurrentChromosome->getSeq();
      FarEndSearchPerRegion* CurrentRegion_output=new FarEndSearchPerRegion (CurrentRegion_input->CurrentChromosome, read.getTOTAL_SNP_ERROR_CHECKED(), Max_size);
      for (int i = 0; i <= read.getTOTAL_SNP_ERROR_CHECKED_Minus(); i++) {
         CategorizePositions( CurrentChar, chromosomeSeq, CurrentRegion_input->PD_Plus, CurrentRegion_output->PD_Plus, i, 1, userSettings -> ADDITIONAL_MISMATCH);
         CategorizePositions( CurrentCharRC, chromosomeSeq, CurrentRegion_input->PD_Minus, CurrentRegion_output->PD_Minus, i, -1, userSettings -> ADDITIONAL_MISMATCH);
      }
      if (CountElements(CurrentRegion_output, 1)) {
         AllEmpty = false;
         WholeGenomeSearchResult_output.push_back(CurrentRegion_output);
      } else {
         delete CurrentRegion_output;
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
      CheckBothPerfect(read, readSeq, WholeGenomeSearchResult_output, minimumLengthToReportMatch, BP_End, CurrentLengthOutput, UP);
   } else {
   } // else-if Sum
   for (unsigned int i=0; i<WholeGenomeSearchResult_output.size(); i++) {
      delete WholeGenomeSearchResult_output[ i ];
   }
}

void ExtendMatch(SPLIT_READ & read,
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
         if (Max_size < CurrentRegion_input->PD_Plus[CheckedIndex].size()) {
            Max_size = CurrentRegion_input->PD_Plus[CheckedIndex].size();
         }
         if (Max_size < CurrentRegion_input->PD_Minus[CheckedIndex].size()) {
            Max_size = CurrentRegion_input->PD_Minus[CheckedIndex].size();
         }
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
      } else {
         delete CurrentRegion_output;
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
   } else {
   } // else-if Sum
   for (unsigned int i=0; i<WholeGenomeSearchResult_output.size(); i++) {
      delete WholeGenomeSearchResult_output[ i ];
   }
}

unsigned int minimumNumberOfMismatches( const std::vector <FarEndSearchPerRegion*> & WholeGenomeSearchResult_input, const unsigned int maxNumberMismatches )
{
   unsigned int Sum=0;
   unsigned int numberOfMismatches=0;
   for (; numberOfMismatches<=maxNumberMismatches; numberOfMismatches++ ) {
      for (unsigned RegionIndex = 0; RegionIndex < WholeGenomeSearchResult_input.size(); RegionIndex++) {
         Sum += WholeGenomeSearchResult_input[RegionIndex]->PD_Plus[numberOfMismatches].size() + WholeGenomeSearchResult_input[RegionIndex]->PD_Minus[numberOfMismatches].size();
      }
      if ( Sum != 0 ) {
         break;
      }
   }
   return numberOfMismatches;
}

void CheckBothPerfect(SPLIT_READ & read,
                      const std::string & readSeq,
                      const std::vector <FarEndSearchPerRegion*> & WholeGenomeSearchResult_input,
                      const short minimumLengthToReportMatch,
                      const short BP_End,
                      const short CurrentLength,
                      SortedUniquePoints &UP)
{
   //UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
   unsigned NumberOfMatchPositionsWithLessMismatches = 0;
   int Sum = 0;

   if (CurrentLength >= minimumLengthToReportMatch && CurrentLength <= BP_End) {
      if (minimumNumberOfMismatches( WholeGenomeSearchResult_input, userSettings -> ADDITIONAL_MISMATCH + 1 ) > g_maxMismatch[CurrentLength] ) {
         return;
      }
      for (short numberOfMismatches = 0; numberOfMismatches < 1; numberOfMismatches++) {
         if (NumberOfMatchPositionsWithLessMismatches) {
            break;
         }

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
                  } else {
                     UniquePoint MinMatch(  hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Minus[numberOfMismatches][0], BACKWARD, ANTISENSE, numberOfMismatches );
                     MatchPosition = MinMatch;
                  }

                  if (CheckMismatches(WholeGenomeSearchResult_input[regionWithMatch]->CurrentChromosome->getSeq(), readSeq, MatchPosition, read.FarEndMismatch)) {
                     UP.push_back (MatchPosition);
                     break;
                  } // if CheckMismatches
               } // if Sum==1
            } // if AdditionalMismatches
         } // if sumsize ==1
      } // for-loop
   } // if length of match is sufficient to be reportable

   if (CurrentLength < BP_End) {
      ExtendMatchPerfect( read, readSeq, WholeGenomeSearchResult_input, minimumLengthToReportMatch, BP_End, CurrentLength, UP );
   }
}

void CheckBoth(SPLIT_READ & read,
               const std::string & readSeq,
               const std::vector <FarEndSearchPerRegion*> & WholeGenomeSearchResult_input,
               const short minimumLengthToReportMatch,
               const short BP_End,
               const short CurrentLength,
               SortedUniquePoints &UP)
{
   //UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
   unsigned NumberOfMatchPositionsWithLessMismatches = 0;
   int Sum = 0;

   if (CurrentLength >= minimumLengthToReportMatch && CurrentLength <= BP_End) {
      if (minimumNumberOfMismatches( WholeGenomeSearchResult_input, read.getMAX_SNP_ERROR() ) > g_maxMismatch[CurrentLength] ) {
         return;
      }
      for (short numberOfMismatches = 0; numberOfMismatches <= read.getMAX_SNP_ERROR(); numberOfMismatches++) {
         if (NumberOfMatchPositionsWithLessMismatches) {
            break;
         }

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
                  } else {
                     UniquePoint MinMatch(  hitRegion->CurrentChromosome, CurrentLength, hitRegion->PD_Minus[numberOfMismatches][0], BACKWARD, ANTISENSE, numberOfMismatches );
                     MatchPosition = MinMatch;
                  }

                  if (CheckMismatches(WholeGenomeSearchResult_input[regionWithMatch]->CurrentChromosome->getSeq(), readSeq, MatchPosition, read.FarEndMismatch)) {
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
         if (LastChr != Input_UP[i].chromosome_p) {
            continue;
         }
         if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
            if (Terminal == Input_UP[i].AbsLoc - Input_UP[i].LengthStr) {
               TempUP.push_back(Input_UP[i]);
            }
         }
      }
   } else if (LastDirection == BACKWARD) {
      Terminal = LastUP.AbsLoc + LastUP.LengthStr;
      for (unsigned i = 0; i < Input_UP.size(); i++) {
         if (LastChr != Input_UP[i].chromosome_p) {
            continue;
         }
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



