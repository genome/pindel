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
#include "pindel.h"
#include "bddata.h"
#include "reader.h"
#include "searcher.h"
#include "reporter.h"
#include "parameter.h"
#include "control_state.h"
#include "search_deletions_nt.h"
#include "search_inversions.h"
#include "search_inversions_nt.h"
#include "search_tandem_duplications.h"
#include "search_tandem_duplications_nt.h"
#include "read_buffer.h"
#include "farend_searcher.h"
#include "search_variant.h"
#include "searchshortinsertions.h"
#include "searchdeletions.h"
#include "logdef.h"
#include "assembly.h"

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

const std::string Pindel_Version_str = "Pindel version 0.2.4q, March 27 2012.";
std::ostream* logStream;
std::ofstream g_logFile;
int findParameter(std::string name);

int g_binIndex = -1; // global variable for the bin index, as I cannot easily pass an extra parameter to the diverse functions
int g_maxPos = -1; // to calculate which is the last position in the chromosome, and hence to calculate the number of bins

std::set<std::string> g_sampleNames;

short Before, After;

BDData g_bdData;


ParCollection par;



unsigned int CountIndels = 0;
const int alphs = 4;
const char alphabet[alphs] = { 'A', 'C', 'G', 'T' };



unsigned long long int TheMax = 0;
const short MAX_MISMATCHES = 4;
float ExtraDistanceRate = 0.1;

double Const_Log_T = 0.0;
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

unsigned int g_NumReadScanned = 0;
unsigned int g_NumReadInWindow = 0;
unsigned int g_InWinPlus = 0;
unsigned int g_InWinMinus = 0;
unsigned int g_CloseMappedPlus = 0;
unsigned int g_CloseMappedMinus = 0;

std::vector<Parameter *> parameters;

// #########################################################
int ADDITIONAL_MISMATCH = 3; // user
unsigned int g_minimalAnchorQuality = 20; // true value set in the defineParameters
int Min_Perfect_Match_Around_BP = 5; // user                   //#
int MIN_IndelSize_NT = 50; //user            //#
int MIN_IndelSize_Inversion = 50; //user       //#
double Seq_Error_Rate = 0.05; // user            //#
unsigned int BalanceCutoff = 0; //                    //#
bool Analyze_INV = true; //# user
bool Analyze_TD = true; //# user
bool Analyze_LI = true; //user
bool Analyze_BP = true; // user
unsigned int NumRead2ReportCutOff = 3; //#
int NumRead2ReportCutOff_BP = 2;
int MaxRangeIndex = 9; // 5 or 6 or 7 or maximum 8      //# // user
double MaximumAllowedMismatchRate = 0.1; //#  // user
int Min_Num_Matched_Bases = 30; //# // user
int Max_Length_NT = 30; // user
bool ReportCloseMappedRead = false; // user
const bool ReportSVReads = false;
const bool ReportLargeInterChrSVReads = false;
const unsigned short Indel_SV_cutoff = 50;
double FLOAT_WINDOW_SIZE = 10.0;
int WINDOW_SIZE = 10000000;
const int AROUND_REGION_BUFFER = 10000; // how much earlier reads should be selected if only a region of the chromosome needs be specified.
// #########################################################

unsigned int Distance = 300;
short MinFar_D = 8; //atoi(argv[3]);
const short MaxDI = 30;

unsigned int SPLIT_READ::getLastAbsLocCloseEnd() const
{
    return UP_Close[ UP_Close.size()-1 ].AbsLoc;
}

bool SPLIT_READ::goodFarEndFound() const
{
    return !UP_Far.empty();
}

bool SPLIT_READ::hasCloseEnd() const
{
    return !UP_Close.empty();
}

unsigned int SPLIT_READ::MaxEndSize( const std::vector<UniquePoint>& upVector) const
{
    if (upVector.size() == 0 ) {
        return 0;
    }
    else {
        int lastElementIndex = upVector.size()-1;
        return upVector[ lastElementIndex ].LengthStr;
    }
}

unsigned int SPLIT_READ::MaxLenFarEnd() const
{
    return MaxEndSize( UP_Far );
}

unsigned int SPLIT_READ::MaxLenFarEndBackup() const
{
    return MaxEndSize( UP_Far_backup);
}

bool doOutputBreakdancerEvents()
{
    return ( par.breakDancerOutputFilename != "" && parameters[findParameter("-b")]->isSet());
}

void outputBreakDancerEvent( const std::string& chromosomeName, const int leftPosition, const int rightPosition,
                             const int svSize, const std::string& svType, const int svCounter)
{
    std::ofstream breakDancerOutputFile(par.breakDancerOutputFilename.c_str(), std::ios::app);
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

bool readTransgressesBinBoundaries(SPLIT_READ & read,
                                   const unsigned int &upperBinBorder)
{
    return (read.BPRight > upperBinBorder - 2 * read.InsertSize);
}

/** 'readInSpecifiedRegion' if a region is specified, check if the read is in it. */
bool readInSpecifiedRegion(const SPLIT_READ & read, // in: the read
                           const bool regionStartDefined,
                           const bool regionEndDefined,
                           const int startOfRegion, // in: the first base of the specified region
                           const int endOfRegion // in: the last base of the specified region (-1 if no region has been specified)
                          )
{
    bool passesFilter = true;

    // if a start position has been defined, and the breakpoint is before it
    if (regionStartDefined && (read.BPLeft + 1 < (unsigned int) startOfRegion)) {
        passesFilter = false;
    }

    // if an end of the region has been specified
    if (regionEndDefined && (read.BPLeft + 1 > (unsigned int) endOfRegion)) {
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
    os << "UnmatchedSeq: " << splitRead.UnmatchedSeq << std::endl;
    os << "HalfMapped: " << splitRead.HalfMapped << std::endl;
    os << "HalfUnMapped: " << splitRead.HalfUnmapped << std::endl;
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
    os << "UP_Far_backup: " ;
    for (unsigned int i=0; i<splitRead.UP_Far_backup.size(); i++) {
        os << "[" << i << "]=" << splitRead.UP_Far_backup[i] << " ";
    }
    os << std::endl;
    os << "ReadLength: " << splitRead.ReadLength << std::endl;
    os << "ReadLengthMinus: " << splitRead.ReadLengthMinus << std::endl;
    os << "MAX_SNP_ERROR:" << splitRead.MAX_SNP_ERROR << std::endl;
    os << "TOTAL_SNP_ERROR_CHECKED:" << splitRead.TOTAL_SNP_ERROR_CHECKED << std::endl;
    os << "TOTAL_SNP_ERROR_CHECKED_Minus:" << splitRead.TOTAL_SNP_ERROR_CHECKED_Minus << std::endl;
    os << "MinClose:" << splitRead.MinClose << std::endl;
    os << "BP:" << splitRead.BP << std::endl;
    os << "Left:" << splitRead.Left << std::endl;
    os << "Right:" << splitRead.Right << std::endl;
    os << "BPLeft:" << splitRead.BPLeft << std::endl;
    os << "BPRight:" << splitRead.BPRight << std::endl;
    os << "IndelSize:" << splitRead.IndelSize << std::endl;
    //os << "UniqueAnchor:" << splitRead.UniqueAnchor << std::endl;
    os << "UniqueRead:" << splitRead.UniqueRead << std::endl;
    os << "score:" << splitRead.score << std::endl;
    //os << "InsertedStr:" << splitRead.InsertedStr << std::endl;
    os << "NT_str:" << splitRead.NT_str  << std::endl;
    os << "NT_size:" << splitRead.NT_size << std::endl;
    return os;
//	cout << ":" <<  << std::endl;
}

/* 'defineParameters' defines the parameters to be used by Pindel. Takes the variables from the calling function as argument for those variables which
 do not need to be stored in the par structure. */
void defineParameters()
{
    parameters. push_back(
        new StringParameter(&par.referenceFileName, "-f", "--fasta",
                            "the reference genome sequences in fasta format", true, ""));
    parameters. push_back(
        new StringParameter(
            &par.pindelFilename,
            "-p",
            "--pindel-file",
            "the Pindel input file; either this, a pindel configuration file (consisting of multiple pindel filenames) or a bam configuration file is required",
            false, ""));
    parameters. push_back(
        new StringParameter(
            &par.bamConfigFileName,
            "-i",
            "--config-file",
            "the bam config file; either this, a pindel input file, or a pindel config file is required. Per line: path and file name of bam, insert size and sample tag.     For example: /data/tumour.bam  400  tumour",
            false, ""));
    parameters. push_back(
        new StringParameter(
            &par.pindelConfigFilename,
            "-P",
            "--pindel-config-file",
            "the pindel config file, containing the names of all Pindel files that need to be sampled; either this, a bam config file or a pindel input file is required. Per line: path and file name of pindel input. Example: /data/tumour.txt",
            false, ""));
    parameters. push_back(
        new StringParameter(&par.outputFileName, "-o", "--output-prefix",
                            "Output prefix;", true, ""));
    parameters. push_back(
        new StringParameter(
            &par.SearchRegion,
            "-c",
            "--chromosome",
            "Which chr/fragment. Pindel will process reads for one chromosome each time. ChrName must be the same as in reference sequence and in read file. '-c ALL' will make Pindel loop over all chromosomes. The search for indels and SVs can also be limited to a specific region; -c 20:10,000,000 will only look for indels and SVs after position 10,000,000 = [10M, end], -c 20:5,000,000-15,000,000 will report indels in the range between and including the bases at position 5,000,000 and 15,000,000 = [5M, 15M]. (default ALL)",
            false, "ALL"));

    parameters. push_back(
        new BoolParameter(&par.showHelp, "-h", "--help",
                          "show the command line options of Pindel", false, false));

    parameters. push_back(
        new IntParameter(&par.numThreads, "-T", "--number_of_threads",
                         "the number of threads Pindel will use (default 1).",
                         false, 1));

    parameters. push_back(
        new IntParameter(
            &MaxRangeIndex,
            "-x",
            "--max_range_index",
            "the maximum size of structural variations to be detected; the higher this number, the greater "
            "the number of SVs reported, but the computational cost and memory requirements increase, as "
            "does the rate of false "
            "positives. 1=128, 2=512, 3=2,048, 4=8,092, 5=32,368, 6=129,472, 7=517,888, 8=2,071,552, 9=8,286,208. "
            "(maximum 9, default 5)", false, 5));
    parameters. push_back(
        new FloatParameter(
            &FLOAT_WINDOW_SIZE,
            "-w",
            "--window_size",
            "for saving RAM, divides the reference in bins of X million bases and only analyzes the reads that belong in the current "
            "bin, " "(default 10 (=10 million))", false, 10));

    parameters. push_back(
        new FloatParameter(&Seq_Error_Rate, "-e",
                           "--sequencing_error_rate",
                           "the expected fraction of sequencing errors "
                           "(default 0.03)", false, 0.03));

    parameters. push_back(
        new FloatParameter(
            &MaximumAllowedMismatchRate,
            "-u",
            "--maximum_allowed_mismatch_rate",
            "Only reads with no less than this fraction of mismatches than the reference genome will be considered. "
            "(default 0.05)", false, 0.05));

    parameters. push_back(
        new BoolParameter(&Analyze_INV, "-r", "--report_inversions",
                          "report inversions " "(default true)", false, true));
    parameters. push_back(
        new BoolParameter(&Analyze_TD, "-t", "--report_duplications",
                          "report tandem duplications " "(default true)", false, true));
    parameters. push_back(
        new BoolParameter(
            &Analyze_LI,
            "-l",
            "--report_long_insertions",
            "report insertions of which the full sequence cannot be deduced because of their length "
            "(default true)", false, true));
    parameters. push_back(
        new BoolParameter(&Analyze_BP, "-k", "--report_breakpoints",
                          "report breakpoints " "(default true)", false, true));

    parameters. push_back(
        new BoolParameter(
            &ReportCloseMappedRead,
            "-s",
            "--report_close_mapped_reads",
            "report reads of which only one end (the one closest to the mapped read of the paired-end read) "
            "could be mapped. " "(default false)", false, false));

    parameters. push_back(
        new BoolParameter(
            &par.reportOnlyCloseMappedReads,
            "-S",
            "--report_only_close_mapped_reads",
            "do not search for SVs, only report reads of which only one end (the one closest to the mapped read of the paired-end read) "
            "could be mapped (the output file can then be used as an input file for another run of pindel, which may save size if you need to transfer files). " "(default false)", false, false));

    parameters. push_back(
        new StringParameter(
            &par.breakdancerFileName,
            "-b",
            "--breakdancer",
            "Pindel is able to use calls from other SV methods such as BreakDancer to further increase sensitivity and specificity.                    BreakDancer result or calls from any methods must in the format:   ChrA LocA stringA ChrB LocB stringB other",
            false, ""));

    parameters. push_back(
        new IntParameter(
            &ADDITIONAL_MISMATCH,
            "-a",
            "--additional_mismatch",
            "Pindel will only map part of a read to the reference genome if there are no other candidate positions with no more than the specified number of mismatches position. The bigger tha value, the more accurate but less sensitive. (default value 1)",
            false, 1));

    parameters. push_back(
        new IntParameter(
            &Min_Perfect_Match_Around_BP,
            "-m",
            "--min_perfect_match_around_BP",
            "at the point where the read is split into two, there should at least be "
            "this number of perfectly matching bases between read and reference (default value 3)",
            false, 3));
    parameters. push_back(
        new IntParameter(&MIN_IndelSize_NT, "-n", "--min_NT_size",
                         "only report inserted (NT) sequences in deletions greater than this size "
                         "(default 50)", false, 50));
    // TODO: Make sure MIN_IndelSize_NT is > 0, make into unsigned.
    parameters. push_back(
        new IntParameter(&MIN_IndelSize_Inversion, "-v",
                         "--min_inversion_size",
                         "only report inversions greater than this number of bases "
                         "(default 50)", false, 50));
    parameters. push_back(
        new IntParameter(
            &Min_Num_Matched_Bases,
            "-d",
            "--min_num_matched_bases",
            "only consider reads as evidence if they map with more than X bases to the reference. "
            "(default 30)", false, 30));
    parameters. push_back(
        new UIntParameter(
            &BalanceCutoff,
            "-B",
            "--balance_cutoff",
            "the number of bases of a SV above which a more stringent filter is applied which demands "
            "that both sides of the SV are mapped with sufficiently long strings of bases "
            "(default 100)", false, 100));
    parameters. push_back(
        new UIntParameter(
            &g_minimalAnchorQuality,
            "-A",
            "--anchor_quality",
            "the minimal mapping quality of the reads Pindel uses as anchor "
            "(default 20)", false, 20));
    parameters. push_back(
        new StringParameter(
            &par.inf_AssemblyInputFilename,
            "-z",
            "--input_SV_Calls_for_assembly",
            "A filename of a list of SV calls for assembling breakpoints \n"
            "Types: DEL, INS, DUP, INV, CTX and ITX \n"    
            "File format: Type chrA posA Confidence_Range_A chrB posB Confidence_Range_B \n"
            "Example: DEL chr1 10000 50 chr2 20000 100 "                
            "", false, ""));
    parameters. push_back(
        new StringParameter(
            &par.breakDancerOutputFilename,
            "-Q",
            "--output_of_breakdancer_events",
            "If breakdancer input is used, you can specify a filename here to write the confirmed breakdancer "
            "events with their exact breakpoints to " "The list of BreakDancer calls with Pindel support information. Format: chr   Loc_left   Loc_right   size   type   index. " "            For example, \"1	72766323 	72811840 	45516	D	11970\" means the deletion event chr1:72766323-72811840 of size 45516 is reported as an event with index 11970 in Pindel report of deletion. "
            "", false, ""));
	parameters.push_back( 
		new StringParameter(
			&par.logFilename,
			"-L",
			"--name_of_logfile",
			"Specifies a file to write Pindel's log to (default: no logfile, log is written to the screen/stdout)", 
			false, ""));

}

/* 'findParameter' returns the index of the parameter with name 'name'; -1 if not found.*/
int findParameter(std::string name)
{
    for (unsigned int parameterCounter = 0; parameterCounter
            < parameters.size(); parameterCounter++) {
        if (parameters[parameterCounter]->hasName(name)) {
            return parameterCounter;
        }
    }
    LOG_DEBUG(cout << "Result of FindParameter is -1\n");
    return -1;
}

/* 'readParameters' reads the parameters as entered in the command line. */
void readParameters(int argc, char *argv[])
{
    // TODO: Ask Kai whether this can be removed
    //for (int argumentIndex=1; argumentIndex<argc; argumentIndex++ ) { log << argumentIndex  << ". " << argv[argumentIndex] << endl; }

    for (int argumentIndex = 1; argumentIndex < argc; argumentIndex++) {
        std::string currentArgument = argv[argumentIndex];

        //find argument in parameterlist
        int parameterIndex = findParameter(currentArgument);
        if (parameterIndex == -1) {
            LOG_ERROR(*logStream << "unknown argument: " << currentArgument << std::endl);
            return;
        }

        if (parameters[parameterIndex]->isUnary()) {
            parameters[parameterIndex]->setValue(true); // default
            if ((argumentIndex + 1 < argc) && (argv[argumentIndex + 1][0]
                                               != '-')) { // so there are more arguments, and next one isn't regular -x
                if (tolower(argv[argumentIndex + 1][0]) == 'f'
                        || (argv[argumentIndex + 1][0] == '0')) {
                    parameters[parameterIndex]->setValue(false);
                }
                argumentIndex++; // in any case increase the argument index
            }
        }
        else {   // argument needs a parameter
            argumentIndex++; // move on to next argument in the list
            if (argumentIndex >= argc) {
                LOG_ERROR(*logStream << "argument of " << currentArgument << " lacking.\n");
                return;
            }
            if (argv[argumentIndex][0] == '-') {
                LOG_ERROR(*logStream << "argument of " << currentArgument
                          << " seems erroneous.\n");
                return;
            }
            // but if everything is allright,
            LOG_DEBUG(*logStream << "Giving " << currentArgument << " the value " << argv[ argumentIndex ] << std::endl);
            parameters[parameterIndex]->setValue(
                std::string(argv[argumentIndex]));
        }
    }
}

/* isReadsFileParam returns whether the parameter points to a read-containing file, and is therefore required,
 even though not both are required. */
bool isReadsFileParam(Parameter * param)
{
    return (param->hasName("-i") || param->hasName("-p"));
}

/* 'printHelp' prints all parameters available. */
void printHelp()
{
    *logStream << std::endl
              << "Program:   pindel (detection of indels and structural variations)"
              << std::endl;
    *logStream << Pindel_Version_str << std::endl;
    *logStream << "Contact:   Kai Ye <k.ye@lumc.nl>" << std::endl << std::endl;

    *logStream << "Usage:     pindel -f <reference.fa> -p <pindel_input>"
              << std::endl;
    *logStream << "           [and/or -i bam_configuration_file]" << std::endl;
    *logStream << "           -c <chromosome_name> -o <prefix_for_output_file>"
              << std::endl << std::endl;

    *logStream << "Required parameters:" << std::endl;

    for (unsigned int i = 0; i < parameters.size(); i++) {
        if (parameters[i]->isRequired() || isReadsFileParam(parameters[i])) {
            parameters[i]->describe();
        }
    }
    *logStream << "\nOptional parameters:" << std::endl;

    for (unsigned int parameterIndex = 0; parameterIndex < parameters.size(); parameterIndex++) {
        if (!parameters[parameterIndex]->isRequired() && !isReadsFileParam(
                    parameters[parameterIndex])) {
            parameters[parameterIndex]->describe();
        }
    }
}

/* 'checkParameters' checks whether all required parameters have been set. */
bool checkParameters()
{
    if (parameters[findParameter("-h")]->getBValue()) {
        printHelp();
        return false;
    }

    for (unsigned int parameterIndex = 0; parameterIndex < parameters.size(); parameterIndex++) {
        if (parameters[parameterIndex]->isRequired()
                && !parameters[parameterIndex]->isSet()) {
            LOG_ERROR(*logStream << "Required parameter "
                      << parameters[parameterIndex]-> getShortName() << "/"
                      << parameters[parameterIndex]-> getLongName() << " "
                      << parameters[parameterIndex]-> getDescription()
                      << " needs to be set." << std::endl);
            return false;
        } //if
    }
    // here handle the tricky fact that at least one of -i or -p needs be set; both are not required.
	bool hasBam = parameters[findParameter("-i")]->isSet();
	bool hasPin = parameters[findParameter("-p")]->isSet();
	bool hasPinConfig = parameters[findParameter("-P")]->isSet();
    if (!hasBam && !hasPin && !hasPinConfig) {
        LOG_ERROR(*logStream
                  << "Bam and/or pindel input file required, use -p, -P and/or -i to designate input file(s)."
                  << std::endl);
        return false;
    }
    LOG_DEBUG(*logStream << "chkP4\n");
    return true;
}

/** 'eliminate' eliminates a character from the input string. */
void eliminate(const char ch, // in: character to be eliminated from the string
               std::string & str // modif: string that needs to be modified
              )
{
    size_t eliminateCharPos = str.find(ch);
    while (eliminateCharPos != std::string::npos) {
        str.erase(eliminateCharPos, 1);
        eliminateCharPos = str.find(ch);
    }
}

/** 'parseRegion' interprets the region specified by the user in the -c option. */
void parseRegion(const std::string & region, // in: region
                 bool &regionStartDefined,
                 bool &regionEndDefined,
                 int &startOfRegion, // out: starting position of the region, -1 if not specified
                 int &endOfRegion, // out: ending position of the region, -1 if not specified
                 std::string & chromosomeName, // out: name of the pure chromosome without region information
                 bool & correctParse // out: whether parsing has succeeded.
                )
{
    size_t separatorPos = region.find(":");
    regionStartDefined = false;
    regionEndDefined = false;
    startOfRegion = -1;
    endOfRegion = -1;
    correctParse = false;

    // found a separator
    if (separatorPos != std::string::npos) {
        chromosomeName = region.substr(0, separatorPos);
        std::string coordinates = region.substr(separatorPos + 1);
        eliminate(',', coordinates); // removes the ',' in 1,000 or 1,000,000 that users may add for readability but wreak havoc with atoi
        size_t startEndSeparatorPos = coordinates.find("-");

        // there are two coordinates
        if (startEndSeparatorPos != std::string::npos) {
            std::string secondPositionStr = coordinates.substr(
                                                startEndSeparatorPos + 1);
            endOfRegion = atoi(secondPositionStr.c_str());
            regionEndDefined = true;
        }
        startOfRegion = atoi(coordinates.c_str());
        regionStartDefined = true;

        LOG_DEBUG(*logStream << "sor: " << startOfRegion << "eor: " << endOfRegion << std::endl);
        if (startOfRegion < 0 || (regionEndDefined && ( endOfRegion < startOfRegion) )) {
            correctParse = false;
        }
        else {
            correctParse = true;
        }

    }
    // no separator found
    else {
        chromosomeName = region;
        correctParse = true;
    }
}

/** 'isFinishedPindel' returns true if there are no more reads to be processed. */
bool isFinishedPindel(const int lastPositionAnalyzed, // in: last position analyzed so far
                      const bool regionEndDefined,
                      const int endOfScan // in: the last position to be scanned
                     )
{
    // if an endOfScan-value has been set
    if ( regionEndDefined ) {
        return (lastPositionAnalyzed >= endOfScan);
    }
    else {   // using g_maxPos
        return (lastPositionAnalyzed >= g_maxPos);
    }
}

/** 'isFinishedBAM' returns true if there are no more reads to be processed. */
bool isFinishedBAM(const int lastPositionAnalyzed, // in: last position analyzed so far
                   const bool regionEndDefined,
                   const int endOfScan, // in: the last position to be scanned
                   const int chromosomeSize // in: the size of the chromosome
                  )
{
    // if an endOfScan-value has been set
    if (regionEndDefined) {
        return (lastPositionAnalyzed >= endOfScan);
    }
    else {   // using chromosomeSize
        return (lastPositionAnalyzed >= chromosomeSize);
    }
}

bool fileExists(const std::string& filename )
{
    std::ifstream test(filename.c_str());
    bool exists = (bool)test;
    test.close();
    return exists;
}

void readBamConfigFile(std::string& bamConfigFileName, ControlState& currentState )
{
    int sampleCounter=0;
    currentState.config_file.open(par.bamConfigFileName.c_str());
    if (currentState.config_file) {
        while (currentState.config_file.good()) {
            currentState.config_file >> currentState.info.BamFile
                                     >> currentState.info.InsertSize;
//*logStream << "Bamfilename: " << currentState.info.BamFile << " Insertsize " << currentState.info.InsertSize << std::endl;
            if (!currentState.config_file.good()) break;
            currentState.info.Tag = "";
            currentState.config_file >> currentState.info.Tag;
//*logStream << "Tag is '" << currentState.info.Tag << "'\n";
            if (currentState.info.Tag=="") {
                *logStream << "Missing tag in line '" << currentState.info.BamFile << "\t" << currentState.info.InsertSize << "' in configuration file " << bamConfigFileName << "\n";
                exit(EXIT_FAILURE);
            }
            g_sampleNames.insert( currentState.info.Tag );
            if (! fileExists( currentState.info.BamFile )) {
                *logStream << "I cannot find the file '"<< currentState.info.BamFile << "'. referred to in configuration file '" << bamConfigFileName << "'. Please change the BAM configuration file.\n\n";
                exit(EXIT_FAILURE);
            }
            if (! fileExists( currentState.info.BamFile+".bai" )) {
                *logStream << "I cannot find the bam index-file '"<< currentState.info.BamFile << ".bai' that should accompany the file " << currentState.info.BamFile << " mentioned in the configuration file " << bamConfigFileName << ". Please run samtools index on " <<
                          currentState.info.BamFile << ".\n\n";
                exit(EXIT_FAILURE);
            }

            //copy kai and throw crap into useless variable
				std::string restOfLine;
            std::getline(currentState.config_file, restOfLine);
            currentState.bams_to_parse.push_back(currentState.info);
            sampleCounter++;
        } // while
        if (sampleCounter==0) {
            *logStream << "Could not find any samples in the sample file '" << bamConfigFileName
                      << "'. Please run Pindel again with a config-file of the specified type (format 'A.bam	<insert-size>	sample_label)\n\n";
            exit( EXIT_FAILURE );
        }
    }
    else {
        // no config-file defined
        *logStream << "BAM configuration file '" << par.bamConfigFileName << "' does not exist. Please run Pindel again with an existing config-file (format 'A.bam	insert-size	sample_label')\n\n";
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


std::string uppercase( const std::string& input )
{
    std::string output = input;
    for(unsigned int pos=0; pos<input.length(); pos++ ) {
        output[ pos] = toupper( input[pos] );
    }
    return output;
}

int init(int argc, char *argv[], ControlState& currentState )
{


    if (NumRead2ReportCutOff == 1) {
        BalanceCutoff = 300000000;
    }

    // define all the parameters you have
    defineParameters();

    // now read the parameters from the command line
    readParameters(argc, argv);



	if (par.logFilename != "" ) {
		g_logFile.open( par.logFilename.c_str() );
		logStream = &g_logFile;
	}

    *logStream << Pindel_Version_str << std::endl;

    if (argc <= 1) { // the user has not given any parameters
        printHelp();
        return EXIT_FAILURE;
    }

    // check parameters
    if (!checkParameters()) {
        exit ( EXIT_FAILURE);
    }
    if (FLOAT_WINDOW_SIZE > 5000.0) {
        LOG_ERROR(*logStream << "Window size of " << FLOAT_WINDOW_SIZE
                  << " million bases is too large" << std::endl);
        return 1;
    }
    else if (FLOAT_WINDOW_SIZE > 100.0) {
        LOG_ERROR(*logStream << "Window size of " << FLOAT_WINDOW_SIZE
                  << " million bases is rather large; this may produce bad::allocs or segmentation faults. If that happens, either try to reduce the window size or deactivate the searching for breakpoints and long insertions by adding the command-line options \"-l false -k false\"." << std::endl);
    }
    WINDOW_SIZE = (int)(1000000 * FLOAT_WINDOW_SIZE);

    // if all parameters are okay, open the files
    currentState.inf_Seq.open(par.referenceFileName.c_str());

    currentState.PindelReadDefined = parameters[findParameter("-p")]->isSet();
    if (currentState.PindelReadDefined) {
        currentState.inf_Pindel_Reads.open(par.pindelFilename.c_str());
    }
	currentState.pindelConfigDefined = parameters[findParameter("-P")]->isSet();
	if (currentState.pindelConfigDefined) {
		readPindelConfigFile( par.pindelConfigFilename, currentState.pindelfilesToParse );
	}
    currentState.BAMDefined = parameters[findParameter("-i")]->isSet();
    if (currentState.BAMDefined) {
        readBamConfigFile( par.bamConfigFileName, currentState );
    }

    currentState.OutputFolder = par.outputFileName;

    currentState.BreakDancerDefined = parameters[findParameter("-b")]->isSet();
    if (currentState.BreakDancerDefined) {
        currentState.inf_BP_test.open(par.breakdancerFileName.c_str());
        currentState.inf_BP.open(par.breakdancerFileName.c_str());
        g_bdData.loadBDFile(par.breakdancerFileName);
    }

    par.AssemblyInputDefined = parameters[findParameter("-z")]->isSet();
    if (par.AssemblyInputDefined) {
        currentState.inf_AssemblyInput.open(par.inf_AssemblyInputFilename.c_str());
    }

    omp_set_num_threads(par.numThreads);

    if (MaxRangeIndex > 9) {
       LOG_ERROR(*logStream
                  << "Maximal range index (-x) exceeds the allowed value (9) - resetting to 9."
                  << std::endl);
        MaxRangeIndex = 9;
    }

    bool WithFolder = false;
    int StartOfFileName = 0;
    for (int i = currentState.bam_file.size(); i >= 0; i--) {
        if (currentState.bam_file[i] == '/') {
            StartOfFileName = i;
            WithFolder = true;
            break;
        }
    }

    if (WithFolder) {
        currentState.bam_file = currentState.bam_file.substr(
                                    StartOfFileName + 1,
                                    currentState.bam_file.size() - 1 - StartOfFileName);
    }

    currentState.SIOutputFilename = currentState.OutputFolder + "_SI"; // output file name
    std::ofstream SIoutputfile_test(currentState.SIOutputFilename.c_str());
    if (!SIoutputfile_test) {
        LOG_ERROR(*logStream << "Sorry, cannot write to the file: "
                  << currentState.SIOutputFilename << std::endl);
        return 1;
    }
    SIoutputfile_test.close();

    currentState.DeletionOutputFilename = currentState.OutputFolder + "_D";
    std::ofstream
    DeletionOutf_test(currentState.DeletionOutputFilename.c_str());
    if (!DeletionOutf_test) {
        LOG_ERROR(*logStream << "Sorry, cannot write to the file: "
                  << currentState.DeletionOutputFilename << std::endl);
        return 1;
    }
    DeletionOutf_test.close();

    currentState.TDOutputFilename = currentState.OutputFolder + "_TD";
    std::ofstream TDOutf_test(currentState.TDOutputFilename.c_str());
    if (!TDOutf_test) {
        LOG_ERROR(*logStream << "Sorry, cannot write to the file: "
                  << currentState.TDOutputFilename << std::endl);
        return 1;
    }
    TDOutf_test.close();

    currentState.InversionOutputFilename = currentState.OutputFolder + "_INV";
    std::ofstream InversionOutf_test(
        currentState.InversionOutputFilename.c_str());
    if (!InversionOutf_test) {
        LOG_ERROR(*logStream << "Sorry, cannot write to the file: "
                  << currentState.InversionOutputFilename << std::endl);
        return 1;
    }
    InversionOutf_test.close();

    currentState.LargeInsertionOutputFilename = currentState.OutputFolder
            + "_LI";
    std::ofstream LargeInsertionOutf_test(
        currentState.LargeInsertionOutputFilename.c_str());
    if (!LargeInsertionOutf_test) {
        LOG_ERROR(*logStream << "Sorry, cannot write to the file: "
                  << currentState.LargeInsertionOutputFilename << std::endl);
        return 1;
    }
    LargeInsertionOutf_test.close();

    currentState.RestOutputFilename = currentState.OutputFolder + "_BP";
    std::ofstream RestOutf_test(currentState.RestOutputFilename.c_str());
    if (!RestOutf_test) {
        LOG_ERROR(*logStream << "Sorry, cannot write to the file: "
                  << currentState.RestOutputFilename << std::endl);
        return 1;
    }
    RestOutf_test.close();

    std::string CloseEndMappedOutputFilename = currentState.OutputFolder + "_CloseEndMapped";
    std::ofstream CloseEndMappedOutput( CloseEndMappedOutputFilename. c_str() );

    currentState.breakDancerOutputFilename = par.breakDancerOutputFilename;

    if ( currentState.breakDancerOutputFilename.compare("") != 0 ) {
        std::ofstream bdOutf_test(currentState.breakDancerOutputFilename.c_str());
        if (!bdOutf_test) {
            LOG_ERROR(*logStream << "Sorry, cannot write to the file: "
                      << currentState.breakDancerOutputFilename << std::endl);
            return 1;
        }
        bdOutf_test.close();
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
    DSizeArray[2] = DSizeArray[1] * 4; // 512
    DSizeArray[3] = DSizeArray[2] * 4; // 2048
    DSizeArray[4] = DSizeArray[3] * 4; // 8092
    DSizeArray[5] = DSizeArray[4] * 4; // 32368
    DSizeArray[6] = DSizeArray[5] * 4; // 129472
    DSizeArray[7] = DSizeArray[6] * 4; // 517888
    DSizeArray[8] = DSizeArray[7] * 4;

    DSizeArray[9] = DSizeArray[8] * 4;
    DSizeArray[10] = DSizeArray[9] * 4;
    DSizeArray[11] = DSizeArray[10] * 4;
    DSizeArray[12] = DSizeArray[11] * 4;
    DSizeArray[13] = DSizeArray[12] * 4;
    DSizeArray[14] = DSizeArray[13] * 4;



    std::vector < std::string > chromosomes;

    char FirstCharOfFasta;
    currentState.inf_Seq >> FirstCharOfFasta;
    if (FirstCharOfFasta != '>') {
        *logStream << "The reference genome must be in fasta format!"
                  << std::endl;
        return 1;
    }

    currentState.SpecifiedChrVisited = false;

    currentState.startOfRegion = -1;
    currentState.endOfRegion = -1;
    bool correctParse = false;
    std::string chrName;
    parseRegion(par.SearchRegion, currentState.regionStartDefined, currentState.regionEndDefined, currentState.startOfRegion,
                currentState.endOfRegion, chrName, correctParse);
    if (!correctParse) {
        LOG_ERROR(*logStream << "I cannot parse the region '" << par.SearchRegion
                  << "'. Please give region in the format -c ALL, -c <chromosome_name> "
                  "(for example -c 20) or -c <chromosome_name>:<start_position>[-<end_position>], for example -c II:1,000 or "
                  "-c II:1,000-50,000. If an end position is specified, it must be larger than the start position."
                  << std::endl);
        exit ( EXIT_FAILURE);
    }
    currentState.TargetChrName = chrName; // removes the region from the 'pure' chromosome name

    if (uppercase(currentState.TargetChrName).compare("ALL") == 0 && par.AssemblyInputDefined == false) {
        *logStream << "Looping over all chromosomes." << std::endl;
        currentState.loopOverAllChromosomes = true;
    }
    return EXIT_SUCCESS;
}


void SearchFarEnd( const std::string& chromosome, SPLIT_READ& read)
{
    const int BD_SPAN = 200;
    const int START_SEARCH_SPAN = 128;

    // when using bins, some reads may already have been assigned far ends already if they were members of the previous bins; they
    // can be skipped here
    if (read.goodFarEndFound()) {
        return;
    }

    std::vector<unsigned int> bdEvents;

    g_bdData.getCorrespondingEvents( read, bdEvents );
    for (unsigned int bdEventIndex=0; bdEventIndex<bdEvents.size(); bdEventIndex++ ) {
        SearchFarEndAtPos( chromosome, read, bdEvents[ bdEventIndex ], BD_SPAN );
        if (read.goodFarEndFound()) {
            return;
        }
    }


    // if breakdancer does not find the event, or not find an event we trust, we turn to regular pattern matching
    int searchSpan=START_SEARCH_SPAN;
    int centerOfSearch = read.getLastAbsLocCloseEnd();
    for (int rangeIndex=1; rangeIndex<=MaxRangeIndex; rangeIndex++ ) {
        SearchFarEndAtPos( chromosome, read, centerOfSearch, searchSpan );
        searchSpan *= 4;
        if (read.goodFarEndFound()) {
            return;
        }
    }

    // if no perfect far end has been found
    if (read.MaxLenFarEndBackup() > read.MaxLenFarEnd() ) {
        read.UP_Far.swap(read.UP_Far_backup);
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

    /* 1 init starts */
    /* 1 init ends */
    /* 2 load genome sequences and reads starts */
    /* 2 load genome sequences and reads ends */

    ControlState currentState;

    int returnValue;

	logStream=&std::cout;
    returnValue = init(argc, argv, currentState );
    if (returnValue != EXIT_SUCCESS) {
        return returnValue;
    }

    std::string emptystr;
    
    
    /* Start of shortcut to assembly */ // currentState.inf_AssemblyInput.open(par.inf_AssemblyInputFilename.c_str());
    //std::cout << "here " << par.AssemblyInputDefined << std::endl;
    if (par.AssemblyInputDefined) {
        
        std::cout << "Entering assembly mode ..." << std::endl;
        doAssembly(currentState, par);
        /*
        // step 1: define output file
        std:string AssemblyOutputFilename = inf_AssemblyInputFilename + "_ASM";
        std::ofstream AssemblyOutput(AssemblyOutputFilename.c_str()); 
        
        // step 1: get whole genome as vector of chr
        
        // step 2: get a list of breakpoints 
         */
        
        std::cout << "Leaving assembly mode and terminating this run." << std::endl;
        exit(EXIT_SUCCESS);
    }
    
    /* End of shortcut to assembly */


    /* 3 loop over chromosomes. this is the most outer loop in the main function. */

    // Get a new chromosome again and again until you have visited the specified chromosome or the file ends
    // CurrentChrName stores the name of the chromosome.
    while (currentState.SpecifiedChrVisited == false && currentState.inf_Seq
            >> currentState.CurrentChrName && !currentState.inf_Seq.eof()) {
        /* 3.1 preparation starts */
        (*logStream << "Processing chromosome: " << currentState.CurrentChrName
         << std::endl);
        //TODO: check with Kai what's the use of this line.
        // dangerous, there may be no other elements on the fasta header line
        std::getline(currentState.inf_Seq, emptystr);
        if (currentState.loopOverAllChromosomes) {
            GetOneChrSeq(currentState.inf_Seq, currentState.CurrentChrSeq, true);
            //currentState.WhichChr = currentState.CurrentChrName;
        }
        else if (currentState.CurrentChrName == currentState.TargetChrName) {   // just one chr and this is the correct one
            GetOneChrSeq(currentState.inf_Seq, currentState.CurrentChrSeq, true);
            currentState.SpecifiedChrVisited = true;
        }
        else {   // not build up sequence
            GetOneChrSeq(currentState.inf_Seq, currentState.CurrentChrSeq, false);
            (*logStream << "Skipping chromosome: " << currentState.CurrentChrName
             << std::endl);
            continue;
        }

        CONS_Chr_Size = currentState.CurrentChrSeq.size() - 2 * g_SpacerBeforeAfter; // #################
        g_maxPos = 0; // #################
        (*logStream << "Chromosome Size: " << CONS_Chr_Size << std::endl);
        CurrentChrMask.resize(currentState.CurrentChrSeq.size());
        g_bdData.loadChromosome( currentState.CurrentChrName, currentState.CurrentChrSeq.size() );
        for (unsigned int i = 0; i < currentState.CurrentChrSeq.size(); i++) {
            CurrentChrMask[i] = 'N';
        }
        BoxSize = currentState.CurrentChrSeq.size() / 30000;
        if (BoxSize == 0) BoxSize = 1;
        unsigned NumBoxes = (unsigned) (currentState.CurrentChrSeq.size() * 2
                                        / BoxSize) + 1; // box size
        (*logStream << "NumBoxes: " << NumBoxes << "\tBoxSize: " << BoxSize << std::endl);

        /* 3.1 preparation ends */

        /* 3.2 apply sliding windows to input datasets starts. This is the 2nd level while loop */
        g_binIndex = 0; // to start with 0... 

        int startOffSet = 0;
        // if a region has been specified
        if (currentState.regionStartDefined ) {
            startOffSet = currentState.startOfRegion - AROUND_REGION_BUFFER;
            if (startOffSet < 0) {
                startOffSet = 0;
            }
        }
        currentState.lowerBinBorder = startOffSet;
        currentState.upperBinBorder = currentState.lowerBinBorder + WINDOW_SIZE;


        if (currentState.regionEndDefined ) {
            currentState.endRegionPlusBuffer = currentState.endOfRegion + AROUND_REGION_BUFFER;
            if (currentState.upperBinBorder > currentState.endRegionPlusBuffer) {
                currentState.upperBinBorder = currentState.endRegionPlusBuffer;
            }
        }


        int displayedStartOfRegion =
            ((currentState.regionStartDefined) ? (currentState.startOfRegion) : currentState.lowerBinBorder);
        int displayedEndOfRegion = displayedStartOfRegion + WINDOW_SIZE;
        if ( displayedEndOfRegion > currentState.upperBinBorder ) {
            displayedEndOfRegion = currentState.upperBinBorder;
        }
        if ( currentState.regionEndDefined && displayedEndOfRegion > currentState.endOfRegion ) {
            displayedEndOfRegion = currentState.endOfRegion;
        }
        
        // loop over one chromosome
        do {

            /* 3.2.1 preparation starts */

            g_NumReadInWindow = 0; // #################
            g_InWinPlus = 0; // #################
            g_InWinMinus = 0; // #################
            g_CloseMappedPlus = 0; // #################
            g_CloseMappedMinus = 0; // #################

            if (displayedStartOfRegion < displayedEndOfRegion) {
                (*logStream << "\nLooking at chromosome " << currentState.CurrentChrName
                 << " bases " << displayedStartOfRegion << " to "
                 << displayedEndOfRegion << "." << std::endl);
            }
            else {
                (*logStream
                 << "Checking out reads near the borders of the specified regions for extra evidence."
                 << std::endl);
            }

            if (Time_Load_S == 0) {
                Time_Load_S = time(NULL);
            }
            //short ReturnFromReadingReads;
            getReads(currentState, par); 
            /*
            if (currentState.BAMDefined) {
                ReturnFromReadingReads = 0;
                for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
                    *logStream << "Insertsize in bamreads: " << currentState.bams_to_parse[i].InsertSize << std::endl;
                    ReturnFromReadingReads = ReadInBamReads(
                                                 currentState.bams_to_parse[i].BamFile.c_str(),
                                                 currentState.CurrentChrName, 
                                                 &currentState.CurrentChrSeq,
                                                 currentState.Reads,
                                                 currentState.bams_to_parse[i].InsertSize,
                                                 currentState.bams_to_parse[i].Tag,
                                                 currentState.lowerBinBorder,
                                                 currentState.upperBinBorder, readBuffer );
                    if (ReturnFromReadingReads == 0) {
                        LOG_ERROR(*logStream << "Bam read failed: "
                                  << currentState.bams_to_parse[i].BamFile
                                  << std::endl);
                        return 1;
                    }
                    else if (currentState.Reads.size() == 0) {
                        LOG_ERROR(*logStream << "No currentState.Reads for "
                                  << currentState.CurrentChrName << " found in "
                                  << currentState.bams_to_parse[i].BamFile
                                  << std::endl);
                    }
                    (*logStream << "BAM file index\t" << i << "\t"
                     << currentState.Reads.size() << std::endl);
                }

            }
            if (currentState.pindelConfigDefined) {	
				for (unsigned int fileIndex=0; fileIndex<currentState.pindelfilesToParse.size(); fileIndex++ ) {
					std::ifstream currentPindelfile( currentState.pindelfilesToParse[ fileIndex ].c_str() );
					readInPindelReads( currentPindelfile, currentState.pindelfilesToParse[ fileIndex ].c_str(), currentState );
				}
			}
            //readBuffer.flush();

            if (currentState.PindelReadDefined) {
					readInPindelReads( currentState.inf_Pindel_Reads, par.pindelFilename, currentState );
            }
            */
            Time_Mine_E = time(NULL);

            if (currentState.Reads.size() ) {
                (*logStream << "There are " << currentState.Reads. size()
                 << " reads for this chromosome region." << std::endl); // what region?

                int TotalNumReads = currentState.Reads.size();
                if (ReportCloseMappedRead || par.reportOnlyCloseMappedReads ) {
                    std::string CloseEndMappedOutputFilename =
                        currentState.OutputFolder + "_CloseEndMapped";
                    std::ofstream CloseEndMappedOutput(
                        CloseEndMappedOutputFilename. c_str(), std::ios::app);
                    for (int Index = 0; Index < TotalNumReads; Index++) {
                        CloseEndMappedOutput << currentState.Reads[Index].Name
                                             << "\n" << currentState.Reads[Index].UnmatchedSeq
                                             << "\n" << currentState.Reads[Index].MatchedD
                                             << "\t" << currentState.Reads[Index].FragName
                                             << "\t" << currentState.Reads[Index].MatchedRelPos
                                             << "\t" << currentState.Reads[Index].MS << "\t"
                                             << currentState.Reads[Index].InsertSize << "\t"
                                             << currentState.Reads[Index].Tag << "\n";
                    }
                }
                Time_Load_E = time(NULL);
                if (!par.reportOnlyCloseMappedReads) {


                    unsigned int Num_Left;
                    Num_Left = currentState.Reads.size();
                    Const_Log_T = log10((double) Num_Left);


                    /* 3.2.1 preparation ends */
                    #pragma omp parallel default(shared)
                    {
                        #pragma omp for
                        for (int readIndex= 0; readIndex < (int)currentState.Reads.size(); readIndex++ ) {
                            SearchFarEnd( currentState.CurrentChrSeq, currentState.Reads[readIndex] );
                        }
                    }





                    (*logStream << "Far end searching completed for this window." << std::endl);

                    SearchDeletions searchD;
                    searchD.Search(currentState, NumBoxes);

                    returnValue = searchIndels(currentState, NumBoxes);

                    if (Analyze_TD) {
                        returnValue = searchTandemDuplications(currentState, NumBoxes);
                        returnValue = searchTandemDuplicationsNT(currentState, NumBoxes);
                    }

                    if (Analyze_INV) {
                        returnValue = searchInversions(currentState, NumBoxes);
                        returnValue = searchInversionsNT(currentState, NumBoxes);
                    }

                    SearchShortInsertions searchSI;
                    searchSI.Search(currentState, NumBoxes);
                    /* 3.2.8 report starts */


                    if (ReportSVReads) {
                        std::string SVReadOutputFilename = currentState.OutputFolder
                                                           + "_SVReads";
                        std::ofstream SVReadOutput(SVReadOutputFilename.c_str(),
                                                   std::ios::app);
                        for (int Index = 0; Index < TotalNumReads; Index++) {
                            if (currentState.Reads[Index].IndelSize > Indel_SV_cutoff
                                    || currentState.Reads[Index].IndelSize == 0)
                                SVReadOutput << currentState.Reads[Index].Name << "\n"
                                             << currentState.Reads[Index]. UnmatchedSeq
                                             << "\n" << currentState.Reads[Index]. MatchedD
                                             << "\t" << currentState.Reads[Index]. FragName
                                             << "\t"
                                             << currentState.Reads[Index]. MatchedRelPos
                                             << "\t" << currentState.Reads[Index]. MS
                                             << "\t"
                                             << currentState.Reads[Index]. InsertSize
                                             << "\t" << currentState.Reads[Index].Tag
                                             << "\n";
                        }
                    }

                    if (ReportLargeInterChrSVReads) {
                        std::string LargeInterChrSVReadsOutputFilename =
                            currentState.OutputFolder + "_LargeORInterChrReads";
                        std::ofstream LargeInterChrSVReadsOutput(
                            LargeInterChrSVReadsOutputFilename.c_str(),
                            std::ios::app);
                        for (int Index = 0; Index < TotalNumReads; Index++) {
                            if (currentState.Reads[Index].IndelSize == 0)
                                LargeInterChrSVReadsOutput
                                        << currentState.Reads[Index].Name << "\n"
                                        << currentState.Reads[Index].UnmatchedSeq
                                        << "\n" << currentState.Reads[Index].MatchedD
                                        << "\t" << currentState.Reads[Index].FragName
                                        << "\t"
                                        << currentState.Reads[Index].MatchedRelPos
                                        << "\t" << currentState.Reads[Index].MS << "\t"
                                        << currentState.Reads[Index].InsertSize << "\t"
                                        << currentState.Reads[Index].Tag << "\n";
                        }
                    }

                    unsigned Count_Far = 0;
                    unsigned Count_Used = 0;
                    unsigned Count_Unused = 0;

                    for (int Index = TotalNumReads - 1; Index >= 0; Index--) {
                        if (!currentState.Reads[Index].UP_Far.empty()) {
                            Count_Far++;
                        }
                        if (!currentState.Reads[Index].UP_Far.empty()
                                || currentState.Reads[Index].Found) {

                        }
                        else {
                            Count_Unused++;
                        }
                        if (currentState.Reads[Index].Used) {
                            Count_Used++;
                        }
                    }

                    (*logStream << "Total: " << TotalNumReads << ";\tClose_end_found "
                     << TotalNumReads << ";\tFar_end_found " << Count_Far
                     << ";\tUsed\t" << Count_Used << "." << std::endl
                     << std::endl);
                    (*logStream << "For LI and BP: " << Count_Unused << std::endl
                     << std::endl);

                    if (Analyze_LI) {
                        time_t Time_LI_S, Time_LI_E;
                        Time_LI_S = time(NULL);
                        std::ofstream LargeInsertionOutf(
                            currentState.LargeInsertionOutputFilename. c_str(),
                            std::ios::app);
                        SortOutputLI(currentState.CurrentChrSeq, currentState.Reads,
                                     LargeInsertionOutf, currentState.lowerBinBorder, currentState.upperBinBorder);
                        LargeInsertionOutf.close();
                        Time_LI_E = time(NULL);
                        (*logStream << "Mining, Sorting and output LI results: "
                         << (unsigned
                             int) difftime(Time_LI_E, Time_LI_S) << " seconds."
                         << std::endl << std::endl);
                        ;
                    }

                    std::vector<SPLIT_READ> BP_Reads;
                    BP_Reads.clear();
                    if (Analyze_BP) {
                        time_t Time_BP_S, Time_BP_E;
                        Time_BP_S = time(NULL);
                        std::ofstream RestOutf(currentState.RestOutputFilename.c_str(),
                                               std::ios::app);
                        SortOutputRest(currentState.CurrentChrSeq, currentState.Reads,
                                       BP_Reads, RestOutf, currentState.lowerBinBorder, currentState.upperBinBorder);
                        RestOutf.close();
                        Time_BP_E = time(NULL);
                        (*logStream << "Mining, Sorting and output BP results: "
                         << (unsigned
                             int) difftime(Time_BP_E, Time_BP_S) << " seconds."
                         << std::endl << std::endl);
                    }
                }
                Time_Sort_E = time(NULL);

                AllLoadings += (unsigned int) difftime(Time_Load_E, Time_Load_S);
                AllSortReport += (unsigned int) difftime(Time_Sort_E, Time_Load_E);
                currentState.Reads.clear();
                (*logStream << "There are " << currentState.FutureReads. size()
                 << " reads saved for the next cycle.\n" << std::endl);
                currentState.Reads.swap(currentState.FutureReads);

            }
            else {
                (*logStream << "There are no reads for this bin." << std::endl);
            }
            Time_Load_S = 0;
            currentState.lowerBinBorder += WINDOW_SIZE;
            currentState.upperBinBorder += WINDOW_SIZE;
            displayedStartOfRegion += WINDOW_SIZE;
            displayedEndOfRegion += WINDOW_SIZE;
            if ( currentState.regionEndDefined && displayedEndOfRegion > currentState.endOfRegion ) {
                displayedEndOfRegion = currentState.endOfRegion;
            }
            if ( currentState.regionEndDefined && currentState.upperBinBorder > currentState.endRegionPlusBuffer ) {
                currentState.upperBinBorder = currentState.endRegionPlusBuffer;
            }
            g_binIndex++;
            /* 3.2.8 report ends */

        } // do {
        while ((currentState.PindelReadDefined && !isFinishedPindel(
                    currentState.lowerBinBorder, currentState.regionEndDefined, currentState.endRegionPlusBuffer))
                || (currentState.BAMDefined && !isFinishedBAM(
                        currentState.lowerBinBorder,
                        currentState.regionEndDefined,
                        currentState.endRegionPlusBuffer,
                        currentState.CurrentChrSeq.size())));

        // Pindel: stop loop if the lowerBinBorder is past the last read, or past the endRegionPlusBuffer
        /* 3.2 apply sliding windows to input datasets ends */

    } // while ( loopOverAllChromosomes && chromosomeIndex < chromosomes.size() );

    (*logStream << "Loading genome sequences and reads: " << AllLoadings
     << " seconds." << std::endl);
    (*logStream << "Mining, Sorting and output results: " << AllSortReport
     << " seconds." << std::endl);
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
    short ReadLength = Deletions[S].ReadLength - Deletions[S].NT_size;
    short Min_Length = (short) ((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
    short Max_Length = (short) (ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
    bool LeftMin = false;
    bool LeftMax = false;
    bool RightMin = false;
    bool RightMax = false;
    for (unsigned i = S; i <= E; i++) {
        ReadLength = Deletions[i].ReadLength - Deletions[i].NT_size;
        Min_Length = (short) ((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
        Max_Length = (short) (ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;
        if (Deletions[i].BP <= Min_Length) {
            LeftMin = true;
        }
        if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size
                <= Min_Length) {
            RightMin = true;
        }
        if (Deletions[i].BP >= Max_Length) {
            LeftMax = true;
        }
        if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size
                >= Max_Length) {
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

/* 'updateReadAfterCloseEndMapping' (EWL, 31thAug2011) */
void updateReadAfterCloseEndMapping( SPLIT_READ& Temp_One_Read )
{
    if (g_reportLength < Temp_One_Read.ReadLength) {
        g_reportLength = Temp_One_Read.ReadLength;
    }
    Temp_One_Read.Used = false;
    //Temp_One_Read.UniqueAnchor = true;
    Temp_One_Read.UniqueRead = true;
    LOG_DEBUG(cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t");

    CleanUniquePoints (Temp_One_Read.UP_Close);

    LOG_DEBUG(cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl);

    Temp_One_Read.CloseEndLength = Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size () - 1].LengthStr;

    if (Temp_One_Read.MatchedD == Plus) {
        Temp_One_Read.LeftMostPos = Temp_One_Read.UP_Close[0].AbsLoc + 1 - Temp_One_Read.UP_Close[0].LengthStr;
        g_CloseMappedPlus++;
    }
    else {
        Temp_One_Read.LeftMostPos = Temp_One_Read.UP_Close[0].AbsLoc +	Temp_One_Read.UP_Close[0].LengthStr - Temp_One_Read.ReadLength;
        g_CloseMappedMinus++;
    }
}



void GetCloseEndInner(const std::string & CurrentChrSeq, SPLIT_READ & Temp_One_Read)
{
    Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
    Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;

    std::string CurrentReadSeq;
    std::vector<unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
    if (Temp_One_Read.InsertSize > g_maxInsertSize) {
        g_maxInsertSize = Temp_One_Read.InsertSize;
    }
    for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
        PD[CheckIndex].reserve(3 * Temp_One_Read.InsertSize);
    }
    std::vector<UniquePoint> UP;
    int Start, End;
    short BP_Start; // = MinClose;
    short BP_End; // = ReadLength - MinClose;

    Temp_One_Read.UP_Close.clear();
    BP_Start = Temp_One_Read.MinClose;
    BP_End = Temp_One_Read.ReadLengthMinus;
    if (Temp_One_Read.MatchedD == Plus) {
        CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
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
        //if (Temp_One_Read.Name == "@1-99550/2") {
        //    logStream << Temp_One_Read.Name << " PD[0].size() " << PD[0].size() << std::endl;
        //}
        if (PD[0].size())
            CheckLeft_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD,
                            BP_Start, BP_End, FirstBase, UP); // LengthStr
        if (UP.empty()) {}
        else {
            Temp_One_Read.Used = false;
            Temp_One_Read.UP_Close.swap(UP);
            UP.clear();
        }
    }
    else if (Temp_One_Read.MatchedD == Minus) {
        CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
        End = Temp_One_Read.MatchedRelPos + g_SpacerBeforeAfter;
        Start = End - 3 * Temp_One_Read.InsertSize;
        char RightChar;
        RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
        if (RightChar != 'N') {
            for (int pos = Start; pos < End; pos++) {
                if (CurrentChrSeq[pos] == RightChar) {
                    PD[0].push_back(pos);
                }
            }
        }

        LOG_DEBUG(*logStream << "1\t" << PD[0].size() << "\t" << PD[1].size() << std::endl);
        CheckRight_Close(Temp_One_Read, CurrentChrSeq, CurrentReadSeq, PD,
                         BP_Start, BP_End, FirstBase, UP);
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
        //*logStream << "Trying to match " << Temp_One_Read.UnmatchedSeq << std::endl;
        Temp_One_Read.UnmatchedSeq = ReverseComplement( Temp_One_Read.UnmatchedSeq );

        GetCloseEndInner( CurrentChrSeq, Temp_One_Read );
        //*logStream << "New attempt: Trying to match " << Temp_One_Read.UnmatchedSeq << "\t"
        //          << Temp_One_Read.UP_Close.size() << std::endl;
    }
}

void CheckBoth(const SPLIT_READ & OneRead, const std::string & TheInput,
               const std::string & CurrentReadSeq,
               const std::vector<unsigned int> PD_Plus[],
               const std::vector<unsigned int> PD_Minus[], const short &BP_Start,
               const short &BP_End, const short &CurrentLength,
               std::vector<UniquePoint> &UP)
{   //std::cout << "in " << CurrentLength << std::endl;
    int Sum;

    if (CurrentLength >= BP_Start && CurrentLength <= BP_End) {
        // put it to LeftUP if unique
        for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
            if (PD_Plus[i].size() + PD_Minus[i].size() == 1 && CurrentLength
                    >= BP_Start + i) {
                Sum = 0;
                if (ADDITIONAL_MISMATCH)
                    for (short j = 0; j <= i + ADDITIONAL_MISMATCH; j++) {
                        Sum += PD_Plus[j].size() + PD_Minus[j].size();
                    }

                if (Sum == 1 && i <= (short) (Seq_Error_Rate * CurrentLength
                                              + 1)) {
                    UniquePoint TempOne;
                    TempOne.LengthStr = CurrentLength;
                    if (PD_Plus[i].size() == 1) {
                        TempOne.Direction = FORWARD;
                        TempOne.Strand = SENSE;
                        TempOne.AbsLoc = PD_Plus[i][0];
                        TempOne.Mismatches = i;
                        if (CheckMismatches(TheInput, OneRead.UnmatchedSeq, TempOne)) {
                            UP.push_back (TempOne);
                            break;
                        }
                    }
                    else if (PD_Minus[i].size() == 1) {
                        TempOne.Direction = BACKWARD;
                        TempOne.Strand = ANTISENSE;
                        TempOne.AbsLoc = PD_Minus[i][0];
                        TempOne.Mismatches = i;
                        if (CheckMismatches(TheInput, OneRead.UnmatchedSeq, TempOne)) { // ######################
                            UP.push_back (TempOne);
                            break;
                        } // ######################
                    }
                }
            }
        }
    }
    if (CurrentLength < BP_End) {
        std::vector<unsigned int>
        PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
        std::vector<unsigned int>
        PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
        for (int CheckedIndex = 0; CheckedIndex
                < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
            PD_Plus_Output[CheckedIndex].reserve(PD_Plus[CheckedIndex]. size());
            PD_Minus_Output[CheckedIndex].reserve(
                PD_Minus[CheckedIndex]. size());
        }
        const char CurrentChar = CurrentReadSeq[CurrentLength];
        const char CurrentCharRC = Convert2RC4N[(short) CurrentChar];
        {
            unsigned int pos;
            int SizeOfCurrent;
            for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
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
                else {
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

            SizeOfCurrent
                = PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
            if (CurrentChar == 'N') {
                for (int j = 0; j < SizeOfCurrent; j++) {
                    pos = PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
                    if (Match2N[(short) TheInput[pos]] == 'N')
                        PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus]. push_back(
                            pos);
                }
            }
            else {
                for (int j = 0; j < SizeOfCurrent; j++) {
                    pos = PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
                    if (TheInput[pos] == CurrentChar)
                        PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus]. push_back(
                            pos);
                }
            }
            SizeOfCurrent
                = PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
            if (CurrentCharRC == 'N') {
                for (int j = 0; j < SizeOfCurrent; j++) {
                    pos = PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j]
                          - 1;
                    if (Match2N[(short) TheInput[pos]] == 'N')
                        PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus]. push_back(
                            pos);
                }
            }
            else {
                for (int j = 0; j < SizeOfCurrent; j++) {
                    pos = PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j]
                          - 1;
                    if (TheInput[pos] == CurrentCharRC)
                        PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus]. push_back(
                            pos);
                }
            }

            Sum = 0;
            for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
                Sum += PD_Plus_Output[i].size() + PD_Minus_Output[i].size();
            }
            if (Sum) {
                const short CurrentLengthOutput = CurrentLength + 1;
                CheckBoth(OneRead, TheInput, CurrentReadSeq, PD_Plus_Output, PD_Minus_Output, BP_Start, BP_End, CurrentLengthOutput, UP);
            }
            else {
                return;
            }
        }
    }
    else {
        return;
    }
//std::cout << "out " << CurrentLength << std::endl;
}

void CleanUniquePoints(std::vector<UniquePoint> &Input_UP)
{
    std::vector<UniquePoint> TempUP; //vector <UniquePoint> UP_Close; UP_Far
    UniquePoint LastUP = Input_UP[Input_UP.size() - 1];
    char LastDirection = LastUP.Direction;
    char LastStrand = LastUP.Strand;
    unsigned int Terminal;

    if (LastDirection == FORWARD) {
        Terminal = LastUP.AbsLoc - LastUP.LengthStr;
        for (unsigned i = 0; i < Input_UP.size(); i++) {
            if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand
                    == LastStrand) {
                if (Terminal == Input_UP[i].AbsLoc - Input_UP[i].LengthStr) {
                    TempUP.push_back(Input_UP[i]);
                }
            }
        }
    }
    else if (LastDirection == BACKWARD) {
        Terminal = LastUP.AbsLoc + LastUP.LengthStr;
        for (unsigned i = 0; i < Input_UP.size(); i++) {
            if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand
                    == LastStrand) {
                if (Terminal == Input_UP[i].AbsLoc + Input_UP[i].LengthStr) {
                    TempUP.push_back(Input_UP[i]);
                }
            }
        }
    }
    Input_UP.clear();
    Input_UP = TempUP;
}


