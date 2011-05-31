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
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <omp.h>

// Pindel header files
#include "pindel.h"
#include "reader.h"
#include "searcher.h"
#include "reporter.h"
#include "parameter.h"
#include "ControlState.h"
#include "search_indels.h"
#include "search_inversions.h"
#include "search_inversions_NT.h"
#include "search_short_insertions.h"
#include "search_tandem_duplications.h"
#include "search_tandem_duplications_NT.h"
#include "search_deletions.h"
#include "FarEndSearcher.h"

/*v EWL update 0.0.1, April 8th 2011; can use the -c option with specified regions, and produces LI output that can be read by vcfcreator */

int g_binIndex = -1; // global variable for the bin index, as I cannot easily pass an extra parameter to the diverse functions
int g_maxPos = -1; // to calculate which is the last position in the chromosome, and hence to calculate the number of bins

// TODO: Ask Kai whether this can be removed
//end charris add
//#include <omp.h>

const std::string Pindel_Version_str = "Pindel version 0.2.2, April 8 2011.";

// TODO: Ask Kai whether this can be removed
//unsigned int DSizeArray[15];

short Before, After;

typedef struct {
	std::string referenceFileName;
	std::string pindelFileName;
	std::string bamConfigFileName;
	std::string outputFileName;
	std::string breakdancerFileName;
	int numThreads;
	bool showHelp;
} ParCollection;
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
//double Const_I = 0.0; // TODO: Ask Kai whether this can be removed
unsigned int BoxSize = 10000; // 10k is fine
const double Min_Filter_Ratio = 0.5;
unsigned int SPACERSIZE = 1;
unsigned int OriginalNumRead = 0;
const std::string NonACGT = "$";
short MIN_Len_Match = 4;
unsigned int NumberOfSIsInstances = 0;
unsigned int NumberOfDeletionsInstances = 0;
unsigned int NumberOfDIInstances = 0;
unsigned int NumberOfInvInstances = 0;
unsigned int NumberOfTDInstances = 0;
short ReportLength = 80;
std::vector<std::string> VectorTag;
char Match[256];
char Match2N[256];
char Convert2RC[256];
char Convert2RC4N[256];
char Cap2LowArray[256];
bool FirstChr = true;
const double InsertSizeExtra = 2;
unsigned int CONS_Chr_Size;
unsigned int DSizeArray[15];

std::string BreakDancerMask;
std::string CurrentChrMask;

unsigned int NumReadScanned = 0;
unsigned int NumReadInChr = 0;
unsigned int InChrPlus = 0;
unsigned int InChrMinus = 0;
unsigned int GetPlus = 0;
unsigned int GetMinus = 0;

// TODO: Ask Kai whether this can be removed
//short MAX_SNP_ERROR = 2;

//short TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
//short TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;

//short MAX_ALLOWED_MISMATCHES = TOTAL_SNP_ERROR_CHECKED_Minus + 5;

std::vector<Parameter *> parameters;

// #########################################################
int ADDITIONAL_MISMATCH = 1; // user
int Min_Perfect_Match_Around_BP = 3; // user                   //#
int MIN_IndelSize_NT = 50; //user            //#
int MIN_IndelSize_Inversion = 50; //user       //#
double Seq_Error_Rate = 0.05; // user            //#
// TODO: Ask Kai whether this can be removed
//float Seq_Error_Rate_1 = 0.05;                         //# 
//float Seq_Error_Rate_2 = 0.02;                         //#
//float Seq_Error_Rate_3 = 0.00;                         //#
unsigned int BalanceCutoff = 100; //                    //#
// TODO: Ask Kai whether this can be removed
//short RangeMaxSensivity = 9;       // 3                //#
//short RangeMediumSensivity = 9;    // 5                //# 
//short RangeLowSensivity = 7;                           //#
bool Analyze_INV = true; //# user
bool Analyze_TD = true; //# user
bool Analyze_LI = true; //EWL070111 user
bool Analyze_BP = true; //EWL070111 user
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
// TODO: Ask Kai whether this can be removed
//bool RemoveDuplicates = false;
double FLOAT_WINDOW_SIZE = 10.0;
int WINDOW_SIZE = 10000000;
const int AROUND_REGION_BUFFER = 10000; // how much earlier reads should be selected if only a region of the chromosome needs be specified.
// #########################################################

// TODO: Ask Kai whether this can be removed
//const float Double_Seq_Error_Rate_Per_Side = Seq_Error_Rate_Per_Side * 2;
unsigned int Distance = 300;
// TODO: Ask Kai whether this can be removed
//short MinClose = 8;//short(log((double)Distance)/log(4.0) + 0.8) + 3 + MAX_SNP_ERROR;//atoi(argv[1]);
//short MinFar_I = MinClose + 1;//atoi(argv[2]);
//cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
short MinFar_D = 8; //atoi(argv[3]);
const short MaxDI = 30;

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

//static struct option long_options[] = { { "fasta", required_argument, 0, 'f' },
//		{ "config-file", required_argument, 0, 'i' }, { "pindel-file",
//				required_argument, 0, 'p' }, { "output-prefix",
//				required_argument, 0, 'o' }, { "chromosome", required_argument,
//				0, 'c' }, { "breakdancer", required_argument, 0, 'b' }, };

bool readTransgressesBinBoundaries(SPLIT_READ & read,
		const unsigned int &upperBinBorder) {
	return (read.BPRight > upperBinBorder - 2 * read.InsertSize);
}

/** 'readInSpecifiedRegion' if a region is specified, check if the read is in it. */
bool readInSpecifiedRegion(const SPLIT_READ & read, // in: the read
		const int startOfRegion, // in: the first base of the specified region
		const int endOfRegion // in: the last base of the specified region (-1 if no region has been specified)
) {
	bool passesFilter = true;

	// if a start position has been defined, and the breakpoint is before it
	if ((startOfRegion != -1) && (read.BPLeft + 1
			< (unsigned int) startOfRegion)) {
		passesFilter = false;
	}

	// if an end of the region has been specified
	if ((endOfRegion != -1) && (read.BPLeft + 1 > (unsigned int) endOfRegion)) {
		passesFilter = false;

	}
	// no region specified, so all positions are okay
	return passesFilter;
}

void saveReadForNextCycle(SPLIT_READ & read,
		std::vector<SPLIT_READ> &futureReads) {
	futureReads.push_back(read);
	read.Used = true; // as it cannot be used for this round of analyses anymore
}

/* 'defineParameters' defines the parameters to be used by Pindel. Takes the variables from the calling function as argument for those variables which
 do not need to be stored in the par structure. */
void defineParameters(std::string & WhichChr) {
	parameters. push_back(
			new StringParameter(&par.referenceFileName, "-f", "--fasta",
					"the reference genome sequences in fasta format", true, ""));
	parameters. push_back(
			new StringParameter(
					&par.pindelFileName,
					"-p",
					"--pindel-file",
					"the Pindel input file; \neither this or a bam configure file is required",
					false, ""));
	parameters. push_back(
			new StringParameter(
					&par.bamConfigFileName,
					"-i",
					"--config-file",
					"the bam config file; either this or a pindel input file is required. Per line: path and file name of bam, insert size and sample tag.     For example: /data/tumour.bam  400  tumour",
					false, ""));
	parameters. push_back(
			new StringParameter(&par.outputFileName, "-o", "--output-prefix",
					"Output prefix;", true, ""));
	parameters. push_back(
			new StringParameter(
					&WhichChr,
					"-c",
					"--chromosome",
					"Which chr/fragment. Pindel will process reads for one chromosome each time. ChrName must be the same as in reference sequence and in read file. '-c ALL' will make Pindel loop over all chromosomes. The search for indels and SVs can also be limited to a specific region; -c 20:10,000,000 will only look for indels and SVs after position 10,000,000 = [10M, end], -c 20:5,000,000-15,000,000 will report indels in the range between and including the bases at position 5,000,000 and 15,000,000 = [5M, 15M]",
					true, ""));

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
						"(default 0.05)", false, 0.05));

	parameters. push_back(
			new FloatParameter(
					&MaximumAllowedMismatchRate,
					"-u",
					"--maximum_allowed_mismatch_rate",
					"Only reads with no less than this fraction of mismatches than the reference genome will be considered. "
						"(default 0.1)", false, 0.1));

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

}

/* 'findParameter' returns the index of the parameter with name 'name'; -1 if not found.*/
int findParameter(std::string name) {
	for (unsigned int parameterCounter = 0; parameterCounter
			< parameters.size(); parameterCounter++) {
		if (parameters[parameterCounter]->hasName(name)) {
			return parameterCounter;
		}
	}
	// TODO: Ask Kai whether this can be removed
	//cout << "Result of FindParameter is -1\n";
	return -1;
}

/* 'readParameters' reads the parameters as entered in the command line. */
void readParameters(int argc, char *argv[]) {
	// TODO: Ask Kai whether this can be removed
	//for (int argumentIndex=1; argumentIndex<argc; argumentIndex++ ) { cout << argumentIndex  << ". " << argv[argumentIndex] << endl; }

	for (int argumentIndex = 1; argumentIndex < argc; argumentIndex++) {
		std::string currentArgument = argv[argumentIndex];

		//find argument in parameterlist
		int parameterIndex = findParameter(currentArgument);
		if (parameterIndex == -1) {
			std::cout << "unknown argument: " << currentArgument << std::endl;
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
		} else { // argument needs a parameter
			argumentIndex++; // move on to next argument in the list
			if (argumentIndex >= argc) {
				std::cout << "argument of " << currentArgument << " lacking.\n";
				return;
			}
			if (argv[argumentIndex][0] == '-') {
				std::cout << "argument of " << currentArgument
						<< " seems erroneous.\n";
				return;
			}
			// but if everything is allright,
			// TODO: Ask Kai whether this can be removed
			//cout << "Giving " << currentArgument << " the value " << argv[ argumentIndex ] << endl;
			parameters[parameterIndex]->setValue(
					std::string(argv[argumentIndex]));
		}
	}
}

/* isReadsFileParam returns whether the parameter points to a read-containing file, and is therefore required,
 even though not both are required. */
bool isReadsFileParam(Parameter * param) {
	return (param->hasName("-i") || param->hasName("-p"));
}

/* 'printHelp' prints all parameters available. */
void printHelp() {
	// TODO: Ask Kai whether this can be removed
	//cout << "Please specify input, either bam configure file and/or pindel input format" << endl;
	//cout.width(0);
	std::cout << std::endl
			<< "Program:   pindel (detection of indels and structural variations)"
			<< std::endl;
	std::cout << "Version:   0.2.2" << std::endl;
	std::cout << "Contact:   Kai Ye <k.ye@lumc.nl>" << std::endl << std::endl;

	std::cout << "Usage:     pindel -f <reference.fa> -p <pindel_input>"
			<< std::endl;
	std::cout << "           [and/or -i bam_configuration_file]" << std::endl;
	std::cout << "           -c <chromosome_name> -o <prefix_for_output_file>"
			<< std::endl << std::endl;

	std::cout << "Required parameters:" << std::endl;
	// TODO: Ask Kai whether this can be removed
	//parameters[1]->describe();
	for (unsigned int i = 0; i < parameters.size(); i++) {
		if (parameters[i]->isRequired() || isReadsFileParam(parameters[i])) {
			parameters[i]->describe();
		}
	}
	std::cout << "\nOptional parameters:" << std::endl;

	for (unsigned int parameterIndex = 0; parameterIndex < parameters.size(); parameterIndex++) {
		if (!parameters[parameterIndex]->isRequired() && !isReadsFileParam(
				parameters[parameterIndex])) {
			parameters[parameterIndex]->describe();
		}
	}
}

/* 'checkParameters' checks whether all required parameters have been set. */
bool checkParameters() {
	if (parameters[findParameter("-h")]->getBValue()) {
		printHelp();
		return false;
	}

	for (unsigned int parameterIndex = 0; parameterIndex < parameters.size(); parameterIndex++) {
		if (parameters[parameterIndex]->isRequired()
				&& !parameters[parameterIndex]->isSet()) {
			std::cout << "Required parameter "
					<< parameters[parameterIndex]-> getShortName() << "/"
					<< parameters[parameterIndex]-> getLongName() << " "
					<< parameters[parameterIndex]-> getDescription()
					<< " needs to be set." << std::endl;
			return false;
		} //if
	}
	// here handle the tricky fact that at least one of -i or -p needs be set; both are not required.
	bool hasBam = parameters[findParameter("-i")]->isSet();
	bool hasPin = parameters[findParameter("-p")]->isSet();
	if (!hasBam && !hasPin) {
		std::cout
				<< "Bam and/or pindel input file required, use -p and/or -i to designate input file(s)."
				<< std::endl;
		return false;
	}
	// TODO: Ask Kai whether this can be removed
	//cout << "chkP4\n";
	return true;
}

/** 'eliminate' eliminates a character from the input string. */
void eliminate(const char ch, // in: character to be eliminated from the string
		std::string & str // modif: string that needs to be modified
) {
	size_t eliminateCharPos = str.find(ch);
	while (eliminateCharPos != std::string::npos) {
		str.erase(eliminateCharPos, 1);
		eliminateCharPos = str.find(ch);
	}
}

/** 'parseRegion' interprets the region specified by the user in the -c option. */
void parseRegion(const std::string & region, // in: region
		int &startOfRegion, // out: starting position of the region, -1 if not specified
		int &endOfRegion, // out: ending position of the region, -1 if not specified
		std::string & chromosomeName, // out: name of the pure chromosome without region information
		bool & correctParse // out: whether parsing has succeeded.
) {
	size_t separatorPos = region.find(":");
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
		}
		startOfRegion = atoi(coordinates.c_str());
		// TODO: Ask Kai whether this can be removed
		//cout << "sor: " << startOfRegion << "eor: " << endOfRegion << endl;
		if (startOfRegion < 0 || (endOfRegion < startOfRegion && endOfRegion
				!= -1)) {
			correctParse = false;
		} else {
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
bool isFinishedPindel(const int upperBinBorder, // in: last position analyzed so far
		const int endOfScan // in: the last position to be scanned
) {
	// if an endOfScan-value has been set
	if (endOfScan != -1) {
		return (upperBinBorder >= endOfScan);
	} else { // using g_maxPos
		return (upperBinBorder >= g_maxPos);
	}
}

/** 'isFinishedBAM' returns true if there are no more reads to be processed. */
bool isFinishedBAM(const int upperBinBorder, // in: last position analyzed so far
		const int endOfScan, // in: the last position to be scanned
		const int chromosomeSize // in: the size of the chromosome
) {
	// if an endOfScan-value has been set
	if (endOfScan != -1) {
		return (upperBinBorder >= endOfScan);
	} else { // using chromosomeSize
		return (upperBinBorder >= chromosomeSize);
	}
}

int init(int argc, char *argv[], ControlState& currentState) {
	std::cout << Pindel_Version_str << std::endl;

	if (NumRead2ReportCutOff == 1)
		BalanceCutoff = 3000000000u;

	// define all the parameters you have
	defineParameters(currentState.WhichChr);

	// now read the parameters from the command line
	readParameters(argc, argv);

	if (argc <= 1) { // the user has not given any parameters
		printHelp();
		return EXIT_FAILURE;
	}

	// check parameters
	if (!checkParameters())
		exit ( EXIT_FAILURE);
	if (FLOAT_WINDOW_SIZE > 5000.0) {
		std::cout << "Window size " << FLOAT_WINDOW_SIZE
				<< " million bases is too large" << std::endl;
		return 1;
	} else if (FLOAT_WINDOW_SIZE > 500.0) {
		std::cout << "Window size " << FLOAT_WINDOW_SIZE
				<< " million bases is rather large" << std::endl;
	}
	WINDOW_SIZE = 1000000 * FLOAT_WINDOW_SIZE;

	// if all parameters are okay, open the files
	currentState.inf_Seq.open(par.referenceFileName.c_str());

	currentState.PindelReadDefined = parameters[findParameter("-p")]->isSet();
	if (currentState.PindelReadDefined) {
		currentState.inf_Pindel_Reads.open(par.pindelFileName.c_str());
	}

	currentState.BAMDefined = parameters[findParameter("-i")]->isSet();
	if (currentState.BAMDefined) {
		currentState.config_file.open(par.bamConfigFileName.c_str());
		while (currentState.config_file.good()) {
			currentState.config_file >> currentState.info.BamFile
					>> currentState.info.InsertSize >> currentState.info.Tag;
			//copy kai and throw crap into useless variable
			std::getline(currentState.config_file, currentState.line);
			if (currentState.config_file.good()) {
				currentState.bams_to_parse.push_back(currentState.info);
			}
		}
	}

	currentState.OutputFolder = par.outputFileName;

	currentState.BreakDancerDefined = parameters[findParameter("-b")]->isSet();
	if (currentState.BreakDancerDefined) {
		currentState.inf_BP_test.open(par.breakdancerFileName.c_str());
		currentState.inf_BP.open(par.breakdancerFileName.c_str());
	}

	omp_set_num_threads(par.numThreads);

	if (MaxRangeIndex > 9) {
		std::cout
				<< "Maximal range index (-x) exceeds the allowed value (9) - resetting to 9."
				<< std::endl;
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
		std::cout << "Sorry, cannot write to the file: "
				<< currentState.SIOutputFilename << std::endl;
		return 1;
	}
	SIoutputfile_test.close();

	currentState.DeletionOutputFilename = currentState.OutputFolder + "_D";
	std::ofstream
			DeletionOutf_test(currentState.DeletionOutputFilename.c_str());
	if (!DeletionOutf_test) {
		std::cout << "Sorry, cannot write to the file: "
				<< currentState.DeletionOutputFilename << std::endl;
		return 1;
	}
	DeletionOutf_test.close();

	currentState.TDOutputFilename = currentState.OutputFolder + "_TD";
	std::ofstream TDOutf_test(currentState.TDOutputFilename.c_str());
	if (!TDOutf_test) {
		std::cout << "Sorry, cannot write to the file: "
				<< currentState.TDOutputFilename << std::endl;
		return 1;
	}
	TDOutf_test.close();

	currentState.InversionOutputFilename = currentState.OutputFolder + "_INV";
	std::ofstream InversionOutf_test(
			currentState.InversionOutputFilename.c_str());
	if (!InversionOutf_test) {
		std::cout << "Sorry, cannot write to the file: "
				<< currentState.InversionOutputFilename << std::endl;
		return 1;
	}
	InversionOutf_test.close();

	currentState.LargeInsertionOutputFilename = currentState.OutputFolder
			+ "_LI";
	std::ofstream LargeInsertionOutf_test(
			currentState.LargeInsertionOutputFilename.c_str());
	if (!LargeInsertionOutf_test) {
		std::cout << "Sorry, cannot write to the file: "
				<< currentState.LargeInsertionOutputFilename << std::endl;
		return 1;
	}
	LargeInsertionOutf_test.close();

	currentState.RestOutputFilename = currentState.OutputFolder + "_BP";
	std::ofstream RestOutf_test(currentState.RestOutputFilename.c_str());
	if (!RestOutf_test) {
		std::cout << "Sorry, cannot write to the file: "
				<< currentState.RestOutputFilename << std::endl;
		return 1;
	}
	RestOutf_test.close();

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
	for (unsigned int i = 0; i < SpacerBeforeAfter; i++)
		Spacer += "N";

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

	std::string TempLie_BD;

	char FirstSharpChar;

	std::string TempLine_BD;
	BreakDancer Temp_BD_event;

	currentState.All_BD_events_WG.push_back(Temp_BD_event);

	while (currentState.inf_BP_test >> FirstSharpChar) {
		if (FirstSharpChar == '#') {
			std::getline(currentState.inf_BP_test, TempLine_BD);
			std::getline(currentState.inf_BP, TempLine_BD);
		} else {
			std::getline(currentState.inf_BP_test, TempLine_BD);
			currentState.inf_BP >> Temp_BD_event.ChrName_A >> Temp_BD_event.S1
					>> TempLine_BD >> Temp_BD_event.ChrName_B
					>> Temp_BD_event.S3 >> TempLine_BD;
			std::getline(currentState.inf_BP, TempLine_BD);

			Temp_BD_event.S2 = Temp_BD_event.S1 + 200;
			if (Temp_BD_event.S1 > 200)
				Temp_BD_event.S1 = Temp_BD_event.S1 - 200;
			else
				Temp_BD_event.S1 = 1;

			Temp_BD_event.S4 = Temp_BD_event.S3 + 200;
			if (Temp_BD_event.S3 > 200)
				Temp_BD_event.S3 = Temp_BD_event.S3 - 200;
			else
				Temp_BD_event.S3 = 1;

			Temp_BD_event.S1 += SpacerBeforeAfter;
			Temp_BD_event.S2 += SpacerBeforeAfter;
			Temp_BD_event.S3 += SpacerBeforeAfter;
			Temp_BD_event.S4 += SpacerBeforeAfter;

			currentState.All_BD_events_WG.push_back(Temp_BD_event);
		}
	}
	std::cout << "BreakDancer events: " << currentState.All_BD_events_WG.size()
			- 1 << std::endl;

	std::vector < std::string > chromosomes;

	char FirstCharOfFasta;
	currentState.inf_Seq >> FirstCharOfFasta;
	if (FirstCharOfFasta != '>') {
		std::cout << "The reference genome must be in fasta format!"
				<< std::endl;
		return 1;
	}

	currentState.SpecifiedChrVisited = false;

	currentState.startOfRegion = -1;
	currentState.endOfRegion = -1;
	bool correctParse = false;
	std::string chrName;
	parseRegion(currentState.WhichChr, currentState.startOfRegion,
			currentState.endOfRegion, chrName, correctParse);
	if (!correctParse) {
		std::cout << "I cannot parse the region '" << currentState.WhichChr
				<< "'. Please give region in the format -c ALL, -c <chromosome_name> "
					"(for example -c 20) or -c <chromosome_name>:<start_position>[-<end_position>], for example -c II:1,000 or "
					"-c II:1,000-50,000. If an end position is specified, it must be larger than the start position."
				<< std::endl;
		exit ( EXIT_FAILURE);
	}
	currentState.WhichChr = chrName; // removes the region from the 'pure' chromosome name

	int startOffSet = 0;
	// if a region has been specified
	if (currentState.startOfRegion >= 0) {
		startOffSet = currentState.startOfRegion - AROUND_REGION_BUFFER;
		if (startOffSet < 0) {
			startOffSet = 0;
		}
	}
	currentState.lowerBinBorder = startOffSet - WINDOW_SIZE;
	currentState.upperBinBorder = currentState.lowerBinBorder + WINDOW_SIZE;

	currentState.endRegionPlusBuffer = -1; // -1 indicates that the chromosome must be read to the end
	if (currentState.endOfRegion > -1) {
		currentState.endRegionPlusBuffer = currentState.endOfRegion
				+ AROUND_REGION_BUFFER;
	}

	if (currentState.WhichChr.compare("ALL") == 0) {
		std::cout << "Looping over ALL chromosomes." << std::endl;
		currentState.loopOverAllChromosomes = true;
	}
	return EXIT_SUCCESS;
}

int searchBreakPoints(ControlState& currentState) {
	/* 3.2.2.1 prepare break dancer events */
	currentState.All_BD_events.clear();
	for (unsigned int All_BD_events_WG_index = 0; All_BD_events_WG_index
			< currentState.All_BD_events_WG.size(); All_BD_events_WG_index++) {
		if (currentState.All_BD_events_WG[All_BD_events_WG_index].ChrName_A
				== currentState.CurrentChrName
				&& currentState.All_BD_events_WG[All_BD_events_WG_index].ChrName_B
						== currentState.CurrentChrName)
			currentState.All_BD_events.push_back(
					currentState.All_BD_events_WG[All_BD_events_WG_index]);
	}

	if (currentState.All_BD_events.size() > 1) {
		std::cout
				<< "Searching additional breakpoints by adding BreakDancer results"
				<< std::endl;
		int *BD_INDEX = new int[currentState.CurrentChr.size()];
		for (unsigned i = 0; i < currentState.CurrentChr.size(); i++)
			BD_INDEX[i] = 0;
		for (unsigned i = 1; i < currentState.All_BD_events.size(); i++) {
			for (unsigned j = currentState.All_BD_events[i].S1; j
					< currentState.All_BD_events[i].S2; j++)
				BD_INDEX[j] = i;
			for (unsigned j = currentState.All_BD_events[i].S3; j
					< currentState.All_BD_events[i].S4; j++)
				BD_INDEX[j] = i * (-1);
		}
		int BD_Plus = 0;
		int BD_Minus = 0;
		for (unsigned i = 0; i < currentState.CurrentChr.size(); i++) {
			if (BD_INDEX[i] > 0)
				BD_Plus++;
			else if (BD_INDEX[i] < 0)
				BD_Minus++;
		}
		std::cout << BD_Plus << "\t" << BD_Minus << std::endl;

		currentState.CountFarEnd = 0;
		currentState.CountFarEndMinus = 0;
		currentState.CountFarEndPlus = 0;
		int Start_pos, End_pos;
		for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads.size(); ReadIndex++) {
			if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
				continue;
			}
			int BD_event_Index =
					BD_INDEX[currentState.Reads[ReadIndex].UP_Close[0].AbsLoc];
			if (BD_event_Index == 0)
				continue;
			if (BD_event_Index > 0) {
				Start_pos = currentState.All_BD_events[BD_event_Index].S3
						- currentState.Reads[ReadIndex].ReadLength;
				End_pos = currentState.All_BD_events[BD_event_Index].S4
						+ currentState.Reads[ReadIndex].ReadLength;
				GetFarEnd(currentState.CurrentChr,
						currentState.Reads[ReadIndex], Start_pos, End_pos);
			} else {
				Start_pos
						= currentState.All_BD_events[BD_event_Index * (-1)].S1
								- currentState.Reads[ReadIndex].ReadLength;
				End_pos = currentState.All_BD_events[BD_event_Index * (-1)].S2
						+ currentState.Reads[ReadIndex].ReadLength;
				GetFarEnd(currentState.CurrentChr,
						currentState.Reads[ReadIndex], Start_pos, End_pos);
			}

			if (!currentState.Reads[ReadIndex].UP_Far.empty()) { // if a far end has been found


				// if there is a non-template sequence present between close and far end
				if (currentState.Reads[ReadIndex]. UP_Far[currentState.Reads[ReadIndex].UP_Far.size()
						- 1].LengthStr
						+ currentState.Reads[ReadIndex].CloseEndLength
						< currentState.Reads[ReadIndex].ReadLength) {

					// if there are backup reads
					if (currentState.Reads[ReadIndex].UP_Far_backup.size()) {

						// if the backup reads are worse
						if (currentState.Reads[ReadIndex]. UP_Far_backup[currentState.Reads[ReadIndex]. UP_Far_backup.size()
								- 1].LengthStr
								< currentState.Reads[ReadIndex].UP_Far[currentState.Reads[ReadIndex]. UP_Far.size()
										- 1].LengthStr) {

							// put current reads in backup
							currentState.Reads[ReadIndex].UP_Far_backup
									= currentState.Reads[ReadIndex].UP_Far;
							currentState.Reads[ReadIndex].UP_Far.clear();
						} else
							currentState.Reads[ReadIndex].UP_Far.clear(); // otherwise keep your backup
					} else { // there are no backup reads, prepare for next cycle
						currentState.Reads[ReadIndex].UP_Far_backup
								= currentState.Reads[ReadIndex].UP_Far;
						currentState.Reads[ReadIndex].UP_Far.clear();
					}
				} else { // no non-template bases present: good match found
					currentState.CountFarEnd++;
					if (currentState.Reads[ReadIndex].MatchedD == Plus)
						currentState.CountFarEndPlus++;
					else
						currentState.CountFarEndMinus++;
					// note:UPFar remains filled, so read will be skipped in next cycle
				}
			}
		}
		std::cout << "\tNumber of reads with far end mapped: "
				<< currentState.CountFarEnd << "\t" << "Far+: "
				<< currentState.CountFarEndPlus << "\tFar-: "
				<< currentState.CountFarEndMinus << std::endl << std::endl; //endl;
		delete[] BD_INDEX;
	}

	/* 3.2.2.2 search breakpoints of deletion */
	std::cout << "Searching breakpoints of deletion" << std::endl;
	for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {
		currentState.CountFarEnd = 0;
		currentState.CountFarEndMinus = 0;
		currentState.CountFarEndPlus = 0;

#pragma omp parallel default(shared)
		{
#pragma omp for
			for (int ReadIndex = 0; ReadIndex
			       < (int)currentState.Reads.size(); ReadIndex++) { // openMP 2.5 requires signed loop
				if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
					continue;
				}
				GetFarEnd_SingleStrandDownStream(currentState.CurrentChr,
						currentState.Reads[ReadIndex], RangeIndex);
				if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
					if (currentState.Reads[ReadIndex]. UP_Far[currentState.Reads[ReadIndex].UP_Far.size()
							- 1].LengthStr
							+ currentState.Reads[ReadIndex].CloseEndLength
							< currentState.Reads[ReadIndex].ReadLength) {
						if (currentState.Reads[ReadIndex].UP_Far_backup.size()) {
							if (currentState.Reads[ReadIndex]. UP_Far_backup[currentState.Reads[ReadIndex]. UP_Far_backup.size()
									- 1].LengthStr
									< currentState.Reads[ReadIndex].UP_Far[currentState.Reads[ReadIndex]. UP_Far.size()
											- 1].LengthStr) {
								currentState.Reads[ReadIndex].UP_Far_backup
										= currentState.Reads[ReadIndex].UP_Far;
								currentState.Reads[ReadIndex].UP_Far.clear();
							} else
								currentState.Reads[ReadIndex].UP_Far.clear();
						} else {
							currentState.Reads[ReadIndex].UP_Far_backup
									= currentState.Reads[ReadIndex].UP_Far;
							currentState.Reads[ReadIndex].UP_Far.clear();
						}
					} else {
#pragma omp critical
						{
							currentState.CountFarEnd++;
							if (currentState.Reads[ReadIndex].MatchedD == Plus)
								currentState.CountFarEndPlus++;
							else
								currentState.CountFarEndMinus++;
						}
					}
				}
			}
		} // #pragma omp parallel default(shared)

		std::cout << RangeIndex << "\tNumber of reads with far end mapped: "
				<< currentState.CountFarEnd << "\t" << "Far+: "
				<< currentState.CountFarEndPlus << "\tFar-: "
				<< currentState.CountFarEndMinus << std::endl;
	}

	/* 3.2.2.3 search breakpoints of short insertion */
	std::cout << "Searching breakpoints of short insertions" << std::endl;
	for (short RangeIndex = 1; RangeIndex < 2; RangeIndex++) {
		currentState.CountFarEnd = 0;
		currentState.CountFarEndMinus = 0;
		currentState.CountFarEndPlus = 0;

#pragma omp parallel default(shared)
		{
#pragma omp for
			for (int ReadIndex = 0; ReadIndex
			       < (int)currentState.Reads.size(); ReadIndex++) { // OpenMP 2.5 requires signed looop
				if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
					continue;
				}

				GetFarEnd_SingleStrandDownStreamInsertions(
						currentState.CurrentChr, currentState.Reads[ReadIndex],
						RangeIndex);
#pragma omp critical
				{
					if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
						currentState.CountFarEnd++;
						if (currentState.Reads[ReadIndex].MatchedD == Plus)
							currentState.CountFarEndPlus++;
						else
							currentState.CountFarEndMinus++;
					}
				}
			}
		}

		std::cout << RangeIndex << "\tNumber of reads with far end mapped: "
				<< currentState.CountFarEnd << "\t" << "Far+: "
				<< currentState.CountFarEndPlus << "\tFar-: "
				<< currentState.CountFarEndMinus << std::endl;
	}

	/* 3.2.2.4 search breakpoints of tandem duplication */
	if (Analyze_TD) {

		std::cout << "Searching breakpoints of tandem duplications"
				<< std::endl;
		for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {

			currentState.CountFarEnd = 0;
			currentState.CountFarEndMinus = 0;
			currentState.CountFarEndPlus = 0;
#pragma omp parallel default(shared)
			{
#pragma omp for
				for (int ReadIndex = 0; ReadIndex
				       < (int)currentState.Reads.size(); ReadIndex++) { // OpenMP 2.5 requires signed loop
					if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
						continue;
					}

					GetFarEnd_SingleStrandUpStream(currentState.CurrentChr,
							currentState.Reads[ReadIndex], RangeIndex);
					if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
						if (currentState.Reads[ReadIndex]. UP_Far[currentState.Reads[ReadIndex].UP_Far.size()
								- 1].LengthStr
								+ currentState.Reads[ReadIndex].CloseEndLength
								>= currentState.Reads[ReadIndex].ReadLength) {
							if (currentState.Reads[ReadIndex].MatchedD == Plus) {
								if (currentState.Reads[ReadIndex].UP_Close[0].AbsLoc
										< currentState.Reads[ReadIndex].ReadLength
												+ currentState.Reads[ReadIndex].UP_Far[0].AbsLoc)
									currentState.Reads[ReadIndex].UP_Far.clear();
							} else { // if (currentState.Reads[ReadIndex].MatchedD == Minus)
								if (currentState.Reads[ReadIndex].UP_Far[0].AbsLoc
										< currentState.Reads[ReadIndex].ReadLength
												+ currentState.Reads[ReadIndex].UP_Close[0].AbsLoc)
									currentState.Reads[ReadIndex].UP_Far.clear();
							}
						} else {
							if (currentState.Reads[ReadIndex].UP_Far_backup.size()) {
								if (currentState.Reads[ReadIndex]. UP_Far_backup[currentState.Reads[ReadIndex]. UP_Far_backup.size()
										- 1].LengthStr
										< currentState.Reads[ReadIndex]. UP_Far[currentState.Reads[ReadIndex].UP_Far. size()
												- 1].LengthStr) {
									currentState.Reads[ReadIndex].UP_Far_backup
											= currentState.Reads[ReadIndex].UP_Far;
									currentState.Reads[ReadIndex].UP_Far.clear();
								} else
									currentState.Reads[ReadIndex].UP_Far.clear();
							} else {
								currentState.Reads[ReadIndex].UP_Far_backup
										= currentState.Reads[ReadIndex].UP_Far;
								currentState.Reads[ReadIndex].UP_Far.clear();
							}
						}
#pragma omp critical
						{
							if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
								currentState.CountFarEnd++;
								if (currentState.Reads[ReadIndex].MatchedD
										== Plus)
									currentState.CountFarEndPlus++;
								else
									currentState.CountFarEndMinus++;
							}
						}
					}
				}
			}
			std::cout << RangeIndex
					<< "\tNumber of reads with far end mapped: "
					<< currentState.CountFarEnd << "\t" << "Far+: "
					<< currentState.CountFarEndPlus << "\tFar-: "
					<< currentState.CountFarEndMinus << std::endl;
		}
	} // if (Analyze_TD)

	/* 3.2.2.5 search breakpoints of inversion */

	if (Analyze_INV) {
		std::cout << "Searching breakpoints of inversion" << std::endl;
		for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {

			currentState.CountFarEnd = 0;
			currentState.CountFarEndMinus = 0;
			currentState.CountFarEndPlus = 0;
#pragma omp parallel default(shared)
			{
#pragma omp for
				for (int ReadIndex = 0; ReadIndex
				       < (int)currentState.Reads.size(); ReadIndex++) { // OpenMP 2.5 requires signed loop
					if (currentState.Reads[ReadIndex].Used
							|| currentState.Reads[ReadIndex].UP_Far.size()) {
						continue;
					}

					GetFarEnd_OtherStrand(currentState.CurrentChr,
							currentState.Reads[ReadIndex], RangeIndex);

					if (!currentState.Reads[ReadIndex].UP_Far.empty()) {
						if (currentState.Reads[ReadIndex]. UP_Far[currentState.Reads[ReadIndex].UP_Far.size()
								- 1].LengthStr
								+ currentState.Reads[ReadIndex].CloseEndLength
								< currentState.Reads[ReadIndex].ReadLength) {
							if (currentState.Reads[ReadIndex].UP_Far_backup.size()) {
								if (currentState.Reads[ReadIndex]. UP_Far_backup[currentState.Reads[ReadIndex]. UP_Far_backup.size()
										- 1].LengthStr
										< currentState.Reads[ReadIndex]. UP_Far[currentState.Reads[ReadIndex].UP_Far. size()
												- 1].LengthStr) {
									currentState.Reads[ReadIndex].UP_Far_backup
											= currentState.Reads[ReadIndex].UP_Far;
									currentState.Reads[ReadIndex].UP_Far.clear();
								} else
									currentState.Reads[ReadIndex].UP_Far.clear();
							} else {
								currentState.Reads[ReadIndex].UP_Far_backup
										= currentState.Reads[ReadIndex].UP_Far;
								currentState.Reads[ReadIndex].UP_Far.clear();
							}
						} else {
#pragma omp critical
							{
								currentState.CountFarEnd++;
								if (currentState.Reads[ReadIndex].MatchedD
										== Plus)
									currentState.CountFarEndPlus++;
								else
									currentState.CountFarEndMinus++;
							}
						}
					}
				}
			}
			std::cout << RangeIndex
					<< "\tNumber of reads with far end mapped: "
					<< currentState.CountFarEnd << "\t" << "Far+: "
					<< currentState.CountFarEndPlus << "\tFar-: "
					<< currentState.CountFarEndMinus << std::endl;
		}
	} // if (Analyze_INV)

	/* 3.2.2.6 compare backup with current value */
	// TODO: check with Kai what is the use of following
	// compare backup with current value
	std::cout << "revisit all breakpoints identified ...";
	for (unsigned ReadIndex = 0; ReadIndex < currentState.Reads.size(); ReadIndex++) {
		if (currentState.Reads[ReadIndex].UP_Far.empty()) {
			if (!currentState.Reads[ReadIndex].UP_Far_backup.empty()) {
				currentState.Reads[ReadIndex].UP_Far
						= currentState.Reads[ReadIndex].UP_Far_backup;
			}
		} else if (!currentState.Reads[ReadIndex].UP_Far_backup.empty()) {
			if (currentState.Reads[ReadIndex]. UP_Far_backup[currentState.Reads[ReadIndex].UP_Far_backup.size()
					- 1].LengthStr
					> currentState.Reads[ReadIndex].UP_Far[currentState.Reads[ReadIndex].UP_Far. size()
							- 1].LengthStr) {
				currentState.Reads[ReadIndex].UP_Far
						= currentState.Reads[ReadIndex].UP_Far_backup;
			}
		}
	}
	std::cout << " done." << std::endl;

	return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {

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

	returnValue = init(argc, argv, currentState);

	if (returnValue != EXIT_SUCCESS) {
		return returnValue;
	}

	std::string emptystr;

	/* 3 loop over chromosomes. this is the most outer loop in the main function. */
	while (currentState.SpecifiedChrVisited == false && currentState.inf_Seq
			>> currentState.CurrentChrName && !currentState.inf_Seq.eof()) {

		/* 3.1 preparation starts */
		std::cout << "Processing chromosome: " << currentState.CurrentChrName
				<< std::endl;

		//TODO: check with Kai what's the use of this line.
		std::getline(currentState.inf_Seq, emptystr);

		if (currentState.loopOverAllChromosomes) {
			GetOneChrSeq(currentState.inf_Seq, currentState.CurrentChr, true);
			currentState.WhichChr = currentState.CurrentChrName;
		} else if (currentState.CurrentChrName == currentState.WhichChr) { // just one chr and this is the correct one
			GetOneChrSeq(currentState.inf_Seq, currentState.CurrentChr, true);
			currentState.SpecifiedChrVisited = true;
		} else { // not build up sequence
			GetOneChrSeq(currentState.inf_Seq, currentState.CurrentChr, false);
			std::cout << "Skipping chromosome: " << currentState.CurrentChrName
					<< std::endl;
			continue;
		}

		CONS_Chr_Size = currentState.CurrentChr.size() - 2 * SpacerBeforeAfter;
		std::cout << "Chromosome Size: " << CONS_Chr_Size << std::endl;
		CurrentChrMask.resize(currentState.CurrentChr.size());
		for (unsigned int i = 0; i < currentState.CurrentChr.size(); i++) {
			CurrentChrMask[i] = 'N';
		}
		unsigned NumBoxes = (unsigned) (currentState.CurrentChr.size()
				/ BoxSize) + 1; // box size
		std::cout << NumBoxes << "\t" << BoxSize << std::endl;

		/* 3.1 preparation ends */

		/* 3.2 apply sliding windows to input datasets starts. This is the 2nd level while loop */
		g_binIndex = -1; // to start with 0...
		int
				displayedStartOfRegion =
						((currentState.startOfRegion >= 0) ? (currentState.startOfRegion
								- WINDOW_SIZE)
								: currentState.lowerBinBorder);
		int displayedEndOfRegion = displayedStartOfRegion + WINDOW_SIZE;

		do {

			/* 3.2.1 preparation starts */
			g_binIndex++;
			currentState.lowerBinBorder += WINDOW_SIZE;
			currentState.upperBinBorder += WINDOW_SIZE;
			displayedStartOfRegion += WINDOW_SIZE;
			displayedEndOfRegion += WINDOW_SIZE;
			if (displayedEndOfRegion > currentState.endOfRegion) {
				displayedEndOfRegion = currentState.endOfRegion;
			}

			// if the region end is specified, and it is before the regular upper border of the bin
			if (currentState.endRegionPlusBuffer > -1
					&& currentState.upperBinBorder
							> currentState.endRegionPlusBuffer) {
				currentState.upperBinBorder = currentState.endRegionPlusBuffer;
			}

			if (displayedStartOfRegion < displayedEndOfRegion) {
				std::cout << "Looking at chromosome " << currentState.WhichChr
						<< " bases " << displayedStartOfRegion << " to "
						<< displayedEndOfRegion << "." << std::endl;
			} else {
				std::cout
						<< "Checking out reads near the borders of the specified regions for extra evidence."
						<< std::endl;
			}

			if (Time_Load_S == 0)
				Time_Load_S = time(NULL);
			short ReturnFromReadingReads;

			VectorTag.clear();
			if (currentState.BAMDefined) {
				ReturnFromReadingReads = 0;
				for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
					ReturnFromReadingReads = ReadInBamReads(
							currentState.bams_to_parse[i].BamFile.c_str(),
							currentState.WhichChr, &currentState.CurrentChr,
							currentState.Reads,
							currentState.bams_to_parse[i].InsertSize,
							currentState.bams_to_parse[i].Tag,
							currentState.lowerBinBorder,
							currentState.upperBinBorder);
					if (ReturnFromReadingReads == 0) {
						std::cout << "Bam read failed: "
								<< currentState.bams_to_parse[i].BamFile
								<< std::endl;
						return 1;
					} else if (currentState.Reads.size() == 0) {
						std::cout << "No currentState.Reads for "
								<< currentState.WhichChr << " found in "
								<< currentState.bams_to_parse[i].BamFile
								<< std::endl;
					}
					std::cout << "BAM file index\t" << i << "\t"
							<< currentState.Reads.size() << std::endl;
				}

			}

			if (currentState.PindelReadDefined) {
				ReturnFromReadingReads = ReadInRead(
						currentState.inf_Pindel_Reads, currentState.WhichChr,
						currentState.CurrentChr, currentState.Reads,
						currentState.lowerBinBorder,
						currentState.upperBinBorder);
				if (ReturnFromReadingReads == 1) {
					std::cout << "malformed record detected!" << std::endl;
					return 1;
				} else if (currentState.Reads.size() == 0) {
					std::cout << "No reads found!?" << std::endl;
				}
			}
			Time_Mine_E = time(NULL);

			if (currentState.Reads.size())
				std::cout << "There are " << currentState.Reads. size()
						<< " reads for this chromosome." << std::endl;
			else {
				std::cout << "There are no reads for this bin." << std::endl;
				continue;
			}

			//TODO: to remove the following 3 lines?
			unsigned int Num_Left;
			Num_Left = currentState.Reads.size();
			Const_Log_T = log10((double) Num_Left);

			Time_Load_E = time(NULL);
			/* 3.2.1 preparation ends */

			returnValue = searchBreakPoints(currentState);

			returnValue = searchDeletions(currentState, NumBoxes);

			returnValue = searchIndels(currentState, NumBoxes);

			if (Analyze_TD) {
				returnValue = searchTandemDuplications(currentState, NumBoxes);
				returnValue = searchTandemDuplicationsNT(currentState, NumBoxes);
			}

			if (Analyze_INV) {
				returnValue = searchInversions(currentState, NumBoxes);
				returnValue = searchInversionsNT(currentState, NumBoxes);
			}

			returnValue = searchShortInsertions(currentState, NumBoxes);

			/* 3.2.8 report starts */
			int TotalNumReads = currentState.Reads.size();
			if (ReportCloseMappedRead) {
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
				if (!currentState.Reads[Index].UP_Far.empty())
					Count_Far++;
				if (!currentState.Reads[Index].UP_Far.empty()
						|| currentState.Reads[Index].Found) {

				} else {
					Count_Unused++;
				}
				if (currentState.Reads[Index].Used)
					Count_Used++;
			}

			std::cout << "Total: " << TotalNumReads << ";\tClose_end_found "
					<< TotalNumReads << ";\tFar_end_found " << Count_Far
					<< ";\tUsed\t" << Count_Used << "." << std::endl
					<< std::endl;
			std::cout << "For LI and BP: " << Count_Unused << std::endl
					<< std::endl;

			if (Analyze_LI) {
				time_t Time_LI_S, Time_LI_E;
				Time_LI_S = time(NULL);
				std::ofstream LargeInsertionOutf(
						currentState.LargeInsertionOutputFilename. c_str(),
						std::ios::app);
				SortOutputLI(currentState.CurrentChr, currentState.Reads,
						LargeInsertionOutf);
				LargeInsertionOutf.close();
				Time_LI_E = time(NULL);
				std::cout << "Mining, Sorting and output LI results: "
						<< (unsigned
						int) difftime(Time_LI_E, Time_LI_S) << " seconds."
						<< std::endl << std::endl;
				;
			}

			std::vector<SPLIT_READ> BP_Reads;
			BP_Reads.clear();
			if (Analyze_BP) {
				time_t Time_BP_S, Time_BP_E;
				Time_BP_S = time(NULL);
				std::ofstream RestOutf(currentState.RestOutputFilename.c_str(),
						std::ios::app);
				SortOutputRest(currentState.CurrentChr, currentState.Reads,
						BP_Reads, RestOutf);
				RestOutf.close();
				Time_BP_E = time(NULL);
				std::cout << "Mining, Sorting and output BP results: "
						<< (unsigned
						int) difftime(Time_BP_E, Time_BP_S) << " seconds."
						<< std::endl << std::endl;
			}

			Time_Sort_E = time(NULL);

			AllLoadings += (unsigned int) difftime(Time_Load_E, Time_Load_S);
			AllSortReport += (unsigned int) difftime(Time_Sort_E, Time_Load_E);
			currentState.Reads.clear();
			std::cout << "I have " << currentState.FutureReads. size()
					<< " reads saved for the next cycle." << std::endl;
			currentState.Reads.swap(currentState.FutureReads);
			Time_Load_S = 0;
			/* 3.2.8 report ends */

		} // do {
		while ((currentState.PindelReadDefined && !isFinishedPindel(
				currentState.upperBinBorder, currentState.endRegionPlusBuffer))
				|| (currentState.BAMDefined && !isFinishedBAM(
						currentState.upperBinBorder,
						currentState.endRegionPlusBuffer,
						currentState.CurrentChr.size())));
		/* 3.2 apply sliding windows to input datasets ends */

	} // while ( loopOverAllChromosomes && chromosomeIndex < chromosomes.size() );

	std::cout << "Loading genome sequences and reads: " << AllLoadings
			<< " seconds." << std::endl;
	std::cout << "Mining, Sorting and output results: " << AllSortReport
			<< " seconds." << std::endl;
	return 0;
} //main

std::vector<std::string> ReverseComplement(
		const std::vector<std::string> &InputPatterns) {
	std::vector < std::string > OutputPatterns; // = InputPatterns;
	unsigned int NumPattern = InputPatterns.size();
	OutputPatterns.resize(NumPattern);
	for (unsigned int i = 0; i < NumPattern; i++) {
		OutputPatterns[i] = ReverseComplement(InputPatterns[i]);
	}
	return OutputPatterns;
}

std::string Reverse(const std::string & InputPattern) {
	std::string OutputPattern = InputPattern;
	unsigned int LenPattern = InputPattern.size();
	for (unsigned int j = 0; j < LenPattern; j++)
		OutputPattern[j] = InputPattern[j];
	return OutputPattern;
}

std::string ReverseComplement(const std::string & InputPattern) {
	std::string OutputPattern = InputPattern;
	unsigned int LenPattern = InputPattern.size();
	for (unsigned int j = 0; j < LenPattern; j++)
		OutputPattern[j] = Convert2RC4N[(unsigned int) InputPattern[LenPattern
				- j - 1]];
	return OutputPattern;
}

// TODO: Ask Kai whether this can be removed
//const short MAX_SNP_ERROR = 2;
//const short ADDITIONAL_MISMATCH = 2;
//const short TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;


std::string Cap2Low(const std::string & input) {
	std::string output = input;
	for (unsigned int i = 0; i < output.size(); i++) {
		output[i] = Cap2LowArray[(short) input[i]];
	}
	return output;
}

bool NotInVector(const std::string & OneTag,
		const std::vector<std::string> &VectorTag) {
	for (unsigned i = 0; i < VectorTag.size(); i++) {
		if (OneTag == VectorTag[i])
			return false;
	}
	return true;
}

bool ReportEvent(const std::vector<SPLIT_READ> &Deletions,
		const unsigned int &S, const unsigned int &E) {
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
			//break; // TODO: Ask Kai whether this can be removed
		}
		if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size
				<= Min_Length) {
			RightMin = true;
			//break; // TODO: Ask Kai whether this can be removed
		}
		if (Deletions[i].BP >= Max_Length) {
			LeftMax = true;
			//break; // TODO: Ask Kai whether this can be removed
		}
		if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size
				>= Max_Length) {
			RightMax = true;
			//break; // TODO: Ask Kai whether this can be removed
		}
	}

	if (LeftMin && LeftMax && RightMin && RightMax)
		return true;
	else
		return false;
}

void GetRealStart4Deletion(const std::string & TheInput,
		unsigned int &RealStart, unsigned int &RealEnd) {
	unsigned int PosIndex = RealStart + SpacerBeforeAfter;
	unsigned int Start = PosIndex + 1;
	unsigned int End = RealEnd + SpacerBeforeAfter - 1;
	while (TheInput[PosIndex] == TheInput[End]) {
		--PosIndex;
		--End;
	}
	RealStart = PosIndex - SpacerBeforeAfter;
	PosIndex = RealEnd + SpacerBeforeAfter;
	while (TheInput[PosIndex] == TheInput[Start]) {
		++PosIndex;
		++Start;
	}
	RealEnd = PosIndex - SpacerBeforeAfter;
}

void GetRealStart4Insertion(const std::string & TheInput,
		const std::string & InsertedStr, unsigned int &RealStart,
		unsigned int &RealEnd) {
	unsigned int IndelSize = InsertedStr.size();
	unsigned int PosIndex = RealStart + SpacerBeforeAfter;
	// TODO: Ask Kai whether this can be removed
	//unsigned int Start = PosIndex + 1;

	//unsigned int End = RealEnd + SpacerBeforeAfter - 1;
	for (int i = IndelSize - 1; i >= 0; i--) {
		if (TheInput[PosIndex] == InsertedStr[i])
			PosIndex--;
		else
			break;
	}
	if (PosIndex == RealStart + SpacerBeforeAfter - IndelSize) {
		while (TheInput[PosIndex] == TheInput[PosIndex + IndelSize]) {
			PosIndex--;
		}
	}
	RealStart = PosIndex - SpacerBeforeAfter;
	PosIndex = RealEnd + SpacerBeforeAfter;
	for (unsigned int i = 0; i < IndelSize; i++) {
		if (TheInput[PosIndex] == InsertedStr[i])
			PosIndex++;
		else
			break;
	}
	if (PosIndex == RealEnd + SpacerBeforeAfter + IndelSize) {
		while (TheInput[PosIndex] == TheInput[PosIndex - IndelSize]) {
			PosIndex++;
		}
	}
	RealEnd = PosIndex - SpacerBeforeAfter;
}

std::vector<Region> Merge(const std::vector<Region> &AllRegions) {
	return AllRegions;
}

void GetCloseEnd(const std::string & CurrentChr, SPLIT_READ & Temp_One_Read) {

	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	std::string CurrentReadSeq;
	std::vector<unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(3 * Temp_One_Read.InsertSize);
	}
	std::vector<UniquePoint> UP;
	// TODO: Ask Kai whether this can be removed
	//char Direction;
	int Start, End;
	short BP_Start; // = MinClose;
	short BP_End; // = ReadLength - MinClose;

	// TODO: Ask Kai whether this can be removed
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//  PD[i].clear();
	//}
	//UP.clear();
	Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose;
	BP_End = Temp_One_Read.ReadLengthMinus;
	// TODO: Ask Kai whether this can be removed
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (Temp_One_Read.MatchedD == Plus) {
			CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
			Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			End = Start + 3 * Temp_One_Read.InsertSize;
			LeftChar = CurrentReadSeq[0];
			if (LeftChar != 'N') {
				// TODO: Ask Kai whether this can be removed
				//#pragma omp parallel default(shared)
				{
					// #pragma omp for
					for (int pos = Start; pos < End; pos++) {
						if (CurrentChr[pos] == LeftChar) {
							PD[0].push_back(pos);
						} else
							PD[1].push_back(pos);
					}
				}
			} else { //Match2N[(short)'A'] = 'N';
				// TODO: Ask Kai whether this can be removed
				//#pragma omp parallel default(shared)
				{
					//#pragma omp for
					for (int pos = Start; pos < End; pos++) {
						if (Match2N[(short) CurrentChr[pos]] == 'N') {
							PD[0].push_back(pos);
						}
						//else Left_PD[1].push_back(pos);
					}
				}
			}
			CheckLeft_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD,
					BP_Start, BP_End, FirstBase, UP); // LengthStr
			// TODO: Ask Kai whether this can be removed
			//Direction = Minus;
			if (UP.empty()) {
			} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
					>= Temp_One_Read.ReadLength) {
			} else {
				for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
					if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
							UP[LeftUP_index])) {
						Temp_One_Read.Used = false;
						Temp_One_Read.UP_Close.push_back(UP[LeftUP_index]);
					}
				}
			}

			UP.clear();
		} else if (Temp_One_Read.MatchedD == Minus) {
			CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
			End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			Start = End - 3 * Temp_One_Read.InsertSize;
			RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
			// TODO: Ask Kai whether this can be removed
			//cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << endl;
			CheckRight_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD,
					BP_Start, BP_End, FirstBase, UP);
			// TODO: Ask Kai whether this can be removed
			//cout << UP.size() << endl;
			//Direction = '+';
			if (UP.empty()) {
			} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
					>= Temp_One_Read.ReadLength) {
			} else {
				for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
					if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
							UP[RightUP_index])) {
						Temp_One_Read.Used = false;
						Temp_One_Read.UP_Close.push_back(UP[RightUP_index]);
					}
				}
			}
			UP.clear();
		}
	} else { // TOTAL_SNP_ERROR_CHECKED_Minus
		if (Temp_One_Read.MatchedD == Plus) {
			CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
			Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			End = Start + 3 * Temp_One_Read.InsertSize;
			LeftChar = CurrentReadSeq[0];
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
			CheckLeft_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD,
					BP_Start, BP_End, FirstBase, UP);
			// TODO: Ask Kai whether this can be removed
			//Direction = Minus;
			if (UP.empty()) {
			} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
					>= Temp_One_Read.ReadLength) {
			} else {
				for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
					if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
							UP[LeftUP_index])) {
						Temp_One_Read.Used = false;
						Temp_One_Read.UP_Close.push_back(UP[LeftUP_index]);
					}
				}
			}

			UP.clear();
		} else if (Temp_One_Read.MatchedD == Minus) {
			CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
			End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			Start = End - 3 * Temp_One_Read.InsertSize;
			RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
			// TODO: Ask Kai whether this can be removed
			//cout << "2" << PD[0].size() << "\t" << PD[1].size() << endl;
			CheckRight_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD,
					BP_Start, BP_End, FirstBase, UP);
			// TODO: Ask Kai whether this can be removed
			//cout << UP.size() << endl;
			//Direction = '+';
			if (UP.empty()) {
			} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
					>= Temp_One_Read.ReadLength) {
			} else {
				for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
					if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
							UP[RightUP_index])) {
						Temp_One_Read.Used = false;
						Temp_One_Read.UP_Close.push_back(UP[RightUP_index]);
					}
				}
			}
			UP.clear();
		}
	}
	// TODO: Ask Kai whether this can be removed
	//if (!Temp_One_Read.UP_Close.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Close);
	return;
}

void GetFarEnd_SingleStrandDownStreamInsertions(const std::string & CurrentChr,
		SPLIT_READ & Temp_One_Read, const short &RangeIndex) {
	// TODO: Ask Kai whether this can be removed
	//Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	//if (Temp_One_Read.MatchedRelPos > 160000)
	//cout << Temp_One_Read.UnmatchedSeq << Temp_One_Read.UnmatchedSeq.size() << endl;
	std::string CurrentReadSeq;
	std::vector<unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(
				Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		// TODO: Ask Kai whether this can be removed
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	std::vector<UniquePoint> UP;
	// TODO: Ask Kai whether this can be removed
	//char Direction;
	unsigned int Start, End;
	short BP_Start; // = MinClose;
	short BP_End; // = ReadLength - MinClose;

	// TODO: Ask Kai whether this can be removed
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//  PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + Range))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	// TODO: Ask Kai whether this can be removed
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.MatchedD == Minus) {
		char LeftChar;
		// TODO: Ask Kai whether this can be removed
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;

		End = Temp_One_Read.UP_Close[0].AbsLoc
				+ Temp_One_Read.UP_Close[0].LengthStr;
		if (End > SpacerBeforeAfter + Temp_One_Read.InsertSize * 2
				+ DSizeArray[RangeIndex])
			Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else
			Start = SpacerBeforeAfter;

		if (End > CurrentChr.size() - SpacerBeforeAfter)
			End = CurrentChr.size() - SpacerBeforeAfter;
		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		} else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start,
				BP_End, FirstBase, UP);
		if (UP.empty()) {
		} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
				>= Temp_One_Read.ReadLength) {
		} else {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[LeftUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
				}
			}
		}
		UP.clear();
	} else if (Temp_One_Read.MatchedD == Plus) {
		char RightChar;
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);

		Start = Temp_One_Read.UP_Close[0].AbsLoc
				- Temp_One_Read.UP_Close[0].LengthStr;
		// TODO: Ask Kai whether this can be removed
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter)
			End = CurrentChr.size() - SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		} else {
			if (RightChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		}

		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start,
				BP_End, FirstBase, UP);
		// TODO: Ask Kai whether this can be removed
		//Direction = '+';
		if (UP.empty()) {
		} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
				>= Temp_One_Read.ReadLength) {
		} else {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[RightUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
				}
			}
		}

		UP.clear();
	}
	// TODO: Ask Kai whether this can be removed
	//if (!Temp_One_Read.UP_Far.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Far);
	return;
}

void GetFarEnd_SingleStrandDownStream(const std::string & CurrentChr,
		SPLIT_READ & Temp_One_Read, const short &RangeIndex) {
	// TODO: Ask Kai whether this can be removed
	//cout << "in GetFarEnd_SingleStrandDownStream" << endl;
	//Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	//if (Temp_One_Read.MatchedRelPos > 160000)
	//cout << "RangeIndex: " << RangeIndex << "\t" << DSizeArray[RangeIndex] << endl;
	//if (RangeIndex == 4)
	//  cout << "\nin\t" << DSizeArray[RangeIndex] << "\t" << Temp_One_Read.UnmatchedSeq << "\t" << Temp_One_Read.UnmatchedSeq.size() << endl;
	std::string CurrentReadSeq;
	std::vector<unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(
				Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	std::vector<UniquePoint> UP;
	// TODO: Ask Kai whether this can be removed
	//char Direction;
	unsigned int Start, End;
	short BP_Start; // = MinClose;
	short BP_End; // = ReadLength - MinClose;

	// TODO: Ask Kai whether this can be removed
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//  PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + Range))/log(4.0) + 0.8) + 3;
	//cout << "Temp_One_Read.MinClose: " << Temp_One_Read.MinClose << endl;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	// TODO: Ask Kai whether this can be removed
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.MatchedD == Minus) {
		char LeftChar;
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;

		End = Temp_One_Read.UP_Close[0].AbsLoc
				+ Temp_One_Read.UP_Close[0].LengthStr
				- Temp_One_Read.ReadLength;
		if (End > SpacerBeforeAfter + Temp_One_Read.InsertSize * 2
				+ DSizeArray[RangeIndex])
			Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else
			Start = SpacerBeforeAfter;

		if (End > CurrentChr.size() - SpacerBeforeAfter)
			End = CurrentChr.size() - SpacerBeforeAfter;
		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		} else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		// TODO: Ask Kai whether this can be removed
		//if (Temp_One_Read.MatchedRelPos > 160000)
		//if (RangeIndex == 4)
		//cout << "- " << PD[0].size() << "\t" << PD[1].size() << endl;
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start,
				BP_End, FirstBase, UP);
		if (UP.empty()) {
		} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
				>= Temp_One_Read.ReadLength) {
		} else {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[LeftUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
				}
			}
		}

		UP.clear();
	} else if (Temp_One_Read.MatchedD == Plus) {
		char RightChar;
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);

		Start = Temp_One_Read.UP_Close[0].AbsLoc
				- Temp_One_Read.UP_Close[0].LengthStr
				+ Temp_One_Read.ReadLength;
		// TODO: Ask Kai whether this can be removed
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter)
			End = CurrentChr.size() - SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		} else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (RightChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		// TODO: Ask Kai whether this can be removed
		//if (Temp_One_Read.MatchedRelPos > 160000)
		//if (RangeIndex == 4)
		//cout << "+ " << PD[0].size() << "\t" << PD[1].size() << endl;
		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start,
				BP_End, FirstBase, UP);
		if (UP.empty()) {
		} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
				>= Temp_One_Read.ReadLength) {
		} else {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[RightUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
				}
			}
		}

		UP.clear();
	}
	// TODO: Ask Kai whether this can be removed
	//if (!Temp_One_Read.UP_Far.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Far);
	return;
}

void GetFarEnd_SingleStrandUpStream(const std::string & CurrentChr,
		SPLIT_READ & Temp_One_Read, const short &RangeIndex) {
	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	std::string CurrentReadSeq;
	std::vector<unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(
				Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		// TODO: Ask Kai whether this can be removed
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	std::vector<UniquePoint> UP;
	// TODO: Ask Kai whether this can be removed
	//char Direction;
	int Start, End;
	short BP_Start; // = MinClose;
	short BP_End; // = ReadLength - MinClose;

	// TODO: Ask Kai whether this can be removed
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//  PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.MatchedD == Minus) {
		char LeftChar;
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;

		Start = Temp_One_Read.UP_Close[0].AbsLoc
				+ Temp_One_Read.UP_Close[0].LengthStr;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > (int)CurrentChr.size() - (int)SpacerBeforeAfter)
			End = (int)CurrentChr.size() - (int)SpacerBeforeAfter;

		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		} else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start,
				BP_End, FirstBase, UP);
		// TODO: Ask Kai whether this can be removed
		//Direction = Minus;
		if (UP.empty()) {
		} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
				>= Temp_One_Read.ReadLength) {
		} else {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[LeftUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
				}
			}
		}

		UP.clear();
	} else if (Temp_One_Read.MatchedD == Plus) {
		char RightChar;
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);

		End = Temp_One_Read.UP_Close[0].AbsLoc
				- Temp_One_Read.UP_Close[0].LengthStr;
		// TODO: Ask Kai whether this can be removed
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;

		if (End > (int)DSizeArray[RangeIndex] + (int)Temp_One_Read.InsertSize * 2
				+ (int)SpacerBeforeAfter)
			Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else
			Start = SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		} else {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else PD[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start,
				BP_End, FirstBase, UP);
		// TODO: Ask Kai whether this can be removed
		//Direction = '+';
		if (UP.empty()) {
		} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
				>= Temp_One_Read.ReadLength) {
		} else {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[RightUP_index])) {
					Temp_One_Read.Used = false;
					Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
				}
			}
		}

		UP.clear();
	}
	// TODO: Ask Kai whether this can be removed
	//if (!Temp_One_Read.UP_Far.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Far);
	return;
}

void GetFarEnd_OtherStrand(const std::string & CurrentChr,
		SPLIT_READ & Temp_One_Read, const short &RangeIndex) {
	// TODO: Ask Kai whether this can be removed
	//short ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//short ReadLengthMinus = ReadLength - 1;
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 5 + Range * 2))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;
	int Start, End;
	short BP_Start = Temp_One_Read.MinClose + RangeIndex;
	short BP_End = Temp_One_Read.ReadLengthMinus;
	// TODO: Ask Kai whether this can be removed
	//char Direction;

	std::vector<UniquePoint> UP;
	std::vector<unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	std::vector<unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];

	if (Temp_One_Read.MatchedD == Plus) {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter
				- Temp_One_Read.ReadLength - 2 * Temp_One_Read.InsertSize
				- DSizeArray[RangeIndex];
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter
				- Temp_One_Read.ReadLength + 3 * Temp_One_Read.InsertSize
				+ DSizeArray[RangeIndex];
	} else {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - 3
				* Temp_One_Read.InsertSize + Temp_One_Read.ReadLength
				- DSizeArray[RangeIndex];
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter + 2
				* Temp_One_Read.InsertSize + Temp_One_Read.ReadLength
				+ DSizeArray[RangeIndex];
	}

	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(End - Start + 1);
		PD_Minus[CheckIndex].reserve(End - Start + 1);
	}

	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (Temp_One_Read.MatchedD == Plus) {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == CurrentBase)
						PD_Plus[0].push_back(pos);
					else
						PD_Plus[1].push_back(pos);
					// TODO: Ask Kai whether this can be removed
					//if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD_Plus[0].push_back(pos);
						// TODO: Ask Kai whether this can be removed
						//PD_Minus[0].push_back(pos);
					} else {
						PD_Plus[1].push_back(pos);
						// TODO: Ask Kai whether this can be removed
						//PD_Minus[1].push_back(pos);
					}
				}
			}
		} else { // -
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					// TODO: Ask Kai whether this can be removed
					//if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					if (CurrentChr[pos] == CurrentBaseRC)
						PD_Minus[0].push_back(pos);
					else
						PD_Minus[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						// TODO: Ask Kai whether this can be removed
						//PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					} else {
						// TODO: Ask Kai whether this can be removed
						//PD_Plus[1].push_back(pos);
						PD_Minus[1].push_back(pos);
					}
				}
			}
		}
	} else { // TOTAL_SNP_ERROR_CHECKED_Minus
		if (Temp_One_Read.MatchedD == Plus) {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == CurrentBase)
						PD_Plus[0].push_back(pos);
					// TODO: Ask Kai whether this can be removed
					//else PD_Plus[1].push_back(pos);
					//if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						PD_Plus[0].push_back(pos);
						// TODO: Ask Kai whether this can be removed
						//PD_Minus[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else {
					//  PD_Plus[1].push_back(pos);
					//PD_Minus[1].push_back(pos);
					//}
				}
			}
		} else { // -
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					// TODO: Ask Kai whether this can be removed
					//if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					if (CurrentChr[pos] == CurrentBaseRC)
						PD_Minus[0].push_back(pos);
					// TODO: Ask Kai whether this can be removed
					//else PD_Minus[1].push_back(pos);
				}
			} else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr[pos]] == 'N') {
						// TODO: Ask Kai whether this can be removed
						//PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					}
					// TODO: Ask Kai whether this can be removed
					//else {
					//PD_Plus[1].push_back(pos);
					//  PD_Minus[1].push_back(pos);
					//}
				}
			}
		}
	}

	CheckBoth(Temp_One_Read, CurrentChr, Temp_One_Read.UnmatchedSeq, PD_Plus,
			PD_Minus, BP_Start, BP_End, FirstBase, UP);
	if (UP.empty()) {
	} else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose
			>= Temp_One_Read.ReadLength) {
	} else {
		for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
			if (UP[UP_index].Direction == Plus) {
				// TODO: Ask Kai whether this can be removed
				//Direction = Minus;
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[UP_index]))
					Temp_One_Read.UP_Far.push_back(UP[UP_index]);
			} else {
				// TODO: Ask Kai whether this can be removed
				//Direction = '+';
				if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq,
						UP[UP_index]))
					Temp_One_Read.UP_Far.push_back(UP[UP_index]);
			}
		}
	}

	UP.clear();
	// TODO: Ask Kai whether this can be removed
	//if (Temp_One_Read.UP_Far.size()) {
	//Before = Temp_One_Read.UP_Far.size();
	//cout << "1: " << Temp_One_Read.UP_Far.size() << "\tt\t";
	//CleanUniquePoints(Temp_One_Read.UP_Far);
	//After = Temp_One_Read.UP_Far.size();
	//if (Before != After)
	//   cout << Before << "\t" << After << endl;
	//}
	return;
}

void GetFarEnd(const std::string & CurrentChr, SPLIT_READ & Temp_One_Read,
		const int &in_start, const int &in_end) {
	FarEndSearcher fes(CurrentChr, Temp_One_Read);
	fes.Get(in_start, in_end);
}

void CheckBoth(const SPLIT_READ & OneRead, const std::string & TheInput,
		const std::string & CurrentReadSeq,
		const std::vector<unsigned int> PD_Plus[],
		const std::vector<unsigned int> PD_Minus[], const short &BP_Start,
		const short &BP_End, const short &CurrentLength,
		std::vector<UniquePoint> &UP) {
	int Sum;
	if (CurrentLength >= BP_Start && CurrentLength <= BP_End) {
		// put it to LeftUP if unique
		for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			if (PD_Plus[i].size() + PD_Minus[i].size() == 1 && CurrentLength
					>= BP_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH)
					for (short j = 1; j <= ADDITIONAL_MISMATCH; j++)
						Sum += PD_Plus[i + j].size() + PD_Minus[i + j].size();

				if (Sum == 0 && i <= (short) (Seq_Error_Rate * CurrentLength
						+ 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					if (PD_Plus[i].size() == 1) {
						TempOne.Direction = FORWARD;
						TempOne.Strand = SENSE;
						TempOne.AbsLoc = PD_Plus[i][0];
					} else if (PD_Minus[i].size() == 1) {
						TempOne.Direction = BACKWARD;
						TempOne.Strand = ANTISENSE;
						TempOne.AbsLoc = PD_Minus[i][0];
					}
					TempOne.Mismatches = i;
					UP.push_back(TempOne);
					break;
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
		// TODO: Ask Kai whether this can be removed
		//const int SizeOfCurrent = Left_PD.size();
		// TODO: Ask Kai whether this can be removed
		//if (TOTAL_SNP_ERROR_CHECKED_Minus)
		{
			unsigned int pos;
			int SizeOfCurrent;
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
				SizeOfCurrent = PD_Plus[i].size();
				if (CurrentChar == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos = PD_Plus[i][j] + 1;
						if (Match2N[(short) TheInput[pos]] == 'N')
							PD_Plus_Output[i].push_back(pos);
						else
							PD_Plus_Output[i + 1].push_back(pos);
					}
				} else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos = PD_Plus[i][j] + 1;
						if (TheInput[pos] == CurrentChar)
							PD_Plus_Output[i].push_back(pos);
						else
							PD_Plus_Output[i + 1].push_back(pos);
					}
				}
				SizeOfCurrent = PD_Minus[i].size();
				if (CurrentCharRC == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos = PD_Minus[i][j] - 1;
						if (Match2N[(short) TheInput[pos]] == 'N')
							PD_Minus_Output[i].push_back(pos);
						else
							PD_Minus_Output[i + 1].push_back(pos);
					}
				} else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos = PD_Minus[i][j] - 1;
						if (TheInput[pos] == CurrentCharRC)
							PD_Minus_Output[i].push_back(pos);
						else
							PD_Minus_Output[i + 1].push_back(pos);
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
			} else {
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
			} else {
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
				CheckBoth(OneRead, TheInput, CurrentReadSeq, PD_Plus_Output,
						PD_Minus_Output, BP_Start, BP_End, CurrentLengthOutput,
						UP);
			} else
				return;
		}
	} else
		return;

}

void CleanUniquePoints(std::vector<UniquePoint> &Input_UP) {
	std::vector<UniquePoint> TempUP; //vector <UniquePoint> UP_Close; UP_Far
	UniquePoint LastUP = Input_UP[Input_UP.size() - 1];
	// TODO: Ask Kai whether this can be removed
	//TempUP.push_back(LastUP);
	char LastDirection = LastUP.Direction;
	char LastStrand = LastUP.Strand;
	unsigned int Terminal;

	if (LastDirection == FORWARD) {
		Terminal = LastUP.AbsLoc - LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand
					== LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc - Input_UP[i].LengthStr)
					TempUP.push_back(Input_UP[i]);
			}
		}
	} else if (LastDirection == BACKWARD) {
		Terminal = LastUP.AbsLoc + LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand
					== LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc + Input_UP[i].LengthStr)
					TempUP.push_back(Input_UP[i]);
			}
		}
	}
	Input_UP.clear();
	Input_UP = TempUP;
}
