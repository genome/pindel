#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <list>
#include <set>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <time.h>
#include <bits/basic_string.h>
//added by charris
#include <getopt.h>
#include "bam.h"
#include "sam.h"
#include "kstring.h"
#include "kseq.h"
#include "khash.h"
#include "ksort.h"
#include <faidx.h>
#include <assert.h>
#include <omp.h>

using namespace std;

/*v EWL update 0.0.1, April 8th 2011; can use the -c option with specified regions, and produces LI output that can be read by vcfcreator */

//init hash/maps for read pairing on the fly
KSORT_INIT_GENERIC(uint32_t)
KHASH_MAP_INIT_STR(read_name, bam1_t*)    

int g_binIndex = -1; // global variable for the bin index, as I cannot easily pass an extra parameter to the diverse functions
int g_maxPos = -1; // to calculate which is the last position in the chromosome, and hence to calculate the number of bins


    //end charris add
    //#include <omp.h>

    const string Pindel_Version_str = "Pindel version 0.2.2, April 8 2011.";

    const unsigned SpacerBeforeAfter = 10000000;
    //unsigned int DSizeArray[15];

    short Before, After;

/* 'Parameter' stores an individual parameter; it is set by the command line parameters, and used by the program. */
class Parameter {
public:
	bool isRequired() const { return d_required; }
	void describe() const;
	string getDescription() const { return d_description; }	
	string getShortName() const { return d_shortName; }
	string getLongName() const { return d_longName; }
	bool hasName( const string& name ) const { return ( d_shortName.compare(name) == 0 || d_longName.compare( name ) == 0 ); }
	bool isSet() const {return d_isSet; }
	virtual void setValue(const string& value ) {cout << "WHAT!" << endl ;};
	virtual void setValue(const int value) {};
	virtual void setValue(const double value) {};
	virtual void setValue(const bool value) {};
	virtual int getIValue() const {};
	virtual bool getBValue() const {};
	virtual string getSValue() const {};
	virtual double getFValue() const {};

	virtual bool isUnary() const { return false; }

	Parameter( const string& shortName, const string& longName, const string& description, const bool required);
protected:
	void set() { d_isSet = true; } //cout << "setting " << d_shortName << endl; }


private:
	bool d_required;
	string d_shortName;
	string d_longName;
	string d_description;
	bool d_isSet;
	static const int d_DESCRIBE_WIDTH = 11;
	static const int d_MAX_LINE_WIDTH = 80;
	const string makeNiceLine( const string& rawDescription ) const;
	const string describeTab() const { //return string(' ',d_DESCRIBE_WIDTH); 
		string descriptionTab=""; for (int i=0; i<d_DESCRIBE_WIDTH; i++) { descriptionTab+=" ";} return descriptionTab; }
};

/* 'getWord' removes the first word from the line, a bit like perl's unshift. */
string getWord(string& line)
{
	string head, tail;
	int endPos = line.find(" ");
	if (endPos == line.npos ) {
		head = line;
		tail = "";
	} else {
		head = line.substr(0, endPos);
		tail = line.substr(endPos + 1 );
	}

	line = tail;
	return head;
}

/* 'makeNiceLine' adds \n's on the right position of the description to get it fitting nicely in a window. */

const string Parameter::makeNiceLine( const string& rawDescription ) const
{
	string neatLine = describeTab();
	string words = rawDescription;
	const int LIMIT = d_MAX_LINE_WIDTH - d_DESCRIBE_WIDTH;
	int lineSize = 0;
	while (words.size()>0) {
		string newWord = getWord( words );
		if (newWord.size() + lineSize > LIMIT ) {
			neatLine += '\n';
			neatLine += describeTab();
			lineSize = 0;
		}
		neatLine += newWord + " ";
		lineSize += newWord.size()+1; // plus space!
	}

	return neatLine;
}

class IntParameter: public Parameter {
public:
	IntParameter( int* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const int value );

	virtual int getIValue() const { return *d_data_ptr; }
	virtual void setValue(const string& value);
	virtual void setValue(const int value);
private:
	int* d_data_ptr;
};

IntParameter::IntParameter( int* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const int value ) :
	Parameter( shortName, longName, description, required ), d_data_ptr( par_ptr )
{
	*d_data_ptr = value;
}

void IntParameter::setValue(const string& value)
{
	setValue( atoi( value.c_str() ));
}

void IntParameter::setValue(const int value)
{
	*d_data_ptr = value;
	set();
}

class BoolParameter: public Parameter {
public:
	BoolParameter( bool* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const bool value );

	virtual bool getBValue() const { return *d_data_ptr; }
	virtual void setValue(const string& value);
	virtual void setValue(const bool value);
	virtual bool isUnary() const { return true; }
private:
	bool* d_data_ptr;
};

BoolParameter::BoolParameter( bool* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const bool value ) :
	Parameter( shortName, longName, description, required )
{
	d_data_ptr = par_ptr;
	*d_data_ptr = value;
}

void BoolParameter::setValue(const string& value)
{
	char firstChar = tolower( value[0] );
	setValue( ( firstChar == 'f' || firstChar == '0' ) ? false : true );
}

void BoolParameter::setValue(const bool value)
{
	*d_data_ptr = value;
	set();
}

class FloatParameter: public Parameter {
public:
	FloatParameter( double* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const double value );

	double getFValue() const { return *d_data_ptr; }
	void setValue(const string& value);
	void setValue(const double value);
private:
	double* d_data_ptr;
};

FloatParameter::FloatParameter( double* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const double value ) :
	Parameter( shortName, longName, description, required )
{
	d_data_ptr = par_ptr;
	*d_data_ptr = value;
}

void FloatParameter::setValue(const string& value)
{
	setValue( atof( value.c_str() ));	
}

void FloatParameter::setValue(const double value)
{
	*d_data_ptr = value;
	set();
}

class StringParameter: public Parameter {
public:
	StringParameter( string* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const string& value );

	string getSValue() const { return *d_data_ptr; }
	void setValue(const string& value);
private:
	string* d_data_ptr;
};

StringParameter::StringParameter( string* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const string& value ) :
	Parameter( shortName, longName, description, required )
{
	d_data_ptr = par_ptr;
	*d_data_ptr = value;
}

void StringParameter::setValue(const string& value)
{
	*d_data_ptr = value;
	set();
}



Parameter::Parameter( const string& shortName, const string& longName, const string& description, const bool required )
{
	d_required = required;
	d_shortName = shortName;
	d_longName = longName;
	d_description = description;
	d_isSet = false;
}

void Parameter::describe() const 
{ 
	for (int i=0; i<d_DESCRIBE_WIDTH; i++)	{ cout << " ";}
	cout << d_shortName << "/";
	cout << d_longName << endl;
	
	//for (int i=0; i<d_DESCRIBE_WIDTH; i++)	{ cout << " ";}
	cout << makeNiceLine( d_description );
	//if ( d_required ) { cout << " required parameter" ; }
	cout << endl << endl; 
}


typedef struct {
	string referenceFileName;
	string pindelFileName;
	string bamConfigFileName;
	string outputFileName;
	string breakdancerFileName;
	int numThreads;
	bool showHelp;
} ParCollection;
ParCollection par;


    unsigned int CountIndels = 0;
    const int alphs = 4;
    const char alphabet[alphs] = {'A','C','G','T'};
unsigned long long int TheMax = 0;
const short MAX_MISMATCHES = 4;
float ExtraDistanceRate = 0.1;

double Const_Log_T = 0.0;
double Const_S = 0.0;
double LOG14 = log10(0.25);
//double Const_I = 0.0;
unsigned int BoxSize = 10000; // 10k is fine
const double Min_Filter_Ratio = 0.5;
unsigned int SPACERSIZE = 1;
unsigned int OriginalNumRead = 0; 
const string NonACGT = "$";
short MIN_Len_Match = 4;
unsigned int NumberOfSIsInstances = 0;
unsigned int NumberOfDeletionsInstances = 0;
unsigned int NumberOfDIInstances = 0;
unsigned int NumberOfInvInstances = 0;
unsigned int NumberOfTDInstances = 0;
short ReportLength = 80;
vector <string> VectorTag;
char Match[256];
char Match2N[256];
char Convert2RC[256];
char Convert2RC4N[256];
char Cap2LowArray[256];
bool FirstChr = true;
const double InsertSizeExtra = 2;
unsigned int CONS_Chr_Size;
const string N_str = "N";
const char N_char = 'N';
const char BD_char = 'A';
unsigned int DSizeArray[15];
const short Max_uint8_t = 100;

string BreakDancerMask;
string CurrentChrMask;

const char Plus = '+';
const char Minus = '-';
const char FirstCharReadName = '@';

unsigned int NumReadScanned = 0;
unsigned int NumReadInChr = 0;
unsigned int InChrPlus = 0;
unsigned int InChrMinus = 0;
unsigned int GetPlus = 0;
unsigned int GetMinus = 0;

//short MAX_SNP_ERROR = 2;

const short Max_short = 128;

//short TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
//short TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;

//short MAX_ALLOWED_MISMATCHES = TOTAL_SNP_ERROR_CHECKED_Minus + 5;

vector <Parameter*> parameters;

// #########################################################
int ADDITIONAL_MISMATCH = 1; // user
int Min_Perfect_Match_Around_BP = 3; // user                   //#
int MIN_IndelSize_NT = 50;            //user            //#  
int MIN_IndelSize_Inversion = 50;          //user       //#
double Seq_Error_Rate = 0.05;           // user            //# 
//float Seq_Error_Rate_1 = 0.05;                         //# 
//float Seq_Error_Rate_2 = 0.02;                         //#
//float Seq_Error_Rate_3 = 0.00;                         //#
unsigned int BalanceCutoff = 100;      //                    //#
//short RangeMaxSensivity = 9;       // 3                //#
//short RangeMediumSensivity = 9;    // 5                //# 
//short RangeLowSensivity = 7;                           //#
bool Analyze_INV = true;             //# user
bool Analyze_TD = true;             //# user
bool Analyze_LI = true; //EWL070111 user
bool Analyze_BP = true; //EWL070111 user
int NumRead2ReportCutOff = 3;           //#
int NumRead2ReportCutOff_BP = 2;
int MaxRangeIndex = 9;// 5 or 6 or 7 or maximum 8      //# // user
double MaximumAllowedMismatchRate = 0.1;            //#  // user
int Min_Num_Matched_Bases = 30;                  //# // user
int Max_Length_NT = 30;  // user
const int NumberOfReadsPerBuffer = 1000;  // estimate later
bool ReportCloseMappedRead = false; // user
const bool ReportSVReads = false; 
const bool ReportLargeInterChrSVReads = false; 
const short Indel_SV_cutoff = 50;
//bool RemoveDuplicates = false;
double FLOAT_WINDOW_SIZE = 10.0;
int WINDOW_SIZE = 10000000;
const int AROUND_REGION_BUFFER = 10000; // how much earlier reads should be selected if only a region of the chromosome needs be specified.
// #########################################################

//const float Double_Seq_Error_Rate_Per_Side = Seq_Error_Rate_Per_Side * 2;
unsigned int Distance = 300;
//short MinClose = 8;//short(log((double)Distance)/log(4.0) + 0.8) + 3 + MAX_SNP_ERROR;//atoi(argv[1]);
//short MinFar_I = MinClose + 1;//atoi(argv[2]);
//cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
short MinFar_D = 8;//atoi(argv[3]);
const short MaxDI = 30; 
const short FirstBase = 1;

//struct Fragment {
//   unsigned int Start; 
//   unsigned int End;  
//};

struct SupportPerSample {
	int NumPlus;
	int NumMinus;
	int NumUPlus;
	int NumUMinus;	
};

struct UniquePoint {
    short LengthStr;
    unsigned int AbsLoc; 
    char Direction; // forward reverse
    char Strand; // sense antisense
    short Mismatches;
};

const char SENSE = '+';
const char ANTISENSE = '-';
const char FORWARD = '+';
const char BACKWARD = '-';

// UP_Close[0].AbsLoc

struct SPLIT_READ {
	SPLIT_READ() { FragName=""; Name=""; UnmatchedSeq=""; HalfMapped="";HalfUnmapped=""; MatchedD=0; MatchedRelPos=0; MS=0;
					  InsertSize=0; Tag=""; ReadLength=0; ReadLengthMinus=0; MAX_SNP_ERROR=0; TOTAL_SNP_ERROR_CHECKED=0;
					  TOTAL_SNP_ERROR_CHECKED_Minus=0; MinClose=0; BP=0; Left=0; Right=0; BPLeft=0; BPRight=0; IndelSize=0; Unique=false;
					  score=0.0; InsertedStr="";  NT_str=""; NT_size=0; Used=false; CloseEndLength=0; Found=false; LeftMostPos = 0; }
	string FragName;
	string Name;
	string UnmatchedSeq;
	string HalfMapped;
	string HalfUnmapped;
	char MatchedD;
	int MatchedRelPos;
	short MS;
	short InsertSize;	
	string Tag;
	vector <UniquePoint> UP_Close; // partial alignment of the unmapped reads close to the mapped read
	vector <UniquePoint> UP_Far;
	vector <UniquePoint> UP_Far_backup;
	short ReadLength;
	short ReadLengthMinus;
	short MAX_SNP_ERROR;// = (short)(Temp_One_Read.UnmatchedSeq.size() * Seq_Error_Rate);
	
	short TOTAL_SNP_ERROR_CHECKED;// = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
	short TOTAL_SNP_ERROR_CHECKED_Minus;// = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
	short MinClose;
	short BP;
	int Left;
	int Right;
	int BPLeft;
	int BPRight;
	int IndelSize;
	bool Unique;
	double score;
	string InsertedStr;
	string NT_str;	
	short NT_size;
	//string NT_2str;
	//short NT_2size;
	bool Used;
	short CloseEndLength;
	bool Found;
	int LeftMostPos;
};

//add charris structs
struct flagshit {
	 flagshit() { mapped=false; unique=false; sw=false; edits=0; suboptimal=false; }
    bool mapped;
    bool unique;
    bool sw;
    int edits;
    bool suboptimal;
};


struct fetch_func_data { 
	 fetch_func_data() { read_to_map_qual=NULL; header=NULL; b1_flags=NULL; b2_flags=NULL; CurrentChr=NULL; Tag=""; InsertSize=0; }
    vector <SPLIT_READ>* LeftReads;
    khash_t(read_name) *read_to_map_qual;
    bam_header_t* header;
    flagshit* b1_flags;
    flagshit* b2_flags;
    string* CurrentChr;
    string Tag;
    int InsertSize;
};

struct bam_info {
	 bam_info() { BamFile=""; InsertSize=0; Tag="";}
    string BamFile;
    int InsertSize;
    string Tag;
};


struct LI_Pos {
	 LI_Pos() { Plus_Pos=0; Minus_Pos=0; WhetherReport=false; }
    unsigned Plus_Pos;
    unsigned Minus_Pos;
    vector <unsigned> Plus_Reads; // put index here
    vector <unsigned> Minus_Reads;
    bool WhetherReport;
};

struct Rest_Pos {
	 Rest_Pos() { Strand='X'; Pos=0; }
    char Strand;
    unsigned Pos;
    vector <unsigned> Pos_Reads; // put index here
};

struct Indel4output {
	 Indel4output() { BPLeft=0; BPRight=0; IndelSize=0; Start=0; End=0; RealStart=0; RealEnd=0; NT_size=0; WhetherReport=false; 
							IndelStr=""; Support=0; }
    unsigned int BPLeft;
    unsigned int BPRight;
    int IndelSize;
    unsigned int Start;
    unsigned int End;
    unsigned int RealStart;
    unsigned int RealEnd;
    short NT_size;
    bool WhetherReport;
    string IndelStr;
    short Support;
};

void ReadInOneChr(ifstream & inf_Seq, string & TheInput, const string & ChrName);
void GetOneChrSeq(ifstream & inf_Seq, string & CurrentChr, bool WhetherBuildUp);
short ReadInRead(ifstream & inf_Seq, const string & CurrentFragName, const string & CurrentFrag, vector <SPLIT_READ> & Reads, const int lowerBinBorder, const int upperBinBorder);
///////ADDED BY CHARRIS
void parse_flags_and_tags(const bam1_t *b, flagshit* flags); 
bool ReadInBamReads(const char* bam_path, const string & FragName, string * CurrentChr, vector <SPLIT_READ> & LeftReads, int InsertSize, string Tag, int binStart, int binEnd);
int32_t bam_cigar2len(const bam1_core_t *c, const uint32_t *cigar);
static int fetch_func(const bam1_t *b1, void *data); 
void build_record(const bam1_t* mapped_read, const bam1_t* unmapped_read, void* data); 
#ifdef __cplusplus
extern "C" {
    int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);
    int32_t bam_aux2i(const uint8_t *s);
    void bam_init_header_hash(bam_header_t *header);
}
#endif

///END CHARRIS ADD

void GetCloseEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read);
void GetFarEnd_BothStrands(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_OtherStrand(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_SingleStrandDownStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_SingleStrandUpStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
void GetFarEnd_SingleStrandDownStreamInsertions(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex);
//void GetFarEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const int & Range, const string & BreakDancerMask);
void GetFarEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const int & start, const int & end);
//void CheckConsistancy(vector <SPLIT_READ> & input);
struct BreakDancer {
	BreakDancer() { ChrName_A=""; ChrName_B=""; Size=0; Score=0; S1=0; S2=0; S3=0; S4=0; }
	string ChrName_A;
	string ChrName_B;
	string Type;
	int Size;
	int Score;
	unsigned S1;
	unsigned S2;
	unsigned S3;
	unsigned S4;
};

struct Region {
	 Region() { start=0; end=0; }
    unsigned start;
    unsigned end;
};


void OutputDeletions(const vector <SPLIT_READ> & Deletions, 
                     const string & TheInput, 
                     const unsigned int & C_S, 
                     const unsigned int & C_E,
		               const unsigned int & RealStart,
		               const unsigned int & RealEnd,
                     ofstream & DeletionOutf);

void OutputTDs(const vector <SPLIT_READ> & TDs, 
					const string & TheInput, 
					const unsigned int & C_S, 
					const unsigned int & C_E,
					const unsigned int & RealStart,
					const unsigned int & RealEnd,
					ofstream & TDOutf);


void OutputInversions(const vector <SPLIT_READ> & Inv, 
							 const string & TheInput, 
							 const unsigned int & C_S, 
							 const unsigned int & C_E,
							 const unsigned int & RealStart,
							 const unsigned int & RealEnd,
							 ofstream & InvOutf);

void OutputSIs(const vector <SPLIT_READ> & SIs,
               const string & TheInput,
               const unsigned int & C_S,
               const unsigned int & C_E,
	            const unsigned int & RealStart,
	            const unsigned int & RealEnd,
               ofstream & SIsOutf);



//void SortAndPutLI(vector <SPLIT_READ> & Input, vector <int> & Output, const string & CurrentFragName, const string & chr, const char & Direction, ofstream & LI_SR_Outf);
//void ProcessLIs(vector <SPLIT_READ> & LIs, ofstream & LIoutputfile);

short CompareTwoReads(const SPLIT_READ & First, const SPLIT_READ & Second);
vector <string> ReverseComplement(const vector <string> & input);
string Reverse(const string & InputPattern);
string ReverseComplement(const string & InputPattern);
string Cap2Low(const string & input);

void CheckLeft_Close(const SPLIT_READ & OneRead,
							const string & TheInput, 
							const string & CurrentReadSeq, 
							const vector <unsigned int> Left_PD[], 
							const short & BP_Left_Start, 
							const short & BP_Left_End, 
							const short & CurrentLength,
							vector <UniquePoint> & LeftUP);

void CheckRight_Close(const SPLIT_READ & OneRead,
							 const string & TheInput, 
							 const string & CurrentReadSeq, 
							 const vector <unsigned int> Right_PD[], 
							 const short & BP_Right_Start, 
							 const short & BP_Right_End, 
							 const short & CurrentPos, 
							 vector <UniquePoint> & RightUP);

void CheckLeft_Far(SPLIT_READ & OneRead,
						 const string & TheInput, 
						 const string & CurrentReadSeq, 
						 const vector <unsigned int> Left_PD[], 
						 const short & BP_Left_Start, 
						 const short & BP_Left_End, 
						 const short & CurrentLength,
						 vector <UniquePoint> & LeftUP);

void CheckRight_Far(SPLIT_READ & OneRead,
						  const string & TheInput, 
						  const string & CurrentReadSeq, 
						  const vector <unsigned int> Right_PD[], 
						  const short & BP_Right_Start, 
						  const short & BP_Right_End, 
						  const short & CurrentPos, 
						  vector <UniquePoint> & RightUP);


void CheckBoth(const SPLIT_READ & OneRead,
					const string & TheInput, 
               const string & CurrentReadSeq,
               const vector <unsigned int> PD_Plus[],
					const vector <unsigned int> PD_Minus[],
               const short & BP_Start,
               const short & BP_End,
               const short & CurrentLength,
               vector <UniquePoint> & UP);
bool CheckMismatches(const string & TheInput, 
                     const string & CurrentReadSeq, 
                     //const unsigned int & Start,
							const UniquePoint & UP);

//void ProcessLIs(vector <SPLIT_READ> & LIs, ofstream & LIoutputfile);
void SortOutputSI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> SIs[], ofstream & SIsOutf);
void SortOutputD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> Deletions[], ofstream & DeletionOutf);
void SortOutputTD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf);
void SortOutputTD_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf);
void SortOutputInv(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> Inv[], ofstream & InvOutf);
void SortOutputInv_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> Inv[], ofstream & InvOutf);
void SortOutputDI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> DI[], ofstream & DIOutf);
void SortOutputLI(const string & CurrentChr, vector <SPLIT_READ> & AllReads, ofstream & Outf_LI);
void SortOutputRest(const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <SPLIT_READ> & BP_Reads, ofstream & Outf_Rest);
bool NotInVector(const string & OneTag, const vector <string> & VectorTag);
void GetIndelTypeAndRealStart(const string & TheInput, const unsigned int & BPLeft, 
                              const unsigned int & IndelSize, const string & IndelStr, 
                              string & IndelType, unsigned int & RealStart, const bool & WhetherD);
void GetRealStart4Deletion(const string & TheInput, unsigned int & RealStart, unsigned int & RealEnd);

void GetRealStart4Insertion(const string & TheInput, const string & InsertedStr, unsigned int & RealStart, unsigned int & RealEnd);
bool ReportEvent(const vector <SPLIT_READ> & Deletions, const unsigned int & Pre_S, const unsigned int & Pre_E);

void CleanUniquePoints (vector <UniquePoint> & Input_UP);

short WhetherRemoveDuplicates;

string TempLine_DB_Unique;

vector <Region> Merge(const vector <Region> & AllRegions);

short CompareTwoString(const string & Str_A, const string & Str_B);

static struct option long_options[] =
{
    {"fasta", required_argument, 0, 'f'},
    {"config-file", required_argument, 0, 'i'},
	 {"pindel-file", required_argument, 0, 'p'},
    {"output-prefix", required_argument, 0, 'o'},
    {"chromosome", required_argument, 0, 'c'},
    {"breakdancer", required_argument, 0, 'b'},
};

bool readTransgressesBinBoundaries(SPLIT_READ &read, const int & upperBinBorder) 
{
   return (read.BPRight > upperBinBorder - 2 * read.InsertSize) ; 
}

/** 'readInSpecifiedRegion' if a region is specified, check if the read is in it. */
bool readInSpecifiedRegion(
	const SPLIT_READ & read,	// in: the read
	const int startOfRegion,	// in: the first base of the specified region
	const int endOfRegion 	   // in: the last base of the specified region (-1 if no region has been specified)
)
{
	bool passesFilter = true; 

	// if a start position has been defined, and the breakpoint is before it
	if ( ( startOfRegion != -1 ) && ( read.BPLeft+1<startOfRegion )) {
		passesFilter = false;
	}

	// if an end of the region has been specified
	if ( ( endOfRegion != -1 ) && (read.BPLeft+1 > endOfRegion ) ) {
		passesFilter = false;
		
	} 
	// no region specified, so all positions are okay
	return passesFilter;
}

void saveReadForNextCycle(SPLIT_READ &read, vector<SPLIT_READ> &futureReads) 
{
   futureReads.push_back(read);
   read.Used = true; // as it cannot be used for this round of analyses anymore
}

/* 'defineParameters' defines the parameters to be used by Pindel. Takes the variables from the calling function as argument for those variables which
	do not need to be stored in the par structure. */
void defineParameters(string& WhichChr)
{
	parameters.push_back(
		new StringParameter( &par.referenceFileName, "-f", "--fasta", "the reference genome sequences in fasta format", true, "" ) );
	parameters.push_back(
		new StringParameter( &par.pindelFileName, "-p", "--pindel-file", "the Pindel input file; \neither this or a bam configure file is required", false, "" ) );
	parameters.push_back(
		new StringParameter( &par.bamConfigFileName, "-i", "--config-file", "the bam config file; either this or a pindel input file is required. Per line: path and file name of bam, insert size and sample tag.     For example: /data/tumour.bam  400  tumour", false, "" ) );
	parameters.push_back(
		new StringParameter( &par.outputFileName, "-o", "--output-prefix", "Output prefix;", true, "" ) );
	parameters.push_back(
		new StringParameter( &WhichChr, "-c", "--chromosome", "Which chr/fragment. Pindel will process reads for one chromosome each time. ChrName must be the same as in reference sequence and in read file. '-c ALL' will make Pindel loop over all chromosomes. The search for indels and SVs can also be limited to a specific region; -c 20:10,000,000 will only look for indels and SVs after position 10,000,000 = [10M, end], -c 20:5,000,000-15,000,000 will report indels in the range between and including the bases at position 5,000,000 and 15,000,000 = [5M, 15M]", true, "" ) );	

	parameters.push_back(
		new BoolParameter( &par.showHelp, "-h", "--help", "show the command line options of Pindel", false, false ) );

	parameters.push_back( 
								new IntParameter( &par.numThreads, "-T", "--number_of_threads", "the number of threads Pindel will use (default 1).", false, 1 ) );
	
	parameters.push_back(
								new IntParameter( &MaxRangeIndex, "-x", "--max_range_index", "the maximum size of structural variations to be detected; the higher this number, the greater "
													  "the number of SVs reported, but the computational cost and memory requirements increase, as "
													  "does the rate of false "
													  "positives. 1=128, 2=512, 3=2,048, 4=8,092, 5=32,368, 6=129,472, 7=517,888, 8=2,071,552, 9=8,286,208. "
													  "(maximum 9, default 5)", false, 5 ) );
	parameters.push_back(
								new FloatParameter( &FLOAT_WINDOW_SIZE, "-w", "--window_size", "for saving RAM, divides the reference in bins of X million bases and only analyzes the reads that belong in the current "
														 "bin, "	
														 "(default 10 (=10 million))", false, 10 ) );
	

	parameters.push_back(
								new FloatParameter( &Seq_Error_Rate, "-e", "--sequencing_error_rate", "the expected fraction of sequencing errors "
														 "(default 0.05)", false, 0.05 ) );
	
	parameters.push_back(
								new FloatParameter( &MaximumAllowedMismatchRate, "-u", "--maximum_allowed_mismatch_rate", "Only reads with no less than this fraction of mismatches than the reference genome will be considered. "
														 "(default 0.1)", false, 0.1 ) );	
	
	parameters.push_back(
								new BoolParameter( &Analyze_INV, "-r", "--report_inversions", "report inversions "
														"(default true)", false, true ) );
	parameters.push_back(
								new BoolParameter( &Analyze_TD, "-t", "--report_duplications", "report tandem duplications "
														"(default true)", false, true ) );
	parameters.push_back(
								new BoolParameter( &Analyze_LI, "-l", "--report_long_insertions", "report insertions of which the full sequence cannot be deduced because of their length "
														"(default true)", false, true ) );
	parameters.push_back(
								new BoolParameter( &Analyze_BP, "-k", "--report_breakpoints", "report breakpoints "
														"(default true)", false, true ) );
	
	parameters.push_back(
								new BoolParameter( &ReportCloseMappedRead, "-s", "--report_close_mapped_reads", "report reads of which only one end (the one closest to the mapped read of the paired-end read) "
														"could be mapped. "
														"(default false)", false, false ) );
	
	parameters.push_back(
								new StringParameter( &par.breakdancerFileName, "-b", "--breakdancer", "Pindel is able to use calls from other SV methods such as BreakDancer to further increase sensitivity and specificity.                    BreakDancer result or calls from any methods must in the format:   ChrA LocA stringA ChrB LocB stringB other", false, "" ) );

	parameters.push_back(
		new IntParameter( &ADDITIONAL_MISMATCH, "-a", "--additional_mismatch", "Pindel will only map part of a read to the reference genome if there are no other candidate positions with no more than the specified number of mismatches position. The bigger tha value, the more accurate but less sensitive. (default value 1)", false, 1 ) );

	parameters.push_back(
		new IntParameter( &Min_Perfect_Match_Around_BP, "-m", "--min_perfect_match_around_BP", "at the point where the read is split into two, there should at least be "
																									  "this number of perfectly matching bases between read and reference (default value 3)"
																									  , false, 3 ) );
	parameters.push_back(
		new IntParameter( &MIN_IndelSize_NT, "-n", "--min_NT_size", "only report inserted (NT) sequences in deletions greater than this size "
																					   "(default 50)", false, 50 ) );
	parameters.push_back(
		new IntParameter( &MIN_IndelSize_Inversion, "-v", "--min_inversion_size", "only report inversions greater than this number of bases "
																					   "(default 50)", false, 50 ) );
	parameters.push_back(
								new IntParameter( &Min_Num_Matched_Bases, "-d", "--min_num_matched_bases", "only consider reads as evidence if they map with more than X bases to the reference. "
													  "(default 30)", false, 30 ) );
	
}

/* 'findParameter' returns the index of the parameter with name 'name'; -1 if not found.*/
int findParameter(string name)
{
	for (int parameterCounter=0; parameterCounter<parameters.size(); parameterCounter++ ) {
		if (parameters[ parameterCounter ]->hasName( name ) ) { return parameterCounter; }
	}	
	//cout << "Result of FindParameter is -1\n";
	return -1;
}

/* 'readParameters' reads the parameters as entered in the command line. */
void readParameters(int argc, char* argv[])
{
	//for (int argumentIndex=1; argumentIndex<argc; argumentIndex++ ) { cout << argumentIndex  << ". " << argv[argumentIndex] << endl; }

	for (int argumentIndex=1; argumentIndex<argc; argumentIndex++ ) {
		string currentArgument = argv[ argumentIndex ];

		//find argument in parameterlist
		int parameterIndex = findParameter( currentArgument );
		if ( parameterIndex == -1 ) { cout << "unknown argument: " << currentArgument << endl; return; }

		if ( parameters[ parameterIndex ]->isUnary() ) {
			 parameters[ parameterIndex ]->setValue(true); // default
			if ( (argumentIndex+1 < argc) && ( argv[ argumentIndex+1][0]!='-' ) ) { // so there are more arguments, and next one isn't regular -x
				if ( tolower(argv[ argumentIndex+1][0]) == 'f' || ( argv[ argumentIndex+1][0] =='0' ) ) {
					parameters[ parameterIndex ]->setValue(false); 
				}
				argumentIndex++; // in any case increase the argument index
			} 
		} else { // argument needs a parameter
			argumentIndex++; // move on to next argument in the list
			if (argumentIndex >= argc ) { cout << "argument of " << currentArgument << " lacking.\n"; return ; }
			if (argv[ argumentIndex ][0]=='-' ) { cout << "argument of " << currentArgument << " seems erroneous.\n"; return ; }
			// but if everything is allright, 
			//cout << "Giving " << currentArgument << " the value " << argv[ argumentIndex ] << endl;
			parameters[ parameterIndex ]->setValue( string(argv[ argumentIndex ]) );
		}
	}
}

/* isReadsFileParam returns whether the parameter points to a read-containing file, and is therefore required,
	even though not both are required. */
bool isReadsFileParam(Parameter* param)
{
	return ( param->hasName("-i") || param->hasName("-p") );
}

/* 'printHelp' prints all parameters available. */
void printHelp()
{
	//cout << "Please specify input, either bam configure file and/or pindel input format" << endl;
	//cout.width(0);	
	cout << "\nProgram:   pindel (detection of indels and structural variations)\n";
	cout << "Version:   0.2.2\n";
	cout << "Contact:   Kai Ye <k.ye@lumc.nl>\n\n";

	cout << "Usage:     pindel -f <reference.fa> -p <pindel_input>\n";
	cout << "           [and/or -i bam_configuration_file]\n";
	cout << "           -c <chromosome_name> -o <prefix_for_output_file>\n\n";

	cout << "Required parameters:\n";
	//parameters[1]->describe();
	for (int i=0; i<parameters.size(); i++) {
		if (parameters[i]->isRequired() || isReadsFileParam(parameters[i]) ) {
			parameters[i]->describe();
		}
	} 
	cout << "\nOptional parameters:\n";

   for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
		if (!parameters[parameterIndex]->isRequired() && !isReadsFileParam(parameters[parameterIndex])) {
			parameters [ parameterIndex ]->describe();
		}
	}
}

/* 'checkParameters' checks whether all required parameters have been set. */
bool checkParameters()
{
   if (parameters[ findParameter("-h") ]->getBValue() ) {
      printHelp();
      return false;
   }

	for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
		if (parameters[ parameterIndex ]->isRequired() && !parameters[ parameterIndex ]->isSet()) {
			cout << "Required parameter " << parameters[ parameterIndex ]->getShortName() << "/" << parameters[ parameterIndex ]->getLongName()
				  << " " << parameters[ parameterIndex ]->getDescription() << " needs to be set.\n";
			return false;
		}  //if
	}
	// here handle the tricky fact that at least one of -i or -p needs be set; both are not required.
	bool hasBam = parameters[ findParameter("-i") ]->isSet();
	bool hasPin = parameters[ findParameter("-p") ]->isSet();
	if (!hasBam && !hasPin) { cout << "Bam and/or pindel input file required, use -p and/or -i to designate input file(s).\n"; return false; }
			//cout << "chkP4\n";
	return true;
}

/* 'gatherChromosomeNames' gets the names of the chromosomes that occur in the file "referenceFile" and returns it in "chromosomeNames". */
void gatherChromosomeNames(const string& referenceFileName, vector<string>& chromosomeNames )
{
	chromosomeNames.clear(); 

	/*string faiName = referenceFileName + ".fam";
	ifstream faiFile(faiName.c_str());
	if (faiFile.is_open()) { cout << "FAI file exists"; }
	else { cout << "FAI does not exist"; }
	
	exit(EXIT_SUCCESS); // for fast debug*/

	// copying Kai's ReadInOneChr-code here
	char TempChar;
   string TempChrName, tempLine;
	ifstream referenceFile( referenceFileName.c_str());
   referenceFile >> TempChar;
   if (TempChar != '>') {
   	cout << "Please use fasta format for the reference file." << endl;
      exit(EXIT_FAILURE);
   }
   while (referenceFile >> TempChrName) { // for every chromosome do
		getline(referenceFile,tempLine); // throws away rest of comments on reference line.
      cout << "Adding chromosome " << TempChrName << " to the chromosome list.\n";
		chromosomeNames.push_back( TempChrName );
		referenceFile >> TempChar;
		while (!referenceFile.eof() && TempChar!='>' ) {
			referenceFile >> TempChar;
		}
		if (referenceFile.eof()) break;
	}

	// now reset the file pointer to the beginning
	referenceFile.clear();
	referenceFile.seekg(0);


}

/** 'eliminate' eliminates a character from the input string. */
void eliminate (
	const char ch,	// in: character to be eliminated from the string
	string& str		// modif: string that needs to be modified
)
{
	int eliminateCharPos = str.find(ch);
	while ( eliminateCharPos != string::npos ) {
		str.erase( eliminateCharPos, 1 );
		eliminateCharPos = str.find(ch);
	} 
}


/** 'parseRegion' interprets the region specified by the user in the -c option. */
void parseRegion( 
	const string& region, // in: region
	int& startOfRegion,    // out: starting position of the region, -1 if not specified 
	int& endOfRegion, 		// out: ending position of the region, -1 if not specified
	string& chromosomeName, // out: name of the pure chromosome without region information
	bool& correctParse 	 // out: whether parsing has succeeded.
) 
{
	int separatorPos = region.find(":");
	startOfRegion = -1;
	endOfRegion = -1;
	correctParse = false;

	// found a separator
	if (separatorPos != string::npos ) {
		chromosomeName = region.substr(0,separatorPos);
		string coordinates = region.substr( separatorPos+1 );
		eliminate(',',coordinates); // removes the ',' in 1,000 or 1,000,000 that users may add for readability but wreak havoc with atoi
		int startEndSeparatorPos = coordinates.find("-");
		
		// there are two coordinates		
		if ( startEndSeparatorPos != string::npos ) {
			string secondPositionStr = coordinates.substr( startEndSeparatorPos+1 );
			endOfRegion = atoi( secondPositionStr.c_str() );
		}
		startOfRegion = atoi( coordinates.c_str() );
//cout << "sor: " << startOfRegion << "eor: " << endOfRegion << endl;
		if (startOfRegion<0 || ( endOfRegion<startOfRegion && endOfRegion!=-1) ) {
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
bool isFinishedPindel(
	const int upperBinBorder, // in: last position analyzed so far
	const int endOfScan		  // in: the last position to be scanned
)
{
	// if an endOfScan-value has been set
	if ( endOfScan != -1 ) {
		return ( upperBinBorder >= endOfScan );
	}
	else { // using g_maxPos
		return ( upperBinBorder >= g_maxPos );
	}
}

/** 'isFinishedBAM' returns true if there are no more reads to be processed. */
bool isFinishedBAM(
	const int upperBinBorder, // in: last position analyzed so far
	const int endOfScan,		  // in: the last position to be scanned
	const int chromosomeSize  // in: the size of the chromosome
)
{
	// if an endOfScan-value has been set
	if ( endOfScan != -1 ) {
		return ( upperBinBorder >= endOfScan );
	}
	else { // using chromosomeSize
		return ( upperBinBorder >= chromosomeSize );
	}
}


int main (int argc, char *argv[]) {
	 cout << Pindel_Version_str << endl; 
	
    if (NumRead2ReportCutOff == 1) BalanceCutoff = 3000000000;

    //begin charris option parsing replacement
    int c;
    int options_set=0;
    ifstream inf_Seq;
	 ifstream inf_Pindel_Reads;
    string bam_file;
    string OutputFolder;
    string WhichChr;
    string line;
    vector <bam_info> bams_to_parse;
    ifstream config_file;
    bam_info info;

    ifstream  inf_ReadsSeq;   // input file name
    ifstream  inf_BP_test;   // input file name
    ifstream  inf_BP;   // input file name
	bool RefDefined = false;
	bool BAMDefined = false; 
	bool PindelReadDefined = false;
	bool OutputDefined = false;
	bool BreakDancerDefined = false;
	bool ChrDefined = false;
	// define all the parameters you have
	defineParameters(WhichChr);

   	// now read the parameters from the command line
   readParameters( argc, argv);	

	if (argc<=1) { // the user has not given any parameters
	   printHelp();
	   return EXIT_FAILURE;            
   }
	   
    // check parameters
   if (!checkParameters()) exit(EXIT_FAILURE);
	if (FLOAT_WINDOW_SIZE > 5000.0) {
		cout << "Window size " << FLOAT_WINDOW_SIZE << " million bases is too large" << endl;
		return 1;
	}
	else if (FLOAT_WINDOW_SIZE > 500.0) {
		cout << "Window size " << FLOAT_WINDOW_SIZE << " million bases is rather large" << endl;
	}
	WINDOW_SIZE = 1000000 * FLOAT_WINDOW_SIZE;

	// if all parameters are okay, open the files
	inf_Seq.open( par.referenceFileName.c_str() );

	PindelReadDefined = parameters[ findParameter("-p") ]->isSet();
	if (PindelReadDefined) { inf_Pindel_Reads.open(par.pindelFileName.c_str()); }

	BAMDefined = parameters[ findParameter("-i") ]->isSet();
   if (BAMDefined) { 
		config_file.open(par.bamConfigFileName.c_str()); 
   	while(config_file.good()) {
     		config_file >> info.BamFile >> info.InsertSize >> info.Tag;
         //copy kai and throw crap into useless variable
         getline(config_file, line);
         if(config_file.good()) {
            bams_to_parse.push_back(info);
         }
      }
	}

	OutputFolder = par.outputFileName;

	BreakDancerDefined = parameters[ findParameter("-b") ]->isSet();
	if ( BreakDancerDefined ) {
      inf_BP_test.open(par.breakdancerFileName.c_str());
      inf_BP.open(par.breakdancerFileName.c_str());
	}

	omp_set_num_threads(par.numThreads);

	if (MaxRangeIndex>9) { cout << "Maximal range index (-x) exceeds the allowed value (9) - resetting to 9.\n"; MaxRangeIndex=9; }

    // #################################################################



    bool WithFolder = false;
    int StartOfFileName = 0;
    for (int i = bam_file.size(); i >= 0; i--) {
        if (bam_file[i] == '/') {
            StartOfFileName = i;
            WithFolder = true;
            break;
        }
    }

    if (WithFolder) {
        bam_file = bam_file.substr(StartOfFileName + 1, bam_file.size() - 1 - StartOfFileName);
    }
    //cout << bam_file << endl;




    string SIOutputFilename = OutputFolder + "_SI";      // output file name
    //strcpy(SIOutputFilename, (OutputFolder + bam_file + "_SI").c_str());
    ofstream SIoutputfile_test(SIOutputFilename.c_str());
    if (!SIoutputfile_test) {
        cout << "Sorry, cannot write to the file: " << SIOutputFilename << endl;
        return 1;
    }
    SIoutputfile_test.close();

    //char DeletionOutputFilename[10000];      // output file name
    string DeletionOutputFilename = OutputFolder + "_D";
    //strcpy(DeletionOutputFilename, (OutputFolder + bam_file + "_D").c_str());
    ofstream DeletionOutf_test(DeletionOutputFilename.c_str());
    if (!DeletionOutf_test) {
        cout << "Sorry, cannot write to the file: " << DeletionOutputFilename << endl;
        return 1;
    }
    DeletionOutf_test.close();

    string TDOutputFilename = OutputFolder + "_TD";
    //strcpy(DeletionInsertinOutputFilename, (OutputFolder + bam_file + "_DI").c_str());
    ofstream TDOutf_test(TDOutputFilename.c_str());
    if (!TDOutf_test) {
        cout << "Sorry, cannot write to the file: " << TDOutputFilename << endl;
        return 1;
    }
    TDOutf_test.close();

    //char InversionOutputFilename[10000];      // output file name
    string InversionOutputFilename = OutputFolder + "_INV";
    //strcpy(InversionOutputFilename, (OutputFolder + bam_file + "_INV").c_str());
    ofstream InversionOutf_test(InversionOutputFilename.c_str());
    if (!InversionOutf_test) {
        cout << "Sorry, cannot write to the file: " << InversionOutputFilename << endl;
        return 1;
    }
    InversionOutf_test.close();

    //char LargeInsertionOutputFilename[10000];      // output file name
    string LargeInsertionOutputFilename = OutputFolder + "_LI";
    //strcpy(LargeInsertionOutputFilename, (OutputFolder + bam_file + "_LI").c_str());
    ofstream LargeInsertionOutf_test(LargeInsertionOutputFilename.c_str());
    if (!LargeInsertionOutf_test) {
        cout << "Sorry, cannot write to the file: " << LargeInsertionOutputFilename << endl;
        return 1;
    }
    LargeInsertionOutf_test.close();

    //char RestOutputFilename[10000];      // output file name
    string RestOutputFilename = OutputFolder + "_BP";
    //strcpy(RestOutputFilename, (OutputFolder + bam_file + "_BP").c_str());
    ofstream RestOutf_test(RestOutputFilename.c_str());
    if (!RestOutf_test) {
        cout << "Sorry, cannot write to the file: " << RestOutputFilename << endl;
        return 1;
    }
    RestOutf_test.close();


    //WhetherRemoveDuplicates = atoi(argv[6]);
    int Count_SI = 0;
    int Count_D = 0;
    int Count_DI = 0;
    int Count_TD = 0;
    int Count_TD_NT = 0;
    int Count_Inv = 0;
    int Count_Inv_NT = 0;
    int Count_D_Plus = 0;
    int Count_D_Minus = 0;
    int Count_DI_Plus = 0;
    int Count_DI_Minus = 0;
    int Count_TD_Plus = 0;
    int Count_TD_Minus = 0;
    int Count_TD_NT_Plus = 0;
    int Count_TD_NT_Minus = 0;
    int Count_Inv_Plus = 0;
    int Count_Inv_Minus = 0;
    int Count_Inv_NT_Plus = 0;
    int Count_Inv_NT_Minus = 0;
    int Count_SI_Plus = 0;
    int Count_SI_Minus = 0;
    //int Entering_D_Plus = 0;
    //int Entering_D_Minus = 0;
    //int Plus_Sum_Left = 0;
    //int Plus_Sum_Right = 0;
    //int Minus_Sum_Left = 0;
    //int Minus_Sum_Right = 0;

    Match[(short)'A'] = 'A';
    Match[(short)'C'] = 'C';
    Match[(short)'G'] = 'G';
    Match[(short)'T'] = 'T';
    Match[(short)'N'] = 'X';
    Match[(short)'$'] = '$';
    Match2N[(short)'A'] = 'N';
    Match2N[(short)'C'] = 'N';
    Match2N[(short)'G'] = 'N';
    Match2N[(short)'T'] = 'N';
    Match2N[(short)'N'] = 'X';
    Match2N[(short)'$'] = '$';
    Convert2RC[(short)'A'] = 'T';
    Convert2RC[(short)'C'] = 'G';
    Convert2RC[(short)'G'] = 'C';
    Convert2RC[(short)'T'] = 'A';
    Convert2RC[(short)'N'] = 'X';
    Convert2RC[(short)'$'] = '$';
    Convert2RC4N[(short)'A'] = 'T';
    Convert2RC4N[(short)'C'] = 'G';
    Convert2RC4N[(short)'G'] = 'C';
    Convert2RC4N[(short)'T'] = 'A';
    Convert2RC4N[(short)'N'] = 'N';
    Cap2LowArray[(short)'A'] = 'a';
    Cap2LowArray[(short)'C'] = 'c';
    Cap2LowArray[(short)'G'] = 'g';
    Cap2LowArray[(short)'T'] = 't';
    Cap2LowArray[(short)'N'] = 'n';
    Cap2LowArray[(short)'$'] = 'n';

    time_t Time_Load_S, Time_Load_E, Time_Mine_E, Time_Sort_E;//, Time_End;
    Time_Load_S = time(NULL);
    unsigned int AllLoadings = 0;
    unsigned int AllMinings = 0;
    unsigned int AllSortReport = 0;


    string Spacer = "";
    for (int i = 0; i < SpacerBeforeAfter; i++)
        Spacer += "N";
    //cout << Distance << endl;
    Distance = 300;
    //MinClose = short(log((double)Distance)/log(4.0) + 0.8) + 3 + MAX_SNP_ERROR;//atoi(argv[1]);
    //MinFar_I = MinClose + 1;//atoi(argv[2]);
    //cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
    //MinFar_D = 6;//atoi(argv[3]);

    DSizeArray[0]  = 0;
    //	DSizeArray[1]  =         25;
    DSizeArray[1]  = 128;
    DSizeArray[2]  = DSizeArray[1] * 4; // 512
    DSizeArray[3]  = DSizeArray[2] * 4; // 2048
    DSizeArray[4]  = DSizeArray[3] * 4; // 8092
    DSizeArray[5]  = DSizeArray[4] * 4; // 32368
    DSizeArray[6]  = DSizeArray[5] * 4; // 129472
    DSizeArray[7]  = DSizeArray[6] * 4; // 517888
    DSizeArray[8]  = DSizeArray[7] * 4;

    DSizeArray[9]  = DSizeArray[8] * 4;
    DSizeArray[10] = DSizeArray[9] * 4;
    DSizeArray[11] = DSizeArray[10] * 4;
    DSizeArray[12] = DSizeArray[11] * 4;
    DSizeArray[13] = DSizeArray[12] * 4;
    DSizeArray[14] = DSizeArray[13] * 4;

    //unsigned int DSizeExtra = 100;   
    //unsigned int D_SIZE = 100;

    string TempLie_BD;

    string CurrentChr;
    char TempStartChrChar;
    //short FragID;

    char FirstSharpChar;

    unsigned int EndOfFragment;// = CurrentChr.size() - Spacer;
    unsigned int StartOfFragment;// = CurrentChr.size() - Spacer;

    //short BP_Left_Start;
    //short BP_Left_End;
    //short BP_Right_Start;
    //short BP_Right_End;
    //unsigned int LeftStart, LeftEnd, RightStart, RightEnd;
    //bool ShortInsertionOK;
    //bool DeletionOK;
    //unsigned int DISTANCE;
    //short ReadLengthMinus;
    //short ReadLength;

    //char LeftChar, RightChar;
    unsigned int Num_Left;
    //short CurrentChrIndex = 0;



    string TempLine_BD;

    //while (inf_Seq >> CurrentFragName) 

   
    //cout << StartOfFragment << "\t" << EndOfFragment << endl;
    vector <SPLIT_READ> InputReads, Reads, BP_Reads, FutureReads;

    vector <BreakDancer> All_BD_events_WG, All_BD_events;
    BreakDancer Temp_BD_event;
    All_BD_events_WG.push_back(Temp_BD_event);

    //short * BD_INDEX;// = new short[CurrentChr.size()];


    //if (WhichChr == CurrentFragName) 
    //{
    while (inf_BP_test >> FirstSharpChar) {
        if (FirstSharpChar == '#') {
            getline(inf_BP_test, TempLine_BD);
            getline(inf_BP, TempLine_BD);
        }
        else {
            getline(inf_BP_test, TempLine_BD);
            inf_BP >> Temp_BD_event.ChrName_A >> Temp_BD_event.S1 >> TempLine_BD 
                >> Temp_BD_event.ChrName_B >> Temp_BD_event.S3 >> TempLine_BD;
            //>> Temp_BD_event.Type >> Temp_BD_event.Size;
            getline(inf_BP, TempLine_BD);

            //<< Temp_BD_event.Size << endl;
            //if (Temp_BD_event.Type == "DEL" ) 
            //if (WhichChr == Temp_BD_event.ChrName_A)
            {
                //if (Temp_BD_event.S1 + 200 > CONS_Chr_Size) Temp_BD_event.S2 = CONS_Chr_Size - 1;
                //else 
					 Temp_BD_event.S2 = Temp_BD_event.S1 + 200;
                if (Temp_BD_event.S1 > 200) Temp_BD_event.S1 = Temp_BD_event.S1 - 200;
                else Temp_BD_event.S1 = 1;
            }
            //if (WhichChr == Temp_BD_event.ChrName_B)
            {
                //if (Temp_BD_event.S3 + 200 > CONS_Chr_Size) Temp_BD_event.S4 = CONS_Chr_Size - 1;
                //else 
					 Temp_BD_event.S4 = Temp_BD_event.S3 + 200;
                if (Temp_BD_event.S3 > 200) Temp_BD_event.S3 = Temp_BD_event.S3 - 200;
                else Temp_BD_event.S3 = 1;
            }					 
            Temp_BD_event.S1 += SpacerBeforeAfter;
            Temp_BD_event.S2 += SpacerBeforeAfter;
            Temp_BD_event.S3 += SpacerBeforeAfter;
            Temp_BD_event.S4 += SpacerBeforeAfter;
            //if (WhichChr == Temp_BD_event.ChrName_A && WhichChr == Temp_BD_event.ChrName_B) {
                //cout << WhichChr << " " << Temp_BD_event.S1 << " " << Temp_BD_event.S2 << " " << Temp_BD_event.S3 << " " << Temp_BD_event.S4 << endl;  
                All_BD_events_WG.push_back(Temp_BD_event);
            //}
        }
    }
    cout << "BreakDancer events: " << All_BD_events_WG.size() - 1 << endl;

    //}


	vector<string> chromosomes; 
	int chromosomeIndex = 0;

	
	string CurrentChrName, emptystr;
	char FirstCharOfFasta;
	inf_Seq >> FirstCharOfFasta;
	if (FirstCharOfFasta != '>') {
		cout << "The reference genome must be in fasta format!" << endl;
		return 1;
	}
	
	bool SpecifiedChrVisited = false;

	int startOfRegion = -1;
	int endOfRegion = -1;
	bool correctParse = false;
	string chrName;
	parseRegion( WhichChr, startOfRegion, endOfRegion, chrName, correctParse );
	if (!correctParse) {
		cout << "I cannot parse the region '" << WhichChr << "'. Please give region in the format -c ALL, -c <chromosome_name> "
				  "(for example -c 20) or -c <chromosome_name>:<start_position>[-<end_position>], for example -c II:1,000 or "
				  "-c II:1,000-50,000. If an end position is specified, it must be larger than the start position.\n";
		exit( EXIT_FAILURE );
	}
	WhichChr = chrName; // removes the region from the 'pure' chromosome name

	int startOffSet = 0;
	// if a region has been specified
	if ( startOfRegion>=0 ) {
		startOffSet = startOfRegion - AROUND_REGION_BUFFER;
		if (startOffSet<0) {
			startOffSet = 0;
		}
	}
	int endRegionPlusBuffer = -1; // -1 indicates that the chromosome must be read to the end
	if ( endOfRegion > -1 ) {
		endRegionPlusBuffer = endOfRegion + AROUND_REGION_BUFFER;
	}
	bool loopOverAllChromosomes = false;
	if (WhichChr.compare("ALL") == 0 ) {
cout << "Looping over ALL chromosomes.\n";
		loopOverAllChromosomes = true;
		//gatherChromosomeNames(par.referenceFileName, chromosomes);
		//WhichChr = chromosomes[ 0 ];
	}

//cout << "Chr: " << chrName << "region: " << startOfRegion << " till " << endOfRegion << endl;

while (SpecifiedChrVisited == false && inf_Seq >> CurrentChrName && !inf_Seq.eof()) { // loop over chromosomes
	cout << "Processing chromosome: " << CurrentChrName << endl;
//cout << "WhichChr is: '" << WhichChr << "', CurrentChrName is '" << CurrentChrName << "'" << endl;
	getline(inf_Seq, emptystr);
	if (loopOverAllChromosomes) {
		GetOneChrSeq(inf_Seq, CurrentChr, true);
		WhichChr = CurrentChrName;
	}
	else if (CurrentChrName == WhichChr) { // just one chr and this is the correct one
		GetOneChrSeq(inf_Seq, CurrentChr, true);
		SpecifiedChrVisited = true;
	}
	else { // not build up sequence
		GetOneChrSeq(inf_Seq, CurrentChr, false);
		cout << "Skipping chromosome: " << CurrentChrName << endl;
		continue;
	}
	
   CONS_Chr_Size = CurrentChr.size() - 2 * SpacerBeforeAfter;
	cout << "Chromosome Size: " << CONS_Chr_Size << endl;
	CurrentChrMask.resize(CurrentChr.size());
	for (int i = 0; i < CurrentChr.size(); i++) {
		CurrentChrMask[i] = 'N';
	}
   unsigned NumBoxes =(unsigned) (CurrentChr.size() / BoxSize) + 1; // box size 
   cout << NumBoxes << "\t" << BoxSize << endl;

    // cout << "NumBoxes: " << NumBoxes << endl;
   vector <unsigned> SIs[NumBoxes];
    //vector <SPLIT_READ> LIs[NumBoxes];
   vector <unsigned> Deletions[NumBoxes];
   vector <unsigned> TD[NumBoxes];
   vector <unsigned> TD_NT[NumBoxes];
   vector <unsigned> DI[NumBoxes];
   vector <unsigned> Inv[NumBoxes];
   vector <unsigned> Inv_NT[NumBoxes];

   EndOfFragment = CurrentChr.size() - SpacerBeforeAfter;
   StartOfFragment = SpacerBeforeAfter;


/* Starting the loop to read the subfiles one by one (EWL070111) -> */
g_binIndex = -1; // to start with 0...
char ch;
   int lowerBinBorder = startOffSet - WINDOW_SIZE ; 
   int upperBinBorder = lowerBinBorder + WINDOW_SIZE;
	int displayedStartOfRegion = ((startOfRegion>=0) ? ( startOfRegion - WINDOW_SIZE ) : lowerBinBorder );
	int displayedEndOfRegion = displayedStartOfRegion + WINDOW_SIZE;
do { 
	g_binIndex++;
	lowerBinBorder += WINDOW_SIZE;
	upperBinBorder += WINDOW_SIZE;
	displayedStartOfRegion += WINDOW_SIZE;
	displayedEndOfRegion += WINDOW_SIZE;
	if (displayedEndOfRegion > endOfRegion) {
		displayedEndOfRegion = endOfRegion; 
	}

	// if the region end is specified, and it is before the regular upper border of the bin
	if (endRegionPlusBuffer>-1 && upperBinBorder > endRegionPlusBuffer ) {
		upperBinBorder = endRegionPlusBuffer;
	}
	
	if (displayedStartOfRegion < displayedEndOfRegion ) {
	   cout << "Looking at chromosome " << WhichChr << " bases " << displayedStartOfRegion << " to " << displayedEndOfRegion << ".\n";
	}
	else {
		cout << "Checking out reads near the borders of the specified regions for extra evidence.\n";
	}

	//cout << "BinBorder " << lowerBinBorder << " " << upperBinBorder << endl; 


//cout << "Step 2" << endl;
//cin >> ch;
   if (Time_Load_S==0) Time_Load_S=time(NULL); // EWL070111
    short ReturnFromReadingReads;
    int GetPlus = 0;
    int GetMinus = 0;
    //if (WhichChr == CurrentFragName) 
	// read Pindel format reads

	 // read bam format reads
	//cout << "here" << endl;
	 VectorTag.clear();
	 if (BAMDefined)
    {
        ReturnFromReadingReads = 0;
        for(int i=0; i < bams_to_parse.size(); i++) {
			  //cout << i << "\t" << bams_to_parse[i].BamFile.c_str() << "\t" << Reads.size()  << endl;
            ReturnFromReadingReads = ReadInBamReads(bams_to_parse[i].BamFile.c_str(), WhichChr, &CurrentChr, Reads, bams_to_parse[i].InsertSize, bams_to_parse[i].Tag, lowerBinBorder, upperBinBorder);
            if (ReturnFromReadingReads == 0) {
                cout << "Bam read failed: " << bams_to_parse[i].BamFile <<  endl;
                return 1;
            }
            else if (Reads.size() == 0) {
               cout << "No Reads for " << WhichChr << " found in " << bams_to_parse[i].BamFile << endl;
            } 
			  cout << "BAM file index\t" << i << "\t" << Reads.size() << endl;
        }

    }
	
	if (PindelReadDefined)
	{
		ReturnFromReadingReads = 0;
		//ReturnFromReadingReads = ReadInRead(inf_Pindel_Reads, WhichChr, CurrentChr, Reads);
		ReturnFromReadingReads = ReadInRead(inf_Pindel_Reads, WhichChr, CurrentChr, Reads, lowerBinBorder, upperBinBorder);
		if (ReturnFromReadingReads == 1) {
			cout << "malformed record detected!" << endl;
			return 1;
		}
		else if (Reads.size() == 0) { cout << "No reads found!?\n";}
	}
	//cout << "There are " << Reads.size() << endl;
//cout << "Step 3" << endl;
//cin >> ch;
    Time_Mine_E = time(NULL);
	
	//for (int i = 0; i < Reads.size(); i++) {
	//	if (Reads[i].Tag == "normal") cout << "normal\t" << i << endl;
	//}

    if (Reads.size())
        cout << "There are " << Reads.size() << " reads for this chromosome." <<endl;
    else {
        cout << "There are no reads for this bin." <<endl;  
        continue;
    }
    Num_Left = Reads.size();
    Const_Log_T = log10((double)Num_Left);
    Time_Load_E = time(NULL);
    int CountFarEnd, CountFarEndPlus, CountFarEndMinus;


	// ################### module 3: search breakpoints #####################
	All_BD_events.clear();
	for (int All_BD_events_WG_index = 0; All_BD_events_WG_index < All_BD_events_WG.size(); All_BD_events_WG_index++) {
		if (All_BD_events_WG[All_BD_events_WG_index].ChrName_A == CurrentChrName && All_BD_events_WG[All_BD_events_WG_index].ChrName_B == CurrentChrName)
			All_BD_events.push_back(All_BD_events_WG[All_BD_events_WG_index]);
			//CurrentChrName ChrName_A 
	} 
	//All_BD_events_WG
	if (All_BD_events.size() > 1) {
		cout << "Searching additional breakpoints by adding BreakDancer results" << endl;
		int * BD_INDEX = new int[CurrentChr.size()];
		//cout << "1" << endl;
		for (unsigned i = 0; i < CurrentChr.size(); i++) BD_INDEX[i] = 0; 
		for (unsigned i = 1; i < All_BD_events.size(); i++) {
			//cout << i << endl;
			for (unsigned j = All_BD_events[i].S1; j < All_BD_events[i].S2; j++) BD_INDEX[j] = i;
			for (unsigned j = All_BD_events[i].S3; j < All_BD_events[i].S4; j++) BD_INDEX[j] = i * (-1);					 
		}
		int BD_Plus = 0;
		int BD_Minus = 0;
		for (unsigned i = 0; i < CurrentChr.size(); i++) {
			if (BD_INDEX[i] > 0) BD_Plus++;
			else if (BD_INDEX[i] < 0) BD_Minus++; 
		}
		cout << BD_Plus << "\t" << BD_Minus << endl;				 
		
		//ADDITIONAL_MISMATCH = 2;
		//Seq_Error_Rate = 0.05;
		//int BD_event_Index = 0;
		CountFarEnd = 0;
		CountFarEndMinus = 0;
		CountFarEndPlus = 0;
		int Start_pos, End_pos;
		for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			if (!Reads[ReadIndex].UP_Far.empty()) {
				//CountFarEnd++;
				continue;
			}
			int BD_event_Index = BD_INDEX[Reads[ReadIndex].UP_Close[0].AbsLoc];
			if (BD_event_Index == 0) continue;
			//All_BD_events[i].S4; j++) BD_INDEX[j]
			//MAX_SNP_ERROR = (short)(Reads[ReadIndex].UnmatchedSeq.size() * Seq_Error_Rate + 1.0);
			
			//TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
			//TOTAL_SNP_ERROR_CHECKED_Minus = MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
			//MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize)/log(4.0) + 0.8)) + 3;//atoi(argv[1]);
			//MinFar_I = MinClose + 1;//atoi(argv[2]);
			//cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
			//MinFar_D = 8;//atoi(argv[3]);	
			if (BD_event_Index > 0) {
				//cout << ReadIndex << "\t" << BD_event_Index << "\t" << All_BD_events[BD_event_Index].S1 << "\t" << All_BD_events[BD_event_Index].S2 << "\t" << All_BD_events[BD_event_Index].S3 << "\t" << All_BD_events[BD_event_Index].S4 << endl;
				Start_pos = All_BD_events[BD_event_Index].S3 - Reads[ReadIndex].ReadLength;
				End_pos = All_BD_events[BD_event_Index].S4 + Reads[ReadIndex].ReadLength;
				//cout << Start_pos << "\t" << End_pos << endl;
				GetFarEnd(CurrentChr, Reads[ReadIndex], Start_pos, End_pos);  
			}
			else { // < 0
				//cout << ReadIndex << "\t" << BD_event_Index << "\t" << All_BD_events[BD_event_Index * (-1)].S1 << "\t" << All_BD_events[BD_event_Index * (-1)].S2 << "\t" << All_BD_events[BD_event_Index * (-1)].S3 << "\t" << All_BD_events[BD_event_Index * (-1)].S4 << endl;
				Start_pos = All_BD_events[BD_event_Index * (-1)].S1 - Reads[ReadIndex].ReadLength;
				End_pos = All_BD_events[BD_event_Index * (-1)].S2 + Reads[ReadIndex].ReadLength;
				//cout << Start_pos << "\t" << End_pos << endl;
				GetFarEnd(CurrentChr, Reads[ReadIndex], Start_pos, End_pos);  
			}
			
			if (!Reads[ReadIndex].UP_Far.empty()) { // if a far end has been found


				// if there is a non-template sequence present between close and far end
				if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength < Reads[ReadIndex].ReadLength) {

					// if there are backup reads
					if (Reads[ReadIndex].UP_Far_backup.size()) {

						// if the backup reads are worse
						if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
							
							// put current reads in backup
							Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
							Reads[ReadIndex].UP_Far.clear();
						}
						else Reads[ReadIndex].UP_Far.clear(); // otherwise keep your backup
					}
					else { // there are no backup reads, prepare for next cycle
						Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
						Reads[ReadIndex].UP_Far.clear();
					}
				}
				else { // no non-template bases present: good match found
					CountFarEnd++; 
					if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
					else CountFarEndMinus++; 
					// note:UPFar remains filled, so read will be skipped in next cycle
				}
			}
		}
		cout << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t" 
		//<< "\tTotal number of reads:" << Reads.size() << "\n" 
		//<< CountFarEnd * 100.0 / Reads.size() << " %\n"
		<< "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << "\n\n";//endl;
		delete[] BD_INDEX;
	}
	
	cout << "Searching breakpoints of deletion events" << endl;
	//int SkippReads = 0;
	for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) { 
		//cout << "RangeIndex: " << RangeIndex << endl;
		//SkippReads = 0;
		CountFarEnd = 0;
		CountFarEndMinus = 0;
		CountFarEndPlus = 0;
		
#pragma omp parallel default(shared)
		{
#pragma omp for
			for (int ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
				//cout << ReadIndex << endl;
				if (!Reads[ReadIndex].UP_Far.empty()) {
					//SkippReads++;
					//CountFarEnd++;
					continue;
				}
				//if (Reads[ReadIndex].MatchedRelPos > 160000) 
				//if (RangeIndex == 4)
				//cout << RangeIndex << "\t" << ReadIndex << "\t" << Reads[ReadIndex].UP_Far.size() << "\t";
				GetFarEnd_SingleStrandDownStream(CurrentChr, Reads[ReadIndex], RangeIndex);
				//if (Reads[ReadIndex].MatchedRelPos > 160000) 
				//if (RangeIndex == 4)
				//cout << "after: " << Reads[ReadIndex].UP_Far.size() << "\t";
				if (!Reads[ReadIndex].UP_Far.empty()) {
					if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength < Reads[ReadIndex].ReadLength) {
						if (Reads[ReadIndex].UP_Far_backup.size()) {
							if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
								Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
								Reads[ReadIndex].UP_Far.clear();
							}
							else Reads[ReadIndex].UP_Far.clear();
						}
						else {
							Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
							Reads[ReadIndex].UP_Far.clear();
						}
					}
					else {
#pragma omp critical
						{
							CountFarEnd++; 
							if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
							else CountFarEndMinus++; 							 
						}
					}
				}
				//if (Reads[ReadIndex].MatchedRelPos > 160000) 
				//if (RangeIndex == 4)
				//cout << "final: " << Reads[ReadIndex].UP_Far.size() << endl;
			}
			//cout << "SkippReads " << SkippReads << endl;
		} // #pragma omp parallel default(shared)
		
		cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t" 
		//<< "\tTotal number of reads:" << Reads.size() << "\n" 
		//<< CountFarEnd * 100.0 / Reads.size() << " %\n"
		<< "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
	}
	
	cout << "Searching breakpoints of SI events" << endl;
	for (short RangeIndex = 1; RangeIndex < 2; RangeIndex++) {
		CountFarEnd = 0;
		CountFarEndMinus = 0;
		CountFarEndPlus = 0;
		
#pragma omp parallel default(shared)
		{
#pragma omp for
			for (int ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
				if (!Reads[ReadIndex].UP_Far.empty()) {
					//CountFarEnd++;
					continue;
				}
				
				GetFarEnd_SingleStrandDownStreamInsertions(CurrentChr, Reads[ReadIndex], RangeIndex);
#pragma omp critical 
				{
					if (!Reads[ReadIndex].UP_Far.empty()) {
						CountFarEnd++; 
						if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
						else CountFarEndMinus++;
					}					 
				}
			}
		}
		
		cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t" 
		//<< "\tTotal number of reads:" << Reads.size() << "\n" 
		//<< CountFarEnd * 100.0 / Reads.size() << " %\n"
		<< "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
	}
	
	if (Analyze_TD) {
		
		
		cout << "Searching breakpoints of tandem duplication events" << endl;
		//int CountFarEnd, CountFarEndPlus, CountFarEndMinus;
		for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {
			
			CountFarEnd = 0;
			CountFarEndMinus = 0;
			CountFarEndPlus = 0;
#pragma omp parallel default(shared)
			{
#pragma omp for
				for (int ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
					if (!Reads[ReadIndex].UP_Far.empty()) {
						//CountFarEnd++;
						continue;
					}
					
					//MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;//atoi(argv[1]);
					//MinFar_I = MinClose + 1;//atoi(argv[2]);
					//cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
					//MinFar_D = 8;//atoi(argv[3]);				 
					//if (RangeIndex <= 5)
					GetFarEnd_SingleStrandUpStream(CurrentChr, Reads[ReadIndex], RangeIndex);
					if (!Reads[ReadIndex].UP_Far.empty()) {
						if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength >= Reads[ReadIndex].ReadLength) {
							if (Reads[ReadIndex].MatchedD == Plus) {
								if (Reads[ReadIndex].UP_Close[0].AbsLoc < Reads[ReadIndex].ReadLength + Reads[ReadIndex].UP_Far[0].AbsLoc)
									Reads[ReadIndex].UP_Far.clear();
							}
							else { // if (Reads[ReadIndex].MatchedD == Minus)
								if (Reads[ReadIndex].UP_Far[0].AbsLoc < Reads[ReadIndex].ReadLength + Reads[ReadIndex].UP_Close[0].AbsLoc)
									Reads[ReadIndex].UP_Far.clear();
							}
						}
						else {
							if (Reads[ReadIndex].UP_Far_backup.size()) {
								if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
									Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
									Reads[ReadIndex].UP_Far.clear();
								}
								else Reads[ReadIndex].UP_Far.clear();
							}
							else {
								Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
								Reads[ReadIndex].UP_Far.clear();
							}
						} 
#pragma omp critical 
						{
							if (!Reads[ReadIndex].UP_Far.empty()) {
								CountFarEnd++; 
								if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
								else CountFarEndMinus++; 
							} 
						}
					}
				}
			}			 
			cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t" 
			//<< "\tTotal number of reads:" << Reads.size() << "\n" 
			//<< CountFarEnd * 100.0 / Reads.size() << " %\n"
			<< "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
		}
	} // if (Analyze_TD)
	
	if (Analyze_INV) {
		//cout << "here" << endl;
		cout << "Searching breakpoints of inversions" << endl;
		unsigned ReadsUsedForD = 0;
		unsigned ReadsUsedForDI = 0;
		for (short RangeIndex = 1; RangeIndex < MaxRangeIndex; RangeIndex++) {
			
			CountFarEnd = 0;
			CountFarEndMinus = 0;
			CountFarEndPlus = 0;
#pragma omp parallel default(shared)
			{
#pragma omp for
				for (int ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
					if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.size()) {
						continue;
					}
					
					//MinClose = short(log((double)(2 * Reads[ReadIndex].InsertSize + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;//atoi(argv[1]);
					//MinFar_I = MinClose + 1;//atoi(argv[2]);
					//cout << "For short insertion: " << MinClose << "\t" << MinFar_I << endl;
					//MinFar_D = 8;//atoi(argv[3]);				 
					//if (RangeIndex <= 5)
				   GetFarEnd_OtherStrand(CurrentChr, Reads[ReadIndex], RangeIndex);
					
					if (!Reads[ReadIndex].UP_Far.empty()) {
						if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr + Reads[ReadIndex].CloseEndLength < Reads[ReadIndex].ReadLength) {
							if (Reads[ReadIndex].UP_Far_backup.size()) {
								if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr <Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
									Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
									Reads[ReadIndex].UP_Far.clear();
								}
								else Reads[ReadIndex].UP_Far.clear();
							}
							else {
								Reads[ReadIndex].UP_Far_backup = Reads[ReadIndex].UP_Far;
								Reads[ReadIndex].UP_Far.clear();
							}
						}
						else {
#pragma omp critical 
							{
								CountFarEnd++; 
								if (Reads[ReadIndex].MatchedD == Plus) CountFarEndPlus++;
								else CountFarEndMinus++;							 
							}
						}
					}
				}
			}
			cout << RangeIndex << "\tNumber of reads with far end mapped: " << CountFarEnd << "\t" 
			//<< "\tTotal number of reads:" << Reads.size() << "\n" 
			//<< CountFarEnd * 100.0 / Reads.size() << " %\n"
			<< "Far+: " << CountFarEndPlus << "\tFar-: " << CountFarEndMinus << endl;
		}
	} // if (Analyze_INV)
	
	// compare backup with current value
	cout << "revisit all breakpoints identified ...";
	for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
		if (Reads[ReadIndex].UP_Far.empty()) {
			if (!Reads[ReadIndex].UP_Far_backup.empty()) {
				Reads[ReadIndex].UP_Far = Reads[ReadIndex].UP_Far_backup;
			}
		}
		else if (!Reads[ReadIndex].UP_Far_backup.empty()) {
			if (Reads[ReadIndex].UP_Far_backup[Reads[ReadIndex].UP_Far_backup.size() - 1].LengthStr >
				 Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].LengthStr) {
				Reads[ReadIndex].UP_Far = Reads[ReadIndex].UP_Far_backup;
			}
		} 
	}
	cout << " done." << endl;	 
	
	// ################### module 4: search variants and report #####################    
	
	//short MAX_MISMATCHES_Per_Read = 0;;
	cout << "Searching deletion events ... " << endl;
	for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
		if (Reads[ReadIndex].UP_Far.empty()) continue;
		//MAX_MISMATCHES_Per_Read = (short)(Seq_Error_Rate * Reads[ReadIndex].ReadLength + 1);
		//if (Reads[ReadIndex].UP_Far.size()) 
		{
			if (Reads[ReadIndex].MatchedD == Plus) { // MAX_SNP_ERROR
				for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
					for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
						if (Reads[ReadIndex].Used) break;
						if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
						for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
							//cout << "+" << FarIndex << endl;
							if (Reads[ReadIndex].Used) break;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
								continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								//cout << "1" << endl;
								//cout << "+\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr 
								//	<< "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr
								//     << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].ReadLength << "\t";
								//cout << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\t>\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 << endl;
								if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) {
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize =  (Reads[ReadIndex].Right - Reads[ReadIndex].Left) - Reads[ReadIndex].ReadLengthMinus;
									Reads[ReadIndex].NT_str = "";
									Reads[ReadIndex].NT_size = 0;
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{  
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{ 
											/*
											int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
											//cout << Reads[ReadIndex].BPLeft << " " << lowerBinBorder << endl;
											if ((int)Reads[ReadIndex].BPLeft < 0)
												Deletions[0].push_back(ReadIndex);
											else */ 
											if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												Deletions[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												Count_D++;
												Count_D_Plus++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			else if (Reads[ReadIndex].MatchedD == Minus) {
				for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
					for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
						if (Reads[ReadIndex].Used) break;
						if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
						for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) {
							if (Reads[ReadIndex].Used) break;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
								continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
									 && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) {
									
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize = (Reads[ReadIndex].Right - Reads[ReadIndex].Left) - Reads[ReadIndex].ReadLengthMinus;
									Reads[ReadIndex].NT_str = "";
									Reads[ReadIndex].NT_size = 0;
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{
										
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{ 
											//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
											//if ((int)Reads[ReadIndex].BPLeft < 0)
												//Deletions[0].push_back(ReadIndex);
											//else 
											if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												Deletions[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												Count_D++;
												Count_D_Minus++;
											}
											//cout << "- " << Count_D << endl; 
										}
									}
								}
							}
						}
					}
				} // for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++)
			}
		}
	}
	cout << "Total: " << Count_D << "\t+" << Count_D_Plus << "\t-" << Count_D_Minus << endl;
	ofstream DeletionOutf(DeletionOutputFilename.c_str(),ios::app);
	SortOutputD(NumBoxes, CurrentChr, Reads, Deletions, DeletionOutf);
	//DeletionOutf.close();  
	for (unsigned int i = 0; i < NumBoxes; i++) Deletions[i].clear(); 
	
	
	cout << "Searching deletion-insertions ... " << endl;
	unsigned CloseIndex, FarIndex;
	for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
		if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
		//if (Reads[ReadIndex].UP_Far.size()) 
		{
			CloseIndex = Reads[ReadIndex].UP_Close.size() - 1;
			FarIndex = Reads[ReadIndex].UP_Far.size() - 1;
			if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > (short)(1 + Seq_Error_Rate * (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr)))
				continue;
			if (Reads[ReadIndex].MatchedD == Plus) {
				//for (unsigned CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) 
				{
					//for (unsigned FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) 
					{
						
						if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
							if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength
								 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) { 
								//DeletionUnique = true;
								Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
								Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
								Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
								Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
								
								Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
								Reads[ReadIndex].InsertedStr = "";                                 
								
								Reads[ReadIndex].IndelSize =  (Reads[ReadIndex].Right - Reads[ReadIndex].Left) + Reads[ReadIndex].NT_size - Reads[ReadIndex].ReadLengthMinus;
								
								Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
								Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
								//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
								//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
								//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
								
								{  
									if (Reads[ReadIndex].IndelSize >= MIN_IndelSize_NT && Reads[ReadIndex].NT_size <= Max_Length_NT /*Reads[ReadIndex].NT_size*/) 
									{ 
										
									if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{ 
											//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
											//if ((int)Reads[ReadIndex].BPLeft < 0)
												//DI[0].push_back(ReadIndex);
											//else 
											if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												DI[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											//DI[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												Count_DI++;
												Count_DI_Plus++;
											}
										}
									}
								}
								
							}
							
							// ############
							//}
						}
					}
				}
			}
			else if (Reads[ReadIndex].MatchedD == Minus) {
				//for (unsigned CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) 
				{
					//for (unsigned FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) 
					{
						if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
							if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength
								 && Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr >= Min_Num_Matched_Bases && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) { 
								//DeletionUnique = true;
								Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
								Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
								Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
								Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].UP_Far[FarIndex].LengthStr;
								Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
								
								Reads[ReadIndex].IndelSize =  (Reads[ReadIndex].Right - Reads[ReadIndex].Left) - Reads[ReadIndex].ReadLengthMinus + Reads[ReadIndex].NT_size;
								//                                 LeftReads[Left_Index].NT_str = "";
								Reads[ReadIndex].InsertedStr = "";
								//cout << LeftReads[Left_Index].IndelSize << endl;
								Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
								Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
								//LeftReads[Left_Index].score = Const_I + Const_S + LOG14 * LeftReads[Left_Index].ReadLength + Const_Log_T;
								//LeftReads[Left_Index].Unique = true;
								//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
								//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
								//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
								{  
									if (Reads[ReadIndex].IndelSize >= MIN_IndelSize_NT && Reads[ReadIndex].NT_size <= Max_Length_NT)
									{ 
										
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{ 
											//DI[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
											//if ((int)Reads[ReadIndex].BPLeft < 0)
											//	DI[0].push_back(ReadIndex);
											//else *
											if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												DI[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												Count_DI++;
												Count_DI_Minus++;
											}
										}
									}
								}
								
							}
							
							// ######################
						}
					}
				}
			}
		}
	}
	cout << "Total: " << Count_DI << "\t+" << Count_DI_Plus << "\t-" << Count_DI_Minus << endl;
	//ofstream DeletionInsertionOutf(DeletionInsertinOutputFilename.c_str());
	SortOutputDI(NumBoxes, CurrentChr, Reads, DI, DeletionOutf);
	//DeletionInsertionOutf.close();
	DeletionOutf.close();
	for (unsigned int i = 0; i < NumBoxes; i++) DI[i].clear(); 
	
	if (Analyze_TD) {
		
		
		cout << "Searching tandem duplication events ... " << endl;
		for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			//short MAX_MISMATCHES_Per_Read = (short)(Seq_Error_Rate * Reads[ReadIndex].ReadLength + 1);
			//if (Reads[ReadIndex].UP_Far.size()) 
			{
				if (Reads[ReadIndex].MatchedD == Plus) {
					for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
						for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
							if (Reads[ReadIndex].Used) break;
							if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
								if (Reads[ReadIndex].Used) break;
								if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
									continue;
								if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
									if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc) {
										Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
										Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
										Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
										
										Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1;
										Reads[ReadIndex].NT_str = "";
										Reads[ReadIndex].NT_size = 0;
										Reads[ReadIndex].InsertedStr = "";
										Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
										//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
										//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
										//LeftReads[Left_Index].Unique = true;
										
										{  
											//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
											if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
											else 
											{  
												//TD[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
													//TD[0].push_back(ReadIndex);
												//else  
												if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
													TD[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
													Reads[ReadIndex].Used = true;
													Count_TD++;
													Count_TD_Plus++;
												}
											}
										}
									}
								}
							}
						}
					}
					
				}
				else if (Reads[ReadIndex].MatchedD == Minus) {
					for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
						for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
							if (Reads[ReadIndex].Used) break;
							if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) {
								if (Reads[ReadIndex].Used) break;
								if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
									continue;
								if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
									if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
										 && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc) {
										Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
										Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
										Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
										
										Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1;;
										Reads[ReadIndex].NT_str = "";
										Reads[ReadIndex].NT_size = 0;
										Reads[ReadIndex].InsertedStr = "";
										Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
										//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
										//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
										//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
										//LeftReads[Left_Index].Unique = true;
										
										{
											//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
											if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
											else 
											{ 
												//TD[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
													//TD[0].push_back(ReadIndex);
												//else 
												if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
													TD[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
													Reads[ReadIndex].Used = true;
												
													Count_TD++;
													Count_TD_Minus++;
												}
												//cout << "- " << Count_D << endl; 
											}
										}
									}
								}
							}
						}
					}
					
				}
			}
		}
		cout << "Total: " << Count_TD << "\t+" << Count_TD_Plus << "\t-" << Count_TD_Minus << endl;
		ofstream TDOutf(TDOutputFilename.c_str(),ios::app);
		SortOutputTD(NumBoxes, CurrentChr, Reads, TD, TDOutf);
		//TDOutf.close();  
		for (unsigned int i = 0; i < NumBoxes; i++) TD[i].clear(); 
		
		cout << "Searching tandem dupliation events with non-template sequence ... " << endl;
		for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			//MAX_MISMATCHES_Per_Read = (short)(Seq_Error_Rate * Reads[ReadIndex].ReadLength + 1);
			CloseIndex = Reads[ReadIndex].UP_Close.size() - 1;
			FarIndex = Reads[ReadIndex].UP_Far.size() - 1;
			if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Reads[ReadIndex].ReadLength) continue;
			if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > (short)(1 + Seq_Error_Rate * (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr)))
				continue;
			//if (Reads[ReadIndex].UP_Far.size()) 
			{
				if (Reads[ReadIndex].MatchedD == Plus) {
					//for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) 
					{
						//if (Reads[ReadIndex].Used) break;
						//for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) 
						{
							//if (Reads[ReadIndex].Used) break;
							
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								if (Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr > Min_Num_Matched_Bases) {
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1;
									Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].UP_Far[FarIndex].LengthStr;
									Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{  
										//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{
											if (Reads[ReadIndex].NT_size <= Max_Length_NT && readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
												//	TD_NT[0].push_back(ReadIndex);
												//else 
													TD_NT[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//TD_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												Count_TD_NT++;
												Count_TD_NT_Plus++;												 
											}
										}
									}
								}
							}
						}
					}
				}
				else if (Reads[ReadIndex].MatchedD == Minus) {
					//for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) 
					{
						//if (Reads[ReadIndex].Used) break;
						//for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) 
						{
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								if (Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].ReadLength < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
									 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr > Min_Num_Matched_Bases) {
									
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1;
									Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].UP_Far[FarIndex].LengthStr;
									Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{
										//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{
											if (Reads[ReadIndex].NT_size <= Max_Length_NT && readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												//TD_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
													//TD_NT[0].push_back(ReadIndex);
												//else 
												TD_NT[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												
												Count_TD_NT++;
												Count_TD_NT_Minus++;
												//cout << "- " << Count_D << endl; 												 
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cout << "Total: " << Count_TD_NT << "\t+" << Count_TD_NT_Plus << "\t-" << Count_TD_NT_Minus << endl;
		//ofstream TDOutf(TDOutputFilename.c_str());
		SortOutputTD_NT(NumBoxes, CurrentChr, Reads, TD_NT, TDOutf);
		TDOutf.close();  
		for (unsigned int i = 0; i < NumBoxes; i++) TD_NT[i].clear(); 
	} // if (Analyze_TD_INV)
	if (Analyze_INV) {
		cout << "Searching inversions ... " << endl;
		for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			if (Reads[ReadIndex].UP_Close[0].Strand != Reads[ReadIndex].UP_Far[0].Strand
				 && 
				 Reads[ReadIndex].UP_Close[0].Direction == Reads[ReadIndex].UP_Far[0].Direction) {
				
				if (Reads[ReadIndex].MatchedD == Plus) {
					if (Reads[ReadIndex].UP_Far[0].AbsLoc > Reads[ReadIndex].UP_Close[Reads[ReadIndex].UP_Close.size() - 1].AbsLoc + MIN_IndelSize_Inversion) { // normal situation
						for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
							for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
								if (Reads[ReadIndex].Used) break;
								if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) {
									if (Reads[ReadIndex].Used) break;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
										continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
										if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion) {
											//cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc 
											//<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
											Reads[ReadIndex].Left = (Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
											Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].ReadLength;
											Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
											
											Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;
											Reads[ReadIndex].NT_str = "";
											Reads[ReadIndex].NT_size = 0;
											Reads[ReadIndex].InsertedStr = "";
											Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
											//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
											//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
											//LeftReads[Left_Index].Unique = true;
											
											{  
												//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
												if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
												else 
												{
													//Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
													//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
													//if ((int)Reads[ReadIndex].BPLeft < 0)
														//Inv[0].push_back(ReadIndex);
													//else 
													if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {														
														Inv[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
														Reads[ReadIndex].Used = true;
														Count_Inv++;
														Count_Inv_Plus++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
					else if (Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Close[0].AbsLoc) {
						for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
							for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
								if (Reads[ReadIndex].Used) break;
								if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
									if (Reads[ReadIndex].Used) break;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
										continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
										// anchor inside reversed block.
										if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr == Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
											 ) {
											//cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc 
											//<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
											Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].ReadLength;
											Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
											Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
											
											Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;
											Reads[ReadIndex].NT_str = "";
											Reads[ReadIndex].NT_size = 0;
											Reads[ReadIndex].InsertedStr = "";
											Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											Reads[ReadIndex].BPLeft = (Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											// CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
											// CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
											//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
											//LeftReads[Left_Index].Unique = true;
											
											{  
												//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
												if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
												else 
												{
													//Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
													//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
													//if ((int)Reads[ReadIndex].BPLeft < 0)
														//Inv[0].push_back(ReadIndex);
													//else 
													if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
														Inv[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
														Reads[ReadIndex].Used = true;
														Count_Inv++;
														Count_Inv_Plus++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
				else if (Reads[ReadIndex].MatchedD == Minus) {
					if (Reads[ReadIndex].UP_Close[Reads[ReadIndex].UP_Close.size() - 1].AbsLoc > Reads[ReadIndex].UP_Far[0].AbsLoc + MIN_IndelSize_Inversion) {
						for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
							for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
								if (Reads[ReadIndex].Used) break;
								if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								for (int FarIndex = 0; FarIndex < Reads[ReadIndex].UP_Far.size(); FarIndex++) {
									if (Reads[ReadIndex].Used) break;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
										continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
										// ######################
										//cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr << "\t" << Reads[ReadIndex].ReadLength << "\t";
										//cout << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc << "\t>\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1 << endl;
										// anchor outside reversed block.
										if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
											 && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion) {
											//cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc 
											//<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
											Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].ReadLength;
											Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
											Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
											
											Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;
											Reads[ReadIndex].NT_str = "";
											Reads[ReadIndex].NT_size = 0;
											Reads[ReadIndex].InsertedStr = "";
											Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											//cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
											//cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
											Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
											//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
											//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
											//LeftReads[Left_Index].Unique = true;
											
											{
												//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
												{ 
													//cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl; 
													//Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
													//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - //lowerBinBorder;
													//if ((int)Reads[ReadIndex].BPLeft < 0)
														//Inv[0].push_back(ReadIndex);
													//else 
													if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
														Inv[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
														Reads[ReadIndex].Used = true;
													
														Count_Inv++;
														Count_Inv_Minus++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
					else if (Reads[ReadIndex].UP_Close[0].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Far[Reads[ReadIndex].UP_Far.size() - 1].AbsLoc) {
						for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
							for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
								if (Reads[ReadIndex].Used) break;
								if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
								for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
									if (Reads[ReadIndex].Used) break;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
										continue;
									if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
										// anchor inside reversed block.
										if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr == Reads[ReadIndex].ReadLength
											 && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
											 ) {
											Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr -  1;
											Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].ReadLength;
											Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
											
											Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;
											Reads[ReadIndex].NT_str = "";
											Reads[ReadIndex].NT_size = 0;
											Reads[ReadIndex].InsertedStr = "";
											Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											//cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
											//cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
											Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
											//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
											//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
											//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
											//LeftReads[Left_Index].Unique = true;
											
											{
												//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
												{ 
													//cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl; 
													//Inv[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
													//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - //lowerBinBorder;
													//if ((int)Reads[ReadIndex].BPLeft < 0)
														//Inv[0].push_back(ReadIndex);
													//else 
													if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
														Inv[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
														Reads[ReadIndex].Used = true;
													
														Count_Inv++;
														Count_Inv_Minus++;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cout << "Total: " << Count_Inv << "\t+" << Count_Inv_Plus << "\t-" << Count_Inv_Minus << endl;
		ofstream InversionOutf(InversionOutputFilename.c_str(),ios::app);
		SortOutputInv(NumBoxes, CurrentChr, Reads, Inv, InversionOutf);
		//InversionOutf.close();
		for (unsigned int i = 0; i < NumBoxes; i++) Inv[i].clear(); 
		
		cout << "Searching inversions with non-template sequence ... " << endl;
		for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
			if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
			CloseIndex = Reads[ReadIndex].UP_Close.size() - 1;
			FarIndex = Reads[ReadIndex].UP_Far.size() - 1;
			if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > (short)(1 + Seq_Error_Rate * (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr)))
				continue;
			if (Reads[ReadIndex].UP_Close[0].Strand != Reads[ReadIndex].UP_Far[0].Strand
				 && 
				 Reads[ReadIndex].UP_Close[0].Direction == Reads[ReadIndex].UP_Far[0].Direction) {
				if (Reads[ReadIndex].MatchedD == Plus) {
					{
						{
							//if (Reads[ReadIndex].Used) break;
							
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength 
									 && 
									 Reads[ReadIndex].UP_Far[FarIndex].AbsLoc > Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion
									 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases) {
									//cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc 
									//<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									Reads[ReadIndex].Left = (Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1) - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].ReadLength;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize =  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;
									
									Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr; // NT_2str
									//cout << "Po " << Reads[ReadIndex].NT_size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{  
										//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{ 
											if (Reads[ReadIndex].NT_size <= Max_Length_NT) {
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
													//Inv_NT[0].push_back(ReadIndex);
												//else 
												if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
													Inv_NT[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
													Reads[ReadIndex].Used = true;
													Count_Inv_NT++;
													Count_Inv_NT_Plus++;												 
												}
											}
										}
									}
								}
								// anchor inside reversed block.
								if (Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength && Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc
									 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases
									 ) {
									//cout << "+\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc 
									//<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].ReadLength;
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;
									
									Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									//cout << "Pi " << Reads[ReadIndex].NT_size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									Reads[ReadIndex].BPLeft = (Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1) - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{  
										//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{
											if (Reads[ReadIndex].NT_size <= Max_Length_NT && readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												//Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
													//Inv_NT[0].push_back(ReadIndex);
												//else 
												Inv_NT[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												Count_Inv_NT++;
												Count_Inv_NT_Plus++;												 
											}
										}
									}
								}
							}
						}
					}
				}
				else if (Reads[ReadIndex].MatchedD == Minus) {
					//for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) 
					{
						//if (Reads[ReadIndex].Used) break;
						//for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) 
						{
							//if (Reads[ReadIndex].Used) break;
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								// ######################
								//cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr << "\t" << Reads[ReadIndex].ReadLength << "\t";
								//cout << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc << "\t>\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1 << endl;
								// anchor outside reversed block.
								if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength
									 && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc > Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + MIN_IndelSize_Inversion
									 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases) {
									//cout << "-\t" << Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc 
									//<< "\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << endl;
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].ReadLength;
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].AbsLoc;
									
									Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									//cout << "Mo " << Reads[ReadIndex].NT_2size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									Reads[ReadIndex].NT_str = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
									//cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{
										//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{
											//cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl; 
											if (Reads[ReadIndex].NT_size <= Max_Length_NT &&  readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												//Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
													//Inv_NT[0].push_back(ReadIndex);
												//else 
												Inv_NT[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												
												Count_Inv_NT++;
												Count_Inv_NT_Minus++;												 
											}
										}
									}
								}
								// anchor inside reversed block.
								if (Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength
									 && Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + MIN_IndelSize_Inversion < Reads[ReadIndex].UP_Far[FarIndex].AbsLoc
									 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr >= Min_Num_Matched_Bases
									 ) {
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr -  1;
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - Reads[ReadIndex].ReadLength;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc;
									
									Reads[ReadIndex].NT_size = Reads[ReadIndex].ReadLength - Reads[ReadIndex].UP_Far[FarIndex].LengthStr - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr;
									//cout << "Mi " << Reads[ReadIndex].NT_2size << "\t" <<  Reads[ReadIndex].ReadLength << "\t" << Reads[ReadIndex].UP_Close[CloseIndex].LengthStr << "\t" << Reads[ReadIndex].UP_Far[FarIndex].LengthStr << endl;
									Reads[ReadIndex].NT_str = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].NT_size);
									Reads[ReadIndex].InsertedStr = "";
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//cout <<  "far\t" << Reads[ReadIndex].UP_Far[FarIndex].AbsLoc << "\tspacer\t" << SpacerBeforeAfter << endl;
									//cout <<  Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter - 1 << endl;
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - 1 - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + Const_S + LOG14 * Reads[ReadIndex].ReadLength + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{
										//if (LeftReads[Left_Index].IndelSize <= DSizeArray[RangeIndex] && LeftReads[Left_Index].IndelSize > DSizeArray[RangeIndex - 1]) 
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{ 
											//cout << "Inv\t" << Reads[ReadIndex].BPLeft << "\t" << BoxSize << endl; 
											if (Reads[ReadIndex].NT_size <= Max_Length_NT && readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												//Inv_NT[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
												//if ((int)Reads[ReadIndex].BPLeft < 0)
													//Inv_NT[0].push_back(ReadIndex);
												//else 
												Inv_NT[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												
												Count_Inv_NT++;
												Count_Inv_NT_Minus++;												 
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		cout << "Total: " << Count_Inv_NT << "\t+" << Count_Inv_NT_Plus << "\t-" << Count_Inv_NT_Minus << endl;
		//ofstream InversionOutf(InversionOutputFilename.c_str());
		SortOutputInv_NT(NumBoxes, CurrentChr, Reads, Inv_NT, InversionOutf);
		//InversionOutf.close();
		for (unsigned int i = 0; i < NumBoxes; i++) Inv_NT[i].clear(); 		 
		
	} // Analyze_INV
	
	cout << "Searching short insertions ... " << endl;
	for (unsigned ReadIndex = 0; ReadIndex < Reads.size(); ReadIndex++) {
		if (Reads[ReadIndex].Used || Reads[ReadIndex].UP_Far.empty()) continue;
		if (!Reads[ReadIndex].UP_Far.empty()) 
		{     
			if (Reads[ReadIndex].MatchedD == Plus) {
				for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
					for (int CloseIndex = 0; CloseIndex < Reads[ReadIndex].UP_Close.size(); CloseIndex++) {
						//cout << "+" << CloseIndex << endl;
						if (Reads[ReadIndex].Used) break;
						if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
						for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
							//cout << "+" << FarIndex << endl;
							if (Reads[ReadIndex].Used) break;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
								continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Minus) {
								if (Reads[ReadIndex].UP_Far[FarIndex].AbsLoc == Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1 && Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + Reads[ReadIndex].UP_Far[FarIndex].LengthStr < Reads[ReadIndex].ReadLength) {
									
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - Reads[ReadIndex].UP_Close[CloseIndex].LengthStr + 1;
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize = Reads[ReadIndex].ReadLengthMinus - (Reads[ReadIndex].Right - Reads[ReadIndex].Left);
									Reads[ReadIndex].InsertedStr = ReverseComplement(Reads[ReadIndex].UnmatchedSeq).substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].IndelSize);
									Reads[ReadIndex].NT_str = "";
									//                                 LeftReads[Left_Index].InsertedStr = "";
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + LOG14 * Reads[ReadIndex].ReadLength - LOG14 * Reads[ReadIndex].IndelSize + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
									{
										//if (RangeIndex == 0) 
										if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
										else 
										{
											//SIs[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
											//if ((int)Reads[ReadIndex].BPLeft < 0)
												//SIs[0].push_back(ReadIndex);
											//else 
											if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
												SIs[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
												Reads[ReadIndex].Used = true;
												Count_SI_Plus++;
												Count_SI++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			else if (Reads[ReadIndex].MatchedD == Minus) {
				for (short MAX_SNP_ERROR_index = 0; MAX_SNP_ERROR_index <= Reads[ReadIndex].MAX_SNP_ERROR; MAX_SNP_ERROR_index++) {
					for (int CloseIndex = Reads[ReadIndex].UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
						//cout << "-: " << CloseIndex << endl;
						if (Reads[ReadIndex].Used) break;
						if (Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index) continue;
						for (int FarIndex = Reads[ReadIndex].UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
							//cout << "-: " << FarIndex << endl;
							if (Reads[ReadIndex].Used) break;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches > MAX_SNP_ERROR_index) continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Mismatches + Reads[ReadIndex].UP_Close[CloseIndex].Mismatches > MAX_SNP_ERROR_index)
								continue;
							if (Reads[ReadIndex].UP_Far[FarIndex].Direction == Plus) {
								if (Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc == Reads[ReadIndex].UP_Far[FarIndex].AbsLoc + 1 && Reads[ReadIndex].UP_Far[FarIndex].LengthStr + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr < Reads[ReadIndex].ReadLength) {
									
									Reads[ReadIndex].Left = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - Reads[ReadIndex].UP_Far[FarIndex].LengthStr + 1;
									Reads[ReadIndex].Right = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc + Reads[ReadIndex].UP_Close[CloseIndex].LengthStr - 1;
									Reads[ReadIndex].BP = Reads[ReadIndex].UP_Far[FarIndex].LengthStr - 1;
									
									Reads[ReadIndex].IndelSize = Reads[ReadIndex].ReadLengthMinus - (Reads[ReadIndex].Right - Reads[ReadIndex].Left);
									Reads[ReadIndex].InsertedStr = Reads[ReadIndex].UnmatchedSeq.substr(Reads[ReadIndex].BP + 1, Reads[ReadIndex].IndelSize);
									Reads[ReadIndex].NT_str = "";
									//   LeftReads[Left_Index].InsertedStr = "";
									Reads[ReadIndex].BPLeft = Reads[ReadIndex].UP_Far[FarIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									Reads[ReadIndex].BPRight = Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc - SpacerBeforeAfter;// - Fragments[LeftReads[Left_Index].MatchedSeqID].Start;
									//CurrentChrMask[Reads[ReadIndex].UP_Close[CloseIndex].AbsLoc] = 'B';
									//CurrentChrMask[Reads[ReadIndex].UP_Far[FarIndex].AbsLoc] = 'B';
									//EWL070111 Reads[ReadIndex].score = Const_I + LOG14 * Reads[ReadIndex].ReadLength - LOG14 * Reads[ReadIndex].IndelSize + Const_Log_T;
									//LeftReads[Left_Index].Unique = true;
									
  									if (readTransgressesBinBoundaries(Reads[ReadIndex],upperBinBorder)) { saveReadForNextCycle(Reads[ReadIndex], FutureReads ); }
									else 
									{
										//SIs[Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
										//int (int)Reads[ReadIndex].BPLeft = (int)Reads[ReadIndex].BPLeft - lowerBinBorder;
										//if ((int)Reads[ReadIndex].BPLeft < 0)
											//SIs[0].push_back(ReadIndex);
										//else 
										if ( readInSpecifiedRegion( Reads[ReadIndex], startOfRegion, endOfRegion ) ) {
											SIs[(int)Reads[ReadIndex].BPLeft / BoxSize].push_back(ReadIndex);
											Reads[ReadIndex].Used = true;
											Count_SI++;
											Count_SI_Minus++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	cout << "Total: " << Count_SI << "\t+" << Count_SI_Plus << "\t-" << Count_SI_Minus << endl;
	ofstream SIoutputfile(SIOutputFilename.c_str(),ios::app);
	SortOutputSI(NumBoxes, CurrentChr, Reads, SIs, SIoutputfile);
	SIoutputfile.close();
	for (unsigned int i = 0; i < NumBoxes; i++) SIs[i].clear();
	
	
	
	int TotalNumReads = Reads.size();
	if (ReportCloseMappedRead) {
		string CloseEndMappedOutputFilename = OutputFolder + "_CloseEndMapped";
		ofstream CloseEndMappedOutput(CloseEndMappedOutputFilename.c_str(),ios::app);
		for (int Index = 0; Index < TotalNumReads; Index++) {
			CloseEndMappedOutput << Reads[Index].Name << "\n"
			<< Reads[Index].UnmatchedSeq << "\n"
			<< Reads[Index].MatchedD << "\t"
			<< Reads[Index].FragName << "\t"
			<< Reads[Index].MatchedRelPos << "\t"
			<< Reads[Index].MS << "\t"
			<< Reads[Index].InsertSize << "\t"
			<< Reads[Index].Tag << "\n";
		}
	} 
	
	if (ReportSVReads) {
		string SVReadOutputFilename = OutputFolder + "_SVReads";
		ofstream SVReadOutput(SVReadOutputFilename.c_str(),ios::app);
		for (int Index = 0; Index < TotalNumReads; Index++) {
			if (Reads[Index].IndelSize > Indel_SV_cutoff || Reads[Index].IndelSize == 0)
				SVReadOutput << Reads[Index].Name << "\n"
				<< Reads[Index].UnmatchedSeq << "\n"
				<< Reads[Index].MatchedD << "\t"
				<< Reads[Index].FragName << "\t"
				<< Reads[Index].MatchedRelPos << "\t"
				<< Reads[Index].MS << "\t"
				<< Reads[Index].InsertSize << "\t"
				<< Reads[Index].Tag << "\n";
		}
	} 
	
	if (ReportLargeInterChrSVReads) {
		string LargeInterChrSVReadsOutputFilename = OutputFolder + "_LargeORInterChrReads";
		ofstream LargeInterChrSVReadsOutput(LargeInterChrSVReadsOutputFilename.c_str(),ios::app);
		for (int Index = 0; Index < TotalNumReads; Index++) {
			if (Reads[Index].IndelSize == 0)
				LargeInterChrSVReadsOutput << Reads[Index].Name << "\n"
				<< Reads[Index].UnmatchedSeq << "\n"
				<< Reads[Index].MatchedD << "\t"
				<< Reads[Index].FragName << "\t"
				<< Reads[Index].MatchedRelPos << "\t"
				<< Reads[Index].MS << "\t"
				<< Reads[Index].InsertSize << "\t"
				<< Reads[Index].Tag << "\n";
		}
	} 
	
	
	unsigned Count_Close = 0;
	unsigned Count_Far = 0;
	unsigned Count_Used = 0;
	unsigned Count_Unused = 0;
	//vector <SPLIT_READ> UnusedClosedReads;
	
	for (int Index = TotalNumReads - 1; Index >= 0; Index--) {
		//cout << Index << endl;
		//if (!Reads[Index].UP_Close.empty()) Count_Close++;
		if (!Reads[Index].UP_Far.empty()) Count_Far++;
		if (!Reads[Index].UP_Far.empty() || Reads[Index].Found) {
			
		} 
		else Count_Unused++;//UnusedClosedReads.push_back(Reads[Index]);
		//cout << "1" << endl;
		if (Reads[Index].Used) Count_Used++;
		//cout << "2" << endl;
		//Reads.pop_back();
	}
	
	cout << "Total: " << TotalNumReads << ";\tClose_end_found " << TotalNumReads << ";\tFar_end_found " << Count_Far << ";\tUsed\t" << Count_Used << ".\n\n";
	cout << "For LI and BP: " << Count_Unused << endl << endl;
	
	if (Analyze_LI) {
		time_t Time_LI_S, Time_LI_E;
		Time_LI_S = time(NULL); 
		ofstream LargeInsertionOutf(LargeInsertionOutputFilename.c_str(),ios::app); 
		SortOutputLI(CurrentChr, Reads, LargeInsertionOutf);
		LargeInsertionOutf.close();
		Time_LI_E = time(NULL);
		cout << "Mining, Sorting and output LI results: " << (unsigned int)difftime(Time_LI_E, Time_LI_S) << " seconds." << endl << endl;;
	}
	
	BP_Reads.clear();
	if (Analyze_BP) {	
		//cout << "BP" << endl;
		time_t Time_BP_S, Time_BP_E;
		Time_BP_S = time(NULL); 
		ofstream RestOutf(RestOutputFilename.c_str(),ios::app);
		SortOutputRest(CurrentChr, Reads, BP_Reads, RestOutf);
		RestOutf.close();	
		Time_BP_E = time(NULL);
		cout << "Mining, Sorting and output BP results: " << (unsigned int)difftime(Time_BP_E, Time_BP_S) << " seconds." << endl << endl;
	}
	
    Time_Sort_E = time(NULL);   
    //Time_Load_E = time(NULL);
    //cout << "SI: " << Count_SI_Plus << " " << Count_SI_Minus << "\n"
    //<< "D: " << Count_D_Plus << " " << Count_D_Minus << endl;
    //cout << Entering_D_Plus << " " << Entering_D_Minus << endl;
    //cout << Plus_Sum_Left << " " << Plus_Sum_Right << "\n" << Minus_Sum_Right << " " << Minus_Sum_Left << endl;

    AllLoadings += (unsigned int)difftime(Time_Load_E, Time_Load_S);
    //AllMinings += (unsigned int)difftime(Time_Mine_E, Time_Load_E);
    AllSortReport += (unsigned int)difftime(Time_Sort_E, Time_Load_E);
    //Time_Load_S = time(NULL); 
   Reads.clear(); 
	cout << "I have " << FutureReads.size() << " reads saved for the next cycle.\n";
   Reads.swap(FutureReads);   
   Time_Load_S=0; 

} while ( 
				( PindelReadDefined && !isFinishedPindel( upperBinBorder, endRegionPlusBuffer ) )
			|| ( BAMDefined && !isFinishedBAM( upperBinBorder, endRegionPlusBuffer, CurrentChr.size() ) )
		  );

}// while ( loopOverAllChromosomes && chromosomeIndex < chromosomes.size() );
    cout << "Loading genome sequences and reads: " << AllLoadings << " seconds." << endl;
    //cout << "Mining indels: " << AllMinings << " seconds." << endl;
    cout << "Mining, Sorting and output results: " << AllSortReport << " seconds." << endl;
    return 0;
    } //main


void ReadInOneChr(ifstream & inf_Seq, string & TheInput, const string & ChrName) {
    TheInput.clear();
    char TempChar;
    string TempLine, TempChrName;
    inf_Seq >> TempChar;
    if (TempChar != '>') {
        cout << "Please use fasta format for the reference file." << endl;
        TheInput.clear();
        return;
    }
    while (inf_Seq >> TempChrName) {
        cout << "Processing chromosome " << TempChrName << endl;
        if (!TheInput.empty()) {
            cout << "Skip the rest of chromosomes.\n";
            break;	
        }
        if (TempChrName == ChrName) {
            getline(inf_Seq, TempLine);
            while (inf_Seq >> TempChar) {
                if (TempChar != '\n' && TempChar != '\r') {
                    if (TempChar == '>') {
                        break;
                    }
                    else {
                        if ( 'a' <= TempChar && TempChar <= 'z' )
                            TempChar = TempChar + ( 'A' - 'a' );
                        switch (TempChar) {
                            case 'A': TheInput += 'A'; break;  // 00000000
                            case 'C': TheInput += 'C'; break; // 00010000
                            case 'G': TheInput += 'G'; break; // 00100000
                            case 'T': TheInput += 'T'; break; // 00110000
                            default:  TheInput += 'N';        // 01000000
                        }
                    } // else TempChar
                }
            }
        }
        else {
            getline(inf_Seq, TempLine);
            while (inf_Seq >> TempChar) {
                if (TempChar != '\n' && TempChar != '\r') {
                    if (TempChar == '>') {
                        break;
                    }
                    else {

                    } // else TempChar
                }
            }
        }
    }
    cout << ChrName << "\t" << TheInput.size() << "\t";
    if (!TheInput.empty()) {
        string Spacer = "";
        for (unsigned i = 0; i < SpacerBeforeAfter; i++) Spacer += "N";
        TheInput = Spacer + TheInput + Spacer;
    }
    cout << TheInput.size() << endl;

	// EWL280311: reset file for to allow re-reading for the next chromosome
	inf_Seq.clear();
	inf_Seq.seekg(0);

   return;
}

void GetOneChrSeq(ifstream & inf_Seq, string & TheInput, bool WhetherBuildUp) {
	TheInput.clear();
	char TempChar;
	string TempLine, TempChrName;

	if (WhetherBuildUp) {
		//getline(inf_Seq, TempLine);
		while (inf_Seq >> TempChar) {
			if (TempChar != '\n' && TempChar != '\r') {
				if (TempChar == '>') {
					break;
				}
				else {
					if ( 'a' <= TempChar && TempChar <= 'z' )
						TempChar = TempChar + ( 'A' - 'a' );
					switch (TempChar) {
						case 'A': TheInput += 'A'; break;  // 00000000
						case 'C': TheInput += 'C'; break; // 00010000
						case 'G': TheInput += 'G'; break; // 00100000
						case 'T': TheInput += 'T'; break; // 00110000
						default:  TheInput += 'N';        // 01000000
					}
				} // else TempChar
			}
		}
	}
	else {
		//getline(inf_Seq, TempLine);
		while (inf_Seq >> TempChar) {
			if (TempChar != '\n' && TempChar != '\r') {
				if (TempChar == '>') {
					break;
				}
				else {
					
				} // else TempChar
			}
		}
	}

	if (!TheInput.empty()) {
		string Spacer = "";
		for (unsigned i = 0; i < SpacerBeforeAfter; i++) Spacer += "N";
		TheInput = Spacer + TheInput + Spacer;
	}
   return;
}


short ReadInRead(ifstream & inf_ReadSeq, const string & FragName, const string & CurrentChr, vector <SPLIT_READ> & Reads, 
	const int lowerBinBorder, const int upperBinBorder)
{
	cout << "Scanning and processing reads anchored in " << FragName << endl; 
	//short ADDITIONAL_MISMATCH = 1;
	SPLIT_READ Temp_One_Read;
	unsigned int NumReadScanned = 0;
	unsigned int NumReadInChr = 0;
	unsigned int InChrPlus = 0;
	unsigned int InChrMinus = 0;
	unsigned int GetPlus = 0;
	unsigned int GetMinus = 0;
	//NumberOfReadsPerBuffer;
	vector <SPLIT_READ> BufferReads;
	ReportLength = 0;
	//cout << LeftReads.size() << endl;
	string TempQC, TempLine, TempStr, TempFragName;
	//int TempInt;
 
   inf_ReadSeq.clear();
   inf_ReadSeq.seekg( 0 );
	VectorTag.clear(); 
//cout << "MINUSEXTRA!\n";
int UPCLOSE_COUNTER =0;		
	//loop over all reads in the file
	while (inf_ReadSeq >> Temp_One_Read.Name) {
		if (Temp_One_Read.Name[0] != FirstCharReadName) { // !='@'
			cout << "Something wrong with the read name: " << Temp_One_Read.Name << endl;
			Reads.clear();
			return 1;
		}
		NumReadScanned++;
      // get (useless) rest of first line
		getline(inf_ReadSeq, TempLine);
		inf_ReadSeq >> Temp_One_Read.UnmatchedSeq;
		getline(inf_ReadSeq, TempLine);

		inf_ReadSeq >> Temp_One_Read.MatchedD;
		if (Temp_One_Read.MatchedD != Minus && Temp_One_Read.MatchedD != Plus) {
			cout << Temp_One_Read.Name << "\n"
			<< Temp_One_Read.UnmatchedSeq << "\n"
			<< Temp_One_Read.MatchedD << " ...\n";
			cout << "+/-" << endl;
			return 1;
		}
		//   >> TempInt 
		inf_ReadSeq >> Temp_One_Read.FragName
		>> Temp_One_Read.MatchedRelPos
		>> Temp_One_Read.MS
		>> Temp_One_Read.InsertSize
		>> Temp_One_Read.Tag;
		getline(inf_ReadSeq, TempLine);

		
		if (Temp_One_Read.MatchedRelPos > g_maxPos ) { g_maxPos = Temp_One_Read.MatchedRelPos; }

		if (Temp_One_Read.FragName == FragName
            && Temp_One_Read.MatchedRelPos>=lowerBinBorder && Temp_One_Read.MatchedRelPos<upperBinBorder) {
			Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
			Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
			NumReadInChr++; 

			Temp_One_Read.MAX_SNP_ERROR = (short)(Temp_One_Read.UnmatchedSeq.size() * Seq_Error_Rate);

			
			Temp_One_Read.TOTAL_SNP_ERROR_CHECKED = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
			Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
			Temp_One_Read.MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;
			//Temp_One_Read.IndelSize = 0;
			Temp_One_Read.Found = false;
			if (Temp_One_Read.MatchedD == Plus)  {
				InChrPlus++;
				//Temp_One_Read.UnmatchedSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
			}
			else InChrMinus++; // this seems to be going correctly
			if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size) Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
			if (Temp_One_Read.MatchedRelPos < 1) Temp_One_Read.MatchedRelPos = 0;
			BufferReads.push_back(Temp_One_Read);
			if (BufferReads.size() >= NumberOfReadsPerBuffer) {
				//cout << "here" << endl;
				#pragma omp parallel default(shared)
				{
					#pragma omp for
				   for (int BufferReadsIndex = 0; BufferReadsIndex < NumberOfReadsPerBuffer; BufferReadsIndex++) {
                  GetCloseEnd(CurrentChr, BufferReads[BufferReadsIndex]);	
					}
				}
				for (int BufferReadsIndex = 0; BufferReadsIndex < NumberOfReadsPerBuffer; BufferReadsIndex++) {
					if (BufferReads[BufferReadsIndex].UP_Close.size()) {
						UPCLOSE_COUNTER++;
						if (ReportLength < BufferReads[BufferReadsIndex].ReadLength) ReportLength = BufferReads[BufferReadsIndex].ReadLength;
						BufferReads[BufferReadsIndex].Used = false;
						BufferReads[BufferReadsIndex].Unique = true;
						//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t";
						CleanUniquePoints(BufferReads[BufferReadsIndex].UP_Close);
						//cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl;
						BufferReads[BufferReadsIndex].CloseEndLength = 
						BufferReads[BufferReadsIndex].UP_Close[BufferReads[BufferReadsIndex].UP_Close.size() - 1].LengthStr;
						if (BufferReads[BufferReadsIndex].MatchedD == Plus) 
						   BufferReads[BufferReadsIndex].LeftMostPos = BufferReads[BufferReadsIndex].UP_Close[BufferReads[BufferReadsIndex].UP_Close.size() - 1].AbsLoc + 1 - BufferReads[BufferReadsIndex].CloseEndLength;
						else BufferReads[BufferReadsIndex].LeftMostPos = BufferReads[BufferReadsIndex].UP_Close[BufferReads[BufferReadsIndex].UP_Close.size() - 1].AbsLoc + BufferReads[BufferReadsIndex].CloseEndLength - BufferReads[BufferReadsIndex].ReadLength;
                  //cout << "Pushing back!\n";
						Reads.push_back(BufferReads[BufferReadsIndex]);
						if (BufferReads[BufferReadsIndex].MatchedD == Plus) GetPlus++;
						else GetMinus++;
						if (NotInVector(BufferReads[BufferReadsIndex].Tag, VectorTag)) {
							VectorTag.push_back(BufferReads[BufferReadsIndex].Tag);
						}
					}
				}
				BufferReads.clear();
			} // if buffer-reads threatens to overflow
		} // if the read is in the correct bin
	}  // while loop over each read
	cout << "last one: " << BufferReads.size() << " and UPCLOSE= " << UPCLOSE_COUNTER <<endl;
#pragma omp parallel default(shared)
	{
      #pragma omp for
		for (int BufferReadsIndex = 0; BufferReadsIndex < BufferReads.size(); BufferReadsIndex++) {
			GetCloseEnd(CurrentChr, BufferReads[BufferReadsIndex]);	
		}
	}
	for (int BufferReadsIndex = 0; BufferReadsIndex < BufferReads.size(); BufferReadsIndex++) {
		if (BufferReads[BufferReadsIndex].UP_Close.size()) {
						UPCLOSE_COUNTER++;
			if (ReportLength < BufferReads[BufferReadsIndex].ReadLength) ReportLength = BufferReads[BufferReadsIndex].ReadLength;
			BufferReads[BufferReadsIndex].Used = false;
			BufferReads[BufferReadsIndex].Unique = true;
			//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t";
			CleanUniquePoints(BufferReads[BufferReadsIndex].UP_Close);
			//cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl;
			BufferReads[BufferReadsIndex].CloseEndLength = 
			BufferReads[BufferReadsIndex].UP_Close[BufferReads[BufferReadsIndex].UP_Close.size() - 1].LengthStr;
			if (BufferReads[BufferReadsIndex].MatchedD == Plus) 
				BufferReads[BufferReadsIndex].LeftMostPos = BufferReads[BufferReadsIndex].UP_Close[BufferReads[BufferReadsIndex].UP_Close.size() - 1].AbsLoc + 1 - BufferReads[BufferReadsIndex].CloseEndLength;
			else BufferReads[BufferReadsIndex].LeftMostPos = BufferReads[BufferReadsIndex].UP_Close[BufferReads[BufferReadsIndex].UP_Close.size() - 1].AbsLoc + BufferReads[BufferReadsIndex].CloseEndLength - BufferReads[BufferReadsIndex].ReadLength;
			Reads.push_back(BufferReads[BufferReadsIndex]);
			if (BufferReads[BufferReadsIndex].MatchedD == Plus) GetPlus++;
			else GetMinus++;
			if (NotInVector(BufferReads[BufferReadsIndex].Tag, VectorTag)) {
				VectorTag.push_back(BufferReads[BufferReadsIndex].Tag);
			}
		}
	}
	BufferReads.clear();
	
	if (FirstChr) {
		cout << "\nThe last read Pindel scanned: \n";
		cout << Temp_One_Read.Name << "\n"
		<< Temp_One_Read.UnmatchedSeq << "\n"
		<< Temp_One_Read.MatchedD << "\t"
		<< Temp_One_Read.FragName << "\t"
		<< Temp_One_Read.MatchedRelPos << "\t"
		<< Temp_One_Read.MS << "\t"
		<< Temp_One_Read.InsertSize << "\t"
		<< Temp_One_Read.Tag << "\n\n";
		FirstChr = false;
      //ReportLength = Temp_One_Read.UnmatchedSeq.size();
	} 
	cout << "NumReadScanned:\t" << NumReadScanned << endl;
	cout << "NumReadInChr:\t" << NumReadInChr << endl;
	cout << "NumReadStored:\t" << Reads.size() << endl;
	cout << "NumReadStored / NumReadInChr = " << Reads.size() * 100.0 / NumReadInChr << " %\n" 
	<< "InChrPlus \t" << InChrPlus << "\tGetPlus \t" << GetPlus << "\t" << GetPlus * 100.0 / InChrPlus 
	<< " %\n" << "InChrMinus\t" << InChrMinus  << "\tGetMinus\t" << GetMinus 
	<< "\t" << GetMinus * 100.0 / InChrMinus << " %\n" << endl;
	//inf_ReadSeq.close();
	if (Reads.size() == 0) return 0;
	//cout << LeftReads.size() << endl;
	cout << "sorting tags ... ";
	string Str4Exchange;
	for (short i = 0; i < VectorTag.size() - 1; i++) {
		for (short j = 1; j < VectorTag.size(); j++) {
			if (VectorTag[i].size() > VectorTag[j].size()) {
				Str4Exchange = VectorTag[i];
				VectorTag[i] = VectorTag[j];
				VectorTag[j] = Str4Exchange; 
			}
			else if (VectorTag[i].size() == VectorTag[j].size()) {
				for (short k = 0; k < VectorTag[i].size(); k++) {
					if ((short)VectorTag[i][k] > (short)VectorTag[j][k]) {
       	         Str4Exchange = VectorTag[i];
						VectorTag[i] = VectorTag[j];
						VectorTag[j] = Str4Exchange;
						break;
					}
				}
			}
		}
	}
	cout << " finished!" << endl;
	return 0;
}



vector <string> ReverseComplement(const vector <string> & InputPatterns) {
    vector <string> OutputPatterns;// = InputPatterns;
    unsigned int NumPattern = InputPatterns.size();
    OutputPatterns.resize(NumPattern);
    for (unsigned int i = 0; i < NumPattern; i++) {
        OutputPatterns[i] = ReverseComplement(InputPatterns[i]);    
    }    
    return OutputPatterns;
}


string Reverse(const string & InputPattern) {
   string OutputPattern = InputPattern;
   unsigned int LenPattern = InputPattern.size();
   for (unsigned int j = 0; j < LenPattern; j++)
      OutputPattern[j] = InputPattern[j];
   return OutputPattern;
}

string ReverseComplement(const string & InputPattern) {
   string OutputPattern = InputPattern;
   unsigned int LenPattern = InputPattern.size();
   for (unsigned int j = 0; j < LenPattern; j++)
      OutputPattern[j] = Convert2RC4N[(unsigned int)InputPattern[LenPattern - j - 1]];
   return OutputPattern;
}

//const short MAX_SNP_ERROR = 2;
//const short ADDITIONAL_MISMATCH = 2;
//const short TOTAL_SNP_ERROR_CHECKED = MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;

void CheckLeft_Close(const SPLIT_READ & OneRead,
							const string & TheInput, 
							const string & CurrentReadSeq,
							const vector <unsigned int> Left_PD[],
							const short & BP_Left_Start,
							const short & BP_Left_End,
							const short & CurrentLength,
							vector <UniquePoint> & LeftUP) {
	int Sum;
   if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
      // put it to LeftUP if unique
		for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			if (Left_PD[i].size() == 1 && CurrentLength >= BP_Left_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH) 
		   		for (short j = 1; j <= ADDITIONAL_MISMATCH; j++) 
			   		Sum += Left_PD[i + j].size();
				
				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Left_PD[i][0]; 
					TempOne.Direction = FORWARD;
					TempOne.Strand = ANTISENSE;
					TempOne.Mismatches = i;
					LeftUP.push_back(TempOne);	
					break;
				}
			}
		}
   }
   if (CurrentLength < BP_Left_End) {
      vector <unsigned int> Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Left_PD_Output[CheckedIndex].reserve(Left_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[CurrentLength];
      //const int SizeOfCurrent = Left_PD.size();
      unsigned int pos;
		for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
			int SizeOfCurrent = Left_PD[i].size();
			if (CurrentChar == 'N') {
            //#pragma omp parallel default(shared)
				{
					// #pragma omp for
				   for (int j = 0; j < SizeOfCurrent; j++) {
					   pos =  Left_PD[i][j] + 1;
					   if (Match2N[(short)TheInput[pos]] == 'N')
						   Left_PD_Output[i].push_back(pos);
				   }     
				}
			}
			else {
            //#pragma omp parallel default(shared)
				{
               //#pragma omp for
				   for (int j = 0; j < SizeOfCurrent; j++) {
					   pos =  Left_PD[i][j] + 1;
					   if (TheInput[pos] == CurrentChar)
						   Left_PD_Output[i].push_back(pos);
					   else Left_PD_Output[i + 1].push_back(pos); 
				   }  
				}
			}			
		}
		
		int SizeOfCurrent = Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
		if (CurrentChar == 'N') {
         //#pragma omp parallel default(shared)
			{
				// #pragma omp for
			   for (int j = 0; j < SizeOfCurrent; j++) {
				   pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
				   if (Match2N[(short)TheInput[pos]] == 'N')
					   Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
			   }
			}
		}
		else {
         //#pragma omp parallel default(shared)
			{
				// #pragma omp for
			   for (int j = 0; j < SizeOfCurrent; j++) {
				   pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
				   if (TheInput[pos] == CurrentChar)
					   Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
					//else Left_PD_Output[i + 1].push_back(pos); 
			   }   
			}
		}	
		Sum = 0;
		for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			Sum += Left_PD_Output[i].size();
		}
      if (Sum) {
         const short CurrentLengthOutput = CurrentLength + 1;
         CheckLeft_Close(OneRead, TheInput, CurrentReadSeq, Left_PD_Output,
								 BP_Left_Start, BP_Left_End,
								 CurrentLengthOutput, LeftUP);         
      }
      else return;
   }   
   else return;
}

void CheckRight_Close(const SPLIT_READ & OneRead,
							 const string & TheInput,
							 const string & CurrentReadSeq,
							 const vector <unsigned int> Right_PD[],
							 const short & BP_Right_Start,
							 const short & BP_Right_End,
							 const short & CurrentLength,
							 vector <UniquePoint> & RightUP) {
	//cout << CurrentLength << "\t" << RightUP.size() << "\t" << Right_PD[0].size() << "\t" << Right_PD[1].size() << endl;
   short ReadLengthMinus = CurrentReadSeq.size() - 1;      
	int Sum;
   if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End) {
	   for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
		   if (Right_PD[i].size() == 1 && CurrentLength >= BP_Right_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH) 
					for (short j = 1; j <= ADDITIONAL_MISMATCH; j++) 
						Sum += Right_PD[i + j].size();
				
				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Right_PD[i][0];   
					TempOne.Direction = BACKWARD;
					TempOne.Strand = SENSE;
					TempOne.Mismatches = i;
					RightUP.push_back(TempOne);	
					break;
				}
			}
		}
   }
	
   if (CurrentLength < BP_Right_End) {
      vector <unsigned int> Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Right_PD_Output[CheckedIndex].reserve(Right_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[ReadLengthMinus - CurrentLength];
      unsigned int pos;
		
		for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
			int SizeOfCurrent = Right_PD[i].size();
			if (CurrentChar == 'N') {
            //#pragma omp parallel default(shared)
				{
					// #pragma omp for
				   for (int j = 0; j < SizeOfCurrent; j++) {
					   pos =  Right_PD[i][j] - 1;
					   if (Match2N[(short)TheInput[pos]] == 'N')
						   Right_PD_Output[i].push_back(pos);
				   }  
				}
			}
			else {
            //#pragma omp parallel default(shared)
				{
					// #pragma omp for
				   for (int j = 0; j < SizeOfCurrent; j++) {
					   pos =  Right_PD[i][j] - 1;
					   if (TheInput[pos] == CurrentChar)
						   Right_PD_Output[i].push_back(pos);
					   else Right_PD_Output[i + 1].push_back(pos); 
				   }
				}
			}			
		}
		
		int SizeOfCurrent = Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
		if (CurrentChar == 'N') {
			for (int j = 0; j < SizeOfCurrent; j++) {
				pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
				if (Match2N[(short)TheInput[pos]] == 'N')
					Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
			}         
		}
		else {
			for (int j = 0; j < SizeOfCurrent; j++) {
				pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
				if (TheInput[pos] == CurrentChar)
					Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				//else Left_PD_Output[i + 1].push_back(pos); 
			}         
		}	
		
		Sum = 0;
		for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			Sum += Right_PD_Output[i].size();
		}
      if (Sum) {
         short CurrentLength_output = CurrentLength + 1;
         CheckRight_Close(OneRead, TheInput, CurrentReadSeq, Right_PD_Output,
								  BP_Right_Start, BP_Right_End,
								  CurrentLength_output, RightUP);
      }
      else return;
   }
   else return;
}

void CheckLeft_Far(SPLIT_READ & OneRead,
						 const string & TheInput, 
						 const string & CurrentReadSeq,
						 const vector <unsigned int> Left_PD[],
						 const short & BP_Left_Start,
						 const short & BP_Left_End,
						 const short & CurrentLength,
						 vector <UniquePoint> & LeftUP) {
	//if (OneRead.MatchedRelPos > 160000 && LeftUP.size())
	//cout << "+ " << TheInput.size() << "\t" << CurrentLength << "\t" << CurrentReadSeq << "\t" << Left_PD[0].size() << "\t" << Left_PD[1].size() << "\t" << Left_PD[2].size() << "\t" << BP_Left_Start << "\t" << BP_Left_End << "\t" << LeftUP.size() << endl;
	
	if (CurrentLength == 20) OneRead.Found = true;
	int Sum;
	//cout << "1" << endl;
   if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End) {
      // put it to LeftUP if unique
		for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			if (Left_PD[i].size() == 1 && CurrentLength >= BP_Left_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH) 
		   		for (short j = 1; j <= ADDITIONAL_MISMATCH; j++) 
			   		Sum += Left_PD[i + j].size();
				
				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Left_PD[i][0]; 
					TempOne.Direction = FORWARD;
					TempOne.Strand = SENSE;
					TempOne.Mismatches = i;
					LeftUP.push_back(TempOne);	
					break;
				}
			}
		}
   }
	//cout << "2" << endl;
   if (CurrentLength < BP_Left_End) {
      vector <unsigned int> Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Left_PD_Output[CheckedIndex].reserve(Left_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[CurrentLength];
      //const int SizeOfCurrent = Left_PD.size();
      unsigned int pos;
		//if (TOTAL_SNP_ERROR_CHECKED_Minus) 
		{
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
				int SizeOfCurrent = Left_PD[i].size();
				if (CurrentChar == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Left_PD[i][j] + 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							Left_PD_Output[i].push_back(pos);
					}         
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Left_PD[i][j] + 1;
						if (TheInput[pos] == CurrentChar)
							Left_PD_Output[i].push_back(pos);
						else Left_PD_Output[i + 1].push_back(pos); 
					}         
				}			
			}
			
			int SizeOfCurrent = Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}         
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (TheInput[pos] == CurrentChar)
						Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
					//else Left_PD_Output[i + 1].push_back(pos); 
				}         
			}	
			
			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
				Sum += Left_PD_Output[i].size();
			}
			//cout << "Sum: " << Sum << endl;
			
			if (Sum) {
				const short CurrentLengthOutput = CurrentLength + 1;
				CheckLeft_Far(OneRead,TheInput, CurrentReadSeq, Left_PD_Output,
								  BP_Left_Start, BP_Left_End,
								  CurrentLengthOutput, LeftUP); 
				//cout << "end0" << endl;
				return;
			}
			else {
				//cout << "end1" << endl;
			   return;	
			}
		}

   }   
   else {
		//cout << "end2" << endl;
		return;
	}
}

void CheckRight_Far(SPLIT_READ & OneRead,
						  const string & TheInput,
						  const string & CurrentReadSeq,
						  const vector <unsigned int> Right_PD[],
						  const short & BP_Right_Start,
						  const short & BP_Right_End,
						  const short & CurrentLength,
						  vector <UniquePoint> & RightUP) {
	//if (OneRead.MatchedRelPos > 160000 && RightUP.size())
		//cout << "- " << TheInput.size() << "\t" << CurrentLength << "\t" << CurrentReadSeq << "\t" << Right_PD[0].size() << "\t" << Right_PD[1].size() << "\t" << Right_PD[2].size() << "\t" << BP_Right_Start << "\t" << BP_Right_End << "\t" << RightUP.size() << endl;
	if (CurrentLength == 20) OneRead.Found = true;
   short ReadLengthMinus = CurrentReadSeq.size() - 1;      
	int Sum;
   if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End) {
	   for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
		   if (Right_PD[i].size() == 1 && CurrentLength >= BP_Right_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH) 
					for (short j = 1; j <= ADDITIONAL_MISMATCH; j++) 
						Sum += Right_PD[i + j].size();
				
				if (Sum == 0 && i <= (short)(CurrentLength * Seq_Error_Rate + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					TempOne.AbsLoc = Right_PD[i][0];   
					TempOne.Direction = BACKWARD;
					TempOne.Strand = ANTISENSE;
					TempOne.Mismatches = i;
					RightUP.push_back(TempOne);	
					break;
				}
			}
		}
   }
	
   if (CurrentLength < BP_Right_End) {
      vector <unsigned int> Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			Right_PD_Output[CheckedIndex].reserve(Right_PD[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[ReadLengthMinus - CurrentLength];
      unsigned int pos;
		
		//if (TOTAL_SNP_ERROR_CHECKED_Minus) 
		{
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
				int SizeOfCurrent = Right_PD[i].size();
				if (CurrentChar == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Right_PD[i][j] - 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							Right_PD_Output[i].push_back(pos);
					}         
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  Right_PD[i][j] - 1;
						if (TheInput[pos] == CurrentChar)
							Right_PD_Output[i].push_back(pos);
						else Right_PD_Output[i + 1].push_back(pos); 
					}         
				}			
			}
			
			int SizeOfCurrent = Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}         
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (TheInput[pos] == CurrentChar)
						Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
					//else Left_PD_Output[i + 1].push_back(pos); 
				}         
			}	
			
			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
				Sum += Right_PD_Output[i].size();
			}
			if (Sum) {
				short CurrentLength_output = CurrentLength + 1;
				CheckRight_Far(OneRead, TheInput, CurrentReadSeq, Right_PD_Output,
									BP_Right_Start, BP_Right_End,
									CurrentLength_output, RightUP);
			}
			else return;
		}
		/*
		 else { // TOTAL_SNP_ERROR_CHECKED_Minus
		 
		 int SizeOfCurrent = Right_PD[TOTAL_SNP_ERROR_CHECKED_Minus].size();
		 if (CurrentChar == 'N') {
		 for (int j = 0; j < SizeOfCurrent; j++) {
		 pos =  Right_PD[TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
		 if (Match2N[(short)TheInput[pos]] == 'N')
		 Right_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
		 }         
		 }
		 else {
		 for (int j = 0; j < SizeOfCurrent; j++) {
		 pos =  Right_PD[TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
		 if (TheInput[pos] == CurrentChar)
		 Right_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
		 //else Left_PD_Output[i + 1].push_back(pos); 
		 }         
		 }	
		 
		 if (!Right_PD_Output[TOTAL_SNP_ERROR_CHECKED_Minus].empty()) {
		 short CurrentLength_output = CurrentLength + 1;
		 CheckRight_Far(TheInput, CurrentReadSeq, Right_PD_Output,
		 BP_Right_Start, BP_Right_End,
		 CurrentLength_output, RightUP);
		 }
		 else return;
		 }
		 */
   }
   else return;
}


short CompareTwoReads(const SPLIT_READ & First, const SPLIT_READ & Second) {
   //if (First.MatchedSeqID < Second.MatchedSeqID) return 0;
   //else if (First.MatchedSeqID > Second.MatchedSeqID) return 1;
   short FirstDISize = First.NT_str.size();
   short SecondDISize = Second.NT_str.size();
   short FirstSISize = First.InsertedStr.size();
   short SecondSISize = Second.InsertedStr.size();
	
   if (First.BPLeft > Second.BPLeft) return 1;
   else if (First.BPLeft < Second.BPLeft) return 0;
   else {
      {
         if (First.IndelSize > Second.IndelSize) return 1; 
         else if (First.IndelSize < Second.IndelSize) return 0;
         else
         {
            if (First.Tag.size() < Second.Tag.size()) return 0;
            else if (First.Tag.size() > Second.Tag.size()) return 1;
            else if (First.Tag.size() == Second.Tag.size()) {
               for (unsigned posindex = 0; posindex < First.Tag.size(); posindex++) {
                  if ((short)First.Tag[posindex] < (short)Second.Tag[posindex]) return 0;
                  else if ((short)First.Tag[posindex] > (short)Second.Tag[posindex]) return 1; 
               }
               if (First.MatchedRelPos < Second.MatchedRelPos)  return 0;
               else if (First.MatchedRelPos > Second.MatchedRelPos) return 1;
               else {
                  if (FirstDISize > SecondDISize) return 1;
						else if (FirstDISize < SecondDISize) return 0;
						else if (FirstDISize != 0) {
                     for (int i = 0; i < FirstDISize; i++) 
                        if ((short)First.NT_str[i] > (short)Second.NT_str[i]) return 1;
                        else if ((short)First.NT_str[i] < (short)Second.NT_str[i]) return 0;
                  }
                  else {
                     if (FirstSISize > SecondSISize) return 1;
                     else if (FirstSISize < SecondSISize) return 0;
                     else if (FirstSISize != 0) {
                        for (int i = 0; i < FirstSISize; i++) 
                           if ((short)First.InsertedStr[i] > (short)Second.InsertedStr[i]) return 1;
                           else if ((short)First.InsertedStr[i] < (short)Second.InsertedStr[i]) return 0;
                     } 
                     else return 2; 
                  } 
               }
               //return 2;                     
            }   
         }    
      } 
   } 
	return 0;
}

string Cap2Low(const string & input) {
   string output = input;
   for (unsigned int i = 0; i < output.size(); i++) {
      output[i] = Cap2LowArray[(short)input[i]];   
   }
   return output;
}

bool CheckMismatches(const string & TheInput, 
                     const string & InputReadSeq, 
                     //const unsigned int & Start,
							const UniquePoint & UP
							) {
	//return true; 	short LengthStr;
	//unsigned int AbsLoc; 
	//cout << "CheckMismatches1" << endl;
	string CurrentReadSeq;
	if (UP.Strand == SENSE) CurrentReadSeq = InputReadSeq;
	else CurrentReadSeq = ReverseComplement(InputReadSeq);
	short CurrentReadLength = CurrentReadSeq.size();
	unsigned int Start = 0;
	//cout << "CheckMismatches2" << endl;
	string BP_On_Read, BP_On_Ref;
	if (UP.Direction == FORWARD) {
		//cout << "+s" << endl;
		//cout << "UP.AbsLoc: " << UP.AbsLoc << "\t" << "UP.LengthStr: " << UP.LengthStr << endl;
		
		Start = UP.AbsLoc - UP.LengthStr + 1; 
		//cout << "Start: " << Start << endl;
		//cout << "Min_Perfect_Match_Around_BP: " << Min_Perfect_Match_Around_BP << endl;
		if (UP.LengthStr <= Min_Perfect_Match_Around_BP) return false;
		BP_On_Read = CurrentReadSeq.substr(UP.LengthStr - Min_Perfect_Match_Around_BP, Min_Perfect_Match_Around_BP);
		//cout << "BP_On_Read: " << BP_On_Read << endl;
		//if (UP.AbsLoc < Min_Perfect_Match_Around_BP) return false;
		BP_On_Ref = TheInput.substr(UP.AbsLoc - Min_Perfect_Match_Around_BP + 1, Min_Perfect_Match_Around_BP);
		//cout << "BP_On_Ref: " << BP_On_Ref << endl;
		if (BP_On_Read != BP_On_Ref) return false;
		//cout << "+e" << endl;
	}
	else if (UP.Direction == BACKWARD) {
		//cout << "-s" << endl;
		Start = UP.AbsLoc + UP.LengthStr - CurrentReadLength;
		if (CurrentReadLength < UP.LengthStr) return false;
		BP_On_Read = CurrentReadSeq.substr(CurrentReadLength - UP.LengthStr, Min_Perfect_Match_Around_BP);
		BP_On_Ref = TheInput.substr(UP.AbsLoc, Min_Perfect_Match_Around_BP);
		if (BP_On_Read != BP_On_Ref) return false;
		//cout << "-e" << endl;
	}
	//cout << "CheckMismatches3" << endl;
	short MAX_ALLOWED_MISMATCHES = (short)(CurrentReadSeq.size() * MaximumAllowedMismatchRate + 1); // 
	
   short NumMismatches = 0;                 // Match2N[(short)'A'] = 'N';    
   
   for (short i = 0; i < CurrentReadLength; i++) {
      if (CurrentReadSeq[i] == N_char) {
         if (Match2N[(short)TheInput[Start + i]] != N_char)  
            NumMismatches++;
      }
      else {
         if (TheInput[Start + i] != CurrentReadSeq[i])
            NumMismatches++;
      }   
   }   
	//cout << "CheckMismatches4" << endl;
   if (NumMismatches > MAX_ALLOWED_MISMATCHES) return true;
   else return false;                     
}

void OutputTDs(const vector <SPLIT_READ> & TDs,
					const string & TheInput, 
					const unsigned int & C_S,
					const unsigned int & C_E,
					const unsigned int & RealStart,
					const unsigned int & RealEnd,
					ofstream & TDOutf) {
   //short ReadLength = Deletions[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (TDs[i].Tag == VectorTag[j]) {
				if (TDs[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (TDs[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (TDs[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (TDs[i].MatchedD == Plus) {
         LeftScore += TDs[i].score;
         LeftS++;
			if (TDs[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += TDs[i].score;
         RightS++;
			if (TDs[i].Unique)
			   RightUNum++;
      }
   }
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples
	
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize =0;
   if (TDs[C_S].IndelSize < 14) GapSize = TDs[C_S].IndelSize;
   else GapSize = 13 + (int)log10(TDs[C_S].IndelSize - 10);
	CurrentChrMask[TDs[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[TDs[C_S].BPRight + SpacerBeforeAfter] = 'B';
	
   TDOutf << "####################################################################################################" << endl;
   TDOutf << NumberOfTDInstances << "\tTD " << TDs[C_S].IndelSize// << " bases " 
	<< "\tNT " << TDs[C_S].NT_size << " \"" << TDs[C_S].NT_str << "\""
	<< "\tChrID " << TDs[C_S].FragName 
	<< "\tBP " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2
	<< "\tBP_range " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum
	<< "\tS1 " << EasyScore ; //EWL070111 << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += TDs[i].MS;
   TDOutf << "\tSUM_MS " << SUM_MS;
	
   TDOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		TDOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   TDOutf << endl;
	
   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
   TDOutf << TheInput.substr(TDs[C_S].BPRight + SpacerBeforeAfter - ReportLength + 1, ReportLength);// << endl;// 
	if (TDs[C_S].NT_size) {
		for (short i = 0; i < TDs[C_S].NT_size; i++) TDOutf << " ";
	}
	//ReportLength 
   TDOutf << Cap2Low(TheInput.substr(TDs[C_S].BPLeft + SpacerBeforeAfter, ReportLength)) << endl;
	
   short SpaceBeforeReadSeq;    
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      SpaceBeforeReadSeq = ReportLength - TDs[GoodIndex].BP - 1;
      
      for (int i = 0; i < SpaceBeforeReadSeq; i++) TDOutf << " ";
      short SpaceBeforeD = ReportLength + ReportLength - SpacerBeforeAfter - TDs[GoodIndex].ReadLength;
      if (TDs[GoodIndex].MatchedD == Minus) {
         TDOutf << TDs[GoodIndex].UnmatchedSeq << endl;
         //for (int i = 0; i < GapSize; i++) TDOutf << " ";    
         //TDOutf << TDs[GoodIndex].UnmatchedSeq.substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
      }
      else {
         TDOutf << ReverseComplement(TDs[GoodIndex].UnmatchedSeq) << endl; 
         //for (int i = 0; i < GapSize; i++) TDOutf << " ";  
         //TDOutf << ReverseComplement(TDs[GoodIndex].UnmatchedSeq).substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
      } 
      //for (int i = 0; i < SpaceBeforeD; i++) TDOutf << " ";
      TDOutf << "\t" << TDs[GoodIndex].MatchedD << "\t" 
		<< TDs[GoodIndex].MatchedRelPos 
		<< "\t" << TDs[GoodIndex].MS
		<< "\t" << TDs[GoodIndex].Tag 
		<< "\t" <<  TDs[GoodIndex].Name << endl;
		//<< "\t" << TDs[C_S].BPLeft
		//<< "\t" << TDs[C_S].BPRight << endl; 
   }     
}

void OutputDeletions(const vector <SPLIT_READ> & Deletions,
                     const string & TheInput, 
                     const unsigned int & C_S,
                     const unsigned int & C_E,
							const unsigned int & RealStart,
							const unsigned int & RealEnd,
                     ofstream & DeletionOutf) {
   //short ReadLength = Deletions[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
	//cout << "d_1" << endl;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	//cout << "d_2" << endl;
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	//cout << "d_3" << endl;
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (Deletions[i].Tag == VectorTag[j]) {
				if (Deletions[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (Deletions[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (Deletions[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (Deletions[i].MatchedD == Plus) {
         LeftScore += Deletions[i].score;
         LeftS++;
			if (Deletions[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += Deletions[i].score;
         RightS++;
			if (Deletions[i].Unique)
			   RightUNum++;
      }
   }
	//cout << "d_4" << endl;
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//cout << "d_5" << endl;
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize = 0;
   if (Deletions[C_S].IndelSize < 14) GapSize = Deletions[C_S].IndelSize;
   else GapSize = 13 + (int)log10(Deletions[C_S].IndelSize - 10);
	//cout << "d_5a" << endl;
	CurrentChrMask[Deletions[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	//cout << "d_5b" << endl;
	CurrentChrMask[Deletions[C_S].BPRight + SpacerBeforeAfter] = 'B';
	//cout << "d_5c" << endl;
	CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	//cout << "d_5d" << endl;
	CurrentChrMask[RealEnd + SpacerBeforeAfter] = 'B';
	//cout << "d_5e" << endl;
   DeletionOutf << "####################################################################################################" << endl;
   DeletionOutf << NumberOfDeletionsInstances << "\tD " << Deletions[C_S].IndelSize// << " bases " 
	<< "\tNT " << Deletions[C_S].NT_size << " \"" << Deletions[C_S].NT_str << "\""
	<< "\tChrID " << Deletions[C_S].FragName 
	<< "\tBP " << Deletions[C_S].BPLeft + 1 << "\t" << Deletions[C_S].BPRight + 1 
	<< "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1  
	<< "\tSupports " << NumberOfReads << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1 << "\t" << LeftUNum
	<< "\t- " << RightS - 1 << "\t" << RightUNum 
	<< "\tS1 " << EasyScore; //EWL070111  << "\tS2 " << PreciseScore; 
	//cout << "d_6" << endl;
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += Deletions[i].MS;
   DeletionOutf << "\tSUM_MS " << SUM_MS;
	
   DeletionOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   DeletionOutf << endl;
	//cout << "d_7" << endl;
   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
   DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength);// << endl;// ReportLength    
   if (Deletions[C_S].IndelSize >= 14) {
      DeletionOutf << Cap2Low(TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1, 5))
		<< "<" << Deletions[C_S].IndelSize - 10 << ">" 
		<< Cap2Low(TheInput.substr(Deletions[C_S].Right - Deletions[C_S].ReadLength + Deletions[C_S].BP - 3, 5));  
   }   
   else DeletionOutf << Cap2Low(TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1, GapSize));             
   DeletionOutf << TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1 + Deletions[C_S].IndelSize, ReportLength - GapSize) << endl;// ReportLength
   short SpaceBeforeReadSeq;    
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      SpaceBeforeReadSeq = ReportLength - Deletions[GoodIndex].BP - 1;
      
      for (int i = 0; i < SpaceBeforeReadSeq; i++) DeletionOutf << " ";
      short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - Deletions[GoodIndex].ReadLength;
      if (Deletions[GoodIndex].MatchedD == Minus) {
         DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(0, Deletions[GoodIndex].BP + 1);// << endl;
         for (int i = 0; i < GapSize; i++) DeletionOutf << " ";    
         DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      else {
         DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(0, Deletions[GoodIndex].BP + 1);// << endl; 
         for (int i = 0; i < GapSize; i++) DeletionOutf << " ";  
         DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      } 
      for (int i = 0; i < SpaceBeforeD; i++) DeletionOutf << " ";
      DeletionOutf << "\t" << Deletions[GoodIndex].MatchedD << "\t" 
		<< Deletions[GoodIndex].MatchedRelPos 
		<< "\t" << Deletions[GoodIndex].MS
		<< "\t" << Deletions[GoodIndex].Tag 
		<< "\t" <<  Deletions[GoodIndex].Name << endl;
		//<< "\t" << Deletions[C_S].BPLeft
		//<< "\t" << Deletions[C_S].BPRight << endl; 
   }     
}

void OutputInversions(const vector <SPLIT_READ> & Inv, 
							 const string & TheInput, 
							 const unsigned int & C_S, 
							 const unsigned int & C_E,
							 const unsigned int & RealStart,
							 const unsigned int & RealEnd,
							 ofstream & InvOutf) {
	int LeftNT_index = -1;
	int RightNT_index = -1;
	for (unsigned Index = C_S; Index <= C_E; Index++) {
		if (Inv[Index].MatchedD == Plus) {
			LeftNT_index = Index;
			break;
		}
	}
	for (unsigned Index = C_S; Index <= C_E; Index++) {
		if (Inv[Index].MatchedD == Minus) {
			RightNT_index = Index;
			break;
		}
	}
	short LeftNT_size = 0;
	short RightNT_size = 0;
	string LeftNT_str = "";
	string RightNT_str = "";
	if (LeftNT_index != -1) {
		LeftNT_size = Inv[LeftNT_index].NT_size;
		LeftNT_str = Inv[LeftNT_index].NT_str;
	}
	if (RightNT_index != -1) {
		RightNT_size = Inv[RightNT_index].NT_size;
		RightNT_str = Inv[RightNT_index].NT_str;
	}
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (Inv[i].Tag == VectorTag[j]) {
				if (Inv[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (Inv[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (Inv[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (Inv[i].MatchedD == Plus) {
         LeftScore += Inv[i].score;
         LeftS++;
			if (Inv[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += Inv[i].score;
         RightS++;
			if (Inv[i].Unique)
			   RightUNum++;
      }
   }
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples
	
   unsigned int EasyScore = LeftS * RightS;
   //double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize =0;
   if (Inv[C_S].IndelSize < 14) GapSize = Inv[C_S].IndelSize;
   else GapSize = 13 + (int)log10(Inv[C_S].IndelSize - 10);
	CurrentChrMask[Inv[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[Inv[C_S].BPRight + SpacerBeforeAfter] = 'B';
   InvOutf << "####################################################################################################" << endl;
   InvOutf << NumberOfInvInstances << "\tINV " << Inv[C_S].IndelSize// << " bases " 
	<< "\tNT " <<  LeftNT_size << ":" << RightNT_size << " \"" << LeftNT_str << "\":\"" << RightNT_str << "\""
	<< "\tChrID " << Inv[C_S].FragName 
	<< "\tBP " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 
	<< "\tBP_range " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum
	<< "\tS1 " << EasyScore; //EWL070111  << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += Inv[i].MS;
   InvOutf << "\tSUM_MS " << SUM_MS;
	
   InvOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		InvOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   InvOutf << endl;
	
	short SpaceBeforeReadSeq;  
   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
   InvOutf << TheInput.substr(Inv[C_S].BPLeft + SpacerBeforeAfter - ReportLength, ReportLength);//;// << endl;// ReportLength  
	//cout << Inv[C_S].NT_size << "\t" << Inv[C_S].NT_2size << endl;
	if (LeftNT_size) {
		for (int i = 0; i < LeftNT_size; i++) {
			InvOutf << " ";
		}
	}
	InvOutf << Cap2Low( ReverseComplement(TheInput.substr(Inv[C_S].BPRight + 1 + SpacerBeforeAfter - ReportLength, ReportLength))) << endl; 
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      if (Inv[GoodIndex].MatchedD == Plus && Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPLeft) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << ReverseComplement(Inv[GoodIndex].UnmatchedSeq);
			for (int i = 0; i < Inv[GoodIndex].BP; i++) InvOutf << " ";
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
			//<< "\t" << Deletions[C_S].BPLeft
			//<< "\t" << Deletions[C_S].BPRight << endl; 
      }
      else if (Inv[GoodIndex].MatchedD == Plus && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPLeft) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << Inv[GoodIndex].UnmatchedSeq;
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
      } 
   }  
	
	InvOutf << "----------------------------------------------------------------------------------------------------" << endl;
	
	
	InvOutf << Cap2Low( ReverseComplement(TheInput.substr(Inv[C_S].BPLeft + SpacerBeforeAfter, ReportLength)));
	if (RightNT_size) {
		for (int i = 0; i < RightNT_size; i++) {
			InvOutf << " ";
		}
	}
	InvOutf << TheInput.substr(Inv[C_S].BPRight + 1 + SpacerBeforeAfter, ReportLength) << endl; 
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      if (Inv[GoodIndex].MatchedD == Minus && Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPRight) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << ReverseComplement(Inv[GoodIndex].UnmatchedSeq);
			for (int i = 0; i < Inv[GoodIndex].BP; i++) InvOutf << " ";
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str 
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
			//<< "\t" << Deletions[C_S].BPLeft
			//<< "\t" << Deletions[C_S].BPRight << endl; 
      }
      else if (Inv[GoodIndex].MatchedD == Minus && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPRight) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << Inv[GoodIndex].UnmatchedSeq;
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
      } 
   } 
	
	
}

void OutputSIs(const vector <SPLIT_READ> & SIs,
               const string & TheInput,
               const unsigned int & C_S,
               const unsigned int & C_E,
					const unsigned int & RealStart,
					const unsigned int & RealEnd,
               ofstream & SIsOutf) {
   //short ReadLength = SIs[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;               
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (SIs[i].Tag == VectorTag[j]) {
				if (SIs[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (SIs[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (SIs[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (SIs[i].MatchedD == Plus) {
         LeftScore += SIs[i].score;
         LeftS++;
			if (SIs[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += SIs[i].score;
         RightS++;
			if (SIs[i].Unique)
			   RightUNum++;
      }
   }
	
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples
	
   unsigned int EasyScore = LeftS * RightS;
   //double PreciseScore = (LeftScore + RightScore) * (-1);
   string CurrentReadSeq;
	
	CurrentChrMask[SIs[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[SIs[C_S].BPRight + SpacerBeforeAfter] = 'B';
	CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	CurrentChrMask[RealEnd + SpacerBeforeAfter] = 'B';
	
   SIsOutf << "####################################################################################################" << endl;
   SIsOutf << NumberOfSIsInstances << "\tI " << SIs[C_S].IndelSize << "\tNT " << SIs[C_S].IndelSize << " \"" << SIs[C_S].InsertedStr << "\""
	<< "\tChrID " << SIs[C_S].FragName 
	<< "\tBP " << SIs[C_S].BPLeft + 1 << "\t" << SIs[C_S].BPRight + 1 
	<< "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum
	<< "\tS1 " << EasyScore; //EWL070111  << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += SIs[i].MS;
   SIsOutf << "\tSUM_MS " << SUM_MS;
	
   SIsOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples  << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		SIsOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   SIsOutf << endl;
	
   SIsOutf << TheInput.substr(SIs[C_S].Left - ReportLength + SIs[C_S].BP + 1, ReportLength);// ReportLength
   for (int i = 0; i < SIs[C_S].IndelSize; i++) SIsOutf << " ";
   SIsOutf << TheInput.substr(SIs[C_S].Left + SIs[C_S].BP + 1, ReportLength) << endl;// ReportLength
   
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      short SpaceBeforeReadSeq = ReportLength - SIs[GoodIndex].BP - 1;
      for (short i = 0; i < SpaceBeforeReadSeq; i++) SIsOutf << " "; 
      if (SIs[GoodIndex].MatchedD == Minus)
         SIsOutf << SIs[GoodIndex].UnmatchedSeq;
      else SIsOutf << ReverseComplement(SIs[GoodIndex].UnmatchedSeq);
      short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - SIs[GoodIndex].ReadLength;
      for (short i = 0; i < SpaceBeforeD; i++) SIsOutf << " "; 
      SIsOutf << "\t" << SIs[GoodIndex].MatchedD 
		<< "\t" << SIs[GoodIndex].MatchedRelPos 
		<< "\t" << SIs[GoodIndex].MS 
		<< "\t" << SIs[GoodIndex].Tag << "\t" << SIs[GoodIndex].Name << endl;  
   }
}

void OutputDI(const vector <SPLIT_READ> & DI,
				  const string & TheInput, 
				  const unsigned int & C_S,
				  const unsigned int & C_E,
				  const unsigned int & RealStart,
				  const unsigned int & RealEnd,
				  ofstream & DeletionOutf) {
   //short ReadLength = DI[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (DI[i].Tag == VectorTag[j]) {
				if (DI[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (DI[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (DI[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (DI[i].MatchedD == Plus) {
         LeftScore += DI[i].score;
         LeftS++;
			if (DI[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += DI[i].score;
         RightS++;
			if (DI[i].Unique)
			   RightUNum++;
      }
   }
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	
   unsigned int EasyScore = LeftS * RightS;
   //double PreciseScore = (LeftScore + RightScore) * (-1);
   //short GapSize =0;
   //if (DI[C_S].IndelSize < 14) GapSize = DI[C_S].IndelSize;
   //else GapSize = 13 + (int)log10(Deletions[C_S].IndelSize - 10);
	CurrentChrMask[DI[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[DI[C_S].BPRight + SpacerBeforeAfter] = 'B';
	//CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	//CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B'; 
   DeletionOutf << "####################################################################################################" << endl;
   DeletionOutf << NumberOfDIInstances + NumberOfDeletionsInstances + 1 << "\tD " << DI[C_S].IndelSize 
	<< "\tNT " << DI[C_S].NT_size << " \"" << DI[C_S].NT_str << "\"" 
	<< "\tChrID " << DI[C_S].FragName 
	<< "\tBP " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 
	<< "\tBP_range " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum  
	<< "\tS1 " << EasyScore; // << "\tS2 0.0";// << PreciseScore << "\t"; 
	//EWL070111 << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += DI[i].MS;
   DeletionOutf << "\tSUM_MS " << SUM_MS;
	
	
   DeletionOutf << "\t"<< VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   DeletionOutf << endl;
	
   //DeletionOutf << TheInput.substr(DI[C_S].Left - ReportLength + DI[C_S].BP + 1, 2 * ReportLength) << endl;       
   DeletionOutf << TheInput.substr(DI[C_S].Left - ReportLength + DI[C_S].BP + 1, ReportLength);// << endl;// ReportLength    
	
   for (short i = 0; i < DI[C_S].NT_size; i++) DeletionOutf << " ";           
   DeletionOutf << TheInput.substr(DI[C_S].Left + DI[C_S].BP + 1 + DI[C_S].IndelSize, ReportLength) << endl;// ReportLength
   short SpaceBeforeReadSeq;    
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      SpaceBeforeReadSeq = ReportLength - DI[GoodIndex].BP - 1;
      
      for (int i = 0; i < SpaceBeforeReadSeq; i++) DeletionOutf << " ";
      //short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - Deletions[GoodIndex].ReadLength;
      if (DI[GoodIndex].MatchedD == Minus) {
         DeletionOutf << DI[GoodIndex].UnmatchedSeq << "\t";
         //for (int i = 0; i < GapSize; i++) DeletionOutf << " ";    
         //DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      else {
         DeletionOutf << ReverseComplement(DI[GoodIndex].UnmatchedSeq) << "\t"; 
         //for (int i = 0; i < GapSize; i++) DeletionOutf << " ";  
         //DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      } 
      //for (int i = 0; i < SpaceBeforeD; i++) DeletionOutf << " ";
      DeletionOutf << "\t" << DI[GoodIndex].MatchedD << "\t" << DI[GoodIndex].MatchedRelPos
		<< "\t" << DI[GoodIndex].MS
		<< "\t" << DI[GoodIndex].Tag << "\t" <<  DI[GoodIndex].Name << endl;
   }     
}


void SortOutputSI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> SIs[], ofstream & SIsOutf) {
   cout << "Sorting and outputing short insertions ..." << endl;
   unsigned int SIsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
	//vector <SPLIT_READ> InputIndels; 
   vector <SPLIT_READ> GoodIndels;
   unsigned int GoodNum;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (SIs[Box_index].size() >= NumRead2ReportCutOff) {
			//InputIndels.clear();
			SIsNum = SIs[Box_index].size();
			//cout << "SIsNum " << SIsNum << endl;
			//for (int i = 0; i < SIsNum; i++) {
			//	InputIndels.push_back(Reads[SIs[Box_index][i]]);
			//}
			for (unsigned int First = 0; First < SIsNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < SIsNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{  
							if (Reads[SIs[Box_index][First]].ReadLength == Reads[SIs[Box_index][Second]].ReadLength) {
								if (Reads[SIs[Box_index][First]].LeftMostPos == Reads[SIs[Box_index][Second]].LeftMostPos)
									Reads[SIs[Box_index][Second]].Unique = false;
							}
							
							if (Reads[SIs[Box_index][First]].BPLeft < Reads[SIs[Box_index][Second]].BPLeft) continue;
							else if (Reads[SIs[Box_index][First]].BPLeft > Reads[SIs[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else {
								if (Reads[SIs[Box_index][First]].IndelSize < Reads[SIs[Box_index][Second]].IndelSize) continue;
								else if (Reads[SIs[Box_index][First]].IndelSize > Reads[SIs[Box_index][Second]].IndelSize) {
									CompareResult = 1;
								}
								//else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
								//		short Compare2Str = CompareTwoString(InputIndels[First].InsertedStr, InputIndels[Second].InsertedStr );
								//if (Compare2Str > 0) CompareResult = 1;
								//else if (Compare2Str == 0) CompareResult = 2;
								//else continue;
								
								//}
								//else if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//	if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//		InputIndels[Second].Unique = false;
								//	}
								//}
							}
                     //CompareResult = CompareTwoReads(InputIndels[First], InputIndels[Second]);
                     if (CompareResult == 1) {
  	                     Temp4Exchange = SIs[Box_index][First];
  	                     SIs[Box_index][First] = SIs[Box_index][Second];
  	                     SIs[Box_index][Second] = Temp4Exchange;
                     }
                     //else if (CompareResult == 2) InputIndels[Second].Unique = false;
                  }
               }
            }
         }
         GoodIndels.clear();
	      IndelEvents.clear();
			//cout << GoodIndels.size() << endl;
			
			for (unsigned int First = 0; First < SIsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[SIs[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
			OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
         //OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].IndelSize == OneIndelEvent.IndelSize)
					//&& OneIndelEvent.IndelStr == GoodIndels[GoodIndex].InsertedStr )
					
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Insertion(CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
					OneIndelEvent.IndelSize =  GoodIndels[GoodIndex].IndelSize;
					OneIndelEvent.IndelStr =  GoodIndels[GoodIndex].InsertedStr;
					//OneIndelEvent.IndelStr = GoodIndels[GoodIndex].InsertedStr;
            }
         }
			
			//if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff) 
			OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
			OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
			
			OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			GetRealStart4Insertion(CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
			IndelEvents.push_back(OneIndelEvent);
			
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			unsigned int IndelSize;
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  IndelSize = IndelEvents[EventIndex].IndelSize;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
						
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else if (IndelEvents[EventIndex_left].IndelSize != IndelSize) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
                  //cout << Max_Support << endl;
                  if (Max_Support >= NumRead2ReportCutOff) {
							OutputSIs(GoodIndels, CurrentChr,
										 IndelEvents[Max_Support_Index].Start,
										 IndelEvents[Max_Support_Index].End,
										 RealStart, RealEnd, SIsOutf);
							NumberOfSIsInstances++;
						}
               }
            }
         }
      }   // if (!insertion[Box_index].empty())
   }
   cout << "Short insertions: " << NumberOfSIsInstances << endl << endl;
}

void SortOutputTD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf) {
   cout << "Sorting and outputing tandem duplications ..." << endl;
   unsigned int TDNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (TDs[Box_index].size() >= NumRead2ReportCutOff) {
			//InputIndels.clear();
			TDNum = TDs[Box_index].size();
			//for (int i = 0 ; i < TDNum; i ++) {
			//	InputIndels.push_back(AllReads[TDs[Box_index][i]]);
			//}
         
         for (unsigned int First = 0; First < TDNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < TDNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (AllReads[TDs[Box_index][First]].ReadLength == AllReads[TDs[Box_index][Second]].ReadLength) {
								if (AllReads[TDs[Box_index][First]].LeftMostPos == AllReads[TDs[Box_index][Second]].LeftMostPos)
									AllReads[TDs[Box_index][Second]].Unique = false;
							}
							if (AllReads[TDs[Box_index][First]].BPLeft < AllReads[TDs[Box_index][Second]].BPLeft) continue;
                     else if (AllReads[TDs[Box_index][First]].BPLeft > AllReads[TDs[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (AllReads[TDs[Box_index][First]].BPLeft == AllReads[TDs[Box_index][Second]].BPLeft) {
								if (AllReads[TDs[Box_index][First]].BPRight < AllReads[TDs[Box_index][Second]].BPRight) continue;
								else if (AllReads[TDs[Box_index][First]].BPRight > AllReads[TDs[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//		
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = TDs[Box_index][First];
  	                     TDs[Box_index][First] = TDs[Box_index][Second];
  	                     TDs[Box_index][Second] = Temp4Exchange;
							}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < TDNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(AllReads[TDs[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size(); 
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Tandem duplications: " << NumberOfTDInstances << endl << endl;
}

void SortOutputTD_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf) {
   cout << "Sorting and outputing tandem duplications with non-template sequence ..." << endl;
   unsigned int TDNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
	int Count_TD_NT_output = 0;
	
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (TDs[Box_index].size() >= NumRead2ReportCutOff) {
			//InputIndels.clear();
			TDNum = TDs[Box_index].size();
			//for (int i = 0 ; i < TDNum; i ++) {
			//	InputIndels.push_back(AllReads[TDs[Box_index][i]]);
			//}
         
         for (unsigned int First = 0; First < TDNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < TDNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (AllReads[TDs[Box_index][First]].ReadLength == AllReads[TDs[Box_index][Second]].ReadLength) {
								if (AllReads[TDs[Box_index][First]].LeftMostPos == AllReads[TDs[Box_index][Second]].LeftMostPos)
									AllReads[TDs[Box_index][Second]].Unique = false;
							}
							if (AllReads[TDs[Box_index][First]].BPLeft < AllReads[TDs[Box_index][Second]].BPLeft) continue;
                     else if (AllReads[TDs[Box_index][First]].BPLeft > AllReads[TDs[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (AllReads[TDs[Box_index][First]].BPLeft == AllReads[TDs[Box_index][Second]].BPLeft) {
								if (AllReads[TDs[Box_index][First]].BPRight < AllReads[TDs[Box_index][Second]].BPRight) continue;
								else if (AllReads[TDs[Box_index][First]].BPRight > AllReads[TDs[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
									if (AllReads[TDs[Box_index][First]].NT_size < AllReads[TDs[Box_index][Second]].NT_size) continue;
									else if (AllReads[TDs[Box_index][First]].NT_size > AllReads[TDs[Box_index][Second]].NT_size) CompareResult = 1;
									//else { // InputIndels[First].NT_size == InputIndels[Second].NT_size
									//	short Compare2Str = CompareTwoString(InputIndels[First].NT_str, InputIndels[Second].NT_str );
									//	if (Compare2Str > 0) CompareResult = 1;
									//	else if (Compare2Str == 0) CompareResult = 2;
									//	else continue;
									//}
								}
							}
							if (CompareResult == 1) {
								Temp4Exchange = TDs[Box_index][First];
  	                     TDs[Box_index][First] = TDs[Box_index][Second];
  	                     TDs[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
							//  InputIndels[First + 1] = InputIndels[Second];
							// InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < TDNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(AllReads[TDs[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size(); 
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
								Count_TD_NT_output++;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
								Count_TD_NT_output++;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Tandem duplications with non-template sequence (TD_NT): " << Count_TD_NT_output << "\n\n";
}

void SortOutputD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Deletions[], ofstream & DeletionOutf) {
   cout << "Sorting and outputing deletions ..." << endl;
   unsigned int DeletionsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
		//cout << Box_index << "\t" << NumBoxes << "\t" << Deletions[Box_index].size() << endl;
		//if (Deletions[Box_index].size() >= NumRead2ReportCutOff)
		//   cout << Box_index << "\t" << Deletions[Box_index].size() << endl;
      if (Deletions[Box_index].size() >= NumRead2ReportCutOff) {
			DeletionsNum = Deletions[Box_index].size();
         
         for (unsigned int First = 0; First < DeletionsNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				//cout << DeletionsNum << " First " << First << endl;
				{
               for (unsigned int Second = First + 1; Second < DeletionsNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						//cout << DeletionsNum << " Second " << Second << endl;
						{
							if (Reads[Deletions[Box_index][First]].ReadLength == Reads[Deletions[Box_index][Second]].ReadLength) {
								if (Reads[Deletions[Box_index][First]].LeftMostPos == Reads[Deletions[Box_index][Second]].LeftMostPos)
									Reads[Deletions[Box_index][Second]].Unique = false;
							}
							if (Reads[Deletions[Box_index][First]].BPLeft < Reads[Deletions[Box_index][Second]].BPLeft) continue;
                     else if (Reads[Deletions[Box_index][First]].BPLeft > Reads[Deletions[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[Deletions[Box_index][First]].BPLeft == Reads[Deletions[Box_index][Second]].BPLeft) {
								if (Reads[Deletions[Box_index][First]].BPRight < Reads[Deletions[Box_index][Second]].BPRight) continue;
								else if (Reads[Deletions[Box_index][First]].BPRight > Reads[Deletions[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								//else CompareResult = 2;
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//			
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = Deletions[Box_index][First];
  	                     Deletions[Box_index][First] = Deletions[Box_index][Second];
  	                     Deletions[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			//for (int i = 0 ; i < DeletionsNum; i ++) {
			//	InputIndels.push_back(Reads[Deletions[Box_index][i]]);
			//}Reads[Deletions[Box_index][First]]
         //cout << "GoodIndels" << endl;
			for (unsigned int First = 0; First < DeletionsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[Deletions[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();
			//cout << Box_index << " box read size " << GoodNum << endl;
         if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
			//cout << "here" << endl;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
				     //cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					//cout << EventIndex << " EventIndex" << endl;
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
						//cout << "max" << endl;
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							//cout << "aa" << endl;
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								//cout << "ba" << endl;
								OutputDeletions(GoodIndels, CurrentChr, 
													 IndelEvents[Max_Support_Index].Start, 
													 IndelEvents[Max_Support_Index].End, 
													 RealStart, RealEnd, DeletionOutf);
								NumberOfDeletionsInstances++;
								//cout << "bb" << endl;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								//cout << "ca" << endl;
								OutputDeletions(GoodIndels, CurrentChr, 
													 IndelEvents[Max_Support_Index].Start, 
													 IndelEvents[Max_Support_Index].End, 
													 RealStart, RealEnd, DeletionOutf);
								NumberOfDeletionsInstances++;
								//cout << "cb" << endl;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Deletions: " << NumberOfDeletionsInstances << endl << endl;
}

void SortOutputInv(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Inv[], ofstream & InvOutf) {
   cout << "Sorting and outputing Inversions ..." << endl;
   unsigned int InversionsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	/*
	 unsigned int C_S = 0;
	 unsigned int C_E = 0;
	 unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	 unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	 unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
		
      if (Inv[Box_index].size() >= NumRead2ReportCutOff) {
			//cout << Box_index << "\t" << Inv[Box_index].size() << endl;
         InversionsNum = Inv[Box_index].size();
			//InputIndels.clear();
			//for (int i = 0; i < InversionsNum; i++) InputIndels.push_back(Reads[Inv[Box_index][i]]);
         for (unsigned int First = 0; First < InversionsNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < InversionsNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (Reads[Inv[Box_index][First]].ReadLength == Reads[Inv[Box_index][Second]].ReadLength) {
								if (Reads[Inv[Box_index][First]].LeftMostPos == Reads[Inv[Box_index][Second]].LeftMostPos)
									Reads[Inv[Box_index][Second]].Unique = false;
							}
							if (Reads[Inv[Box_index][First]].BPLeft < Reads[Inv[Box_index][Second]].BPLeft) continue;
                     else if (Reads[Inv[Box_index][First]].BPLeft > Reads[Inv[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[Inv[Box_index][First]].BPLeft == Reads[Inv[Box_index][Second]].BPLeft) {
								if (Reads[Inv[Box_index][First]].BPRight < Reads[Inv[Box_index][Second]].BPRight) continue;
								else if (Reads[Inv[Box_index][First]].BPRight > Reads[Inv[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								//else CompareResult = 2;
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//			
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = Inv[Box_index][First];
  	                     Inv[Box_index][First] = Inv[Box_index][Second];
  	                     Inv[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < InversionsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[Inv[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();  
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					//GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         //GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					RealStart = IndelEvents[EventIndex].RealStart;
					RealEnd = IndelEvents[EventIndex].RealEnd;
					if (IndelEvents[EventIndex].Support < NumRead2ReportCutOff) continue;
					// report max one
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;						
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
					}
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Inversions (INV): " << NumberOfInvInstances << endl << endl;
}

void SortOutputInv_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Inv[], ofstream & InvOutf) {
   cout << "Sorting and outputing Inversions with non-template sequence ..." << endl;
   unsigned int InversionsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
	int Count_INV_NT_output = 0;
	/*
	 unsigned int C_S = 0;
	 unsigned int C_E = 0;
	 unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	 unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	 unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
		
      if (Inv[Box_index].size() >= NumRead2ReportCutOff) {
         InversionsNum = Inv[Box_index].size();
			//cout << Box_index << "\t" << Inv[Box_index].size() << endl;
			//InputIndels.clear();
			//for (int i = 0; i < InversionsNum; i++) InputIndels.push_back(Reads[Inv[Box_index][i]]);
         for (unsigned int First = 0; First < InversionsNum - 1; First++) {
				
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < InversionsNum; Second++) {
						//cout << InputIndels[First].BPLeft << "\t" << InputIndels[First].BPRight << "\t" 
						//     << InputIndels[Second].BPLeft << "\t" << InputIndels[Second].BPRight<< endl;
                  //if (InputIndels[Second].Unique) 
						{
							if (Reads[Inv[Box_index][First]].ReadLength == Reads[Inv[Box_index][Second]].ReadLength) {
								if (Reads[Inv[Box_index][First]].LeftMostPos == Reads[Inv[Box_index][Second]].LeftMostPos)
									Reads[Inv[Box_index][Second]].Unique = false;
							}
							if (Reads[Inv[Box_index][First]].BPLeft < Reads[Inv[Box_index][Second]].BPLeft) continue;
                     else if (Reads[Inv[Box_index][First]].BPLeft > Reads[Inv[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[Inv[Box_index][First]].BPLeft == Reads[Inv[Box_index][Second]].BPLeft) {
								if (Reads[Inv[Box_index][First]].BPRight < Reads[Inv[Box_index][Second]].BPRight) continue;
								else if (Reads[Inv[Box_index][First]].BPRight > Reads[Inv[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
									if (Reads[Inv[Box_index][First]].NT_size < Reads[Inv[Box_index][Second]].NT_size) continue;
									else if (Reads[Inv[Box_index][First]].NT_size > Reads[Inv[Box_index][Second]].NT_size) CompareResult = 1;
									//else { // InputIndels[First].NT_size == InputIndels[Second].NT_size
									//	short Compare2Str = CompareTwoString(InputIndels[First].NT_str, InputIndels[Second].NT_str );
									//	if (Compare2Str > 0) CompareResult = 1;
									//	else if (Compare2Str == 0) CompareResult = 2;
									//	else continue;
									//}
								}
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//		
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = Inv[Box_index][First];
  	                     Inv[Box_index][First] = Inv[Box_index][Second];
  	                     Inv[Box_index][Second] = Temp4Exchange;
							}
							//if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < InversionsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[Inv[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();  
			//cout << "GoodNum " << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
				//cout << GoodIndex << "\t" << GoodIndels[GoodIndex].BPLeft << "\t" << GoodIndels[GoodIndex].BPRight << "\t" << OneIndelEvent.BPLeft << "\t" << OneIndelEvent.BPRight << endl;
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					//GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			//cout << OneIndelEvent.Support << "\t" << OneIndelEvent.Start << "\t" << OneIndelEvent.End << endl;
         //GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					//cout << IndelEvents[EventIndex].Start << "\t" << IndelEvents[EventIndex].End << "\t" << IndelEvents[EventIndex].Support << endl;
					RealStart = IndelEvents[EventIndex].RealStart;
					RealEnd = IndelEvents[EventIndex].RealEnd;
					if (IndelEvents[EventIndex].Support < NumRead2ReportCutOff) continue;
					// report max one
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
						Count_INV_NT_output++;
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
						Count_INV_NT_output++;
					}
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Inversions with non-template sequence (INV_NT): " << Count_INV_NT_output << endl << endl;
}

void SortOutputDI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> DI[], ofstream & DIOutf) {
   cout << "Sorting and outputing deletions with non-template sequences ..." << endl;
   unsigned int DINum;
   short CompareResult;
   unsigned Temp4Exchange;
	/*
	 unsigned int C_S = 0;
	 unsigned int C_E = 0;
	 unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	 unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	 unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      //cout << "Box_index: "   << Box_index << endl;
      if (DI[Box_index].size() >= NumRead2ReportCutOff) {
         DINum = DI[Box_index].size();
			//InputIndels.clear();
			//for (int i = 0; i < DINum; i++) InputIndels.push_back(Reads[DI[Box_index][i]]);
         for (unsigned int First = 0; First < DINum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < DINum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (Reads[DI[Box_index][First]].ReadLength == Reads[DI[Box_index][Second]].ReadLength) {
								if (Reads[DI[Box_index][First]].LeftMostPos == Reads[DI[Box_index][Second]].LeftMostPos)
									Reads[DI[Box_index][Second]].Unique = false;
							}
							if (Reads[DI[Box_index][First]].BPLeft < Reads[DI[Box_index][Second]].BPLeft) continue;
                     else if (Reads[DI[Box_index][First]].BPLeft > Reads[DI[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[DI[Box_index][First]].BPLeft == Reads[DI[Box_index][Second]].BPLeft) {
								if (Reads[DI[Box_index][First]].BPRight 
									 < Reads[DI[Box_index][Second]].BPRight) continue;
								else if (Reads[DI[Box_index][First]].BPRight > Reads[DI[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								else {
									if (Reads[DI[Box_index][First]].NT_size < Reads[DI[Box_index][Second]].NT_size) continue;
									else if (Reads[DI[Box_index][First]].NT_size > Reads[DI[Box_index][Second]].NT_size) CompareResult = 1;
									//else CompareResult = 2;
									//else {
									//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
									//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
									//			InputIndels[Second].Unique = false;
									//		}
									//			
									//	}
									//}
							   }
                     }
							if (CompareResult == 1) {
								Temp4Exchange = DI[Box_index][First];
  	                     DI[Box_index][First] = DI[Box_index][Second];
  	                     DI[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
						}
               }
            }
         }
         GoodIndels.clear();
	      IndelEvents.clear();
			
			for (unsigned int First = 0; First < DINum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[DI[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size(); 
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
         //cout << "GoodNum: " << GoodNum << endl;   string InsertedStr;   string NT_str;  short NT_size;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
			OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.NT_size = GoodIndels[0].NT_size;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
			OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			//OneIndelEvent.IndelStr = GoodIndels[0].NT_str;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].IndelSize == OneIndelEvent.IndelSize  
                && GoodIndels[GoodIndex].NT_size == OneIndelEvent.NT_size)
					//&& OneIndelEvent.IndelStr == GoodIndels[GoodIndex].NT_str)
               OneIndelEvent.End = GoodIndex;
				else  {
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.IndelSize =  GoodIndels[GoodIndex].IndelSize;
					OneIndelEvent.NT_size = GoodIndels[GoodIndex].NT_size;
					//OneIndelEvent.IndelStr = GoodIndels[GoodIndex].NT_str;
            }
         }
			
			//if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff) 
			IndelEvents.push_back(OneIndelEvent);
			unsigned int RealStart;
			unsigned int RealEnd;
			for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
				if (IndelEvents[EventIndex].End - IndelEvents[EventIndex].Start + 1 >= NumRead2ReportCutOff) 
				{
					RealStart = GoodIndels[IndelEvents[EventIndex].Start].BPLeft;
					RealEnd = GoodIndels[IndelEvents[EventIndex].Start].BPRight;
					//if (IndelEvents[EventIndex].IndelSize < 100) {}
					//if (IndelSize < BalanceCutoff) {
					//        OutputDI(GoodIndels, CurrentChr, 
					//                 IndelEvents[EventIndex].Start, 
					//       	        IndelEvents[EventIndex].End, 
					// 		RealStart, RealEnd, DIOutf);
					//        NumberOfDIInstances++;
					//}
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputDI(GoodIndels, CurrentChr, 
									IndelEvents[EventIndex].Start, 
									IndelEvents[EventIndex].End, 
									RealStart, RealEnd, DIOutf);
						NumberOfDIInstances++;					
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputDI(GoodIndels, CurrentChr, 
									IndelEvents[EventIndex].Start, 
									IndelEvents[EventIndex].End, 
									RealStart, RealEnd, DIOutf);
						NumberOfDIInstances++;
					}
				}
			} 
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "deletions with non-template sequences: " << NumberOfDIInstances << endl << endl;
}


bool NotInVector(const string & OneTag, const vector <string> & VectorTag) {
	for (unsigned i = 0; i < VectorTag.size(); i++) {
		if (OneTag == VectorTag[i]) return false;
	}
	return true;
}


bool ReportEvent(const vector <SPLIT_READ> & Deletions, const unsigned int & S, const unsigned int & E) {
   short ReadLength = Deletions[S].ReadLength - Deletions[S].NT_size;
   short Min_Length = (short)((ReadLength  * Min_Filter_Ratio) + 0.5) - 1;  
   short Max_Length = (short)(ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1;  
   bool LeftMin = false;
   bool LeftMax = false;
	bool RightMin = false;
   bool RightMax = false;
   for (unsigned i = S; i <= E; i++) {
		ReadLength = Deletions[i].ReadLength - Deletions[i].NT_size;
		Min_Length = (short)((ReadLength * Min_Filter_Ratio) + 0.5) - 1;
		Max_Length = (short)(ReadLength * (1 - Min_Filter_Ratio) - 0.5) - 1; 
      if (Deletions[i].BP <= Min_Length) {
         LeftMin = true;
	      //break;
      }
		if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size <= Min_Length) {
         RightMin = true;
	      //break;
      }
		if (Deletions[i].BP >= Max_Length) {
         LeftMax = true;
	      //break;
      }
		if (Deletions[i].ReadLength - Deletions[i].BP - Deletions[i].NT_size >= Max_Length) {
         RightMax = true;
	      //break;
      }
   }
	
   if (LeftMin && LeftMax && RightMin && RightMax) return true;
   else return false;
}

void GetRealStart4Deletion(const string & TheInput, unsigned int & RealStart, unsigned int & RealEnd) {
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

void GetRealStart4Insertion(const string & TheInput, const string & InsertedStr, unsigned int & RealStart, unsigned int & RealEnd) {
   unsigned int IndelSize = InsertedStr.size();
   unsigned int PosIndex = RealStart + SpacerBeforeAfter;
   //unsigned int Start = PosIndex + 1;
	
   //unsigned int End = RealEnd + SpacerBeforeAfter - 1;
   for (int i = IndelSize - 1; i >=0; i--) {
      if (TheInput[PosIndex] == InsertedStr[i]) PosIndex--;
      else break;
   }
   if (PosIndex == RealStart + SpacerBeforeAfter - IndelSize) {
      while (TheInput[PosIndex] == TheInput[PosIndex + IndelSize]) {
         PosIndex--;
      }
   }
   RealStart = PosIndex - SpacerBeforeAfter;
   PosIndex = RealEnd + SpacerBeforeAfter;
   for (unsigned int i = 0; i < IndelSize; i++) {
      if (TheInput[PosIndex] == InsertedStr[i]) PosIndex++;
      else break;
   }
   if (PosIndex == RealEnd + SpacerBeforeAfter + IndelSize) {
      while (TheInput[PosIndex] == TheInput[PosIndex - IndelSize]) {
         PosIndex++;
      }
   }
   RealEnd = PosIndex - SpacerBeforeAfter;
}

vector <Region> Merge(const vector <Region> & AllRegions) {
	return AllRegions;
}

void GetCloseEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read) {
	
	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(3 * Temp_One_Read.InsertSize);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;
	
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (Temp_One_Read.MatchedD == Plus) {
			CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
			Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			End = Start + 3 * Temp_One_Read.InsertSize;
			LeftChar = CurrentReadSeq[0];
			if (LeftChar != 'N') {
            //#pragma omp parallel default(shared)
				{
					// #pragma omp for
				   for (int pos = Start; pos < End; pos++) {
					   if (CurrentChr[pos] == LeftChar) {
						   PD[0].push_back(pos);
					   }
					   else PD[1].push_back(pos);
				   }
				}
			}
			else { //Match2N[(short)'A'] = 'N';
            //#pragma omp parallel default(shared)
				{
               //#pragma omp for
				   for (int pos = Start; pos < End; pos++) {
					   if (Match2N[(short)CurrentChr[pos]] == 'N') {
						   PD[0].push_back(pos);
					   }
						//else Left_PD[1].push_back(pos);
				   }
				}
			}
			CheckLeft_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP); // LengthStr
			//Direction = Minus;
			if (UP.empty()) {}
			else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
			else {
				for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
					   Temp_One_Read.Used = false;
					   Temp_One_Read.UP_Close.push_back(UP[LeftUP_index]);
				   }
			   }	
			}

			UP.clear();
		}
		else if (Temp_One_Read.MatchedD == Minus) {
			CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
			End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			Start = End - 3 * Temp_One_Read.InsertSize;
			RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
			//cout << "1\t" << PD[0].size() << "\t" << PD[1].size() << endl;
			CheckRight_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
			//cout << UP.size() << endl;
			//Direction = '+';
			if (UP.empty()) {}
			else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
			else {
				for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
					if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
						Temp_One_Read.Used = false;
						Temp_One_Read.UP_Close.push_back(UP[RightUP_index]);
					} 
				}				
			}
			UP.clear();
		}
	}
	else { // TOTAL_SNP_ERROR_CHECKED_Minus
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
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
			CheckLeft_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
			//Direction = Minus;
			if (UP.empty()) {}
			else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
			else {
				for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
					   Temp_One_Read.Used = false;
					   Temp_One_Read.UP_Close.push_back(UP[LeftUP_index]);
				   }
			   }	
			}

			UP.clear();
		}
		else if (Temp_One_Read.MatchedD == Minus) {
			CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
			End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
			Start = End - 3 * Temp_One_Read.InsertSize;
			RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
			//cout << "2" << PD[0].size() << "\t" << PD[1].size() << endl;
			CheckRight_Close(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
			//cout << UP.size() << endl;
			//Direction = '+';
			if (UP.empty()) {}
			else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
			else {
				for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
					if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
						Temp_One_Read.Used = false;
						Temp_One_Read.UP_Close.push_back(UP[RightUP_index]);
					} 
				}
			}
			UP.clear();
		}
	}
	//if (!Temp_One_Read.UP_Close.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Close);
	return;
}

void GetFarEnd_SingleStrandDownStreamInsertions(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	//Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	//if (Temp_One_Read.MatchedRelPos > 160000)
	//cout << Temp_One_Read.UnmatchedSeq << Temp_One_Read.UnmatchedSeq.size() << endl;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;
	
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + Range))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.MatchedD == Minus) {
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
		
		End = Temp_One_Read.UP_Close[0].AbsLoc + Temp_One_Read.UP_Close[0].LengthStr;
		if (End > SpacerBeforeAfter + Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex])
		   Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else Start = SpacerBeforeAfter;
		
		
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;
		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		if (UP.empty()) {}
		else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
		else {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
				   Temp_One_Read.Used = false;
				   Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
			   }
		   }
		}
		UP.clear();
	}
	else if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		
		Start = Temp_One_Read.UP_Close[0].AbsLoc - Temp_One_Read.UP_Close[0].LengthStr;
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;
		
		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		else {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		
		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = '+';
		if (UP.empty()) {}
		else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
		else {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
				   Temp_One_Read.Used = false;
				   Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
			   } 
		   }
		}

		UP.clear();
	}		
	//if (!Temp_One_Read.UP_Far.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Far);
	return;
}

void GetFarEnd_SingleStrandDownStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	//cout << "in GetFarEnd_SingleStrandDownStream" << endl;
	//Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	//if (Temp_One_Read.MatchedRelPos > 160000)
	//cout << "RangeIndex: " << RangeIndex << "\t" << DSizeArray[RangeIndex] << endl;
	//if (RangeIndex == 4)
	//	cout << "\nin\t" << DSizeArray[RangeIndex] << "\t" << Temp_One_Read.UnmatchedSeq << "\t" << Temp_One_Read.UnmatchedSeq.size() << endl;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;
	
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + Range))/log(4.0) + 0.8) + 3;
	//cout << "Temp_One_Read.MinClose: " << Temp_One_Read.MinClose << endl;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.MatchedD == Minus) {
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
		
		End = Temp_One_Read.UP_Close[0].AbsLoc + Temp_One_Read.UP_Close[0].LengthStr - Temp_One_Read.ReadLength;
		if (End > SpacerBeforeAfter + Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex])
		   Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else Start = SpacerBeforeAfter;
		
		
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;
		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		//if (Temp_One_Read.MatchedRelPos > 160000)
		//if (RangeIndex == 4)
		//cout << "- " << PD[0].size() << "\t" << PD[1].size() << endl;
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		if (UP.empty()) {}
		else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
		else {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
				   Temp_One_Read.Used = false;
				   Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
			   }
		   }
		}

		UP.clear();
	}
	else if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		
		Start = Temp_One_Read.UP_Close[0].AbsLoc - Temp_One_Read.UP_Close[0].LengthStr + Temp_One_Read.ReadLength;
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;
		
		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		//if (Temp_One_Read.MatchedRelPos > 160000)
		//if (RangeIndex == 4)
		//cout << "+ " << PD[0].size() << "\t" << PD[1].size() << endl;
		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		if (UP.empty()) {}
		else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
		else {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
				   Temp_One_Read.Used = false;
				   Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
			   } 
		   }
		}

		UP.clear();
	}		
	//if (!Temp_One_Read.UP_Far.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Far);
	return;
}

void GetFarEnd_SingleStrandUpStream(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
	Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
	char LeftChar, RightChar;
	string CurrentReadSeq;
	vector <unsigned int> PD[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]);
		//PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	vector <UniquePoint> UP;
	//char Direction;
	int Start, End;
	short BP_Start;// = MinClose;
	short BP_End;// = ReadLength - MinClose;
	
	//for (int i = 0; i < TOTAL_SNP_ERROR_CHECKED; i++) {
	//	PD[i].clear();
	//}
	//UP.clear();
	//Temp_One_Read.UP_Close.clear();
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 2 + DSizeArray[RangeIndex]))/log(4.0) + 0.8) + 3;
	BP_Start = Temp_One_Read.MinClose + RangeIndex;
	BP_End = Temp_One_Read.ReadLengthMinus;
	//Temp_One_Read.Unique = true;
	if (Temp_One_Read.MatchedD == Minus) {
		//CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		CurrentReadSeq = Temp_One_Read.UnmatchedSeq;
		
		Start = Temp_One_Read.UP_Close[0].AbsLoc + Temp_One_Read.UP_Close[0].LengthStr;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2;
		if (End > CurrentChr.size() - SpacerBeforeAfter) End = CurrentChr.size() - SpacerBeforeAfter;
		
		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else { // TOTAL_SNP_ERROR_CHECKED_Minus
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == LeftChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckLeft_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = Minus;
		if (UP.empty()) {}
		else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
		else {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[LeftUP_index])) {
				   Temp_One_Read.Used = false;
				   Temp_One_Read.UP_Far.push_back(UP[LeftUP_index]);
			   }
		   }
		}

		UP.clear();
	}
	else if (Temp_One_Read.MatchedD == Plus) {
		CurrentReadSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		
		End = Temp_One_Read.UP_Close[0].AbsLoc - Temp_One_Read.UP_Close[0].LengthStr;
		//Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter;
		
		if (End > DSizeArray[RangeIndex] + Temp_One_Read.InsertSize * 2 + SpacerBeforeAfter) 
			Start = End - DSizeArray[RangeIndex] - Temp_One_Read.InsertSize * 2;
		else Start = SpacerBeforeAfter;
		
		RightChar = CurrentReadSeq[Temp_One_Read.ReadLengthMinus];
		if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
      else {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == RightChar) {
						PD[0].push_back(pos);
					}
					//else PD[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD[0].push_back(pos);
					}
					//else Left_PD[1].push_back(pos);
				}
			}
		}
		CheckRight_Far(Temp_One_Read, CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);
		//Direction = '+';
		if (UP.empty()) {}
		else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
		else {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[RightUP_index])) {
				   Temp_One_Read.Used = false;
				   Temp_One_Read.UP_Far.push_back(UP[RightUP_index]);
			   }   
		   }
		}

		UP.clear();
	}		
	//if (!Temp_One_Read.UP_Far.empty())
	//   CleanUniquePoints(Temp_One_Read.UP_Far);
	return;
}

void GetFarEnd_OtherStrand(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	//short ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//short ReadLengthMinus = ReadLength - 1;
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 5 + Range * 2))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;
	int Start, End;
	short BP_Start = Temp_One_Read.MinClose + RangeIndex;
	short BP_End = Temp_One_Read.ReadLengthMinus;
	//char Direction;
	
	vector <UniquePoint> UP;
	vector <unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	vector <unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	
	if (Temp_One_Read.MatchedD == Plus) {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength - 2 * Temp_One_Read.InsertSize - DSizeArray[RangeIndex] ;
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength + 3 * Temp_One_Read.InsertSize + DSizeArray[RangeIndex];		
	}
	else {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - 3 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength - DSizeArray[RangeIndex];
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter + 2 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength + DSizeArray[RangeIndex];		
	}
	
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(End - Start + 1);
		PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	
	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short)CurrentBase];
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (Temp_One_Read.MatchedD == Plus) {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					else PD_Plus[1].push_back(pos);
					//if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD_Plus[0].push_back(pos);
						//PD_Minus[0].push_back(pos);
					}
					else {
						PD_Plus[1].push_back(pos);
						//PD_Minus[1].push_back(pos);
					}
				}
			}
		}
		else { // -
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					//if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						//PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					}
					else {
						//PD_Plus[1].push_back(pos);
						PD_Minus[1].push_back(pos);
					}
				}
			}
		}
	}
   else { // TOTAL_SNP_ERROR_CHECKED_Minus
		if (Temp_One_Read.MatchedD == Plus) {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					//if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						PD_Plus[0].push_back(pos);
						//PD_Minus[0].push_back(pos);
					}
					//else {
					//	PD_Plus[1].push_back(pos);
					//PD_Minus[1].push_back(pos);
					//}
				}
			}
		}
		else { // -
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					//if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
					//else PD_Plus[1].push_back(pos);
					if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
					//else PD_Minus[1].push_back(pos);
				}
			}
			else { //Match2N[(short)'A'] = 'N';
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short)CurrentChr[pos]] == 'N') {
						//PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					}
					//else {
					//PD_Plus[1].push_back(pos);
					//	PD_Minus[1].push_back(pos);
					//}
				}
			}
		}
	}
	
	CheckBoth(Temp_One_Read, CurrentChr, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
	if (UP.empty()) {}
	else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
	else {
		for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
		   if (UP[UP_index].Direction == Plus) {
			//Direction = Minus;
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				   Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		   }
		   else {
			   //Direction = '+';
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index])) 
				   Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		   }
	   }
	}

	UP.clear();
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

void GetFarEnd_BothStrands(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const short & RangeIndex) {
	//short ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//short ReadLengthMinus = ReadLength - 1;
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 5 + DSizeArray[RangeIndex] * 2))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;
	int Start, End;
	short BP_Start = Temp_One_Read.MinClose + RangeIndex;
	short BP_End = Temp_One_Read.ReadLengthMinus;
	//char Direction;
	
	vector <UniquePoint> UP;
	vector <unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	vector <unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	
	
	
	if (Temp_One_Read.MatchedD == Plus) {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength - 2 * Temp_One_Read.InsertSize - DSizeArray[RangeIndex] ;
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - Temp_One_Read.ReadLength + 3 * Temp_One_Read.InsertSize + DSizeArray[RangeIndex];		
	}
	else {
		Start = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter - 3 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength - DSizeArray[RangeIndex];
		End = Temp_One_Read.MatchedRelPos + SpacerBeforeAfter + 2 * Temp_One_Read.InsertSize + Temp_One_Read.ReadLength + DSizeArray[RangeIndex];		
	}
	
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(End - Start + 1);
		PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
	
	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short)CurrentBase];
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				else PD_Minus[1].push_back(pos);
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				else {
					PD_Plus[1].push_back(pos);
					PD_Minus[1].push_back(pos);
				}
			}
		}
	}
	else { // TOTAL_SNP_ERROR_CHECKED_Minus
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				//else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				//else PD_Minus[1].push_back(pos);
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				//else {
				//	PD_Plus[1].push_back(pos);
				//	PD_Minus[1].push_back(pos);
				//}
			}
		}
	}
	
	CheckBoth(Temp_One_Read, CurrentChr, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
	if (UP.empty()) {}
	else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
	else {
		for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
		   if (UP[UP_index].Direction == Plus) {
			   //Direction = Minus;
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				   Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		   }
		   else {
			   //Direction = '+';
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index])) 
				   Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		   }
	   }
	}

	UP.clear();
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

void GetFarEnd(const string & CurrentChr, SPLIT_READ & Temp_One_Read, const int & in_start, const int & in_end) {
	//short ReadLength = Temp_One_Read.UnmatchedSeq.size();
	//short ReadLengthMinus = ReadLength - 1;
	//MinClose = short(log((double)(in_end - in_start))/log(4.0) + 0.8) + 3 + MAX_SNP_ERROR;
	int Start, End;
	short BP_Start = Temp_One_Read.MinClose;
	short BP_End = Temp_One_Read.ReadLengthMinus;
	//char Direction;
	vector <UniquePoint> UP;
	vector <unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	vector <unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(in_end - in_start + 1);
		PD_Minus[CheckIndex].reserve(in_end - in_start + 1);
	}
	Start = in_start;
	End = in_end;		
	
	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short)CurrentBase];
	if (Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus) {
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				else PD_Minus[1].push_back(pos);
				//}
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				else {
					PD_Plus[1].push_back(pos);
					PD_Minus[1].push_back(pos);
				}	
				//}
			}
		}
	}
   else {
		if (CurrentBase != 'N') {
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (CurrentChr[pos] == CurrentBase) PD_Plus[0].push_back(pos);
				//else PD_Plus[1].push_back(pos);
				if (CurrentChr[pos] == CurrentBaseRC) PD_Minus[0].push_back(pos);
				//else PD_Minus[1].push_back(pos);
				//}
			}
		}
		else { //Match2N[(short)'A'] = 'N';
			for (int pos = Start; pos < End; pos++) {
				//if (BreakDancerMask[pos] == BD_char) {
				if (Match2N[(short)CurrentChr[pos]] == 'N') {
					PD_Plus[0].push_back(pos);
					PD_Minus[0].push_back(pos);
				}
				//else {
				//	PD_Plus[1].push_back(pos);
				//	PD_Minus[1].push_back(pos);
				//}	
				//}
			}
		}
	}
	//cout << PD_Plus[0].size() << "\t" << PD_Plus[1].size() << endl;
	//cout << PD_Minus[0].size() << "\t" << PD_Minus[1].size() << endl;
	CheckBoth(Temp_One_Read, CurrentChr, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
	//cout << UP.size() << endl;
	if (UP.empty()) {}
	else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {}
	else {
		for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
		   if (UP[UP_index].Direction == Plus) {
			   //Direction = Minus;
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index]))
				   Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		   }
		   else {
			   //Direction = '+';
			   if (CheckMismatches(CurrentChr, Temp_One_Read.UnmatchedSeq, UP[UP_index])) 
				   Temp_One_Read.UP_Far.push_back(UP[UP_index]);
		   }
	   }
	}

	UP.clear();
	
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

void CheckBoth(const SPLIT_READ & OneRead, 
					const string & TheInput, 
               const string & CurrentReadSeq,
               const vector <unsigned int> PD_Plus[],
					const vector <unsigned int> PD_Minus[],
               const short & BP_Start,
               const short & BP_End,
               const short & CurrentLength,
               vector <UniquePoint> & UP) {
	int Sum;
   if (CurrentLength >= BP_Start && CurrentLength <= BP_End) {
      // put it to LeftUP if unique
		for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
			if (PD_Plus[i].size() + PD_Minus[i].size() == 1 && CurrentLength >= BP_Start + i) {
				Sum = 0;
				if (ADDITIONAL_MISMATCH) 
		   		for (short j = 1; j <= ADDITIONAL_MISMATCH; j++) 
			   		Sum += PD_Plus[i + j].size() + PD_Minus[i + j].size();
				
				if (Sum == 0 && i <= (short)(Seq_Error_Rate * CurrentLength + 1)) {
					UniquePoint TempOne;
					TempOne.LengthStr = CurrentLength;
					if (PD_Plus[i].size() == 1) {
						TempOne.Direction = FORWARD;
						TempOne.Strand = SENSE;
						TempOne.AbsLoc = PD_Plus[i][0];
					}
					else if (PD_Minus[i].size() == 1) {
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
      vector <unsigned int> PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		vector <unsigned int> PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
		for (int CheckedIndex = 0; CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++) {
			PD_Plus_Output[CheckedIndex].reserve(PD_Plus[CheckedIndex].size());
			PD_Minus_Output[CheckedIndex].reserve(PD_Minus[CheckedIndex].size());
		}
      const char CurrentChar = CurrentReadSeq[CurrentLength];
		const char CurrentCharRC = Convert2RC4N[(short)CurrentChar];
      //const int SizeOfCurrent = Left_PD.size();
      unsigned int pos;
		int SizeOfCurrent;
		//if (TOTAL_SNP_ERROR_CHECKED_Minus) 
		{
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++) {
				SizeOfCurrent = PD_Plus[i].size();
				if (CurrentChar == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Plus[i][j] + 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							PD_Plus_Output[i].push_back(pos);
						else PD_Plus_Output[i + 1].push_back(pos);
					}         
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Plus[i][j] + 1;
						if (TheInput[pos] == CurrentChar)
							PD_Plus_Output[i].push_back(pos);
						else PD_Plus_Output[i + 1].push_back(pos); 
					}         
				}		
				SizeOfCurrent = PD_Minus[i].size();
				if (CurrentCharRC == 'N') {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Minus[i][j] - 1;
						if (Match2N[(short)TheInput[pos]] == 'N')
							PD_Minus_Output[i].push_back(pos);
						else PD_Minus_Output[i + 1].push_back(pos);
					}         
				}
				else {
					for (int j = 0; j < SizeOfCurrent; j++) {
						pos =  PD_Minus[i][j] - 1;
						if (TheInput[pos] == CurrentCharRC)
							PD_Minus_Output[i].push_back(pos);
						else PD_Minus_Output[i + 1].push_back(pos); 
					}         
				}
			}
			
			SizeOfCurrent = PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentChar == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}         
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Plus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
					if (TheInput[pos] == CurrentChar)
						PD_Plus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}         
			}		
			SizeOfCurrent = PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size();
			if (CurrentCharRC == 'N') {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (Match2N[(short)TheInput[pos]] == 'N')
						PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}         
			}
			else {
				for (int j = 0; j < SizeOfCurrent; j++) {
					pos =  PD_Minus[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
					if (TheInput[pos] == CurrentCharRC)
						PD_Minus_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].push_back(pos);
				}         
			}
			
			
			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++) {
				Sum += PD_Plus_Output[i].size() + PD_Minus_Output[i].size();
			}
			if (Sum) {
				const short CurrentLengthOutput = CurrentLength + 1;
				CheckBoth(OneRead, TheInput, CurrentReadSeq, PD_Plus_Output, PD_Minus_Output,
							 BP_Start, BP_End,
							 CurrentLengthOutput, UP);         
			}
			else return;
		}
   }   
   else return;

}


void CleanUniquePoints (vector <UniquePoint> & Input_UP) {
	vector <UniquePoint> TempUP; //vector <UniquePoint> UP_Close; UP_Far
	UniquePoint LastUP = Input_UP[Input_UP.size() - 1];
	//TempUP.push_back(LastUP);
	char LastDirection = LastUP.Direction;
	char LastStrand = LastUP.Strand;
	unsigned int Terminal;
	
	if (LastDirection == FORWARD) {
		Terminal = LastUP.AbsLoc - LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc - Input_UP[i].LengthStr) 
					TempUP.push_back(Input_UP[i]);
			}
		}
	}
	else if (LastDirection == BACKWARD) {
		Terminal = LastUP.AbsLoc + LastUP.LengthStr;
		for (unsigned i = 0; i < Input_UP.size(); i++) {
			if (Input_UP[i].Direction == LastDirection && Input_UP[i].Strand == LastStrand) {
				if (Terminal == Input_UP[i].AbsLoc + Input_UP[i].LengthStr) 
					TempUP.push_back(Input_UP[i]);
			}
		}			
	}
   Input_UP.clear();
	Input_UP = TempUP;
}

void SortOutputLI(const string & CurrentChr, vector <SPLIT_READ> & Reads, ofstream & LargeInsertionOutf) {
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	// find LI combinations
	uint8_t * plus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	uint8_t * minus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	int32_t * EventIndex_Pos = new int32_t[CurrentChr.size() + 1];
	for (unsigned i = 0; i < CurrentChr.size() + 1; i++) {
		plus_LI_Pos[i] = 0;
		minus_LI_Pos[i] = 0;
		EventIndex_Pos[i] = -1;
	}
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		//UP_Close_index = ;
		if (Reads[Index].Found || Reads[Index].Used || !Reads[Index].UP_Far.empty()) continue;
		temp_AbsLoc = Reads[Index].UP_Close[Reads[Index].UP_Close.size() - 1].AbsLoc;
		if (plus_LI_Pos[temp_AbsLoc] < Max_short) {
			if (Reads[Index].MatchedD == Plus) plus_LI_Pos[temp_AbsLoc]++;	
		}
		if (minus_LI_Pos[temp_AbsLoc] < Max_short) {
			if (Reads[Index].MatchedD == Minus) minus_LI_Pos[temp_AbsLoc]++;
		}
	}
	vector <LI_Pos> LI_Positions;
	LI_Pos temp_LI_pos;
	bool SkipThis;
	int LI_Positions_Size = 0;
	bool SkipPlus;
	for (int Index_Minus = SpacerBeforeAfter; Index_Minus < CurrentChr.size() - SpacerBeforeAfter; Index_Minus++) {
		SkipPlus = false;
		for (int MaskedPosIndexMinus = Index_Minus + 10; MaskedPosIndexMinus >= Index_Minus - 10; MaskedPosIndexMinus--) {
			if (CurrentChrMask[MaskedPosIndexMinus] == 'B') {
				Index_Minus = MaskedPosIndexMinus + 10;
				SkipPlus = true;
				break;
			} 
		}
		if (SkipPlus == false && minus_LI_Pos[Index_Minus] >= NumRead2ReportCutOff) {
			for (int Index_Plus = Index_Minus - 1; Index_Plus <= Index_Minus + 30; Index_Plus++) {
				SkipThis = false;
				for (int MaskedPosIndexPlus = Index_Plus + 10; MaskedPosIndexPlus >= Index_Plus - 10; MaskedPosIndexPlus--) {
					if (CurrentChrMask[MaskedPosIndexPlus] == 'B') {
						if ((MaskedPosIndexPlus + 10) > Index_Minus ) { Index_Minus = MaskedPosIndexPlus + 10; }
						SkipThis = true;
						break;
					} 
				}
				if (SkipThis == false && plus_LI_Pos[Index_Plus] >= NumRead2ReportCutOff) {
					temp_LI_pos.Plus_Pos = Index_Plus;
					temp_LI_pos.Minus_Pos = Index_Minus;
					CurrentChrMask[Index_Plus] == 'B';
					CurrentChrMask[Index_Minus] == 'B';
					//cout << Index_Plus << "\t" << (short)plus_LI_Pos[Index_Plus] << "\t" 
					//     << Index_Minus << "\t" << (short)minus_LI_Pos[Index_Minus] << endl;
					temp_LI_pos.WhetherReport = false;
					LI_Positions.push_back(temp_LI_pos);
					EventIndex_Pos[Index_Plus] = LI_Positions_Size;
					EventIndex_Pos[Index_Minus] = LI_Positions_Size;
					LI_Positions_Size++;
					//Index_Minus += 30;
				}
         }
		}
	}
   //cout << "LI: " << LI_Positions.size() << endl;
	static int Count_LI = 0;
	// find LI supporting reads
	
	
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		//UP_Close_index = Reads[Index].UP_Close.size() - 1;
		if (Reads[Index].Used || !Reads[Index].UP_Far.empty()) continue;
		temp_AbsLoc = Reads[Index].UP_Close[Reads[Index].UP_Close.size() - 1].AbsLoc;
		if (EventIndex_Pos[temp_AbsLoc] == -1) continue;
		Reads[Index].Used = true;
		if (Reads[Index].MatchedD == Plus) {
			LI_Positions[EventIndex_Pos[temp_AbsLoc]].Plus_Reads.push_back(Index);
		}
		else {
			LI_Positions[EventIndex_Pos[temp_AbsLoc]].Minus_Reads.push_back(Index);	
		}
	}
	
	vector <SPLIT_READ> temp_Plus_Reads, temp_Minus_Reads;
	
	bool temp_BalancedPlus_Plus, temp_BalancedPlus_Minus, temp_BalancedMinus_Plus, temp_BalancedMinus_Minus;
	short temp_LengthStr;
	for (unsigned LI_index = 0; LI_index < LI_Positions.size(); LI_index++) {
		if (LI_Positions[LI_index].Minus_Reads.empty() || LI_Positions[LI_index].Plus_Reads.empty()) continue;
		temp_BalancedPlus_Plus = false;
		temp_BalancedPlus_Minus = false;
		temp_BalancedMinus_Plus = false;
		temp_BalancedMinus_Minus = false;
		temp_Plus_Reads.clear();
		temp_Minus_Reads.clear();
		//cout << "Here: " << LI_index << "\t" << LI_Positions[LI_index].Minus_Reads.size() << "\t" << LI_Positions[LI_index].Plus_Reads.size() << endl;
		for (int i = 0; i < LI_Positions[LI_index].Minus_Reads.size(); i++) {
			temp_Minus_Reads.push_back(Reads[LI_Positions[LI_index].Minus_Reads[i]]);
		}
		//CheckConsistancy(temp_Minus_Reads);
		for (int i = 0; i < LI_Positions[LI_index].Minus_Reads.size(); i++) {
			UP_Close_index = temp_Minus_Reads[i].UP_Close.size() - 1;
			temp_LengthStr = temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > temp_Minus_Reads[i].ReadLength * 0.5)
				temp_BalancedMinus_Plus = true;
			else if ((float)temp_LengthStr < temp_Minus_Reads[i].ReadLength * 0.5)
				temp_BalancedMinus_Minus = true;
		}	
		for (int i = 0; i < LI_Positions[LI_index].Plus_Reads.size(); i++) {
			temp_Plus_Reads.push_back(Reads[LI_Positions[LI_index].Plus_Reads[i]]);
		}
		//CheckConsistancy(temp_Plus_Reads);
		for (int i = 0; i < LI_Positions[LI_index].Plus_Reads.size(); i++) {
			UP_Close_index = temp_Plus_Reads[i].UP_Close.size() - 1;
			temp_LengthStr = temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > temp_Plus_Reads[i].ReadLength * 0.5)
				temp_BalancedPlus_Plus = true;
			else if ((float)temp_LengthStr < temp_Plus_Reads[i].ReadLength * 0.5)
				temp_BalancedPlus_Minus = true;
		}
		
		unsigned int NumSupportPerTagPlus[VectorTag.size()];
		unsigned int NumSupportPerTagMinus[VectorTag.size()];
		for (short i = 0; i < VectorTag.size(); i++) {
		   NumSupportPerTagPlus[i] = 0;
			NumSupportPerTagMinus[i] = 0;
		}
		for (int i = 0; i < temp_Minus_Reads.size(); i++) {
			for (short j = 0; j < VectorTag.size(); j++) {
			   if (temp_Minus_Reads[i].Tag == VectorTag[j]) {
					NumSupportPerTagMinus[j]++;
					break;
				}	
			}
		}
		for (int i = 0; i < temp_Plus_Reads.size(); i++) {
			for (short j = 0; j < VectorTag.size(); j++) {
			   if (temp_Plus_Reads[i].Tag == VectorTag[j]) {
					NumSupportPerTagPlus[j]++;
					break;
				}	
			}
		}
		bool SupportedByOneSample = false;
		for (short j = 0; j < VectorTag.size(); j++) {
			//cout << NumSupportPerTagPlus[j] << "\t" << NumSupportPerTagMinus[j] << endl;
			if (NumSupportPerTagPlus[j] > 0 && NumSupportPerTagMinus[j] > 0) {
				SupportedByOneSample = true;
				break;
			}
		}
		
		short PositiveBool = 0;
		if (temp_BalancedPlus_Plus) PositiveBool++;
		if (temp_BalancedPlus_Minus) PositiveBool++;
		if (temp_BalancedMinus_Plus) PositiveBool++;
		if (temp_BalancedMinus_Minus) PositiveBool++;
		
		if (SupportedByOneSample && PositiveBool >= 3) 
		{
			//if (LI_Positions[LI_index].WhetherReport) 
			{
				LargeInsertionOutf << "########################################################" << endl;
				LargeInsertionOutf << Count_LI++ << "\tLI\tChrID " << temp_Plus_Reads[0].FragName << "\t" 
				<< LI_Positions[LI_index].Plus_Pos - SpacerBeforeAfter + 1 << "\t" 
				<< temp_Plus_Reads.size() << "\t"
				<< LI_Positions[LI_index].Minus_Pos - SpacerBeforeAfter  + 1 << "\t" 
				<< temp_Minus_Reads.size() << endl;
				
				LargeInsertionOutf <<  ( CurrentChr.substr(LI_Positions[LI_index].Plus_Pos - ReportLength + 1, ReportLength) )
				<< Cap2Low( CurrentChr.substr(LI_Positions[LI_index].Plus_Pos + 1, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Plus_Reads.size(); i++) {
					UP_Close_index = temp_Plus_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength - temp_LengthStr; j++) {
						LargeInsertionOutf << " ";
					}
					LargeInsertionOutf << ReverseComplement(temp_Plus_Reads[i].UnmatchedSeq) << endl;
				}
				
				LargeInsertionOutf << "--------------------------------------------------------" << endl;
				//LargeInsertionOutf << "-\t" << minus_LI_Pos[Index_Minus].NumReads << endl;
				LargeInsertionOutf << Cap2Low( CurrentChr.substr(LI_Positions[LI_index].Minus_Pos - ReportLength, ReportLength) )
				<< ( CurrentChr.substr(LI_Positions[LI_index].Minus_Pos, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Minus_Reads.size(); i++) {
					UP_Close_index = temp_Minus_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength + temp_LengthStr - temp_Minus_Reads[i].ReadLength; j++) {
						LargeInsertionOutf << " ";
					}
					LargeInsertionOutf << (temp_Minus_Reads[i].UnmatchedSeq) << endl;
				}	
			}
			
		}
		
		//LI_Positions[LI_index].WhetherReport = true; 
	}
	// output
	delete[] plus_LI_Pos ;
	delete[] minus_LI_Pos ;
	delete[] EventIndex_Pos;
	
	cout << "Breakpoints for large insertions (LI): " << Count_LI << "\n\n";
}


void SortOutputRest(const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <SPLIT_READ> & BP_Reads, ofstream & Outf_Rest) {
	SPLIT_READ one_BP_read;
	string HalfMapped, HalfUnmapped;
	int HalfMappedIndex, HalfUnmappedIndex;
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	//cout << "1" << endl;
	// find LI combinations
	uint8_t * plus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	uint8_t * minus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	for (unsigned i = 0; i < CurrentChr.size() + 1; i++) {
		plus_LI_Pos[i] = 0;
		minus_LI_Pos[i] = 0;
	}
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		if (Reads[Index].Found || Reads[Index].Used  || !Reads[Index].UP_Far.empty()) continue;
		UP_Close_index = Reads[Index].UP_Close.size() - 1;
		temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
		if (plus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP) {
			if (Reads[Index].MatchedD == Plus)
				plus_LI_Pos[temp_AbsLoc]++;			
		}
      if (minus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP) {
			if (Reads[Index].MatchedD == Minus)
				minus_LI_Pos[temp_AbsLoc]++;	
		}
	}
	vector <Rest_Pos> Rest_Positions;
	Rest_Pos temp_Rest_pos;
	
	//cout << "2" << endl;
	
	bool SkipThisPos;
	for (int Index = SpacerBeforeAfter; Index < CurrentChr.size() - SpacerBeforeAfter; Index++) {
		SkipThisPos = false;
		//for (int MaskIndex = Index + 10; MaskIndex >= Index - 10; MaskIndex--) {
		//	if (CurrentChrMask[MaskIndex] == 'B') {
		//		Index = MaskIndex + 10;
		//		SkipThisPos = true;
		//		break;
		//	}
		//}
		if (SkipThisPos == true) continue;
		if (plus_LI_Pos[Index] >= NumRead2ReportCutOff_BP) {
			temp_Rest_pos.Strand = Plus;
			temp_Rest_pos.Pos = Index;
			Rest_Positions.push_back(temp_Rest_pos);
		}
		if (minus_LI_Pos[Index] >= NumRead2ReportCutOff_BP) {
			temp_Rest_pos.Strand = Minus;
			temp_Rest_pos.Pos = Index;
			Rest_Positions.push_back(temp_Rest_pos);
		}
	}
	
	// find supporting reads
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		if (Reads[Index].Used || !Reads[Index].UP_Far.empty()) continue;
		UP_Close_index = Reads[Index].UP_Close.size() - 1;
		temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
		for (unsigned Pos_index = 0; Pos_index < Rest_Positions.size(); Pos_index++) {
			if (Reads[Index].MatchedD == Rest_Positions[Pos_index].Strand) {
				if (temp_AbsLoc == Rest_Positions[Pos_index].Pos) {
					Reads[Index].Used = true;
					Rest_Positions[Pos_index].Pos_Reads.push_back(Index);	 // copy index to save memory
				}
			}
		}
	}
	
	//cout << "Other unassigned breakpoints (BP): " << Rest_Positions.size() << "\n\n";
	int Count_BP = 0;
	bool temp_BalancedPlus, temp_BalancedMinus;
	short temp_LengthStr;
	vector <SPLIT_READ> temp_Pos_Reads;
	for (unsigned LI_index = 0; LI_index < Rest_Positions.size(); LI_index++) {
		temp_Pos_Reads.clear();
		temp_BalancedPlus = false;
		temp_BalancedMinus = false;
		for (int i = 0; i < Rest_Positions[LI_index].Pos_Reads.size(); i++) {
			temp_Pos_Reads.push_back(Reads[Rest_Positions[LI_index].Pos_Reads[i]]);
			UP_Close_index = Reads[Rest_Positions[LI_index].Pos_Reads[i]].UP_Close.size() - 1;
			temp_LengthStr = Reads[Rest_Positions[LI_index].Pos_Reads[i]].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength * 0.5)
				temp_BalancedPlus = true;
			else if ((float)temp_LengthStr < Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength * 0.5)
				temp_BalancedMinus = true;
		}	
		if (temp_BalancedPlus && temp_BalancedMinus) {
			Count_BP++;
			if (Rest_Positions[LI_index].Strand == Plus) {
				Outf_Rest << "########################################################" << endl;
				Outf_Rest << "ChrID " << temp_Pos_Reads[0].FragName << "\t" 
				<< Rest_Positions[LI_index].Pos - SpacerBeforeAfter + 1 << "\t" 
				<< Rest_Positions[LI_index].Pos_Reads.size() << "\t+" << endl;
				
				Outf_Rest <<  ( CurrentChr.substr(Rest_Positions[LI_index].Pos - ReportLength + 1, ReportLength) )
				<< Cap2Low( CurrentChr.substr(Rest_Positions[LI_index].Pos + 1, ReportLength ) ) << endl;
				HalfMappedIndex = 0;
				HalfUnmappedIndex = 0;
				for (int i = 0; i < temp_Pos_Reads.size(); i++) {
					UP_Close_index = temp_Pos_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
					if (temp_LengthStr > temp_Pos_Reads[HalfMappedIndex].UP_Close[temp_Pos_Reads[HalfMappedIndex].UP_Close.size() - 1].LengthStr)
						HalfMappedIndex = i;
					if (temp_Pos_Reads[i].ReadLength - temp_LengthStr > temp_Pos_Reads[HalfUnmappedIndex].ReadLength - temp_Pos_Reads[HalfUnmappedIndex].UP_Close[temp_Pos_Reads[HalfUnmappedIndex].UP_Close.size() - 1].LengthStr)
						HalfUnmappedIndex = i;
					for (int j = 0; j < ReportLength - temp_LengthStr; j++) {
						Outf_Rest << " ";
					}
					Outf_Rest << ReverseComplement(temp_Pos_Reads[i].UnmatchedSeq)
					<< "\t" << temp_Pos_Reads[i].MatchedD  
					<< "\t" << temp_Pos_Reads[i].MatchedRelPos 
					<< "\t" << temp_Pos_Reads[i].MS
					<< "\t" << temp_Pos_Reads[i].Tag 
					<< "\t" <<  temp_Pos_Reads[i].Name << endl;
					//temp_Pos_Reads[0].
					//BP_Reads.push_back();
				}
				
			}
			else {
				Outf_Rest << "########################################################" << endl;
				Outf_Rest << "ChrID " << temp_Pos_Reads[0].FragName << "\t" 
				<< Rest_Positions[LI_index].Pos - SpacerBeforeAfter + 1 << "\t" 
				<< Rest_Positions[LI_index].Pos_Reads.size() << "\t-" << endl;
				Outf_Rest << Cap2Low( CurrentChr.substr(Rest_Positions[LI_index].Pos - ReportLength, ReportLength) )
				<< ( CurrentChr.substr(Rest_Positions[LI_index].Pos, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Pos_Reads.size(); i++) {
					UP_Close_index = temp_Pos_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength + temp_LengthStr - temp_Pos_Reads[i].ReadLength; j++) {
						Outf_Rest << " ";
					}
					Outf_Rest << (temp_Pos_Reads[i].UnmatchedSeq)
					<< "\t" << temp_Pos_Reads[i].MatchedD  
					<< "\t" << temp_Pos_Reads[i].MatchedRelPos 
					<< "\t" << temp_Pos_Reads[i].MS
					<< "\t" << temp_Pos_Reads[i].Tag 
					<< "\t" <<  temp_Pos_Reads[i].Name << endl;				}					 
			}
		} 
	}
	delete[] plus_LI_Pos ;
	delete[] minus_LI_Pos ;
	cout << "Other unassigned breakpoints (BP): " << Count_BP << "\n\n";
}

short CompareTwoString(const string & Str_A, const string & Str_B) {
	short Str_Len = Str_A.size();
	short CompareResult;
	for (short i = 0; i < Str_Len; i++) {
		CompareResult = (short)Str_A[i] - (short)Str_B[i];
		if (CompareResult < 0) return -1;
		else if (CompareResult > 0) return 1;
	}
	return 0;
}

bool ReadInBamReads(const char* bam_path, const string & FragName, string * CurrentChr, vector <SPLIT_READ> & LeftReads, int InsertSize, string Tag, int binStart, int binEnd) {
	bamFile fp;
	fp = bam_open(bam_path, "r");
	assert(fp);
	bam_index_t *idx;
	idx = bam_index_load(bam_path); // load BAM index
	assert(idx);  
	bam_header_t* header = bam_header_read(fp);
	bam_init_header_hash(header);
	assert(header);
	//need thing that converts "tid" to "chromosome name" 
	int ret;
	int tid;
	unsigned long int sw_count=0;
	unsigned long int one_end_count=0;      
	tid = bam_get_tid(header, FragName.c_str());
	//kai does the below line in readinreads. dunno why yet
	//VectorTag.clear();
	
	
	fetch_func_data data;
	data.header=header;
	data.CurrentChr=CurrentChr;
	data.LeftReads=&LeftReads;
	data.read_to_map_qual=NULL;
	data.read_to_map_qual=kh_init(read_name); 
	flagshit b1_flags, b2_flags;
	data.b1_flags=&b1_flags;
	data.b2_flags=&b2_flags;
	data.InsertSize=InsertSize;
	data.Tag=Tag;
	//ADDITIONAL_MISMATCH = 1;
	//Seq_Error_Rate = 0.05;
	bam_fetch(fp, idx, tid, binStart , binEnd, &data, fetch_func);
	cout << "Bam:\t" << bam_path << "\tTag:\t" << Tag << endl;
	cout << "NumReadScanned:\t" << NumReadScanned << endl;
	cout << "NumReadInChr:\t" << NumReadInChr << endl;
	cout << "NumReadStored:\t" << LeftReads.size() << endl;
	cout << "NumReadStored / NumReadInChr = " << LeftReads.size() * 100.0 / NumReadInChr << " %\n" 
	<< "InChrPlus \t" << InChrPlus << "\tGetPlus \t" << GetPlus << "\t" << GetPlus * 100.0 / InChrPlus 
	<< " %\n" << "InChrMinus\t" << InChrMinus  << "\tGetMinus\t" << GetMinus 
	<< "\t" << GetMinus * 100.0 / InChrMinus << " %\n" << endl;
	khint_t key; 
	if(kh_size(data.read_to_map_qual) > 0) {
		for (key = kh_begin(data.read_to_map_qual); key != kh_end(data.read_to_map_qual); ++key) {
			if(kh_exist(data.read_to_map_qual,key)) {
				bam_destroy1(kh_value(data.read_to_map_qual, key));
				free((char *) kh_key(data.read_to_map_qual, key));
			}
		}
	}
	kh_clear(read_name, data.read_to_map_qual);
	kh_destroy(read_name, data.read_to_map_qual);
	
	NumReadScanned = 0;
	NumReadInChr = 0;
	InChrPlus = 0;
	InChrMinus = 0;
	GetPlus = 0;
	GetMinus = 0;
	
	
	//kai does the below in read inreads dunno why yet
	string Str4Exchange;
	if(VectorTag.size() > 0) {
		for (short i = 0; i < VectorTag.size() - 1; i++) {
			for (short j = 1; j < VectorTag.size(); j++) {
				if (VectorTag[i].size() > VectorTag[j].size()) {
					Str4Exchange = VectorTag[i];
					VectorTag[i] = VectorTag[j];
					VectorTag[j] = Str4Exchange; 
				}
				else if (VectorTag[i].size() == VectorTag[j].size()) {
					for (short k = 0; k < VectorTag[i].size(); k++) {
						if ((short)VectorTag[i][k] > (short)VectorTag[j][k]) {
							Str4Exchange = VectorTag[i];
							VectorTag[i] = VectorTag[j];
							VectorTag[j] = Str4Exchange;
							break;
						}
					}
				}
			}
		}
	}
	//end kai
	bam_header_destroy(header);
	bam_index_destroy(idx);
	bam_close(fp);
	return true;
}


static int fetch_func(const bam1_t *b1, void *data) {
	
	NumReadScanned++;
	fetch_func_data* data_for_bam=(fetch_func_data*)data;
	vector <SPLIT_READ>* LeftReads =(vector <SPLIT_READ>*) data_for_bam->LeftReads;
	bam_header_t* header = (bam_header_t*)data_for_bam->header;
	khash_t(read_name)* read_to_map_qual = (khash_t(read_name)*)data_for_bam->read_to_map_qual;
	flagshit* b1_flags = data_for_bam->b1_flags;
	flagshit* b2_flags = data_for_bam->b2_flags;
	const string CurrentChr = *(string*)data_for_bam->CurrentChr;
	
	SPLIT_READ Temp_One_Read;
	const bam1_core_t* b1_core;
	bam1_t *b2;
	bam1_core_t* b2_core;
	b1_core = &b1->core;
	string read_name = bam1_qname(b1); 
	//    if(!(b1_core->flag & BAM_FPROPER_PAIR)) { 
	//            return 0;
	//NO BUENO!        }
	khint_t key = kh_get(read_name, read_to_map_qual, bam1_qname(b1));
	if(key == kh_end(read_to_map_qual)) {
		int ret;
		key = kh_put(read_name, read_to_map_qual, strdup(bam1_qname(b1)), &ret);
		kh_value(read_to_map_qual, key) = bam_dup1(b1);
		return 0;
	}
	else {
		b2 = bam_dup1(kh_value(read_to_map_qual, key));
		bam_destroy1(kh_value(read_to_map_qual, key));
		b2_core = &b2->core;
		//this seems stupid, but in order to manage the read names, necessary
		free((char *) kh_key(read_to_map_qual, key));
		kh_del(read_name, read_to_map_qual, key);
		string c_sequence;
	}
	
	parse_flags_and_tags(b1, b1_flags);
	parse_flags_and_tags(b2, b2_flags);
	//read_name = bam1_qname(b1); 
	if(!(b1_flags->mapped) && !(b2_flags->mapped)) { 
		bam_destroy1(b2);
		return 0;
	}
	
	if(!b1_flags->unique && !b1_flags->sw && !b2_flags->unique && !b2_flags->sw) {
		bam_destroy1(b2);
		return 0;
	}
	
	if(b1_flags->unique && b2_flags->unique && b1_flags->edits ==0 && b2_flags->edits ==0 ) {
		
		bam_destroy1(b2);
		
		return 0;
	}
	if(b1_flags->suboptimal && b2_flags->suboptimal) {
		bam_destroy1(b2);
		return 0;
	}
	int max_edits1 = int(b1_core->l_qseq * .05) +1;
	int max_edits2 = int(b2_core->l_qseq * .05) +1;
	
	if(b1_flags->edits > max_edits1 && b2_flags->edits > max_edits2) {
		bam_destroy1(b2);
		return 0;
	}
	
	int which_read_is_mapped=0;
	if(b1_flags->mapped && b2_flags->mapped) {
		if((b1_flags->unique || b1_flags->sw) && b1_flags->edits <= max_edits1 && b2_flags->edits >0) {
			which_read_is_mapped=1;
		}
		else if((b2_flags->unique || b2_flags->sw) && b2_flags->edits <= max_edits1 && b1_flags->edits > 0) {
			which_read_is_mapped=2;
		}
		else {
			bam_destroy1(b2);
			return 0;
		}
	}else if (!b1_flags->mapped && b2_flags->edits <= max_edits2) {
		which_read_is_mapped=2;
	}else if (!b2_flags->mapped && b1_flags->edits <= max_edits1) {
		which_read_is_mapped=1;
	}
	else {
		bam_destroy1(b2);
		return 0;
	}
	if(which_read_is_mapped==1) {
		build_record(b1, b2, data);
		if((b2_flags->unique || b2_flags->sw) && b2_flags->edits <= max_edits2 && b1_flags->edits >0) {
			build_record(b2, b1, data);
		}
	}else if(which_read_is_mapped ==2) {
		build_record (b2,b1, data);
	}else {
		printf("Should not get here\n");
		exit(1);
	}
	bam_destroy1(b2);
	return 0;
}

/* 'isInBin' returns whether the read "read" is in the designated bin. */
bool isInBin(const SPLIT_READ &read)
{
	if (read.MatchedRelPos > g_maxPos ) { g_maxPos = read.MatchedRelPos; } 	
	return ( ( read.MatchedRelPos >= (g_binIndex*WINDOW_SIZE) ) &&
				( read.MatchedRelPos < ((g_binIndex+1)*WINDOW_SIZE) ) );
}

void build_record(const bam1_t* mapped_read, const bam1_t* unmapped_read, void* data ) {
	
	SPLIT_READ Temp_One_Read;
	fetch_func_data* data_for_bam=(fetch_func_data*)data;
	vector <SPLIT_READ>* LeftReads =(vector <SPLIT_READ>*) data_for_bam->LeftReads;
	bam_header_t* header = (bam_header_t*)data_for_bam->header;
	string CurrentChr = *(string*)data_for_bam->CurrentChr;
	string Tag = (string)data_for_bam->Tag;
	int InsertSize = (int)data_for_bam->InsertSize;
	
	const bam1_core_t* mapped_core;
	const bam1_core_t* unmapped_core;
	mapped_core= &mapped_read->core;
	unmapped_core = &unmapped_read->core;
	Temp_One_Read.Name="@";
	Temp_One_Read.Name.append((const char*)bam1_qname(unmapped_read));
	if(unmapped_core->flag & BAM_FREAD1) {
		Temp_One_Read.Name.append("/1");
	}
	else if(unmapped_core->flag & BAM_FREAD2) {
		Temp_One_Read.Name.append("/2");
	}
	string c_sequence;
	int i;
	uint8_t *s = bam1_seq(unmapped_read);
	for (i = 0; i < unmapped_core->l_qseq; ++i) c_sequence.append(1, bam_nt16_rev_table[bam1_seqi(s, i)]);
	//rudimentary n filter
	int length=unmapped_core->l_qseq;
	while(c_sequence[0] == 'N') {
		c_sequence.erase(0, 1);
		length--;
	}
	if(c_sequence.size() > 0) {
		while(c_sequence[length-1] == 'N') {
			c_sequence.erase(length-1, 1);
			length--;
		}
	}
	int n_count=0, found=c_sequence.find('N', 0);
	int max_ns = length * .10;
	while(found != string::npos) {
		n_count++;
		found = c_sequence.find('N', found+1);
	}
	if(n_count > max_ns || length < 22) {
		return;
	}
	//rudimentary n filter end
	Temp_One_Read.ReadLength=length;
	Temp_One_Read.ReadLengthMinus=length - 1;
	if(unmapped_core->flag & BAM_FREVERSE) {
		Temp_One_Read.UnmatchedSeq=ReverseComplement(c_sequence);
	}else {
		Temp_One_Read.UnmatchedSeq=c_sequence;
	}
	Temp_One_Read.MatchedRelPos=mapped_core->pos;
	if(mapped_core->flag & BAM_FREVERSE) {
		Temp_One_Read.MatchedD='-';
		uint32_t* cigar_pointer= bam1_cigar(mapped_read); 
		//reusing length to be something else now. eat it.
		length = bam_cigar2len(mapped_core, cigar_pointer);
		Temp_One_Read.MatchedRelPos+=length;
	}    else {
		Temp_One_Read.MatchedD='+';
	}
	Temp_One_Read.MS=mapped_core->qual;
	//FIXME pass these through from the command line with a struct
	Temp_One_Read.InsertSize=InsertSize;
	Temp_One_Read.Tag=Tag; 
	if(((mapped_core->tid == unmapped_core->tid)&&
		 (mapped_core->pos != unmapped_core->pos)) &&
		(abs(mapped_core->isize) < Temp_One_Read.InsertSize)){
		if(Temp_One_Read.MatchedD=='+'){
			Temp_One_Read.MatchedRelPos-=Temp_One_Read.InsertSize;
		} else {
			Temp_One_Read.MatchedRelPos+=Temp_One_Read.InsertSize;
		}
	}
	
	string FragName=header->target_name[mapped_core->tid];
	Temp_One_Read.FragName=FragName;
	//cout << Temp_One_Read.Name << "\n" << Temp_One_Read.UnmatchedSeq << "\n";
	//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.FragName << "\t" << Temp_One_Read.MatchedRelPos << "\t" << Temp_One_Read.MS << "\t" << Temp_One_Read.InsertSize << "\t" << Temp_One_Read.Tag << "\n";
	NumReadInChr++;
	Temp_One_Read.MAX_SNP_ERROR = (short)(Temp_One_Read.UnmatchedSeq.size() * Seq_Error_Rate);
	Temp_One_Read.TOTAL_SNP_ERROR_CHECKED = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
	Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
	Temp_One_Read.MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;//atoi(argv[1]);
	//MinFar_I = MinClose + 1;//atoi(argv[2]);
	if (Temp_One_Read.MatchedD == Plus)  {
		InChrPlus++;
		//Temp_One_Read.UnmatchedSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
	}
	else InChrMinus++;
	if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size) Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
	if (Temp_One_Read.MatchedRelPos < 1) Temp_One_Read.MatchedRelPos = 0;
	GetCloseEnd(CurrentChr, Temp_One_Read);
	if (Temp_One_Read.UP_Close.size()) {
		if (ReportLength < Temp_One_Read.ReadLength) ReportLength = Temp_One_Read.ReadLength;
		Temp_One_Read.Used = false;
		Temp_One_Read.Unique = true;
		//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t";
		CleanUniquePoints(Temp_One_Read.UP_Close);
		//cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl;
		Temp_One_Read.CloseEndLength = 
		Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size() - 1].LengthStr;
		
		//if (Temp_One_Read.UP_Close.size()) {
		
		if (Temp_One_Read.MatchedD == Plus) {
			Temp_One_Read.LeftMostPos = Temp_One_Read.UP_Close[0].AbsLoc + 1 - Temp_One_Read.UP_Close[0].LengthStr;
			GetPlus++;
		}
		else {
			Temp_One_Read.LeftMostPos = Temp_One_Read.UP_Close[0].AbsLoc + Temp_One_Read.UP_Close[0].LengthStr - Temp_One_Read.ReadLength;
		   GetMinus++;	
		}
		
		//if (isInBin(Temp_One_Read)) {
			LeftReads->push_back(Temp_One_Read);
		//}
		if(LeftReads->size() % 10000 == 0) { cout << LeftReads->size() << endl; }
	}
	if (NotInVector(Temp_One_Read.Tag, VectorTag)) {
		VectorTag.push_back(Temp_One_Read.Tag);
	}
	return;   
}        

void parse_flags_and_tags(const bam1_t *b, flagshit* flags) {
	const bam1_core_t* c = &b->core;
	char xt_code = 0;
	int mf_code=0, nm_code=0, best_hits=0, sub_hits=0;
	flags->unique=0;
	flags->mapped = !(c->flag & BAM_FUNMAP);
	
	uint8_t *s = bam_aux_get(b, "XT");
	if(s!=0) {
		xt_code = bam_aux2A(s);
		if(xt_code == 'U') {
			flags->unique=1;
		}
		else {
			flags->unique=0;
		}
		
	}
	s = NULL;
	s = bam_aux_get(b, "MF");
	if(s!=0) {
		mf_code = bam_aux2i(s);
		//the below paradigm doesn't exist in maq i think we should assume the read is unique so this hack is here
		if(mf_code != 130) {
			flags->unique=1;
		}
		else {
			flags->unique=0;
		}
		flags->suboptimal=0;
	}
	s = NULL;
	s = bam_aux_get(b, "X0");
	if(s!=0) {
		best_hits = bam_aux2i(s);
	}
	s = NULL;
	s = bam_aux_get(b, "X1");
	if(s!=0) {
		sub_hits = bam_aux2i(s);
		
		if(best_hits+sub_hits == 1 ) {
			flags->suboptimal=0;
		}
		else {
			flags->suboptimal=1;
		}
	}    
	if(xt_code == 'M' || mf_code == 130 ) {
		flags->sw=1;
		//short term fix to unset unique if the maq read was s-w mapped. bwa can't set U and M at once.
	}
	else {
		flags->sw=0;
	}
	s = NULL;
	s = bam_aux_get(b, "NM");
	if(s!=0) {
		nm_code = bam_aux2i(s);
		flags->edits=nm_code;
	}
	else {
		nm_code=0;
		flags->edits=nm_code;
	}
	
	return; 
}   


int32_t bam_cigar2len(const bam1_core_t *c, const uint32_t *cigar)
{
	uint32_t k;
	int32_t l = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP)
			l += cigar[k] >> BAM_CIGAR_SHIFT;
		if (op == BAM_CDEL) 
			l -= cigar[k] >> BAM_CIGAR_SHIFT;
	}
	return l;
}

static void show_flag_info(flagshit* flags1, flagshit* flags2) {
	printf("Read One flags: Mapped:%d Unique:%d sw:%d edits:%d suboptimal:%d\n", flags1->mapped, flags1->unique, flags1->sw, flags1->edits, flags1->suboptimal);
	printf("Read Two flags: Mapped:%d Unique:%d sw:%d edits:%d suboptimal:%d\n", flags2->mapped, flags2->unique, flags2->sw, flags2->edits, flags2->suboptimal);
}


