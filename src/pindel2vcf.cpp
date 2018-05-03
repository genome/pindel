/* vcfcreator.cpp

 Transforms indel and SV-files created by Pindel into VCF4.0 format according to 1000-genome SV specifications.
 This version sorts the records before outputting; may need to make version that uses external memory later!

 Created by Eric-Wubbo Lameijer, section of Molecular Epidemiology, Leiden University Medical Center, March 3rd, 2011.
 e.m.w.lameijer@gmail.com
 +31(0)71-5 125 831

 Version 0.6.3 [February 19th, 2014] Clearer text on usage of -P option
 Version 0.6.2 [December 12th, 2014] Now robust against fasta files that have non-standard line lengths (C++'s getline does not work well on lines of over a million characters)
 Version 0.6.1 [December 12th, 2014] Now has special code to recognize lines that contain SV-data, instead of relying on indirect establishment of their identity from context
		(for example: seeing that the previous line starts with '#', which fails for a grepped Pindel output file.)
 Version 0.6.0 [April 21nd, 2014] automatic warning when finding oversized SVs - and set -co by default as a safeguard. Thanks to user Hong Ching Lee for reporting this problem.
 Version 0.5.9 -ho -he for cutoffs
 Version 0.5.8 [April 28th, 2013] Added "-co" option for compact output; -co followed by an integer indicates
 the longest variation displayed with full base sequence. If, for example, the user gives as command line parameter
 "-co 10", then if either the reference or the alternative allele exceeds 10 bases, the output will be transformed
 into the format chrom pos firstrefbase <SVType>.
 Version 0.5.7 [February 1st, 2013] Now updated my contact information.
 Version 0.5.6 [February 1st, 2013]  Quickfix since user detected empty labels. Note: since no input was delivered that enabled me to reproduce the bug, I could not find the true cause. So problems may still crop up in the future.
 Version 0.5.5 [December 10th, 2012] Modified the code so that LI / -P now gives correct GT:AD instead of GT:RD:AD output. Also fixed small bug in creating ALT of NT-inversions
 Version 0.5.4 [December 6th, 2012] Error found by David Hannah on using the -P option debugged
 Version 0.5.3 [November 8th, 2012] Now should genotype newer pindel output as 0/0, 1/1 or 0/1 based on apparent balance between alleles
 Version 0.5.2 [November 8th, 2012] Bugs removed from 0.5.1 and added GT:RD:AD output to show reference coverage
 Version 0.5.1 [October 31st, 2012] Now also compatible with Pindel's new output format
 Version 0.5.0 [July 9th, 2012] removed small error that caused the -c option to process all other chromosomes as well...
 Version 0.4.9 [July 2nd, 2012] added -P option to allow fusing all output files of one pindel run
 Version 0.4.8 [July 2nd, 2012] debugged .:10 genotype (due to 'dual-encoding' of genotype, while Pindel not having actual chromosome data
 Version 0.4.7 [June 21th, 2012] LI now also have proper labels <as per updated Pindel>
 Version 0.4.6 [June 20th, 2012] END-position now proper according to http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41, access date June 20, 2012
 Version 0.4.5 [June 20th, 2012] displays proper warnings when using the -G option
 Version 0.4.4 [June 19th, 2012] now also prints inversions correctly in GATK-format (when using the -G option)
 Version 0.4.3 [June 19th, 2012] adding equilength replacement calls to -G output (2)
 Version 0.4.2 [June 19th, 2012] adding equilength replacement calls to -G output (1)
 Version 0.4.1 [June 18th, 2012] Now automatically sets the end of long insertions to occur after the beginning...
 Version 0.4.0 [June 18th, 2012] Added -G option to make output GATK-compatible
 Version 0.3.9 [June 8th, 2012] Added maximum total coverage to protect from duplicated regions
 Veriosn 0.3.8 [June 8th, 2012] Long insertions support debugged.
 Version 0.3.7 [June 8th, 2012] Next stage in debugging repeats/postindel options
 Version 0.3.6 [June 8th, 2012] Debugged int repeats/postindel options[1]
 Version 0.3.5 [June 1st, 2012] Refined microhomology/microsattelite sequencing to distinguish internal repeats and 'postindel' repeats
 Version 0.3.4 [May 21st, 2012] Alternative homology calling; now checks whether inserted/deleted sequence is part of a longer repetitive sequence
 Version 0.3.3 [May 2nd, 2012] Debugged -sb/-ss option: now works if 1 sample is selected
 Version 0.3.2 [May 1st, 2012] Adds -sb and -ss options to allow users to only count samples with sufficient individual support.
 Version 0.3.1 [March 27th, 2012] -c option now works!
 Version 0.3.0 [March 27th, 2012] now also works correctly for small inversions (as NT is then only "TG", gave wrong code, ACA ATGTGTG instead of ACA ATG
 Version 0.2.9 [March 27th, 2012] Adds an option for counting samples only if they are balanced with at least N reads on each side
 Version 0.2.8 [March 2nd, 2012] Now does not skip chromosome when compiling long insertions
 Version 0.2.7 [March 2nd, 2012] Skip chromosomes in which Pindel has not called any SVs
 Version 0.2.6 [March 2nd, 2012] Fixed problem in which chromosome was not properly freed from memory
 Version 0.2.5 [March 2nd, 2012] Fixed segmentation fault error when there are no reads in a chromosome
 Version 0.2.4 [February 6th, 2012] Can save memory at the cost of time by implementing a window-size option.
 Version 0.2.3 [February 6th, 2012] More memory-efficient by not explicitly storing REF and ALT strings, but by using recalculation instead of storage
 Version 0.2.1 [January 3rd, 2012] Now also correctly reads in fasta files which do not only contain uppercase characters.
 Version 0.2.0 [November 10th, 2011. The support of an indel, previously called "DP" is now more appropriately called "AD".
 Also, genotypes are somewhat more correctly indicated with . and 1/.
 Also, replaces -1 by . to indicate unknown number of fields in the declarations in the header
 Version 0.1.9 [August 19th, 2011. To save memory, now reads in chromosomes only when needed, so doesn't put the entire genome in memory at once.
 */

/* CONTENTS

 PREFACE: include files and global constants
 CHAPTER 1. General utilities for DNA-string manipulation (reverse-complementing a string)
 CHAPTER 2. Defining the parameters and the 'Parameter class' to handle them.
 */


/*** PREFACE: include files and global constants ->***/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <vector>

#include <ctype.h> // tolower
#include <math.h> //sqrt
#include <stdlib.h> // for atoi
#include <time.h>

#include <stdio.h>

#include <stdint.h>
#include <inttypes.h>


const int FIRST_SAMPLE_INDEX = 32; // index of first sample name

using namespace std;

string g_versionString = "0.6.3";
string g_programName = "pindel2vcf";

bool g_normalBaseArray[256];
int g_sizeToWarnFor = 1000000;

// ################################
// log\binom{n}{k}
static double lbinom(int n, int k)
{
   if (k == 0 || n == k) {
      return 0;
   }
   return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
static double hypergeo(int n11, int n1_, int n_1, int n)
{
   return exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
}

typedef struct {
   int n11, n1_, n_1, n;
   double p;
} hgacc_t;

// incremental version of hypergenometric distribution
static double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux)
{
   if (n1_ || n_1 || n) {
      aux->n11 = n11;
      aux->n1_ = n1_;
      aux->n_1 = n_1;
      aux->n = n;
   } else { // then only n11 changed; the rest fixed
      if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
         if (n11 == aux->n11 + 1) { // incremental
            aux->p *= (double)(aux->n1_ - aux->n11) / n11
                      * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
            aux->n11 = n11;
            return aux->p;
         }
         if (n11 == aux->n11 - 1) { // incremental
            aux->p *= (double)aux->n11 / (aux->n1_ - n11)
                      * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
            aux->n11 = n11;
            return aux->p;
         }
      }
      aux->n11 = n11;
   }
   aux->p = hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
   return aux->p;
}

double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two)
{
   int i, j, max, min;
   double p, q, left, right;
   hgacc_t aux;
   int n1_, n_1, n;

   n1_ = n11 + n12;
   n_1 = n11 + n21;
   n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
   max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
   min = n1_ + n_1 - n;
   if (min < 0) {
      min = 0;   // min n11, for left tail
   }
   *two = *_left = *_right = 1.;
   if (min == max) {
      return 1.;   // no need to do test
   }
   q = hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
   // left tail
   p = hypergeo_acc(min, 0, 0, 0, &aux);
   for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) { // loop until underflow
      left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
   }
   --i;
   if (p < 1.00000001 * q) {
      left += p;
   } else {
      --i;
   }
   // right tail
   p = hypergeo_acc(max, 0, 0, 0, &aux);
   for (right = 0., j = max - 1; p < 0.99999999 * q; --j) { // loop until underflow
      right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
   }
   ++j;
   if (p < 1.00000001 * q) {
      right += p;
   } else {
      ++j;
   }
   // two-tail
   *two = left + right;
   if (*two > 1.) {
      *two = 1.;
   }
   // adjust left and right
   if (abs(i - n11) < abs(j - n11)) {
      right = 1. - left + q;
   } else {
      left = 1.0 - right + q;
   }
   *_left = left;
   *_right = right;
   return q;
}



// ################################



/* the global parameter g_par stores the values of all parameters set by the user; the Parameter class, in contrast, is for user-friendly IO. */
struct ParameterSettings {
   string reference;
   string referenceName;
   string referenceDate;
   string pindelfile;
   string pindelroot;
   string vcffile;
   string chromosome;
   int windowSize;
   int MinCoverage;

   double HetCutoff;
   double HomCutoff;
   int minsize;
   int maxsize;
   bool bothstrands;
   int minsuppSamples;
   int minsuppReads;
   int maxSuppReads;
   int regionStart;
   int regionEnd;
   int maxInterRepeatNo;
   int maxInterRepeatLength;
   int maxPostRepeatNo;
   int maxPostRepeatLength;
   bool onlyBalancedSamples;
   int minimumStrandSupport;
   int compactOutput;
   bool showHelp;
   bool somatic;
   bool gatkCompatible;
} g_par;

bool pindel024uOrLater = false;

/*** END OF PREFACE: include files and global constants ***/

/*** CHAPTER 1. General utilities for DNA-string manipulation (reverse-complementing a string) ->***/

/* returns the complementary DNA-base of base 'inputbase' */
char complementBase( char inputBase )
{
   switch (inputBase) {
   case 'A':
      return 'T';
   case 'C':
      return 'G';
   case 'G':
      return 'C';
   case 'T':
      return 'A';
   default:
      return 'N';
   }
}

/* 'createComplement' creates the complement of a DNA string. */ // HEY!!! don't forget to complement CpG's properly!
void createComplement( const string& dna, string& complement )
{
   int dnaLength = dna.length();
   complement = "";
   for (int position=dnaLength-1; position>=0; position-- ) {
      complement += complementBase( dna[ position ] );
   }
}

/** 'InputReader' can house a vector of files, allowing access as if it were one huge file. */
class InputReader
{

public:
   InputReader();

   string getLine();
   bool eof();
   void addFile(const string filename);
   void rewind();


private:
   vector<string> m_filenames;
   int m_nextFileIndex;
   bool m_readable;
   ifstream m_currentFile;

   bool canReadMore();
   void moveToNextFile();
};

string InputReader::getLine()
{
   if (canReadMore()) {
      string line;
      std::getline( m_currentFile, line );
      return line;
   } else {
      return "";
   }
}

bool InputReader::canReadMore()
{
   // default case: current file is okay
   if (m_currentFile && !m_currentFile.eof()) {
      return true;
   }

   while (!m_currentFile || m_currentFile.eof()) {
      moveToNextFile();
      if (!m_readable) {
         break;   // final EOF
      }
   }

   return m_readable;
}

void InputReader::moveToNextFile()
{
   if (m_nextFileIndex<m_filenames.size()) {
      m_currentFile.close();
      m_currentFile.open( m_filenames[ m_nextFileIndex ].c_str() );
      m_nextFileIndex++;
   } else {
      m_readable = false;
   }
}



void InputReader::rewind()
{
   m_currentFile.open("");
   m_nextFileIndex = 0;
   m_readable = true;
}

InputReader::InputReader()
{
   rewind();
}

void InputReader::addFile(const string filename)
{
   m_filenames.push_back( filename );
}

bool InputReader::eof()
{
   return !canReadMore();
}

/* 'Parameter' stores an individual parameter; it is set by the command line parameters, and used by the program. */
class Parameter
{
public:
   bool isRequired() const
   {
      return d_required;
   }
   void describe() const;
   string getDescription() const
   {
      return d_description;
   }
   string getShortName() const
   {
      return d_shortName;
   }
   string getLongName() const
   {
      return d_longName;
   }
   bool hasName( const string& name ) const
   {
      return ( d_shortName.compare(name) == 0 || d_longName.compare( name ) == 0 );
   }
   bool isSet() const
   {
      return d_isSet;
   }
   virtual void setValue(const string& value )
   {
      cout << "WHAT!" << endl ;
   };
   virtual void setValue(const int value) {};
   virtual void setValue(const double value) {};
   virtual void setValue(const bool value) {};
   virtual int getIValue() const {};
   virtual bool getBValue() const {};
   virtual string getSValue() const {};
   virtual double getFValue() const {};

   virtual bool isUnary() const
   {
      return false;
   }

   Parameter( const string& shortName, const string& longName, const string& description, const bool required);
protected:
   void set()
   {
      d_isSet = true;
      //cout << "setting " << d_shortName << endl;
   }


private:
   bool d_required;
   string d_shortName;
   string d_longName;
   string d_description;
   bool d_isSet;
};

class IntParameter: public Parameter
{
public:
   IntParameter( int* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const int value );

   virtual int getIValue() const
   {
      return *d_data_ptr;
   }
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

class BoolParameter: public Parameter
{
public:
   BoolParameter( bool* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const bool value );

   virtual bool getBValue() const
   {
      return *d_data_ptr;
   }
   virtual void setValue(const string& value);
   virtual void setValue(const bool value);
   virtual bool isUnary() const
   {
      return true;
   }
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

class FloatParameter: public Parameter
{
public:
   FloatParameter( double* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const double value );

   double getFValue() const
   {
      return *d_data_ptr;
   }
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

class StringParameter: public Parameter
{
public:
   StringParameter( string* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const string& value );

   string getSValue() const
   {
      return *d_data_ptr;
   }
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
   cout << d_shortName << "/" << d_longName << "  " << d_description;
   if ( d_required ) {
      cout << ": required parameter" ;
   }
   cout << endl;
}

vector<Parameter*> parameters;

/** 'convertToUppercase' returns the input string in full uppercase. */
string convertToUppercase( const string& inputString)
{
   string outputString = inputString;
   for (int i=0; i<outputString.length(); i++ ) {
      outputString[ i ] = toupper( outputString[ i ] );
   }
   return outputString;
}

bool normalBase(char ch )
{
   return g_normalBaseArray[ ch ];
}

void makeStrangeBasesN(string& dna)
{
   int chromLength = dna.size();
   for (int position=0; position<chromLength; position++ ) {
      if (!normalBase(dna[ position ])) {
         dna[ position ] ='N';
      }
   }
}

/* 'Chromosome' contains a string identifier as well as the base sequence of the chromosome itself. */
class Chromosome
{
public:
   Chromosome(const string& identifier, const string& fastaFilename)
   {
      d_identifier = identifier;
      d_sequence=NULL;
      d_fastaFilename=fastaFilename;
   }
   ~Chromosome()
   {
      delete d_sequence;
   }
   const string* getChromPtr();
   const string getID() const
   {
      return d_identifier;
   }

   void removeFromMemory()
   {
      cout << "Removing chromosome " << d_identifier << " from memory.\n";
      delete d_sequence;
      d_sequence=new string("");
   }

private:
   void readFromFile();
   string d_identifier;
   string* d_sequence;
   string d_fastaFilename;
};

/* 'Chromosome::readFromFile' reads in the reference chromosome sequence from the FASTA file. */
void Chromosome::readFromFile()
{
   ifstream referenceFile( d_fastaFilename.c_str() );
   referenceFile.clear(); // reset file for reading from the start again
   referenceFile.seekg(0);
   if (referenceFile.fail()) {
      cout << "Cannot open reference file. Exiting.\n";
      exit(EXIT_FAILURE);
   }
   string refLine, refName, currentLine;
   string tempChromosome = "";

   std::getline(referenceFile,refLine); // FASTA format always has a first line with the name of the reference in it
   // loop over each chromosome
   bool targetChromosomeRead = false;
   do {
      int counter=1;
      refName = "";
      do {
         refName += refLine[ counter++ ];
      } while ( counter<refLine.size() && (refLine[ counter ] != ' ') && (refLine[ counter ] != '\t') && (refLine[ counter ] != '\n') );

      if (refName == d_identifier ) {
         cout << "Reading chromosome " << refName << " into memory." << endl;
         targetChromosomeRead = true;
         tempChromosome+="N"; // for 1-shift
      }
      std::getline(referenceFile,currentLine);
      while (!referenceFile.eof() && currentLine[0]!='>') {
         if (refName == d_identifier ) {
            tempChromosome += convertToUppercase( currentLine );
         }
         std::getline(referenceFile,currentLine);
      }
      makeStrangeBasesN(tempChromosome);
      d_sequence = new string( tempChromosome );
      refLine = currentLine;
   } while (!referenceFile.eof() && !targetChromosomeRead);
}

const string* Chromosome::getChromPtr()
{
   if (d_sequence == NULL) { // sequence not read into memory
      readFromFile();
      //cout << "Have read " << d_sequence;
   }
   return d_sequence;
}





/* 'Genome' contains a collection of chromosomes, and returns a pointer to the requested chromosome (the chromosome with the requested ID) */
class Genome
{
public:
   const string* getChromosome( const string& id );
   void addChromosome( const Chromosome& chr )
   {
      d_chromosomes.push_back( chr );
   }
   string firstChromosomeName();
   string nextChromosomeName();

   vector<Chromosome> d_chromosomes; // working fast here, but ideally I'd keep this private and implement an iterator
};


const string* Genome::getChromosome( const string& id )
{

   for (int chromosomeIndex=0; chromosomeIndex<d_chromosomes.size(); chromosomeIndex++ ) {
      if ( id.compare( d_chromosomes[ chromosomeIndex ].getID() ) == 0 ) {
         return d_chromosomes[ chromosomeIndex ].getChromPtr();
      }
   }
   return NULL; // default value if chromosome not found
}


/* 'createHeader' writes a VCF4.0 (1k genome SV) compatible header to output VCF-file "outFile". */
void createHeader(ofstream &outFile, const string& sourceProgram, const string& reference, const set<string>& samples)
{
   // file format
   outFile << "##fileformat=VCFv4.0\n";

   // date of file creation
   time_t rawtime;
   struct tm * timeinfo;
   time ( &rawtime );
   timeinfo = localtime ( &rawtime );
   outFile << "##fileDate=" << g_par.referenceDate << endl;

   // source
   outFile << "##source=" << sourceProgram << endl;

   // reference
   outFile << "##reference=" << reference << endl;

   // info fields (selected from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20%28Variant%20Call%20Format%29%20version%204.0/encoding-structural-variants)
   outFile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
   outFile << "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">" << endl;
   outFile << "##INFO=<ID=PF,Number=1,Type=Integer,Description=\"The number of samples carry the variant\">" << endl;
   outFile << "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">" << endl;
   outFile << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
   outFile << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
   outFile << "##INFO=<ID=NTLEN,Number=.,Type=Integer,Description=\"Number of bases inserted in place of deleted code\">" << endl;
   outFile << "##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">" << endl;
   //outFile << "##ALT=<ID=DEL,Description=\"Deletion\">" << endl; /*EWL040311: probably not needed, as our calls are precise so we should rather give the exact replacing sequence instead of a label. */
   //outFile << "##ALT=<ID=DUP,Description=\"Duplication\">" << endl;
   //outFile << "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">" << endl;
   //outFile << "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">" << endl;
   //outFile << "##ALT=<ID=INV,Description=\"Inversion\">" << endl;
   //outFile << "##ALT=<ID=CNV,Description=\"Copy number variable region\">" << endl;
   outFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
   if (pindel024uOrLater) {
      outFile << "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference depth, how many reads support the reference\">" << endl;
   }
   outFile << "##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allele depth, how many reads support this allele\">" << endl;

   // headers of columns
   outFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
   if (samples.size()>0) {
      outFile << "\tFORMAT";
      for (set<string>::iterator counter=samples.begin(); counter!=samples.end(); counter++ ) {
         outFile << "\t" << *counter;
      }
   }
   outFile << "\n";
}

/* 'Pair' contains a pair of strings, outputted as "str1:str2". e*/
class Pair
{
   friend ostream& operator<<(ostream& os, const Pair& pair);

public:
   Pair(const string& first, const string& second)
   {
      d_first=first;
      d_second=second;
   }


private:
   string d_first, d_second;

};

ostream& operator<<(ostream& os, const Pair& pair)
{
   os << pair.d_first << "," << pair.d_second;
}


/* 'Genotype' stores the genotype of a mutation and its read depth. So, for example, "0/1:5". */
class Genotype
{

   friend ostream& operator<<(ostream& os, const Genotype& gt);

public:
   Genotype();
   Genotype( const int readDepthPlus, const int readDepthMinus, const int totalRefSupport );
   void fuse( const Genotype& gt );
   void reset()
   {
      d_readDepthPlus=0;
      d_readDepthMinus=0;
      d_totalRefSupport = 0;
   }
   int getReadDepthPlus() const
   {
      return d_readDepthPlus;
   }
   int getReadDepthMinus() const
   {
      return d_readDepthMinus;
   }
   int getTotalReads() const
   {
      return ( d_readDepthPlus + d_readDepthMinus ) ;
   }

   int getAverageRefSupport() const
   {
      return d_totalRefSupport;
   }

   int getTotalRefSupport() const
   {
      return d_totalRefSupport;
   }
   const string getGTRDAD() const;
   const string getGTAD() const;

private:
   int d_readDepthPlus;
   int d_readDepthMinus;
   int d_totalRefSupport;
   const string getGTold() const;
   const string getGTnew() const;
};

Genotype::Genotype()
{
   Genotype( 0, 0, 0 );
}

Genotype::Genotype( const int readDepthPlus, const int readDepthMinus, const int totalRefSupport )
{
   d_readDepthPlus = readDepthPlus;
   d_readDepthMinus = readDepthMinus;
   d_totalRefSupport = totalRefSupport;
}

void Genotype::fuse( const Genotype& gt )
{
   d_readDepthPlus += gt.d_readDepthPlus;
   d_readDepthMinus += gt.d_readDepthMinus;
   d_totalRefSupport = max(d_totalRefSupport, gt.d_totalRefSupport);
}

/** 'balanced' returns whether the a and b seem derived from a normal genomic binomial distribution of a and b, chance should be over 95% that the distribution is believable given heterozygosy (theoretically, as many
 reference strands as event strands. */
bool balanced( const unsigned int a, const unsigned int b)
{
   if (a==b) {
      return true;
   }
   unsigned int smallest = min( a, b );
   unsigned int largest = max( a, b);
   unsigned int sum = a + b;
   if (sum>22) {
      return (smallest >= (sum/2) - sqrt( sum ));
   } else if (sum>19) {
      return ( smallest >= 6 );
   } else if (sum>16) {
      return ( smallest >= 5 );
   } else if (sum>14) {
      return ( smallest >= 4 );
   } else if (sum>11) {
      return ( smallest >= 3 );
   } else if (sum>8) {
      return ( smallest >= 2 );
   } else if (sum>5) {
      return ( smallest >= 1 );
   } else {
      return true;
   }
}

string deriveGenotype( const Genotype& rawGenotype )
// we may want to replace this simple logic by a function that assesses if reference support and event support are balanced, like 'isBalanced( eventSupp, refSupp );
{
   int totalEventReads = rawGenotype.getTotalReads();
   int totalReferenceReads = rawGenotype.getAverageRefSupport();

   if (totalEventReads + totalReferenceReads < g_par.MinCoverage) {
      return "0/0";
   }
   float AF = (float)totalEventReads / (totalEventReads + totalReferenceReads);
   if (AF < g_par.HetCutoff) {
      return "0/0";
   } else if (AF >= g_par.HetCutoff && AF < g_par.HomCutoff) {
      return "0/1";
   } else if (AF >= g_par.HomCutoff) {
      return "1/1";
   }

   /*
     if (totalEventReads == 0 ) {
     if (totalReferenceReads == 0) {
     return ".";
     }
     else {
     return "0/0";
     }
     }
     else {
     if (totalEventReads > totalReferenceReads) {
     if ( balanced( totalReferenceReads, totalEventReads )) {
     return "0/1";
     }
     else {
     return "1/1";
     }
     }
     else if ( totalEventReads < totalReferenceReads ) {
     if ( balanced( totalEventReads, totalReferenceReads )) {
     return "0/1";
     }
     else {
     return "0/0";
     }
     }
     else {
     return "0/1";
     }
     }
     */
}


const string Genotype::getGTold() const
{
   if (g_par.gatkCompatible) {
      if (d_readDepthPlus==0 && d_readDepthMinus==0) {
         return "0/0";
      } else { // at least one read depth is 1
         return "0/1";
      }
   } else {
      if (d_readDepthPlus==0 && d_readDepthMinus==0) {
         return ".";
      } else { // at least one genotype is 1
         return "1/.";
      }
   } // gatk compatibility not required

   return "ERROR"; // this should not happen
}

const string Genotype::getGTnew() const
{
   return deriveGenotype( *this );
}

const string Genotype::getGTRDAD() const
{
   stringstream ss;
   ss << getGTnew() << ":" << getTotalRefSupport() << "," << getTotalReads();
   return ss.str();
}

const string Genotype::getGTAD() const
{
   stringstream ss;
   ss << getGTold() << ":" << getTotalReads();
   return ss.str();
}


/* 'SVData' stores the data of a certain structural variant. */
class SVData
{

   friend ostream& operator<<(ostream& os, const SVData& svd);

public:
   SVData(const int genotypeTotal);

   int getPosition() const
   {
      // attempted workaround for GATK undocumented feature
      if (g_par.gatkCompatible && altSameLengthAsRef()) {
         return d_position+1;
      } else {
         return d_position;
      }
   }
   int getSize() const
   {
      return d_svlen;
   }
   bool bothStrands() const;
   int getNumSupportSamples(const bool onlyBalancedSamples, const int minimumStrandSupport) const;
   int getNumSupportReads() const;

   string getChromosome() const
   {
      return d_chromosome;
   }

   void setChromosome(const string& chromosome)
   {
      d_chromosome = chromosome;
   }
   void setPosition(const int position)
   {
      d_position = position;
   }
   void setID(const string& id)
   {
      d_id = id;
   }

   void setQuality(const double quality)
   {
      ostringstream strs;
      strs << quality;
      d_quality = strs.str();
   }
   void setFilter(const string& filter)
   {
      d_filter = filter;
   }

   void setEnd(const int end)
   {
      //cout << "Setting end to " << end << endl;

      d_end = end;
   }
   int getVCFPrintEnd() const
   {
      if (d_end <= d_position ) {
         //cout << "Warning: end position of the SV (" << d_end << ") appears to be before the startposition (" << d_position << "). Will adjust end to be startposition+reflen-1.\n";
      } else {} // empty else.
      return d_position+getReference().size()-1;
   }

   void setHomlen(const int homlen)
   {
      d_homlen = homlen;
   }
   void setHomseq(const string& homseq)
   {
      d_homseq = homseq;
   }
   void setSVlen(const int svlen)
   {
      d_svlen = svlen;
   }
   void setSVlen(const string& svlen)
   {
      d_svlen = atoi(svlen.c_str());
   }
   void setSVtype(const string& svtype)
   {
      d_svtype = svtype;
   }

   //void setSupportingReads(const int supportingreads) { d_supportingreads = supportingreads; }
   void addGenotype(const int sampleID, const int readDepthPlus, const int readDepthMinus, const int totalRefSupport)
   {
      Genotype gt(readDepthPlus,readDepthMinus, totalRefSupport);
      d_format[ sampleID ] = gt ;
   }

   void setReplace( int replaceLength, int secondReplaceLen=-1 )
   {
      d_replaceLen = replaceLength;
      d_replaceLenTwo = secondReplaceLen;
   }
   void setBPrange( const int bpr_start, const int bpr_end )
   {
      d_bpr_start = bpr_start;
      d_bpr_end = bpr_end ;
   }
   void setGenome(Genome& genome)
   {
      d_genome_ptr = &genome;
   };
   bool operator<(const SVData& otherSV ) const;
   bool operator==(const SVData& otherSV ) const;
   void fuse( SVData& otherSV );
   string getReference() const; // reference sequence, for indels including the base before it
   string getAlternative() const;
   string getOutputFormattedReference() const;
   string getOutputFormattedAlternative() const;
   void setNT(string nt)
   {
      d_nt = nt;
   };
   void setSecondNT(string secondNT)
   {
      d_nt2 = secondNT;
   };
   bool withinAllowedRepeatsPostIndel(const int maxRepeatLen, const int maxNoRepeats) const;
   bool withinAllowedRepeatsInternal(const int maxRepeatLen, const int maxNoRepeats) const;
   bool isEquilengthReplacement() const
   {
      if ((d_svtype=="RPL" && d_svlen == d_replaceLen ) || (d_svtype=="INV" && d_replaceLen==0 && d_replaceLenTwo==0 )) {
         return true;
      } else {
         return false;
      }
   }

private:

   bool altSameLengthAsRef() const
   {
      return ((d_svtype=="RPL" && d_svlen == d_replaceLen) ||
              (d_svtype=="INV" && d_replaceLen==0 && d_replaceLenTwo ==0 ));
   }

   int d_position; // 1-based

   int d_end;
   int d_homlen;
   int d_bpr_start, d_bpr_end;
   int d_svlen;
   int d_replaceLen;
   int d_replaceLenTwo; // for inversions
   double p_somatic;

   string d_nt;
   string d_nt2; // for inversions
   vector<Genotype> d_format;
   Genome* d_genome_ptr;
   string d_chromosome;
   string d_id; // default '.', as we don't mine variant databases yet
   string d_quality; // '.' by default, but can be floating-point number
   string d_filter;  // "PASS" by default
   string d_svtype;
   string d_homseq;

   string getSVSequence() const;
};

string SVData::getAlternative() const
{
   if (d_svtype == "INS" && d_svlen == 0)  { // long insertion
      return "<INS>";
   }
   string altVariant = "";
   const string* reference = d_genome_ptr->getChromosome( d_chromosome );
   if ( d_svtype == "INS" || d_svtype == "DEL" || d_svtype == "RPL" ) {
      if (!(g_par.gatkCompatible && altSameLengthAsRef())) {
         altVariant += (*reference)[ d_position ];
      }
      altVariant += d_nt;
   } else if ( d_svtype == "DUP:TANDEM" ) {
      string refVariant = getReference();
      altVariant += refVariant;
      altVariant += d_nt;
      altVariant += refVariant.substr(1); // skip the position before...
   } else if ( d_svtype == "INV" ) {
      string refVariant = getReference();
      if (g_par.gatkCompatible && altSameLengthAsRef()) {
         string complement = "";
         createComplement( refVariant, complement );
         altVariant = complement;
      } else {
         altVariant += (*reference)[ d_position ];
         altVariant += d_nt;
         string referenceInv = refVariant.substr(1);
         string complement = "";
         createComplement( referenceInv, complement );
         altVariant += complement;
         altVariant += d_nt2;
      }
   }
   return altVariant;
}


// getReference returns the reference sequence, for indels including the base before it
string SVData::getReference() const
{
   const string* reference = d_genome_ptr->getChromosome( d_chromosome );
   if (d_svtype == "INS" && d_svlen == 0)  { // long insertion
      string refVariant = "";
      refVariant+= (*reference)[ d_position ];
      return refVariant;
   } else { // normal insertion/deletion or whatever
      string refVariant="";
      int startPosition = d_position;
      if (g_par.gatkCompatible && altSameLengthAsRef() ) {
         startPosition = d_position + 1; // workaround GATK
      }
      for (int position=startPosition; position<d_end; position++ ) {
         refVariant += (*reference)[ position ];
      }
      return refVariant;
   }
}

void reportProblematicSVSize( int size, const string& chromosome, int position )
{
   if ( size >= g_sizeToWarnFor ) {
      cout << "Warning! The SV at chromosome " << chromosome << ", position " << position << " is of size " << size << ".";
      cout << " It won't be put in the VCF in base sequence detail, but in the format";
      cout << " 'chrom pos firstrefbase <SVType>'.";
      cout << " This behaviour can be overridden by setting the -co parameter to -1 or to a larger" ;
      cout << " value than the current SV size.\n";
   }
}

string SVData::getOutputFormattedReference() const
{
   string defaultRef = getReference();
   string defaultAlt = getAlternative();
   reportProblematicSVSize( defaultRef.size(), getChromosome(), getPosition() );
   if (defaultAlt == "<INS>" ) {
      return defaultRef;
   } else {
      if (g_par.compactOutput>1) {
         if ( defaultRef.size() > g_par.compactOutput || defaultAlt.size() > g_par.compactOutput ) {
            defaultRef.erase(1);
         }
      }
   }
   return defaultRef;
}

string SVData::getOutputFormattedAlternative() const
{
   string defaultRef = getReference();
   string defaultAlt = getAlternative();

   reportProblematicSVSize( defaultAlt.size(), getChromosome(), getPosition() );
   if (defaultAlt == "<INS>" ) {
      return defaultAlt;
   } else {
      if ( g_par.compactOutput > 1 ) {
         if ( defaultRef.size() > g_par.compactOutput || defaultAlt.size() > g_par.compactOutput ) {
            defaultAlt = "<" + d_svtype + ">";
         }
      }
   }
   return defaultAlt;
}

SVData::SVData(const int genotypeTotal) // default settings
{
   d_id=".";
   d_quality=".";
   d_filter=".";
   d_replaceLen=0;
   d_homlen=0;
   d_end = 0;
   d_homseq="";
   d_nt="";
   d_nt2="";
   int numberOfSamples = genotypeTotal;
   if (numberOfSamples<=0) {
      numberOfSamples=1;
   }
   d_format.resize( numberOfSamples, Genotype(0,0,0) );
};


/* 'bothStrands' Is a SV supported by reads on both strands? */
bool SVData::bothStrands() const
{
   bool strandPlus = false;
   bool strandMinus = false;

   for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
      if ( d_format[ sampleIndex ].getReadDepthPlus() > 0 ) {
         strandPlus = true;
      }
      if ( d_format[ sampleIndex ].getReadDepthMinus() > 0 ) {
         strandMinus = true;
      }
   }
   return ( strandPlus && strandMinus );
}


/* 'getNumSupportSamples': how many samples support this SV? */
int SVData::getNumSupportSamples(const bool onlyBalancedSamples, const int minimumStrandSupport) const
{
   int numSupportingSamples = 0;
   for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
      int plusSupport = d_format[ sampleIndex ].getReadDepthPlus();
      int minSupport = d_format[ sampleIndex ].getReadDepthMinus();
      if (onlyBalancedSamples) {
         if (plusSupport>=minimumStrandSupport && minSupport>=minimumStrandSupport ) {
            numSupportingSamples++;
         }
      } else { // Sample does not need to be balanced
         if (plusSupport>=minimumStrandSupport || minSupport>=minimumStrandSupport ) {
            numSupportingSamples++;
         }
      }
   }
   return numSupportingSamples;
}


/* 'getNumSupportReads': how many reads support this SV? */
int SVData::getNumSupportReads() const
{
   int numReads = 0;
   for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
      numReads += d_format[ sampleIndex ].getTotalReads();
   }
   return numReads;
}


/* comparison operator of SVData objects; helps sorting them for output. */
bool SVData::operator<(const SVData& otherSV ) const
{
   if ( d_chromosome.compare( otherSV.d_chromosome ) != 0 ) {
      return ( d_chromosome.compare( otherSV.d_chromosome ) < 0 );
   }
   if ( getPosition() != otherSV.getPosition()) {
      return ( getPosition() < otherSV.getPosition() );
   } else {
      return ( d_svlen < otherSV.d_svlen );   // position equal: then do smallest SV first
   }
}


/* 'testHypothesis' tests whether the DNA string consists of a number of repeating "hypothesis" strings. If so, it
 returns the number of repeats; otherwise, it returns zero. */
int testHypothesis(const string& hypothesis, const string& sequence )
{
   int hypLen = hypothesis.size();
   int sequenceLen = sequence.size();

   for (int testedBaseIndex=0; testedBaseIndex<sequence.size(); testedBaseIndex++ ) {
      char currentHypothesisBase = hypothesis[ testedBaseIndex % hypLen ];
      if ( currentHypothesisBase != sequence[ testedBaseIndex ] ) {
         return 0;
      }
   }
   return sequenceLen/hypLen;
}


/* 'countRepeats' counts the repeats in a sequence. Returns the shortest result which can explain (allowing for rotation, CACAC => CA) the input sequence */
int countRepeats( const string& sequence, const int maxRepeatLength, int &bestSize )
{
   int maximumLen = min( maxRepeatLength, (int)(sequence.size()/2) );
   if (maxRepeatLength<0 ) {
      // no maximum repeat length set; assume it is infinite
      maximumLen = (int)(sequence.size()/2);
   }
   string hypothesis = "";
   int bestRepeatLength = 0;
   int bestRepeatNumber = 0;
   for (int repeatLen=1; repeatLen<=maximumLen; repeatLen++ ) {
      hypothesis += sequence[ repeatLen - 1 ];
      int repeats = testHypothesis( hypothesis, sequence );
      if (repeats>0 && repeats*hypothesis.size() > bestRepeatLength * bestRepeatNumber) {
         bestRepeatLength = hypothesis.size();
         bestRepeatNumber = repeats;
      }
   }
   bestSize = bestRepeatLength;
   return bestRepeatNumber;
}


/* 'getSVSequence' returns the inserted or deleted sequence; if the mutation is complex or weird, it returns the new ('alt') sequence. */
string SVData::getSVSequence() const
{
   string ref = getReference();
   string alt = getAlternative();
   string modifiedSequence = "";
   int pos=0;
   int maxPos=min( ref.size(), alt.size() );
   while (pos < maxPos && ref[pos]==alt[pos]) {
      pos++;
   };
   if (pos == maxPos ) { // simple insertion or deletion...
      modifiedSequence = ( maxPos==ref.size() ? alt.substr(pos) : ref.substr(pos));
   } else { // replacement
      modifiedSequence = alt.substr( pos );
   }
   return modifiedSequence;
}


/* 'withinAllowedRepeatsPostIndel' Is the number of repeats after the detected SV within the maximum allowed number of repeats (maxNoRepeats) of the
 fundamental repetitive unit of the SV? */
bool SVData::withinAllowedRepeatsPostIndel(const int maxRepeatLen, const int maxNoRepeats) const
{
   // 1. get the minimum repeating unit from the inserted sequence
   string sequenceToAnalyze = getSVSequence();
   int actualRepeatLength = 0;
   int repeatCount = countRepeats( sequenceToAnalyze, maxRepeatLen, actualRepeatLength );
   if (actualRepeatLength>0) {
      string hypothesis = sequenceToAnalyze.substr( 0 , actualRepeatLength );
      int extendedRepeatCount = testHypothesis( hypothesis, sequenceToAnalyze + d_homseq );
      return (extendedRepeatCount - repeatCount <= maxNoRepeats );
   } else {
      int bestSize = 0;
      int extendedRepeatCount = countRepeats( sequenceToAnalyze + d_homseq, maxRepeatLen, bestSize );
      int repetitiveLength = bestSize * extendedRepeatCount;
      int repetitivePartPostIndel = repetitiveLength - sequenceToAnalyze.size();
      double noRepPostIndel = (double)repetitivePartPostIndel / bestSize;
      return (int)noRepPostIndel <= maxNoRepeats;
   }
}


/* 'withinAllowedRepeatsInternal' Is the number of repeats of the smallest repetitive unit in the detected SV within the maximum allowed number of repeats (maxNoRepeats)? */
bool SVData::withinAllowedRepeatsInternal(const int maxRepeatLen, const int maxNoRepeats) const
{
   string sequenceToAnalyze = getSVSequence();
   int actualRepeatLength = 0;
   int repeatCount = countRepeats( sequenceToAnalyze, maxRepeatLen, actualRepeatLength );
   if ( repeatCount > maxNoRepeats ) {
      return false;
   } else {
      return true;
   }
}



/* first version of operator==. Are two events the same? We could make this more complicated, but first see if it works. */
bool SVData::operator==(const SVData& otherSV ) const
{
   // for non-NT
   if ( ( d_svtype == "DEL" )
         && ( otherSV.d_svtype == "DEL" )
         && ( d_bpr_start == otherSV.d_bpr_start )
         && ( d_bpr_end == otherSV.d_bpr_end )
         && ( d_svlen == otherSV.d_svlen )
         && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
         //( d_position == otherSV.d_position )
         //&& ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
         //&& ( d_svtype.compare( otherSV.d_svtype ) == 0 )
         //&& ( d_end == otherSV.d_end )
         //&& ( d_replaceLen == otherSV.d_replaceLen )
      ) {
      return true;
   }
   // for NT
   if ( ( d_svtype == "RPL" )
         && ( otherSV.d_svtype == "RPL" )
         && ( ( d_svlen - d_replaceLen) == ( otherSV.d_svlen - otherSV.d_replaceLen ) )
         && ( d_bpr_start == otherSV.d_bpr_start )
         && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
      ) {
      return true;
   }
   if ( ( d_svtype == "INS" )
         && ( otherSV.d_svtype == "INS" )
         //&& ( ( d_svlen - d_replaceLen) == ( otherSV.d_svlen - otherSV.d_replaceLen ) )
         && ( d_bpr_start == otherSV.d_bpr_start )
         && ( d_bpr_end == otherSV.d_bpr_end )
         && ( d_svlen == otherSV.d_svlen )
         && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
      ) {
      return true;
   }
   return false;
}

/* 'fuse' takes over all parameters from the other (earlier-occurring and hence better) SV, but adds the total support. */
void SVData::fuse( SVData& otherSV )
{
   for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
      otherSV.d_format[ sampleIndex ].fuse( d_format[ sampleIndex ] );
   }
   *this = otherSV;
}

int FACT(int n)
{

   if (n == 0 || n == 1) {
      return 1;
   }
   int fact=1;
   int i;

   for (i=1; i<=n; i++) {
      fact*=i;
   }
   std::cout << n << " " << fact << std::endl;
   return fact;
}

double fisher_test(int a, int c, int b, int d)
{
   double p = 0.0;
   int n = a + b + c + d;
   p = (FACT(a+b)*FACT(c+d)*FACT(a+c)*FACT(d+b)) / (double)(FACT(a)*FACT(b)*FACT(c)*FACT(d)*FACT(n));
   std::cout << a << " " << c << " " << b << " " << d << " " << p << std::endl;
   return p;
}


ostream& operator<<(ostream& os, const SVData& svd)
{
   double somatic_p_value = 0.0;



   os << svd.d_chromosome << "\t";
   os << svd.getPosition() << "\t";
   os << svd.d_id << "\t";
   os << svd.getOutputFormattedReference() << "\t";
   os << svd.getOutputFormattedAlternative() << "\t";
   os << svd.d_quality << "\t";
   if (svd.d_format.size() == 2 && g_par.somatic) {

      somatic_p_value = fisher_test(svd.d_format[0].getTotalReads(), svd.d_format[0].getTotalRefSupport(), svd.d_format[1].getTotalReads(), svd.d_format[1].getTotalRefSupport());
      //if (somatic_p_value < 0.05) svd.d_filter = "PASS";
   }
   if (somatic_p_value < 0.05) {
      os << "PASS\t";
   } else {
      os << svd.d_filter << "\t";
   }


   os << "END=" << svd.getVCFPrintEnd() << ";";
   os << "HOMLEN=" << svd.d_homlen << ";";
   if ( svd.d_homlen != 0 ) {
      os << "HOMSEQ=" << svd.d_homseq << ";";
   }

   os	<< "SVLEN=";
   if (svd.d_svtype.compare("RPL")==0 || svd.d_svtype.compare("DEL")==0 ) {
      if (svd.d_svlen > 0 ) {
         os << "-";  // deletions need negative SVLEN
      }
   }
   os << svd.d_svlen << ";";
   os << "SVTYPE=" << svd.d_svtype;
   if ( svd.d_svtype.compare("RPL")==0 || svd.d_svtype.compare("DUP:TANDEM")==0 || svd.d_svtype.compare("INV")==0 ) {
      os << ";NTLEN=" << svd.d_replaceLen;
   }
   if ( svd.d_svtype.compare("INV")==0 ) {
      os << "," << svd.d_replaceLenTwo;
   }

   if (svd.d_format.size() == 2 && g_par.somatic) {
      os << ";" << somatic_p_value;
   }


   if (pindel024uOrLater && svd.getAlternative()!="<INS>") {
      os << "\tGT:AD";
   } else {
      os << "\tGT:AD";
   }

   for (int counter=0; counter<svd.d_format.size(); counter++ ) {
      os << "\t";
      if (pindel024uOrLater && svd.getAlternative()!="<INS>") {
         os << svd.d_format[ counter ].getGTRDAD();
      } else {
         os << svd.d_format[ counter ].getGTAD();
      }
   }

   os << endl;

   return os;
}

/* 'fetchElement' returns the "index"th element of "infile"; so the first element if index is 1, the second if index is 2, etc. */
string fetchElement( istream& instream, const int index )
{
   string element = "";
   for ( int currentCounter=0; currentCounter<index; currentCounter++ ) {
      instream >> element;
   }
   return element;
}

/* 'countElements' counts the number of elements in stream "instream" */
int countElements( istream& instream )
{
   string element = "";
   int counter = 0;
   while (! instream.fail() ) {
      instream >> element;
      counter++;
   }
   return counter;
}

void showSet( set<string> aSet )
{

   set<string>::iterator index;
   int counter=1;
   for (index=aSet.begin(); index!=aSet.end(); index++ ) {
      cout << counter++ << ". " << *index << endl;
   }
}


/** 'isPindelSVIdentifier' returns whether the string passed to it represents a valid SV type (insertion, deletion, etc.)
	 NOTE: this is the SV identifier as used in pindel output files, pindel2vcf uses other SV types internally,
	 like 'INS' or 'RPL' */
bool isPindelSVIdentifier( string identifier )
{
   if (identifier == "D" || identifier == "I" || identifier == "LI" || identifier == "TD" || identifier == "INV" ) {
      return true;
   } else {
      return false;
   }
}


/** 'isSVSummarizingLine' checks whether 'line' is a SV summarizing line and therefore one that pindel2vcf needs to process
	 (other lines can be skipped) */
bool isSVSummarizingLine( string line )
{
   // first check if the line is a true summarizing name by counting if it has sufficient elements
   // (this is also useful to prevent reading of the SV type at the second position of the line to
   // go wrong if there is no second element)
   stringstream streamForCounting;
   streamForCounting << line;
   int elementsInLine = countElements( streamForCounting );
   if (elementsInLine < 2 ) {
      return false;
   }

   // there are at least two elements. Check if the second one indicates a SV-type
   stringstream lineStream;
   lineStream << line;
   string item = fetchElement( lineStream, 2 );
   return isPindelSVIdentifier( item );
}

void getSampleNamesAndChromosomeNames(InputReader& pindelInput, set<string>& sampleNames, set<string>&chromosomeNames)
{
   //cout << "DEBUG:start GetSampleNamesAndChromosomeNames\n";
   string line;
   unsigned int counter=0;

   // loop over the input file, registering each 'interesting' pindel output line (pindel output has the structure (ABCDn)m
   // A is a ##### separator line, B the interesting line (with sample names and coverage), C the line with the reference DNA
   // sequence (starting with ACTG) and the D-lines start with spaces and give the individual reads that support an event
   // HOWEVER, it is also possible that an user greps only the 'useful' output lines.
   while (!pindelInput.eof()) {

      // find the next 'B-line' (the line containing the sample names, coverage etc.)
      do {
         line = pindelInput.getLine();
      } while (!pindelInput.eof() && !isSVSummarizingLine( line ));

      if (pindelInput.eof()) {
         //cout << "DEBUG:end GetSampleNamesAndChromosomeNames\n";
         return;
      }

      // 'line' should contain a B-line now. Analyze it.
      counter++;
      stringstream lineStream;
      lineStream << line;
      stringstream streamForCounting;
      streamForCounting << line;
      int elementsInLine = countElements( streamForCounting );
      string svType = fetchElement( lineStream, 2 );


      // 'LI' types don't have the same format as normal Pindel lines.
      if ( svType.compare("LI")==0 ) {
         string chromosomeName = fetchElement( lineStream, 2 );
         chromosomeNames.insert( chromosomeName );
         string firstSampleName = fetchElement( lineStream, 7);
         sampleNames.insert( firstSampleName );
         string newSampleName = fetchElement( lineStream, 5 );
         while (!lineStream.fail()) {
            sampleNames.insert( newSampleName );
            newSampleName = fetchElement( lineStream, 5 );
         }
         continue;
      }
      string chromosomeName = fetchElement( lineStream, 6 );
      chromosomeNames.insert( chromosomeName );
      //		cout << "Studying chromosome " << chromosome << endl;

      // 8 = 2+6, so corrects for previous reads
      int numberOfSamples = atoi( fetchElement( lineStream, FIRST_SAMPLE_INDEX - 12 ).c_str() );
      string firstSampleName = fetchElement( lineStream, 4 );
      if ( firstSampleName != "" ) {
         sampleNames.insert( firstSampleName );
      }

      //cout << "ElInLine: " << elementsInLine << ", FSINDEX: " << FIRST_SAMPLE_INDEX << ", NoS=" << numberOfSamples << endl;
      if (elementsInLine> FIRST_SAMPLE_INDEX + 5* numberOfSamples ) {
         pindel024uOrLater = true;
      }
      /*      else { //*** This code seems to give trouble with some pindel output; pindel output format not consistent?
       pindel024uOrLater = false;
       }*/
      int numberOfElementsPerSample = ( pindel024uOrLater ? 7 : 5 );
      string newSampleName = fetchElement( lineStream, numberOfElementsPerSample );
      while (!lineStream.fail()) {
         if ( newSampleName != "" ) {
            sampleNames.insert( newSampleName );
         }
         newSampleName = fetchElement( lineStream, numberOfElementsPerSample );
      }
   }

   //cout << "DEBUG:end GetSampleNamesAndChromosomeNames\n";
}

template<class T> void showVector( vector<T> vect )
{
   cout << "Printing vector: " << endl;
   for (int index=0; index<vect.size(); index++ ) {
      cout << index << " " << vect[ index ] << endl;
   }
   cout << "End of printing.\n";
}




/* 'convertIndelToSVdata' converts insertions and deletions to nicely formatted SV-data. */
void convertIndelToSVdata( InputReader& pindelInput, map< string, int>& sampleMap, Genome& genome, SVData& svd, const string& targetChromosomeID)
{
   string line;
   svd.setGenome( genome );
   do {
      line = pindelInput.getLine();
   } while (!pindelInput.eof() && !isSVSummarizingLine( line ));

   if (pindelInput.eof()) {
      return;
   }

   stringstream lineStream;
   lineStream << line;
   string svType = fetchElement( lineStream, 2 ); // to 2
   if ( svType.compare("LI") == 0 ) {
      svd.setSVtype("INS");
      svd.setSVlen( 0 );
      string chromosomeID = fetchElement( lineStream, 2);
      const string* reference = genome.getChromosome( chromosomeID );
      if ( reference== NULL ) {
         cout << "Error! Reference chromosome \"" << chromosomeID << "\" not found!" << endl;
         exit(EXIT_FAILURE);
      }
      svd.setChromosome( chromosomeID );
      if (chromosomeID!=targetChromosomeID) {
         return;
      }
      int beforeStartPos = atoi( fetchElement( lineStream, 1 ).c_str() );
      svd.setPosition( beforeStartPos );
      int totalPlusSupport = atoi( fetchElement( lineStream, 2 ).c_str());
      int rightmostEndPos = atoi (fetchElement( lineStream, 1 ).c_str()); // now at position 14
      //cout << "plusSupport, righmostEndPos is " << plusSupport << ", " << rightmostEndPos << endl;
      svd.setEnd( rightmostEndPos );
      svd.setBPrange( beforeStartPos, rightmostEndPos );
      int totalMinSupport = atoi( fetchElement( lineStream, 2 ).c_str());
      // if the file has been created by a recent version of pindel, read in the extra elements
      string sampleName = fetchElement( lineStream, 1);
      int refSupportAtStartOfEvent = 0;
      int refSupportAtEndOfEvent = 0;

      /*if ( pindel024uOrLater ) {
         refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str());
         refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str());
         }*/
      int totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
      int numberItemsUntilNextSupport = ( pindel024uOrLater ? 2 : 2 );
      int samplePlusSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str());
      int sampleMinSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str()); // now at position 35, total +supports sample 1
      //int count=0;
      while (!lineStream.fail()) {
         if (sampleMap.find( sampleName )==sampleMap.end() ) {
            cout << "Error: could not find sample " << sampleName << endl;
         } else {
            int sampleID = sampleMap[ sampleName ];
            svd.addGenotype( sampleID, samplePlusSupport , sampleMinSupport, totalRefSupport );
         }
         /*if ( pindel024uOrLater ) {
          refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
          refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str());
          }*/
         int totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
         sampleName = fetchElement( lineStream, 1); // for unique support, 2->1,
         samplePlusSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str()); // for unique support, 2->1
         sampleMinSupport = atoi( fetchElement( lineStream, numberItemsUntilNextSupport ).c_str()); // now at position 33, total +supports sample 1
      }
      return;
   }
   svd.setSVlen( fetchElement( lineStream, 1 ) ); // to 3

   // get number(s) of NT bases added (two numbers for inversions!)
   string numNTaddedStr = fetchElement( lineStream, 2 ); // to 5

   int numNTadded = atoi( numNTaddedStr.c_str() ); // should get first number
   bool simpleInversion = false;
   int numNTinvAdded=-1;
   //cout << "Printing " << numNTaddedStr << endl;
   if ( svType.compare("INV") == 0 ) { // two numbers separated by : instead of one
      //cout << "Found ':' at position " << numNTaddedStr.find(":") << endl;
      if (numNTaddedStr.find(":")==string::npos) {
         //cout << "Found simple inversion!\n";
         simpleInversion = true;
      } else {
         int separatorPos = numNTaddedStr.find(":");
         string secondNumber = numNTaddedStr.substr(separatorPos+1);
         numNTinvAdded = atoi( secondNumber.c_str() );
      }
   }

   string ntAdded = fetchElement( lineStream, 1 ); // to 6
   string ntInvAdded = "";

   // basically, there are two type of inversions:
   //	a) The 'alternatively called small deletions' INV 2 NT 2 "TG"
   // b) the regular inversions INV98 NT 0:60 "":"GCT"
   if ( svType.compare("INV") == 0 ) {
      if (ntAdded.find(":")==string::npos ) {
         simpleInversion = true;
      } else {
         int separatorPos = ntAdded.find(":");
         ntInvAdded = ntAdded.substr( separatorPos+2, numNTinvAdded ); // erases ""
         svd.setSecondNT( ntInvAdded );
         ntAdded = ntAdded.substr(0,separatorPos);
      }
   }
   ntAdded.erase(0,1); // erases opening "
   ntAdded.erase(numNTadded); // erases closing "
   if (!simpleInversion) {
      svd.setNT( ntAdded );
   }

   string chromosomeID = fetchElement( lineStream, 2); // now at position 8
   if (chromosomeID!=targetChromosomeID) {
      return;
   }
   const string* reference = genome.getChromosome( chromosomeID );
   //cout << "reference is " << *reference << endl;
   if ( reference== NULL ) {
      cout << "Error! Reference chromosome \"" << chromosomeID << "\" not found!" << endl;
      exit(EXIT_FAILURE);
   }
   svd.setChromosome( chromosomeID );
   int beforeStartPos = atoi( fetchElement( lineStream, 2 ).c_str() ); // pos 10
   svd.setPosition( beforeStartPos ); // now at position 10
   int leftmostEndPos = atoi( fetchElement( lineStream, 1 ).c_str()); // now at position 11
   int leftmostStartPos = atoi (fetchElement( lineStream, 2 ).c_str());  // at position 13
   int rightmostEndPos = atoi (fetchElement( lineStream, 1 ).c_str()); // now at position 14
   svd.setBPrange( leftmostStartPos, rightmostEndPos );
   svd.setEnd( leftmostEndPos );
   svd.setHomlen( rightmostEndPos - leftmostEndPos );
   string homSeq="";
   for (int position=leftmostEndPos; position<rightmostEndPos; position++ ) {
      homSeq += (*reference)[ position ];
   }
   svd.setHomseq( homSeq );
   if ( svType.compare("D")==0 ) {
      if (numNTadded==0 ) {
         svd.setSVtype( "DEL" );
         svd.setReplace( 0 );
      } else { // some NT-bases added
         svd.setSVtype( "RPL" );
         svd.setReplace( numNTadded );
      }
   } else if ( svType.compare("I")==0 ) {
      svd.setSVtype("INS");
      svd.setReplace( 0 );
   } else if ( svType.compare("TD")==0 ) {
      svd.setSVtype("DUP:TANDEM");
      svd.setReplace( numNTadded );
   } else if ( svType.compare("INV") == 0 ) {
      svd.setSVtype("INV");
      if (simpleInversion) {
         svd.setReplace( 0, 0 );
      } else {
         svd.setReplace( numNTadded, numNTinvAdded );
      }
   }
   string sampleName = fetchElement( lineStream, 18);
   int refSupportAtStartOfEvent = 0;
   int refSupportAtEndOfEvent = 0;

   if ( pindel024uOrLater ) {
      refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
      refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
      //std::cout << "outside refSupportAtStartOfEvent: " <<  sampleName << "\t" << refSupportAtStartOfEvent << "\t" << "refSupportAtEndOfEvent: " << refSupportAtEndOfEvent << std::endl;
   }
   int totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
   //std::cout << "outside totalRefSupport: " << totalRefSupport << std::endl;
   int numberOfItemsUntilNextSupport = ( pindel024uOrLater ? 2 : 2 );
   int plusSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport - 1 ).c_str()); // now at position 33, total +supports sample 1; for unique support 1->2
   int minSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport ).c_str()); // now at position 35, total +supports sample 1
   int count=0;
   while (!lineStream.fail()) {
      if (sampleMap.find( sampleName )==sampleMap.end() ) {
         cout << "Error: could not find sample " << sampleName << endl;
      } else {
         int sampleID = sampleMap[ sampleName ];
         //std::cout << "Adding " << sampleID << "\t" << plusSupport << "\t" << minSupport << "\t" << totalRefSupport << std::endl;
         svd.addGenotype( sampleID, plusSupport , minSupport, totalRefSupport );
      }
      sampleName = fetchElement( lineStream, 2); // for unique support, 2->1,
      if ( pindel024uOrLater ) {
         refSupportAtStartOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
         refSupportAtEndOfEvent = atoi( fetchElement( lineStream, 1 ).c_str() );
         //std::cout << "inside refSupportAtStartOfEvent: " << sampleName << "\t" << refSupportAtStartOfEvent << "\t" << "refSupportAtEndOfEvent: " << refSupportAtEndOfEvent << std::endl;
      }
      totalRefSupport = max(refSupportAtStartOfEvent, refSupportAtEndOfEvent);
      //std::cout << "insert totalRefSupport: " << totalRefSupport << std::endl;
      plusSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport - 1 ).c_str()); // for unique support, 2->1
      minSupport = atoi( fetchElement( lineStream, numberOfItemsUntilNextSupport ).c_str()); // now at position 33, total +supports sample 1
   }
}

/* 'readReference' reads in the reference. */
void readReference( const string& referenceName, Genome& genome )
{
   ifstream referenceFile( referenceName.c_str() );
   if (referenceFile.fail()) {
      cout << "Cannot open reference file. Exiting.\n";
      exit(EXIT_FAILURE);
   }
   //reference="N"; // trick to make the reference automatically 1-positioned
   string refLine, refName, currentLine;

   std::getline(referenceFile,refLine); // FASTQ format always has a first line with the name of the reference in it
   // loop over each chromosome
   do {

      int counter=1;
      refName = "";
      do {
         refName += refLine[ counter++ ];
      } while ( counter<refLine.size() && (refLine[ counter ] != ' ') && (refLine[ counter ] != '\t') && (refLine[ counter ] != '\n') && (refLine[ counter ] != '\r'));
      cout << "Scanning chromosome: " << refName << endl;
      Chromosome newChrom( refName, referenceName );
      std::getline(referenceFile,currentLine);
      while (!referenceFile.eof() && currentLine[0]!='>') {
         std::getline(referenceFile,currentLine);
      }
      genome.addChromosome( newChrom );
      refLine = currentLine;
   } while (!referenceFile.eof());
   //cout << "DEBUG:Exiting reference scanning.\n";
}


/* 'createParameters' creates the default parameters that the VCF converter uses. */
void createParameters()
{
   parameters.push_back(
      new StringParameter( &g_par.reference, "-r", "--reference", "The name of the file containing the reference genome", true, "" ) );
   parameters.push_back(
      new StringParameter( &g_par.referenceName, "-R", "--reference_name", "The name and version of the reference genome", true, "" ) );
   parameters.push_back(
      new StringParameter( &g_par.referenceDate, "-d", "--reference_date", "The date of the version of the reference genome used", true, "" ) );
   parameters.push_back(
      new StringParameter( &g_par.pindelfile, "-p", "--pindel_output", "The name of the pindel output file containing the SVs", false, "" ) );
   parameters.push_back(
      new StringParameter( &g_par.pindelroot, "-P", "--pindel_output_root", "The root-name of the pindel output file; this will result in\n"
                           "one big output file containing deletions, short and long insertions, tandem duplications and inversions.\n"
                           "For example, if the pindel output files are called sample1_D, sample1_SI, sample1_TD etc. then -P sample1 would combine the\n"
                           "information in all those sample files into one big vcf file.", false, "" ) );
   parameters.push_back(
      new StringParameter( &g_par.vcffile, "-v", "--vcf", "The name of the output vcf-file (default: name of pindel output file +\".vcf\"", false, "" ) );
   parameters.push_back(
      new StringParameter( &g_par.chromosome, "-c", "--chromosome", "The name of the chromosome (default: SVs on all chromosomes are processed)", false, "" ) );
   parameters.push_back(
      new IntParameter( &g_par.windowSize, "-w", "--window_size", "Memory saving option: the size of the genomic region in a chromosome of which structural variants are calculated separately, in millions of bases (default 300, for memory saving 100 or 50 recommended)", false, 300 ) );
   parameters.push_back(
      new IntParameter( &g_par.MinCoverage, "-mc", "--min_coverage", "The minimum number of reads to provide a genotype (default 10)", false, 10 ) );
   parameters.push_back(
      new FloatParameter( &g_par.HetCutoff, "-he", "--het_cutoff", "The propertion of reads to call het (default 0.2)", false, 0.2 ) );
   parameters.push_back(
      new FloatParameter( &g_par.HomCutoff, "-ho", "--hom_cutoff", "The propertion of reads to call het (default 0.8)", false, 0.8 ) );

   parameters.push_back(
      new IntParameter( &g_par.minsize, "-is", "--min_size", "The minimum size of events to be reported (default 1)", false, 1 ) );
   parameters.push_back(
      new IntParameter( &g_par.maxsize, "-as", "--max_size", "The maximum size of events to be reported (default infinite)", false, -1 ) );
   parameters.push_back(
      new BoolParameter( &g_par.bothstrands, "-b", "--both_strands_supported", "Only report events that are detected on both strands (default false)", false, false ) );
   parameters.push_back(
      new IntParameter( &g_par.minsuppSamples, "-m", "--min_supporting_samples", "The minimum number of samples an event needs to occur in with sufficient support to be reported (default 0)", false, 1 ) );
   parameters.push_back(
      new IntParameter( &g_par.minsuppReads, "-e", "--min_supporting_reads", "The minimum number of supporting reads required for an event to be reported (default 1)", false, 1 ) );
   parameters.push_back(
      new IntParameter( &g_par.maxSuppReads, "-f", "--max_supporting_reads", "The maximum number of supporting reads allowed for an event to be reported, allows protection against miscalls in due to segmental duplications or poorly mapped regions (default infinite)", false, -1 ) );
   parameters.push_back(
      new IntParameter( &g_par.regionStart, "-sr", "--region_start", "The start of the region of which events are to be reported (default 0)", false, 0 ) );
   parameters.push_back(
      new IntParameter( &g_par.regionEnd, "-er", "--region_end", "The end of the region of which events are to be reported (default infinite)", false, -1 ) );
   parameters.push_back(
      new IntParameter( &g_par.maxInterRepeatNo, "-ir", "--max_internal_repeats", "Filters out all indels where the inserted/deleted sequence is a homopolymer/microsatellite of more than X repetitions (default infinite). For example: T->TCACACA has CACACA as insertion, which is a microsattelite of 3 repeats; this would be filtered out by setting -ir to 2", false, -1 ) );
   parameters.push_back(
      new IntParameter( &g_par.compactOutput, "-co", "--compact_output_limit", "Puts all structural variations of which either the ref allele or the alt allele exceeds the specified size (say 10 in '-co 10') in the format 'chrom pos first_base <SVType>'", false, g_sizeToWarnFor ) );
   parameters.push_back(
      new IntParameter( &g_par.maxInterRepeatLength, "-il", "--max_internal_repeatlength", "Filters out all indels where the inserted/deleted sequence is a homopolymers/microsatellite with an unit size of more than Y, combine with the option -ir. Default value of -il is infinite. For example: T->TCAGCAG has CAGCAG as insertion, which has the fundamental repetitive unit CAG of length 3. This would be filtered out if -il has been set to 3 or above, but would be deemed 'sufficiently unrepetitive' if -il is 2", false, -1 ) );
   parameters.push_back(
      new IntParameter( &g_par.maxPostRepeatNo, "-pr", "--max_postindel_repeats", "Filters out all indels where the inserted/deleted sequence is followed by a repetition (of over X times) of the fundamental repeat unit of the inserted/deleted sequence. For example, T->TCACA would usually be a normal insertion, which is not filtered out, but if the real sequence change is TCACACA->TCACACACACA, it will be filtered out by -pr of 1 or above, as the fundamental repeat unit of the inserted sequence (CA) is repeated more than one time in the postindel sequence [indel sequence CACA, postindel sequence CACACA]. Note: when CAC is inserted next to ACACAC, the repeat sequence is recognized as CA, even though the 'postrepeat' sequence is ACACAC", false, -1 ) );
   parameters.push_back(
      new IntParameter( &g_par.maxPostRepeatLength, "-pl", "--max_postindel_repeatlength", "Filters out all indels where the inserted/deleted sequence is followed by a repetition of  the fundamental repeat unit of the inserted/deleted sequence; the maximum size of that 'fundamental unit' given by the value of -pl (default infinite) For example: TCAG->TCAGCAG has insertion CAG and post-insertion sequence CAG. This insertion would be filtered out if -pl has been set to 3 or above, but would be deemed 'sufficiently unrepetitive' if -pl is 2", false, -1 ) );
   parameters.push_back(
      new BoolParameter( &g_par.onlyBalancedSamples, "-sb", "--only_balanced_samples", "Only count a sample as supporting an event if it is supported by reads on both strands, minimum reads per strand given by the -ss parameter. (default false)", false, 0 ) );
//   parameters.push_back(
  //    new BoolParameter( &g_par.somatic, "-so", "--somatic_p", "compute somatic p value when two samples are present, assume the order is normal and tumor. (default false)", false, 0 ) );

   parameters.push_back(
      new IntParameter( &g_par.minimumStrandSupport, "-ss", "--minimum_strand_support", "Only count a sample as supporting an event if at least one of its strands is supported by X reads (default 1)", false, 1 ) );
   parameters.push_back(
      new BoolParameter( &g_par.gatkCompatible, "-G", "--gatk_compatible", "calls genotypes which could either be homozygous or heterozygous not as ./1 but as 0/1, to ensure compatibility with GATK", false, false ) );
   parameters.push_back(
      new BoolParameter( &g_par.showHelp, "-h", "--help", "Print the help of this converter", false, false ) );
}

/* 'findParameter' returns the index of the parameter with name 'name'; -1 if not found.*/
int findParameter(string name)
{
   for (int parameterCounter=0; parameterCounter<parameters.size(); parameterCounter++ ) {
      if (parameters[ parameterCounter ]->hasName( name ) ) {
         return parameterCounter;
      }
   }
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
      if ( parameterIndex == -1 ) {
         cout << "unknown argument: " << currentArgument << endl;
         return;
      }

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
         if (argumentIndex >= argc ) {
            cout << "argument of " << currentArgument << " lacking.\n";
            return ;
         }
         if (argv[ argumentIndex ][0]=='-' ) {
            cout << "argument of " << currentArgument << " seems erroneous.\n";
            return ;
         }
         // but if everything is allright,
         //cout << "Giving " << currentArgument << " the value " << argv[ argumentIndex ] << endl;
         parameters[ parameterIndex ]->setValue( string(argv[ argumentIndex ]) );
      }
   }
}

/* 'printHelp' prints all parameters available. */
void printHelp()
{
   cout << "\nProgram:    " << g_programName <<" (conversion of Pindel output to VCF format)\n";
   cout << "Version:    " << g_versionString << endl;
   cout << "Contact:    Eric-Wubbo Lameijer <e.m.w.lameijer@gmail.com>\n";
   cout << "Usage:      " << g_programName << " -p <pindel_output_file> -r <reference_file>\n";
   cout << "              -R <name_and_version_of_reference_genome> -d <date_of_reference_genome_version>\n";
   cout << "              [-v <vcf_output_file>]\n\n";
   cout << "           the -v parameter is optional; when no output file name is given, output is written\n";
   cout << "           to a file with the name <pindel_output_file>.vcf.\n\n";
   cout << "Example:    " << g_programName << " -p sample3chr20_D -r human_g1k_v36.fasta -R 1000GenomesPilot-NCBI36\n";
   cout << "              -d 20101123 -v sample3chr20_D.vcf\n\n";
   cout << "or (with -P): " << g_programName << " -P sample3chr20 -r human_g1k_v36.fasta -R 1000GenomesPilot-NCBI36\n";
   cout << "              -d 20101123 -v sample3chr20_all.vcf\n\n";
   cout << "Note:      -is only guaranteed to work correctly on output files produced by pindel version 0.2.3 and above.\n";
   cout << "           -LI and BP files (long insertion and break point files) have a different type of header and\n";
   cout << "            are not supported yet.\n\n";

   for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
      parameters[ parameterIndex ]->describe();
   }
   exit( EXIT_SUCCESS );
}

bool isRegularPindelInput()
{
   return parameters[ findParameter("-p") ]->isSet();
}
bool isRootPindelInput()
{
   return parameters[ findParameter("-P") ]->isSet();
}

/* 'checkParameters' checks whether all required parameters have been set. */
bool checkParameters()
{
   if (parameters[ findParameter("-h") ]->getBValue() ) {
      printHelp();
   }

   bool canRun = true;
   for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
      if (parameters[ parameterIndex ]->isRequired() && !parameters[ parameterIndex ]->isSet()) {
         cout << "\nRequired parameter " << parameters[ parameterIndex ]->getShortName() << "/" << parameters[ parameterIndex ]->getLongName()
              << " " << parameters[ parameterIndex ]->getDescription() << " needs to be set.\n\n";
         canRun = false;
      }  //if
   }

   if ( isRegularPindelInput() && isRootPindelInput() ) {
      cout << "Sorry, you can't use -p and -P at the same time, please choose one option.\n\n";
      canRun = false;
   } else if ( !isRegularPindelInput() && !isRootPindelInput() ) {
      cout << "Pindel2vcf needs a pindel input file, either use the -p or the -P option, please.\n\n";
      canRun = false;
   }

   if (!canRun) {
      cout << "For further information, please run " << g_programName <<" without arguments or with option -h/--help.\n\n";
   }
   return canRun;
}

/* 'setParameters' sets the filters to be used in the rest of the program. */
void setParameters()
{
   g_par.vcffile = parameters[ findParameter( "-v" )]->getSValue();
   if (g_par.vcffile.compare("")==0) {
      if (isRegularPindelInput()) {
         g_par.vcffile = g_par.pindelfile + ".vcf";   // default
      } else if (isRootPindelInput()) {
         g_par.vcffile = g_par.pindelroot + ".vcf";
      } else {
         cout << "Error trying to construct output filename!\n";
         exit( EXIT_FAILURE );
      }
   }
}

/* 'throughFilter' checks whether the event is good enough to be written to the output file. */
bool throughFilter(SVData sv)
{
   if (( g_par.minsize > 1 ) && ( abs( sv.getSize()) < g_par.minsize ) ) {
      return false;
   }
   if (( g_par.maxsize > 0 ) && ( abs( sv.getSize()) > g_par.maxsize ) ) {
      return false;
   }
   if ( g_par.bothstrands && !sv.bothStrands() ) {
      return false;
   }
   if ( ( g_par.minsuppSamples >= 1 ) && ( sv.getNumSupportSamples(g_par.onlyBalancedSamples, g_par.minimumStrandSupport) < g_par.minsuppSamples ) ) {
      return false;
   }
   if ( ( g_par.minsuppReads >= 1 ) && ( sv.getNumSupportReads() < g_par.minsuppReads ) ) {
      return false;
   }
   if ( ( g_par.maxSuppReads >= 1 ) && ( sv.getNumSupportReads() > g_par.maxSuppReads ) ) {
      return false;
   }
   if ( ( g_par.regionStart > 0 ) && ( sv.getPosition() < g_par.regionStart ) ) {
      return false;
   }
   if ( ( g_par.regionEnd > 0 ) && ( sv.getPosition() > g_par.regionEnd ) ) {
      return false;
   }
   if ( g_par.maxInterRepeatNo >= 0 && !sv.withinAllowedRepeatsInternal(g_par.maxInterRepeatLength, g_par.maxInterRepeatNo )) {
      return false;
   }
   if ( g_par.maxPostRepeatNo >= 0 && !sv.withinAllowedRepeatsPostIndel(g_par.maxPostRepeatLength, g_par.maxPostRepeatNo )) {
      return false;
   }

   /*if (g_par.gatkCompatible && sv.isEquilengthReplacement() ) {
     return false;
     }*/



   // all filters passed
   return true;
}

/* 'makeSampleMap' converts a set of sample names to a map containing both sample names and the genotype of a sample. */
void makeSampleMap( const set<string>& sampleNames, map<string, int>& sampleMap )
{
   int count=0;
   for (set<string>::iterator setIt=sampleNames.begin(); setIt!=sampleNames.end(); setIt++ ) {
      sampleMap.insert( pair<string,int>( *setIt, count++) );
   }
}

void initBaseArray()
{
   for (int i=0; i<256; i++) {
      g_normalBaseArray[i] = false;
   }
   g_normalBaseArray['A'] = true;
   g_normalBaseArray['C'] = true;
   g_normalBaseArray['G'] = true;
   g_normalBaseArray['T'] = true;
   g_normalBaseArray['N'] = true;
}

void reportSVsInChromosome(
   const string& chromosomeID,
   const set<string>& chromosomeNames,
   const set<string>& sampleNames,
   InputReader& pindelInput,
   map< string, int >& sampleMap,
   Genome& genome,
   ofstream& vcfFile
)
{
   // if no reads have been found for this chromosome, skip it
   if (chromosomeNames.find(chromosomeID) == chromosomeNames.end() ) {
      cout << "No reads for chromosome " << chromosomeID << ", skipping it.\n";
      return;
   }
   cout << "Processing chromosome " << chromosomeID << endl;
   // rewind file to start
   int regionStart = 0;
   int regionEnd = 0;
   SVData backupSV(sampleNames.size() );
   bool backupAvailable = false;
   do {
      cout << "reportSVsInChromosome: start reading region.\n";
      regionEnd = regionStart + g_par.windowSize*1000000;
      cout << "Reading region " << regionStart << "-" << regionEnd << endl;
      pindelInput.rewind();
      int counter=0;
      vector<SVData> svs;
      if (backupAvailable) {
         svs.push_back( backupSV );
      }
      while (!pindelInput.eof()) {
         SVData svd( sampleNames.size() ) ;
         convertIndelToSVdata( pindelInput, sampleMap, genome, svd, chromosomeID);
         if (!pindelInput.eof() && ( chromosomeID=="" || (svd.getChromosome()==chromosomeID && svd.getPosition()>=regionStart && svd.getPosition()<regionEnd)) ) {
            svs.push_back( svd );
         }
         counter++;
         //if (counter%10==0) cout << "At counter " << counter << " pos " << svd.getPosition() << endl;
      }

      cout << "Total reads: " << svs.size() << endl;
      sort ( svs.begin(), svs.end() );
      cout << "Sorting completed" << endl;
      // now output the SVs
      for (int svIndex=0; svIndex<svs.size(); svIndex++ ) {
         //if ( (svIndex+1)<svs.size() && ( svs[ svIndex ] == svs[ svIndex+1 ] ) ) {
         //svs[ svIndex+1 ].fuse( svs[ svIndex ]);
         //}
         //else { // if not fused with the next element, output this element (unless it's the last element, then it must be saved)
         if ( svIndex!=svs.size()-1 && throughFilter( svs[ svIndex ]) ) {
            vcfFile << svs[ svIndex ];
         } else  { //empty

         } // if else: whether the SV passes through the filter
         //} // if else: whether the SV can be fused with the next SV
      }  // for: loop over all SVs
      if (svs.size()>0) {
         backupSV = svs[ svs.size()-1 ];
         backupAvailable = true;
      }
      regionStart += (g_par.windowSize*1000000);
   } while (regionEnd<genome.getChromosome( chromosomeID )->size());
   if ( backupAvailable && throughFilter( backupSV) ) {
      vcfFile << backupSV;
   }
   //cout << "DEBUG:reportSVsInChromosome: exit.\n";
}

int main(int argc, char* argv[])
{
   initBaseArray();
   createParameters();
   readParameters(argc,argv);
   if (argc==1) {
      printHelp();
      exit(EXIT_SUCCESS);
   }
   if (!checkParameters()) {
      exit(EXIT_FAILURE);
   }
   setParameters();
   ofstream vcfFile(g_par.vcffile.c_str());
   set<string> sampleNames;
   set<string> chromosomeNames;

   InputReader pindelInput;
   if (isRegularPindelInput()) {
      pindelInput.addFile( g_par.pindelfile );
   } else if (isRootPindelInput()) {
      string rootFilename = g_par.pindelroot;
      pindelInput.addFile( rootFilename + "_D");
      pindelInput.addFile( rootFilename + "_SI");
      pindelInput.addFile( rootFilename + "_LI");
      pindelInput.addFile( rootFilename + "_INV");
      pindelInput.addFile( rootFilename + "_TD");
   }
   if (pindelInput.eof()) {
      cout << "The pindel file (-p) does not exist.\n";
      exit( EXIT_FAILURE );
   }
   //cout<< "Samples0:\n";
   //showSet( sampleNames );
   getSampleNamesAndChromosomeNames(pindelInput,sampleNames,chromosomeNames);
   cout<< "Samples:\n";
   showSet( sampleNames );
   cout << "Chromosomes in which SVs have been found:\n";
   showSet( chromosomeNames );

   map< string, int > sampleMap;
   makeSampleMap( sampleNames, sampleMap );
   createHeader(vcfFile,"pindel",g_par.referenceName, sampleNames);

   // read in reference
   Genome genome;
   readReference(g_par.reference,genome);

   if (g_par.chromosome != "" ) {
      // a specific chromosome has been specified
      reportSVsInChromosome( g_par.chromosome, chromosomeNames, sampleNames, pindelInput, sampleMap, genome, vcfFile );
   }

   else for (int chromosomeCount=0; chromosomeCount<genome.d_chromosomes.size(); chromosomeCount++ ) {
         reportSVsInChromosome( genome.d_chromosomes[ chromosomeCount ].getID(), chromosomeNames, sampleNames, pindelInput, sampleMap, genome, vcfFile );
         genome.d_chromosomes[ chromosomeCount ].removeFromMemory(); // to prevent memory overload
      }

   if (g_par.gatkCompatible) {
      cout << "\nNote: for this conversion, the -G (GATK) compatibility option was used; it's possible that this format is not compatible with other VCF-reading software. " <<
           "Also note that GATK requires genotypes to be 0/0, 0/1 or 1/1 instead of undefined, like ./1 or . ('not detected'). However, since pindel cannot yet " <<
           "genotype events " <<
           "(distinguish between 0/1 and 1/1) all events are called as 0/0 (not found) or 0/1, even while some may very well be homozygous alternative (1/1).\n\n";
   } else {
      cout << "\nNote: for this conversion, the -G (GATK) compatibility option was not used; while this allows Pindel to indicate the uncertainty in genotypes, and should be " <<
           "compatible with most software, this format " <<
           "will not be compatible with GATK pipelines and tools such as GATK ValidateVariants; if you wish to input the vcf-file into the GATK pipeline, " <<
           "please use the -G option.\n\n";
   }
}

#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

double fisher22_1sided(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t m11_is_greater_alt)
{
   double cur_prob = EXACT_TEST_BIAS;
   double left_prob = cur_prob;
   double right_prob = 0;
   uint32_t uii;
   double cur11;
   double cur12;
   double cur21;
   double cur22;
   double preaddp;
   // Ensure m11 <= m22 and m12 <= m21.
   if (m12 > m21) {
      uii = m12;
      m12 = m21;
      m21 = uii;
   }
   if (m11 > m22) {
      uii = m11;
      m11 = m22;
      m22 = uii;
   }
   // Flipping m11<->m12 and m21<->m22 also flips the direction of the
   // alternative hypothesis.  So we flip on m11-is-greater alternative
   // hypothesis here to allow the rest of the code to assume m11-is-less.
   if (m11_is_greater_alt) {
      uii = m11;
      m11 = m12;
      m12 = uii;
      uii = m21;
      m21 = m22;
      m22 = uii;
   }
   cur11 = m11;
   cur12 = m12;
   cur21 = m21;
   cur22 = m22;
   if ((((uint64_t)m11) * m22) >= (((uint64_t)m12) * m21)) {
      // starting right of (or at) center, p > 0.5
      // 1. left_prob = sum leftward to precision limit
      // 2. total_prob := left_prob
      // 3. total_prob += sum rightward to total_prob precision limit
      // return left_prob / total_prob
      while (cur11 > 0.5) {
         cur12 += 1;
         cur21 += 1;
         cur_prob *= (cur11 * cur22) / (cur12 * cur21);
         cur11 -= 1;
         cur22 -= 1;
         preaddp = left_prob;
         left_prob += cur_prob;
         if (left_prob <= preaddp) {
            break;
         }
         if (left_prob >= 1.0) {
            // Probability mass of our starting table was represented as 2^{-83},
            // so this would mean the left probability mass partial sum is greater
            // than 2^83 times that.  In which case the final p-value will
            // be indistinguishable from 1 at 53-bit precision if our input just
            // had 32-bit integers.  (Yes, the constant can be reduced.)
            return 1;
         }
      }
      cur11 = m11;
      cur12 = m12;
      cur21 = m21;
      cur22 = m22;
      cur_prob = EXACT_TEST_BIAS;
      right_prob = left_prob; // actually total_prob
      while (cur12 > 0.5) {
         cur11 += 1;
         cur22 += 1;
         cur_prob *= (cur12 * cur21) / (cur11 * cur22);
         cur12 -= 1;
         cur21 -= 1;
         preaddp = right_prob;
         right_prob += cur_prob;
         if (right_prob <= preaddp) {
            break;
         }
      }
      return left_prob / right_prob;
   } else {
      // starting left of center, p could be small
      // 1. right_prob = sum rightward to precision limit
      // 2. left_prob = sum leftward to left_prob precision limit
      // return left_prob / (left_prob + right_prob)
      while (cur12 > 0.5) {
         cur11 += 1;
         cur22 += 1;
         cur_prob *= (cur12 * cur21) / (cur11 * cur22);
         cur12 -= 1;
         cur21 -= 1;
         preaddp = right_prob;
         right_prob += cur_prob;
         if (right_prob == INFINITY) {
            return 0;
         }
         if (right_prob <= preaddp) {
            break;
         }
      }
      cur11 = m11;
      cur12 = m12;
      cur21 = m21;
      cur22 = m22;
      cur_prob = EXACT_TEST_BIAS;
      while (cur11 > 0.5) {
         cur12 += 1;
         cur21 += 1;
         cur_prob *= (cur11 * cur22) / (cur12 * cur21);
         cur11 -= 1;
         cur22 -= 1;
         preaddp = left_prob;
         left_prob += cur_prob;
         if (left_prob <= preaddp) {
            break;
         }
      }
      return left_prob / (left_prob + right_prob);
   }
}
