/* vcfcreator.cpp

	Transforms indel and SV-files created by Pindel into VCF4.0 format according to 1000-genome SV specifications.
	This version sorts the records before outputting; may need to make version that uses external memory later!

   Created by Eric-Wubbo Lameijer, section of Molecular Epidemiology, Leiden University Medical Center, March 3rd, 2011.
	e.m.w.lameijer@lumc.nl
	+31(0)71-526 9745

	Version 0.1.7 [July 19th, 2011. Has updated/corrected help text, mentioning -d] */

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
#include <stdlib.h> // for atoi
#include <time.h>

const int FIRST_SAMPLE_INDEX = 32; // index of first sample name

using namespace std;

string g_versionString = "0.1.7";
string g_programName = "pindel2vcf";

bool g_normalBaseArray[256];

struct ParameterCollection {
	string reference;
	string referenceName;
	string referenceDate;
	string pindelfile;
	string vcffile;
	string chromosome;
	int minsize;
	int maxsize;
	bool bothstrands;
	int minsuppSamples;
	int minsuppReads;
	int regionStart;
	int regionEnd;
	int maxHomopolyRepeats;
	int maxHomopolyLength;
	bool showHelp;
} par;

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
	void set() { 
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
	cout << d_shortName << "/" << d_longName << "  " << d_description;
	if ( d_required ) { cout << ": required parameter" ; }
	cout << endl; 
}

vector<Parameter*> parameters;


/* 'Chromosome' contains a string identifier as well as the base sequence of the chromosome itself. */
class Chromosome {
public:
	Chromosome(const string& identifier) { d_identifier = identifier; d_sequence="";}
	const string* getChromPtr() const { return &d_sequence; }
	const string getID() const { return d_identifier; }
	void addBases(const string& bases ) { d_sequence += bases; }
	void makeStrangeBasesN();

private:
	string d_identifier;
	string d_sequence;
};

bool normalBase(char ch )
{
	return g_normalBaseArray[ ch ];
}

void Chromosome::makeStrangeBasesN()
{
	int chromLength = d_sequence.size();
	for (int position=0; position<chromLength; position++ ) {
		if (!normalBase(d_sequence[ position ])) {
			d_sequence[ position ] ='N';
		}
	}

}

/* 'Genome' contains a collection of chromosomes, and returns a pointer to the requested chromosome (the chromosome with the requested ID) */
class Genome {
public:
	const string* getChromosome( const string& id ) const; 
	void addChromosome( const Chromosome& chr ) { d_chromosomes.push_back( chr ); }
	string firstChromosomeName();
	string nextChromosomeName();

	vector<Chromosome> d_chromosomes; // working fast here, but ideally I'd keep this private and implement an iterator
};


const string* Genome::getChromosome( const string& id ) const
{

	for (int chromosomeIndex=0; chromosomeIndex<d_chromosomes.size(); chromosomeIndex++ ) {
		if ( id.compare( d_chromosomes[ chromosomeIndex ].getID() ) == 0 ) { return d_chromosomes[ chromosomeIndex ].getChromPtr(); }
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
   outFile << "##fileDate=" << par.referenceDate << endl;

   // source
   outFile << "##source=" << sourceProgram << endl;

   // reference
   outFile << "##reference=" << reference << endl;

   // info fields (selected from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/VCF%20%28Variant%20Call%20Format%29%20version%204.0/encoding-structural-variants)
   outFile << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
   outFile << "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">" << endl;
   outFile << "##INFO=<ID=HOMSEQ,Number=-1,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">" << endl;
   outFile << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
   outFile << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << endl;
	outFile << "##INFO=<ID=NTLEN,Number=.,Type=Integer,Description=\"Number of bases inserted in place of deleted code\">" << endl;
   //outFile << "##ALT=<ID=DEL,Description=\"Deletion\">" << endl; /*EWL040311: probably not needed, as our calls are precise so we should rather give the exact replacing sequence instead of a label. */
   //outFile << "##ALT=<ID=DUP,Description=\"Duplication\">" << endl;
   //outFile << "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">" << endl;
   //outFile << "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">" << endl;
   //outFile << "##ALT=<ID=INV,Description=\"Inversion\">" << endl;
   //outFile << "##ALT=<ID=CNV,Description=\"Copy number variable region\">" << endl;
   outFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
   outFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
   outFile << "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of supporting reads\">" << endl; // To INFO?

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


/* 'StringCollection' can save various separate strings, outputting them when needed as one ';'-separated string, but also offering access to its individual 	
	 elements. */
class StringCollection {
public:
   StringCollection() {};
   string totalString() const;   
   void addString(const string& inputString) { d_contents.push_back( inputString ); };

private:
   vector< string > d_contents;

};

string StringCollection::totalString() const
{
   string returnString="";
   int numItems = d_contents.size();
   if (numItems>=1) { returnString+=d_contents[ 0 ]; }
   for (int counter=1; counter<=numItems; counter++ ) {
      returnString += ";" + d_contents[counter];
   } //for
   return returnString;
}

/* 'Pair' contains a pair of strings, outputted as "str1:str2". e*/
class Pair {
friend ostream& operator<<(ostream& os, const Pair& pair);

public:
   Pair(const string& first, const string& second) { d_first=first; d_second=second; }


private:
   string d_first, d_second; 

};

ostream& operator<<(ostream& os, const Pair& pair)
{
   os << pair.d_first << ":" << pair.d_second;
}


/* 'Genotype' stores the genotype of a mutation and its read depth. So, for example, "0/1:5". */
class Genotype {

friend ostream& operator<<(ostream& os, const Genotype& gt);

public:
	Genotype();
   Genotype( int genFirstChrom, int genSecondChrom, int readDepthPlus, int readDepthMinus );
	void fuse( const Genotype& gt );
	void reset() { d_genFirstChromosome=0; d_genSecondChromosome=0; d_readDepthPlus=0; d_readDepthMinus=0; }
	int getReadDepthPlus() const { return d_readDepthPlus; }
	int getReadDepthMinus() const { return d_readDepthMinus; }
	int getTotalReads() const { return ( d_readDepthPlus + d_readDepthMinus ) ; }

private:
   int d_genFirstChromosome;
   int d_genSecondChromosome;
   int d_readDepthPlus;
	int d_readDepthMinus;
};

Genotype::Genotype() 
{
   Genotype( 0, 0, 0, 0 );
}

Genotype::Genotype( int genFirstChrom, int genSecondChrom, int readDepthPlus, int readDepthMinus ) 
{ 
   d_genFirstChromosome = genFirstChrom;
	d_genSecondChromosome = genSecondChrom;
	d_readDepthPlus = readDepthPlus;
	d_readDepthMinus = readDepthMinus;
}

void Genotype::fuse( const Genotype& gt )
{
   d_readDepthPlus += gt.d_readDepthPlus;
	d_readDepthMinus += gt.d_readDepthMinus;
}

ostream& operator<<(ostream& os, const Genotype& gt)
{
  os << gt.d_genFirstChromosome << "/" << gt.d_genSecondChromosome << ":" << ( gt.d_readDepthPlus + gt.d_readDepthMinus );
  return os;
}




/* 'SVData' stores the data of a certain structural variant. */
class SVData {

friend ostream& operator<<(ostream& os, const SVData& svd);

public: 
   SVData(const int genotypeTotal);

	int getPosition() const { return d_position; }
	int getSize() const { return d_svlen; }
	bool bothStrands() const;
	int getNumSupportSamples() const;
	int getNumSupportReads() const;
	int detectHomoPolyRepeats(const int maxHPlen) const;
	string getChromosome() const { return d_chromosome; }

   void setChromosome(const string& chromosome) { d_chromosome = chromosome; }
   void setPosition(const int position) { d_position = position; }
   void setID(const string& id) { d_id = id; }
   void setReference(const string& reference) { d_reference = reference; }
	//void addAlternative(const string& alternative) { d_alternatives.addString( alternative ); }
	void setAlternative(const string& alternative) { d_alternative = alternative ; }
   void setQuality(const double quality) { ostringstream strs; strs << quality; d_quality = strs.str(); }
   void setFilter(const string& filter) { d_filter = filter; }

   void setEnd(const int end) { d_end = end; }
   void setHomlen(const int homlen) { d_homlen = homlen; }
   void setHomseq(const string& homseq) { d_homseq = homseq; }
   void setSVlen(const int svlen) { d_svlen = svlen; }
	void setSVlen(const string& svlen) { d_svlen = atoi(svlen.c_str()); }
   void setSVtype(const string& svtype) { d_svtype = svtype; }

   //void setSupportingReads(const int supportingreads) { d_supportingreads = supportingreads; }
	void addGenotype(const int sampleID, const int readDepthPlus, const int readDepthMinus) { 
		int secondGenome = ( ( (readDepthPlus+readDepthMinus) > 0) ? 1 : 0 );
		Genotype gt(0,secondGenome,readDepthPlus,readDepthMinus); 
		d_format[ sampleID ] = gt ; }

	void setReplace( int replaceLength, int secondReplaceLen=-1 ) { d_replaceLen = replaceLength; d_replaceLenTwo = secondReplaceLen; }
	void setBPrange( const int bpr_start, const int bpr_end ) { d_bpr_start = bpr_start; d_bpr_end = bpr_end ; }

	bool operator<(const SVData& otherSV ) const;
	bool operator==(const SVData& otherSV ) const;
	void fuse( SVData& otherSV );
private:
   string d_chromosome;
   int d_position; // 1-based
	string d_id; // default '.', as we don't mine variant databases yet
   string d_reference; // reference sequence, for indels including the base before it
   //StringCollection d_alternatives;
	string d_alternative;
   string d_quality; // '.' by default, but can be floating-point number
   string d_filter;  // "PASS" by default
   int d_end;
   int d_homlen;
	int d_bpr_start, d_bpr_end;
	string d_homseq;
   int d_svlen;
   string d_svtype;
	int d_replaceLen;
	int d_replaceLenTwo; // for inversions
   vector<Genotype> d_format;
   //int d_supportingreads;
};

SVData::SVData(const int genotypeTotal) // default settings
{
	d_id=".";
	d_quality=".";
	d_filter="PASS";
	d_replaceLen=0;
	d_homlen=0;
	d_homseq="";
	int numberOfSamples = genotypeTotal;
	if (numberOfSamples<=0) numberOfSamples=1;
	d_format.resize( numberOfSamples, Genotype(0,0,0,0) );
}; 

bool SVData::bothStrands() const
{
   bool strandPlus = false; 
	bool strandMinus = false;

  	for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
   	if ( d_format[ sampleIndex ].getReadDepthPlus() > 0 ) { strandPlus = true; } 
   	if ( d_format[ sampleIndex ].getReadDepthMinus() > 0 ) { strandMinus = true; } 
  	}
	return ( strandPlus && strandMinus );
}

int SVData::getNumSupportSamples() const
{
	int numSamples = 0;
	for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
   	if ( d_format[ sampleIndex ].getTotalReads() > 0 ) { numSamples++; }
  	}
	return numSamples;
}

int SVData::getNumSupportReads() const
{
	int numReads = 0;
	for (int sampleIndex=0; sampleIndex<d_format.size(); sampleIndex++ ) {
   	numReads += d_format[ sampleIndex ].getTotalReads();
  	}	
	return numReads;
}

bool SVData::operator<(const SVData& otherSV ) const 
{ 
	if ( d_chromosome.compare( otherSV.d_chromosome ) != 0 ) return ( d_chromosome.compare( otherSV.d_chromosome ) < 0 );	
	if ( d_position != otherSV.d_position) return ( d_position < otherSV.d_position ); 
	else return ( d_svlen < otherSV.d_svlen ); // position equal: then do smallest SV first
}

/* 'testHypothesis' tests whether the DNA string consists of a number of repeating "hypothesis" strings. If so, it
	returns the number of repeats; otherwise, it returns zero. */
int testHypothesis(const string& hypothesis, const string& dna )
{
	int hypLen = hypothesis.size();
	int dnaLen = dna.size();
	bool hypothesisCorrect = true;
	int repeatCounter = 0;	
	// this version assumes that there must be a whole number of repeats; this can be corrected/modified later
	if ( dnaLen % hypLen != 0 ) return 0;
	else {
		for (int repeat=0; repeat<dnaLen/hypLen; repeat++ ) {
			string checkedDNA = dna.substr( repeat*hypLen, hypLen );
			if ( checkedDNA.compare(hypothesis) != 0 ) { hypothesisCorrect = false; }
		} // <for
		return ( hypothesisCorrect ? dnaLen/hypLen : 0 );
	} // <else
}

int countHPRepeats( const string& bases, const int maxRepeatLength )
{
	int maximumLen = min( maxRepeatLength, (int)bases.size() );
	string hypothesis = "";
	for (int repeatLen=1; repeatLen<maximumLen; repeatLen++ ) {
      hypothesis += bases[ repeatLen - 1 ];
		int repeats = testHypothesis( hypothesis, bases );	
		if (repeats>0) return repeats;
	}
	return 0;
}

int SVData::detectHomoPolyRepeats(const int maxHPlen) const
{
	int refHP = countHPRepeats( d_reference, maxHPlen );
	int altHP = countHPRepeats( d_alternative, maxHPlen );
	int hseqHP = countHPRepeats( d_homseq, maxHPlen );
	//if ( max( refHP, max ( altHP, hseqHP ) ) > 0) { cout << "REPEATS!" << d_reference << " " << d_alternative << " " << d_homseq << endl; }
	return max( refHP, max ( altHP, hseqHP ) );
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
			) return true;
   // for NT
	if ( ( d_svtype == "RPL" )  
		  && ( otherSV.d_svtype == "RPL" ) 
		  && ( ( d_svlen - d_replaceLen) == ( otherSV.d_svlen - otherSV.d_replaceLen ) ) 
   	  && ( d_bpr_start == otherSV.d_bpr_start )
		  && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
		  ) return true;
   if ( ( d_svtype == "INS" )  
		  && ( otherSV.d_svtype == "INS" ) 
		  //&& ( ( d_svlen - d_replaceLen) == ( otherSV.d_svlen - otherSV.d_replaceLen ) ) 
   	  && ( d_bpr_start == otherSV.d_bpr_start )
		  && ( d_bpr_end == otherSV.d_bpr_end )
		  && ( d_svlen == otherSV.d_svlen )
		  && ( d_chromosome.compare( otherSV.d_chromosome ) == 0 )
		  ) return true;
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


ostream& operator<<(ostream& os, const SVData& svd)
{
   os << svd.d_chromosome << "\t";
   os << svd.d_position << "\t";
   os << svd.d_id << "\t";
   os << svd.d_reference << "\t";
   os << svd.d_alternative << "\t";
   os << svd.d_quality << "\t";
   os << svd.d_filter << "\t";

   os << "END=" << svd.d_end << ";";
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
	

   os << "\tGT:DP";

   for (int counter=0; counter<svd.d_format.size(); counter++ ) {
      os << "\t" << svd.d_format[ counter ];
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

void getSampleNames(ifstream& svfile, set<string>& sampleNames)
{
	string line;
   while (!svfile.eof()) {
	   do {
   	   getline( svfile, line );
	   } while (!svfile.eof() && !isdigit(line[0])); // skip ###-lines
   	if (svfile.eof()) return; 
	   stringstream lineStream;
   	lineStream << line;
		string svType = fetchElement( lineStream, 2 ); 

		if ( svType.compare("LI")==0 ) {	
			sampleNames.insert("TOTAL");
			return;
		}
		string firstSampleName = fetchElement( lineStream, FIRST_SAMPLE_INDEX-2 );
		sampleNames.insert( firstSampleName );
		string newSampleName = fetchElement( lineStream, 5 );
		while (!lineStream.fail()) {
   	   sampleNames.insert( newSampleName );
			newSampleName = fetchElement( lineStream, 5 );
		}   
	}
}

template<class T> void showVector( vector<T> vect )
{
   cout << "Printing vector: " << endl;
   for (int index=0; index<vect.size(); index++ ) {
     cout << index << " " << vect[ index ] << endl;
   }
   cout << "End of printing.\n";
}

void showSet( set<string> aSet )
{
   cout << "Printing set: " << endl;

	set<string>::iterator index;
	int counter=1;
   for (index=aSet.begin(); index!=aSet.end(); index++ ) {
      cout << counter++ << " " << *index << endl;
   }
   cout << "End of printing.\n";
}

/* returns the complementary DNA-base of base 'inputbase' */
char complementBase( char inputBase )
{
   switch (inputBase) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      default: return 'N';
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




/* 'convertIndelToSVdata' converts insertions and deletions to nicely formatted SV-data. */
void convertIndelToSVdata( ifstream& svfile, map< string, int>& sampleMap, const Genome& genome, SVData& svd)
{
   string line;
	do {
   	getline( svfile, line );
	} while (!svfile.eof() && !isdigit(line[0]));

   if (svfile.eof()) return;

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
		if (chromosomeID!=par.chromosome) {
			return;
		}
		int beforeStartPos = atoi( fetchElement( lineStream, 1 ).c_str() ); 
		svd.setPosition( beforeStartPos );
		string refVariant = "";
		refVariant+= (*reference)[ beforeStartPos ];
		svd.setReference( refVariant );
		svd.setAlternative( "<INS>" );
		int plusSupport = atoi( fetchElement( lineStream, 1 ).c_str()); 	
		int rightmostEndPos = atoi (fetchElement( lineStream, 1 ).c_str()); // now at position 14
	   svd.setEnd( rightmostEndPos );	
	   svd.setBPrange( beforeStartPos, rightmostEndPos );
   	int minSupport = atoi( fetchElement( lineStream, 1 ).c_str());
		svd.addGenotype( 0, plusSupport , minSupport );	
		return;
	}
   svd.setSVlen( fetchElement( lineStream, 1 ) ); // to 3

	// get number(s) of NT bases added (two numbers for inversions!)
	string numNTaddedStr = fetchElement( lineStream, 2 ); // to 5
	int numNTadded = atoi( numNTaddedStr.c_str() ); // should get first number
	int numNTinvAdded=-1;
	if ( svType.compare("INV") == 0 ) { // two numbers separated by : instead of one
		int separatorPos = numNTaddedStr.find(":");
		string secondNumber = numNTaddedStr.substr(separatorPos+1);
		numNTinvAdded = atoi( secondNumber.c_str() ); 
	} 

   string ntAdded = fetchElement( lineStream, 1 ); // to 6
	string ntInvAdded = "";
	if ( svType.compare("INV") == 0 ) {
		int separatorPos = ntAdded.find(":");
		ntInvAdded = ntAdded.substr( separatorPos+2, numNTinvAdded ); // erases ""
		ntAdded = ntAdded.substr(0,separatorPos);
	}
   ntAdded.erase(0,1); // erases opening "
	ntAdded.erase(numNTadded); // erases closing "
	string chromosomeID = fetchElement( lineStream, 2); // now at position 8
	if (chromosomeID!=par.chromosome) {
		return;
	}
	const string* reference = genome.getChromosome( chromosomeID );
	if ( reference== NULL ) {
		cout << "Error! Reference chromosome \"" << chromosomeID << "\" not found!" << endl;
		exit(EXIT_FAILURE);
	} 
   svd.setChromosome( chromosomeID ); 
   int beforeStartPos = atoi( fetchElement( lineStream, 2 ).c_str() ); // pos 10
	svd.setPosition( beforeStartPos ); // now at position 10
   int leftmostEndPos = atoi( fetchElement( lineStream, 1 ).c_str()); // now at position 11
   string refVariant="";
	for (int position=beforeStartPos; position<leftmostEndPos; position++ ) {
		refVariant += (*reference)[ position ];
	}
	svd.setReference( refVariant );
   string altVariant = "";
	if ( svType.compare( "I" ) == 0 || svType.compare( "D" ) == 0 ) {
	   altVariant += (*reference)[ beforeStartPos ];
		altVariant += ntAdded;	
	}
	else if ( svType.compare("TD")==0 ) {
		altVariant += refVariant;
		altVariant += ntAdded;
		altVariant += refVariant.substr(1); // skip the position before...
	} else if ( svType.compare("INV") == 0 ) {
		altVariant += (*reference)[ beforeStartPos ];
		altVariant += ntAdded;

		string referenceInv = refVariant.substr(1);
		string complement = "";
		createComplement( referenceInv, complement );

		altVariant += complement;
		altVariant += ntInvAdded;
	}

   svd.setAlternative( altVariant ); 
	int leftmostStartPos = atoi (fetchElement( lineStream, 2 ).c_str());  // at position 13
   int rightmostEndPos = atoi (fetchElement( lineStream, 1 ).c_str()); // now at position 14
   svd.setBPrange( leftmostStartPos, leftmostEndPos );
   svd.setEnd( rightmostEndPos );
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
	} else if ( svType.compare("I")==0 ){ 
		svd.setSVtype("INS");
		svd.setReplace( 0 ); 
	} else if ( svType.compare("TD")==0 ){ 
		svd.setSVtype("DUP:TANDEM");
		svd.setReplace( numNTadded );
	} else if ( svType.compare("INV") == 0 ){
		svd.setSVtype("INV");
		svd.setReplace( numNTadded, numNTinvAdded );
	}
	string sampleName = fetchElement( lineStream, 18);
   int plusSupport = atoi( fetchElement( lineStream, 1 ).c_str()); // now at position 33, total +supports sample 1; for unique support 1->2
   int minSupport = atoi( fetchElement( lineStream, 2 ).c_str()); // now at position 35, total +supports sample 1
//cout << "SN PS MS" << sampleName << plusSupport << minSupport << endl;
int count=0;
   while (!lineStream.fail()) {
		if (sampleMap.find( sampleName )==sampleMap.end() ) {
			cout << "Error: could not find sample " << sampleName << endl;
		} 
		else {
			int sampleID = sampleMap[ sampleName ];
			//cout << "Found sample " << sampleName << " with number " << sampleID << endl;
		   svd.addGenotype( sampleID, plusSupport , minSupport );
		}
		sampleName = fetchElement( lineStream, 2); // for unique support, 2->1, 
		plusSupport = atoi( fetchElement( lineStream, 1 ).c_str()); // for unique support, 2->1
   	minSupport = atoi( fetchElement( lineStream, 2 ).c_str()); // now at position 33, total +supports sample 1
		//cout << "SN PS MS" << sampleName << plusSupport << minSupport << endl;
		//cout << "Adding sample: " << ++count << endl;
	}
	//cout << "Added a " << svType << "!\n";
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

   getline(referenceFile,refLine); // FASTQ format always has a first line with the name of the reference in it
	// loop over each chromosome
	do {

   	int counter=1;
		refName = "";
		do {
			refName += refLine[ counter++ ]; 
		} while ( counter<refLine.size() && (refLine[ counter ] != ' ') && (refLine[ counter ] != '\t') && (refLine[ counter ] != '\n') );
		cout << "Name of current chromosome: " << refName << endl;
		Chromosome newChrom( refName );
		newChrom.addBases("N"); // for 1-shift
   	getline(referenceFile,currentLine);
	   while (!referenceFile.eof() && currentLine[0]!='>') {
	      newChrom.addBases( currentLine );
	      getline(referenceFile,currentLine);  
	   }
		newChrom.makeStrangeBasesN(); 
		genome.addChromosome( newChrom );
		refLine = currentLine;
	} while (!referenceFile.eof());
	cout << "Exiting reference reading.\n";
}

/* 'createParameters' creates the default parameters that the VCF converter uses. */
void createParameters()
{
	parameters.push_back(
		new StringParameter( &par.reference, "-r", "--reference", "The name of the file containing the reference genome", true, "" ) );
	parameters.push_back(
		new StringParameter( &par.referenceName, "-R", "--reference_name", "The name and version of the reference genome", true, "" ) );
	parameters.push_back(
		new StringParameter( &par.referenceDate, "-d", "--reference_date", "The date of the version of the reference genome used", true, "" ) );
	parameters.push_back(
		new StringParameter( &par.pindelfile, "-p", "--pindel_output", "The name of the pindel output file containing the SVs", true, "" ) );
   parameters.push_back(
		new StringParameter( &par.vcffile, "-v", "--vcf", "The name of the output vcf-file (default: name of pindel output file +\".vcf\"", false, "" ) );
   parameters.push_back(
		new StringParameter( &par.chromosome, "-c", "--chromosome", "The name of the chromosome (default: SVs on all chromosomes are processed)", false, "" ) );
	parameters.push_back(
		new IntParameter( &par.minsize, "-is", "--min_size", "The minimum size of events to be reported (default 1)", false, 1 ) );
	parameters.push_back(
		new IntParameter( &par.maxsize, "-as", "--max_size", "The maximum size of events to be reported (default infinite)", false, -1 ) );
	parameters.push_back(
		new BoolParameter( &par.bothstrands, "-b", "--both_strands_supported", "Only report events that are detected on both strands (default false)", false, false ) ); 
	parameters.push_back(
		new IntParameter( &par.minsuppSamples, "-m", "--min_supporting_samples", "The minimum number of samples an event needs to occur in to be reported (default 1)", false, 1 ) );
	parameters.push_back(
		new IntParameter( &par.minsuppReads, "-e", "--min_supporting_reads", "The minimum number of supporting reads required for an event to be reported (default 1)", false, 1 ) );
	parameters.push_back(
		new IntParameter( &par.regionStart, "-sr", "--region_start", "The start of the region of which events are to be reported (default 0)", false, 0 ) );
	parameters.push_back(
		new IntParameter( &par.regionEnd, "-er", "--region_end", "The end of the region of which events are to be reported (default infinite)", false, -1 ) );
	parameters.push_back(
		new IntParameter( &par.maxHomopolyRepeats, "-lr", "--max_homopolymer_repeats", "Filters out all homopolymers of more than X repetitions (default infinite)", false, -1 ) );
	parameters.push_back(
		new IntParameter( &par.maxHomopolyLength, "-ll", "--max_homopolymer_length", "The maximum size of a repeating unit to be considered a homopolymer, use with the option -lr. (default infinite)", false, -1 ) );
	parameters.push_back(
		new BoolParameter( &par.showHelp, "-h", "--help", "Print the help of this converter", false, false ) );
}

/* 'findParameter' returns the index of the parameter with name 'name'; -1 if not found.*/
int findParameter(string name)
{
	for (int parameterCounter=0; parameterCounter<parameters.size(); parameterCounter++ ) {
		if (parameters[ parameterCounter ]->hasName( name ) ) { return parameterCounter; }
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

/* 'printHelp' prints all parameters available. */
void printHelp()
{
	cout << "\nProgram:   " << g_programName <<" (conversion of Pindel output to VCF format)\n";
	cout << "Version:   " << g_versionString << endl;
	cout << "Contact:   Eric-Wubbo Lameijer <e.m.w.lameijer@lumc.nl>\n\n";
	cout << "Usage:     " << g_programName << " -p <pindel_output_file> -r <reference_file>\n";
	cout << "              -R <name_and_version_of_reference_genome> -d <date_of_reference_genome_version>\n"; 
   cout << "              [-v <vcf_output_file>]\n\n";
	cout << "           the -v parameter is optional; when no output file name is given, output is written\n";
	cout << "           to a file with the name <pindel_output_file>.vcf.\n\n";
	cout << "Example:   " << g_programName << " -p sample3chr20_D -r human_g1k_v36.fasta -R 1000GenomesPilot-NCBI36\n";
   cout << "              -d 20101123-v sample3chr20_D.vcf\n\n";
	cout << "Note:      -is only guaranteed to work correctly on output files produced by pindel version 0.2.3 and above.\n";
	cout << "           -LI and BP files (long insertion and break point files) have a different type of header and\n";
	cout << "            are not supported yet.\n\n";				 

   for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
		parameters[ parameterIndex ]->describe();
	}
	exit( EXIT_SUCCESS );
}

/* 'checkParameters' checks whether all required parameters have been set. */
bool checkParameters()
{
   if (parameters[ findParameter("-h") ]->getBValue() ) printHelp();

	bool canRun = true;
	for (int parameterIndex=0; parameterIndex<parameters.size(); parameterIndex++ ) {
		if (parameters[ parameterIndex ]->isRequired() && !parameters[ parameterIndex ]->isSet()) {
			cout << "\nRequired parameter " << parameters[ parameterIndex ]->getShortName() << "/" << parameters[ parameterIndex ]->getLongName()
				  << " " << parameters[ parameterIndex ]->getDescription() << " needs to be set.\n\n";
			canRun = false;
		}  //if
	}
	if (!canRun) {
		cout << "For further information, please run " << g_programName <<" without arguments or with option -h/--help.\n\n";
	}
	return canRun;
}

/* 'setParameters' sets the filters to be used in the rest of the program. */
void setParameters()
{
	par.vcffile = parameters[ findParameter( "-v" )]->getSValue();
	if (par.vcffile.compare("")==0) { par.vcffile = par.pindelfile + ".vcf"; } // default
}

/* 'throughFilter' checks whether the event is good enough to be written to the output file. */
bool throughFilter(SVData sv)
{
	if (( par.minsize > 1 ) && ( abs( sv.getSize()) < par.minsize ) ) return false;
	if (( par.maxsize > 0 ) && ( abs( sv.getSize()) > par.maxsize ) ) return false;
	if ( par.bothstrands && !sv.bothStrands() ) return false; 
	if ( ( par.minsuppSamples > 1 ) && ( sv.getNumSupportSamples() < par.minsuppSamples ) ) return false;
	if ( ( par.minsuppReads > 1 ) && ( sv.getNumSupportReads() < par.minsuppReads ) ) return false;
	if ( ( par.regionStart > 0 ) && ( sv.getPosition() < par.regionStart ) ) return false;
	if ( ( par.regionEnd > 0 ) && ( sv.getPosition() > par.regionEnd ) ) return false;
	if ( par.maxHomopolyRepeats > 0 ) { // filter out repeating homopolymers
	   //cout << "Checking repeats: seeing " << sv.detectHomoPolyRepeats(par.maxHomopolyLength);
		if ( sv.detectHomoPolyRepeats(par.maxHomopolyLength) > par.maxHomopolyRepeats ) return false;
	}
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
	for (int i=0;i<256;i++) {
		g_normalBaseArray[i] = false;
	}
	g_normalBaseArray['A'] = true;
	g_normalBaseArray['C'] = true;
	g_normalBaseArray['G'] = true;
	g_normalBaseArray['T'] = true;
	g_normalBaseArray['N'] = true;
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
	if (!checkParameters()) exit(EXIT_FAILURE);
	setParameters();
   ofstream vcfFile(par.vcffile.c_str());
   set<string> sampleNames;
   ifstream svfile(par.pindelfile.c_str());
	if (!svfile) {
		cout << "The pindel file (-p) does not exist.\n";
		exit( EXIT_FAILURE );
	}
   getSampleNames(svfile,sampleNames);
   showSet( sampleNames );

	map< string, int > sampleMap;
	makeSampleMap( sampleNames, sampleMap );
   createHeader(vcfFile,"pindel",par.referenceName, sampleNames);

   // read in reference
	Genome genome;
   readReference(par.reference,genome);

// loop over chromosomes
	for (int chromosomeCount=0; chromosomeCount<genome.d_chromosomes.size(); chromosomeCount++ ) {
		par.chromosome = genome.d_chromosomes[ chromosomeCount ].getID();	
		cout << "Processing chromosome " << par.chromosome << endl;
	   // rewind file to start
	   svfile.clear();
		svfile.seekg(0); 
		int counter=0;
	   vector<SVData> svs;
		while (!svfile.eof()) {
   	   SVData svd( sampleNames.size() );
   	   convertIndelToSVdata( svfile, sampleMap, genome, svd);
   	   if (!svfile.eof() && ( par.chromosome=="" || svd.getChromosome()==par.chromosome ) ) svs.push_back( svd );
			counter++;
//if (counter%10==0) cout << "At counter " << counter << " pos " << svd.getPosition() << endl;
   	}
		cout << "Total reads: " << svs.size() << endl;
		sort ( svs.begin(), svs.end() );
		cout << "Sorting completed" << endl;
	   for (int svIndex=0; svIndex<svs.size(); svIndex++ ) {
			if ( (svIndex+1)<svs.size() && ( svs[ svIndex ] == svs[ svIndex+1 ] ) ) {
			//cout << "Fusion!\n";
			//cout << svs[ svIndex ];
			//cout << svs[ svIndex+1 ];
				svs[ svIndex+1 ].fuse( svs[ svIndex ]);	
			//cout << "results in: \n";
  			//cout << svs[ svIndex+1 ];
			} else {
				if ( throughFilter( svs[ svIndex ]) ) {  vcfFile << svs[ svIndex ]; }
				else { 
				// cout << "Filtered out: " << svs[ svIndex ]; 
				}
			}
		}
	}
}


