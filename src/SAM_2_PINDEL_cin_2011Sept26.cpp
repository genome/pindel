// Kai Ye, 25 June 2010
// C++ program to convert sam file to pindel input

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <list>
#include <iomanip>
#include <cmath>
#include <math.h>


using namespace std;

//bool Report = false;
// some global variables and constants:
const int alphs = 4;
const char alphabet[alphs] = {'A','C','G','T'};

const char ForwardStrand = '+';
const char ReverseStrand = '-';

char Convert2RC4N[256];

const string AbnormalAlignment = "IDNSH";

string ReverseComplement(const string & InputPattern);

const short Len_Check = AbnormalAlignment.size();
struct READ {
   string QNAME;
   unsigned FLAG;
   string RNAME;
   unsigned POS;
   short MAPQ;
   string CIGAR;
   string MRNM;
   unsigned MPOS;
   int ISIZE;
   string SEQ;
   string QUAL;
   string OPT;
   bool Paired;
   bool ProperPaired;
   bool QueryUnmapped;
   bool MateUnmapped;
   char StrandOfQuery; // 1 for reverse;
   char StrandOfMate;
   bool FirstRead;
   bool SecondRead;
   bool NotPrimaryAlignment;
   bool FailQualityCheck;
   bool Duplication;
   bool Report;
};

void ReadInRead(ifstream & inf_ReadsSeqs,
                READ & Left,
                READ & Right,
                ofstream & OutOneMapped);

void GetFlag(READ & input);

short WhetherReport(const READ & input);

bool WhetherN(const string & input);
int main (int argc, char *argv[])
{
   if (argc != 7) {
      cout << "Welcome to sam2pindel.\n\n"
           << "6 parameters are required here:\n"
           << "1. Input sam file or - for cin.\n"
           << "2. Output for pindel.\n"
           << "3. insert size.\n"
           << "4. tag.\n"
           << "5. number of extra lines (not start with @) in the beginning of the file. \n"
           << "6. Which sequence platform: Illumina-PairEnd or Illumina-MatePair.\n\n"
           << "Usage: ./sam2pindel input.sam output.pindel 500 tumor 0 Illumina-PairEnd\n"
           << "       samtools view input.bam | ./sam2pindel - output.pindel 500 tumor 0 Illumina-PairEnd\n"
           << endl;
      return 1;
   }


   unsigned CountOneEndMapped = 0;
   unsigned CountDifficultMapped = 0;
   int Num_Lines = atoi(argv[5]);
   const string AbnormalAlignment = "NSIDH";
   const short Len_Check = AbnormalAlignment.size();

   // #################################################################

   ifstream  inf_Read;   // input file name
   string argv_str(argv[1]);
   int argv1 = argv_str.compare("-"); // 0 for cin, others for argv[1]
   //for (short i = 1; i < 5; i++) cout << argv[i] << endl;
   if(argv1==0) {
      if(!cin) {
         cout << "Sorry, there is actually no standard input. " << endl;
      }
   }
   else {
      inf_Read.open(argv[1],ifstream::in);
      if(!inf_Read) {
         cout << "Sorry, cannot find the file: " << argv[1] << endl;
      }
   }

   string OutOneMappedFilename = argv[2];
   ofstream OutOneMapped(OutOneMappedFilename.c_str());
   if (!OutOneMapped) {
      cout << "Sorry, cannot write to the file: " << argv[2] << endl;
      return 1;
   }


   const short InsertSize = short(atoi(argv[3]));
   const string Tag = argv[4];

   const string Platform = argv[6];
   if (Platform != "Illumina-PairEnd" && Platform != "Illumina-MatePair") {
      cout << "Please indicate either Illumina-PairEnd or Illumina-MatePair." << endl;
      return 0;
   }
   // ######################### input sequences database ##############



   Convert2RC4N[(short)'A'] = 'T';
   Convert2RC4N[(short)'C'] = 'G';
   Convert2RC4N[(short)'G'] = 'C';
   Convert2RC4N[(short)'T'] = 'A';
   Convert2RC4N[(short)'N'] = 'N';

   READ TempOne;//, TempRight;

   char TempChar;
   string TempStr;
   if (Num_Lines) {
      for (int i = 0; i < Num_Lines; i++) {
         if (argv1 != 0) {
            getline(inf_Read, TempStr);
         }
         else {
            getline(cin, TempStr);
         }
      }
   }

   unsigned Count = 0;
   short ResultWhetherReport;
   unsigned TempPos;
   while (argv1==0 && cin >> TempOne.QNAME || argv1!=0 && inf_Read >> TempOne.QNAME) {
      if (TempOne.QNAME[0] == '@') {
         if (argv1!=0) {
            getline(inf_Read, TempStr);
         }
         else {
            getline(cin, TempStr);
         }
         continue;
      }

      if (argv1 != 0) {
         inf_Read >> TempOne.FLAG >> TempOne.RNAME >> TempOne.POS
                  >> TempOne.MAPQ >> TempOne.CIGAR >> TempOne.MRNM
                  >> TempOne.MPOS >> TempOne.ISIZE >> TempOne.SEQ
                  >> TempOne.QUAL;
         getline(inf_Read, TempOne.OPT);
      }
      else {
         cin >> TempOne.FLAG >> TempOne.RNAME >> TempOne.POS
             >> TempOne.MAPQ >> TempOne.CIGAR >> TempOne.MRNM
             >> TempOne.MPOS >> TempOne.ISIZE >> TempOne.SEQ
             >> TempOne.QUAL;
         getline(cin, TempOne.OPT);
      }
      TempOne.Report = false;

      if (TempOne.MRNM == "=") {
         TempOne.MRNM = TempOne.RNAME;
      }

      GetFlag(TempOne);

      Count++;

      if (Platform == "Illumina-PairEnd") {
         ResultWhetherReport = WhetherReport(TempOne);
         if (TempOne.MateUnmapped || TempOne.MRNM == "*" || TempOne.SEQ == "*") {}
         else {
            if (ResultWhetherReport == 1) { // one-end mapped
               if (TempOne.StrandOfMate == ForwardStrand) {
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << TempOne.SEQ << "\n"
                               //<< TempOne.QUAL << "\n"
                               << TempOne.StrandOfMate << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               else { // anchor read is on reverse strand
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << TempOne.SEQ << "\n"
                               //<< TempOne.QUAL << "\n"
                               << TempOne.StrandOfMate << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS + TempOne.SEQ.size() << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               CountOneEndMapped++;
            }
            else if (ResultWhetherReport == 2) { //  difficult to map or mapped with indels or partially mapped
               if (TempOne.StrandOfMate == ForwardStrand) {
                  //if (TempLeft.POS <= InsertSize) TempPos = 1;
                  //else TempPos = TempLeft.POS - InsertSize;
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << ReverseComplement(TempOne.SEQ) << "\n"
                               //<< TempOne.QUAL << "\n"
                               << TempOne.StrandOfMate << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               else { // reverse strand
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << TempOne.SEQ << "\n"
                               //<< TempLeft.QUAL << "\n"
                               << TempOne.StrandOfMate << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS + TempOne.SEQ.size() << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               CountDifficultMapped++;
            }
         }
      }
      else if (Platform == "Illumina-MatePair") {
         ResultWhetherReport = WhetherReport(TempOne);
         if (TempOne.MateUnmapped || TempOne.MRNM == "*" || TempOne.SEQ == "*") {}
         else {
            if (ResultWhetherReport == 1) { // one-end mapped
               if (TempOne.StrandOfMate == ForwardStrand) {
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << ReverseComplement(TempOne.SEQ) << "\n"
                               //<< TempOne.QUAL << "\n"
                               << ReverseStrand << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               else { // anchor read is on reverse strand
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << ReverseComplement(TempOne.SEQ) << "\n"
                               //<< TempOne.QUAL << "\n"
                               << ForwardStrand << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS + TempOne.SEQ.size() << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               CountOneEndMapped++;
            }
            else if (ResultWhetherReport == 2) { //  difficult to map or mapped with indels or partially mapped
               if (TempOne.StrandOfMate == ForwardStrand) {
                  //if (TempLeft.POS <= InsertSize) TempPos = 1;
                  //else TempPos = TempLeft.POS - InsertSize;
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << TempOne.SEQ << "\n"
                               //<< TempOne.QUAL << "\n"
                               << ReverseStrand << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               else { // reverse strand
                  OutOneMapped << '@' << TempOne.QNAME << "\n" << ReverseComplement(TempOne.SEQ) << "\n"
                               //<< TempLeft.QUAL << "\n"
                               << ForwardStrand << "\t" << TempOne.MRNM << "\t"
                               << TempOne.MPOS + TempOne.SEQ.size() << "\t" << TempOne.MAPQ << "\t"
                               << InsertSize << "\t" << Tag << "\n";
               }
               CountDifficultMapped++;
            }
         }
      }

   }
   cout << Count << "\t" << CountOneEndMapped << "\t" << CountDifficultMapped << endl;

   return 0;
}//main


bool WhetherN(const string & input)
{
   for (unsigned int i = 0; i < input.size(); i++)
      if (input[i] == 'N') {
         return true;
      }

   return false;
}

void GetFlag(READ & input)
{
   if (input.FLAG % 2 == 1) {
      input.Paired = true;
   }
   else {
      input.Paired = false;
   }
   if ((input.FLAG / 2) % 2 == 1) {
      input.ProperPaired = true;
   }
   else {
      input.ProperPaired = false;
   }
   if ((input.FLAG / 4)  % 2 == 1) {
      input.QueryUnmapped = true;
   }
   else {
      input.QueryUnmapped = false;
   }
   if ((input.FLAG / 8)  % 2 == 1) {
      input.MateUnmapped = true;
   }
   else {
      input.MateUnmapped = false;
   }
   if ((input.FLAG / 16)  % 2 == 1) {
      input.StrandOfQuery = ReverseStrand;
   }
   else {
      input.StrandOfQuery = ForwardStrand;
   }
   if ((input.FLAG / 32)  % 2 == 1) {
      input.StrandOfMate = ReverseStrand;
   }
   else {
      input.StrandOfMate = ForwardStrand;
   }
   if ((input.FLAG / 64)  % 2 == 1) {
      input.FirstRead = true;
   }
   else {
      input.FirstRead = false;
   }
   if ((input.FLAG / 128) % 2 == 1) {
      input.SecondRead = true;
   }
   else {
      input.SecondRead = false;
   }
   if ((input.FLAG / 256) % 2 == 1) {
      input.NotPrimaryAlignment = true;
   }
   else {
      input.NotPrimaryAlignment = false;
   }
   if ((input.FLAG / 512) % 2 == 1) {
      input.FailQualityCheck = true;
   }
   else {
      input.FailQualityCheck = false;
   }
   if ((input.FLAG / 1024) % 2 == 1) {
      input.Duplication = true;
   }
   else {
      input.Duplication = false;
   }
}

short WhetherReport(const READ & input)
{
   //if (input.Duplication || input.FailQualityCheck) return 0;
   short NumN = 0;
   for (int i = 0; i < input.SEQ.size(); i++) {
      if (input.SEQ[i] == 'N') {
         NumN++;
      }
   }
   if (NumN * 10 > input.SEQ.size()) {
      return 0;
   }
   if (input.QueryUnmapped == true && input.MateUnmapped == false) {
      return 1;   //  one-end mapped
   }
   for (short i = 0; i < input.CIGAR.size(); i++) {
      if (input.CIGAR[i] >= 'A' && input.CIGAR[i] <= 'Z' && input.CIGAR[i] != 'M') {
         if (input.MateUnmapped == false) {
            return 2;   // mapped with difficulties
         }
      }
   }
   return 0;
}

string ReverseComplement(const string & InputPattern)
{
   string OutputPattern = InputPattern;
   unsigned int LenPattern = InputPattern.size();
   for (unsigned int j = 0; j < LenPattern; j++) {
      OutputPattern[j] = Convert2RC4N[(unsigned int)InputPattern[LenPattern - j - 1]];
   }
   return OutputPattern;
}


