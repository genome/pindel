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
#include <map>
#include <set>
#include <sstream>
#include <algorithm>

// Pindel header files
#include "logstream.h"
#include "pindel.h"
#include "logdef.h"
#include "output_file_data.h"
#include "output_sorter.h"
#include "reporter.h"
#include "shifted_vector.h"
//#include "genotyping.h"
#include "bam2depth.h"



OutputFileData deletionFileData;

std::string GetConsensusInsertedStr(const std::vector <SPLIT_READ> & Reads, const int & StartIndex, const int & EndIndex);

/* 'createMaps' (EWL, Aug31st2011) creates maps to get a sample name from an index, and the other way around. */

void createMaps( 	std::map<std::string,int>& sampleToIndexMap, std::map<int,std::string>& indexToSampleMap )
{
   int counter=0;
   for (std::set<std::string>::iterator it=g_sampleNames.begin(); it!=g_sampleNames.end(); it++ ) {
      sampleToIndexMap[ *it ] = counter;
      indexToSampleMap[ counter ] = *it;
      counter++;
   }
}

/* 'calculateSupportPerTag (EWL, Aug31st2011) calculates the number of reads per tag (sample) supporting the event. */
void calculateSupportPerTag( std::vector< SPLIT_READ >& reads, const unsigned int firstReadIndex, const unsigned int lastReadIndex,
                             std::map<std::string,int>& sampleToIndexMap, SupportPerSample* NumSupportPerTag )
{

	std::map <std::string, unsigned>::iterator it;
	int tagIndex;
	for (unsigned int readIndex = firstReadIndex; readIndex <= lastReadIndex; readIndex++) {
	      //std::string currentTag = reads[readIndex].Tag;
      		//int tagIndex = sampleToIndexMap[ currentTag ];
		std::map <std::string, unsigned> * SampleName2Number = &(reads[readIndex].SampleName2Number);
		if (reads[readIndex].MatchedD == Plus)	{
		// = reads[readIndex].SampleName2Number.find(m_rawreads[i].Tag);
			for (it = SampleName2Number -> begin(); it != SampleName2Number -> end() ; ++it) {
				tagIndex = sampleToIndexMap[ it -> first ];
				//std::cout << "plus\t" << it -> first << "\t" << it -> second << std::endl;
				NumSupportPerTag[tagIndex].NumPlus += it -> second;
				if (reads[readIndex].UniqueRead) {
					NumSupportPerTag[tagIndex].NumUPlus += it -> second;
				}
			}
		}
		else {
			for (it = SampleName2Number -> begin(); it != SampleName2Number -> end() ; ++it) {
				tagIndex = sampleToIndexMap[ it -> first ];
				//std::cout << "minus\t" << it -> first << "\t" << it -> second << std::endl;
				NumSupportPerTag[tagIndex].NumMinus += it -> second;
				if (reads[readIndex].UniqueRead) {
					NumSupportPerTag[tagIndex].NumUMinus += it -> second;
				}
			}
		}
	}
}

/* 'calculateSupportPerStrand (EWL, Aug31st2011) calculates the number of reads per strand (sample) supporting the event. */
void calculateSupportPerStrand( SupportPerSample* NumSupportPerTag, unsigned int& LeftS, unsigned int& LeftUNum, unsigned int& RightS, unsigned int& RightUNum )
{
	LeftS = 0;
	LeftUNum = 0;
	RightS = 0;
	RightUNum = 0;
	//std::cout << "g_sampleNames.size(): " << g_sampleNames.size() << std::endl;
	for (unsigned index = 0; index < g_sampleNames.size(); index++) {
		
		LeftS += NumSupportPerTag[index].NumPlus;
		LeftUNum += NumSupportPerTag[index].NumUPlus;
		RightS += NumSupportPerTag[index].NumMinus;
		RightUNum += NumSupportPerTag[index].NumUMinus;
		//std::cout << LeftS << "\t" << LeftUNum << "\t" << RightS << "\t" << RightUNum << std::endl;
	}
}

/*
void calculateSupportPerStrand( std::vector< SPLIT_READ >& reads, const unsigned int firstReadIndex, const unsigned int lastReadIndex,
                                unsigned int& LeftS, unsigned int& LeftUNum, unsigned int& RightS, unsigned int& RightUNum )
{
   for (unsigned int readIndex = firstReadIndex; readIndex <= lastReadIndex; readIndex++) {
      if (reads[readIndex].MatchedD == Plus) {
         LeftS++;
         if (reads[readIndex].UniqueRead) {
            LeftUNum++;
         }
      }
      else {
         RightS++;
         if (reads[readIndex].UniqueRead) {
            RightUNum++;
         }
      }
   }
}
*/

void UpdateSampleID(ControlState& currentState, std::vector < SPLIT_READ > & GoodIndels, Indel4output & OneIndelEvent, std::vector <unsigned> & SampleIDs) {
    std::set<std::string> SampleNames;
    for (unsigned index = OneIndelEvent.Start; index <= OneIndelEvent.End; index++) {
        if (SampleNames.find(GoodIndels[index].Tag) == SampleNames.end()) { // not in the set
            SampleNames.insert(GoodIndels[index].Tag);
            //std::cout << GoodIndels[index].Tag << std::endl;
        }
    }
    //bams_to_parse Tag
    for (unsigned index = 0; index < currentState.bams_to_parse.size(); index++) {
        if (SampleNames.find(currentState.bams_to_parse[index].Tag) != SampleNames.end()) {
            SampleIDs.push_back(index);
        }
    }
}

void OutputTDs (std::vector < SPLIT_READ > &TDs,
           const std::string & TheInput,
           const unsigned int &C_S,
           const unsigned int &C_E,
           const unsigned int &RealStart,
           const unsigned int &RealEnd, std::ofstream & TDOutf)
{
    
   unsigned int NumberOfReads = 0;// = C_E - C_S + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;

   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];
   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( TDs, C_S, C_E, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( NumSupportPerTag, LeftS, LeftUNum, RightS, RightUNum );

   short NumberSupSamples = 0;
   short NumU_SupSamples = 0;
   int Num_U_Reads = 0;
   for (unsigned short i = 0; i < g_sampleNames.size (); ++i) {
      if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) {
         ++NumberSupSamples;
      }
      if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) {
         ++NumU_SupSamples;
      }
	NumberOfReads += NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus;
      Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = (LeftS + 1) * (RightS + 1 );

   CurrentChrMask[TDs[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[TDs[C_S].BPRight + g_SpacerBeforeAfter] = 'B';

   // reports BreakDancer event if applicable (-Q-option set and both ends of the SV referring to the same event)
   reportBreakDancerEvent(TDs[C_S].FragName, TDs[C_S].BPLeft, TDs[C_S].BPRight, TDs[C_S].IndelSize, "TD", NumberOfTDInstances);
   TDOutf <<
          "####################################################################################################"
          << std::endl;
   TDOutf << NumberOfTDInstances << "\tTD " << TDs[C_S].IndelSize	// << " bases "
          << "\tNT " << TDs[C_S].NT_size << " \"" << TDs[C_S].NT_str << "\"" << "\tChrID " << TDs[C_S].FragName << "\tBP " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2 << "\tBP_range " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS << "\t" << LeftUNum << "\t- " << RightS << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111 << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += TDs[i].MS;
   }
   TDOutf << "\tSUM_MS " << SUM_MS;

   TDOutf << "\t" << g_sampleNames.size() << "\tNumSupSamples " << NumberSupSamples << "\t" <<
          NumU_SupSamples;
   std::vector <int> CoverageStart, CoverageEnd; //if (!(TDs[C_S].BPLeft + 2 >= g_RegionStart &&  TDs[C_S].BPRight + 2 < g_RegionEnd)) return;
   if (TDs[C_S].BPLeft + 1 >= g_RegionStart && TDs[C_S].BPLeft + 1 < g_RegionEnd) {
       for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
           CoverageStart.push_back(g_RefCoverageRegion[TDs[C_S].BPLeft + 1 - g_RegionStart].RefCoveragePerSample[i]);
       }
   }
   else {
       for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
           CoverageStart.push_back(-1);
       }
   }
    if (TDs[C_S].BPRight + 1 > g_RegionStart && TDs[C_S].BPRight + 1 < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(g_RefCoverageRegion[TDs[C_S].BPRight + 1 - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(-1);
        }
    }
   for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
      TDOutf << "\t" << indexToSampleMap[i]
             << " " << CoverageStart[i] << " " << CoverageEnd[i]
             << " " << NumSupportPerTag[i].NumPlus
             << " " << NumSupportPerTag[i].NumUPlus
             << " " << NumSupportPerTag[i].NumMinus
             << " " << NumSupportPerTag[i].NumUMinus;
   }
   TDOutf << std::endl;

   TDOutf << TheInput.substr (TDs[C_S].BPRight + g_SpacerBeforeAfter - g_reportLength + 1, g_reportLength);	// << endl;//
   if (TDs[C_S].NT_size) {
      for (short i = 0; i < TDs[C_S].NT_size; i++) {
         TDOutf << " ";
      }
   }
   TDOutf << Cap2Low (TheInput.
                      substr (TDs[C_S].BPLeft + g_SpacerBeforeAfter,
                              g_reportLength)) << std::endl;

   short SpaceBeforeReadSeq;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      SpaceBeforeReadSeq = g_reportLength - TDs[GoodIndex].BP - 1;

      for (int i = 0; i < SpaceBeforeReadSeq; i++) {
         TDOutf << " ";
      }
      if (TDs[GoodIndex].MatchedD == Minus) {
         TDOutf << TDs[GoodIndex].getUnmatchedSeq() << std::endl;
      }
      else {
         TDOutf << TDs[GoodIndex].getUnmatchedSeqRev() << std::endl;
      }
      TDOutf << "\t" << TDs[GoodIndex].MatchedD << "\t"
             << TDs[GoodIndex].MatchedRelPos
             << "\t" << TDs[GoodIndex].MS
             << "\t" << TDs[GoodIndex].Tag << "\t" << TDs[GoodIndex].Name << std::endl;
   }
}

void OutputDeletions (std::vector < SPLIT_READ > &Deletions,
                 const std::string & TheInput,
                 const unsigned int &C_S,
                 const unsigned int &C_E,
                 const unsigned int &RealStart,
                 const unsigned int &RealEnd, std::ofstream & DeletionOutf)
{  //if (!(Deletions[C_S].BPLeft + 2 >= g_RegionStart &&  Deletions[C_S].BPRight + 2 < g_RegionEnd)) return;
   //std::cout << "Start " << g_RegionStart << "\tEnd " << g_RegionEnd << std::endl;
   LOG_DEBUG(*logStream << "d_1" << std::endl);
   unsigned int NumberOfReads = 0;//C_E - C_S + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;
   LOG_DEBUG(*logStream << "d_2" << std::endl);
   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];

//DeletionOutf << "EWL BREAKING JENKINS! (see if it still works)\n";

   for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
      NumSupportPerTag[i].NumPlus = 0;
      NumSupportPerTag[i].NumMinus = 0;
      NumSupportPerTag[i].NumUPlus = 0;
      NumSupportPerTag[i].NumUMinus = 0;
   }
   LOG_DEBUG(*logStream << "d_3" << std::endl);
   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( Deletions, C_S, C_E, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( NumSupportPerTag, LeftS, LeftUNum, RightS, RightUNum );


   LOG_DEBUG(*logStream << "d_4" << std::endl);
   short NumberSupSamples = 0;
   short NumU_SupSamples = 0;
   int Num_U_Reads = 0;
   for (unsigned short i = 0; i < g_sampleNames.size (); ++i) {
      if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) {
         ++NumberSupSamples;
      }
      if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) {
         ++NumU_SupSamples;
      }
		NumberOfReads += NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus;
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }
   LOG_DEBUG(*logStream << "d_5" << std::endl);
   unsigned int EasyScore = (LeftS + 1) * (RightS + 1 );
   short GapSize = 0;
   if (Deletions[C_S].IndelSize < 14) {
      GapSize = Deletions[C_S].IndelSize;
   }
   else {
      GapSize = 13 + (int) log10 (Deletions[C_S].IndelSize - 10);
   }
   LOG_DEBUG(*logStream << "d_5a" << std::endl);
   CurrentChrMask[Deletions[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   LOG_DEBUG(*logStream << "d_5b" << std::endl);
   CurrentChrMask[Deletions[C_S].BPRight + g_SpacerBeforeAfter] = 'B';
   LOG_DEBUG(*logStream << "d_5c" << std::endl);
   CurrentChrMask[RealStart + g_SpacerBeforeAfter] = 'B';
   LOG_DEBUG(*logStream << "d_5d" << std::endl);
   CurrentChrMask[RealEnd + g_SpacerBeforeAfter] = 'B';
   LOG_DEBUG(*logStream << "d_5e" << std::endl);
   reportBreakDancerEvent(Deletions[C_S].FragName, Deletions[C_S].BPLeft+1, Deletions[C_S].BPRight+1, Deletions[C_S].IndelSize, "D", deletionFileData.getSvIndex());
   DeletionOutf <<
                "####################################################################################################"
                << std::endl;
   DeletionOutf << deletionFileData.getSvIndex() << "\tD " << Deletions[C_S].IndelSize	// << " bases "
                << "\tNT " << Deletions[C_S].NT_size << " \"" << Deletions[C_S].NT_str << "\"" << "\tChrID " << Deletions[C_S].FragName << "\tBP " << Deletions[C_S].BPLeft + 1 << "\t" << Deletions[C_S].BPRight + 1 << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS << "\t" << LeftUNum << "\t- " << RightS << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore;
   LOG_DEBUG(*logStream << "d_6" << std::endl);
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += Deletions[i].MS;
   }
   DeletionOutf << "\tSUM_MS " << SUM_MS;

   DeletionOutf << "\t" << g_sampleNames.
                size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
                NumU_SupSamples;
    std::vector <int> CoverageStart, CoverageEnd; //if (!(TDs[C_S].BPLeft + 2 >= g_RegionStart &&  TDs[C_S].BPRight + 2 < g_RegionEnd)) return;
    if (Deletions[C_S].BPLeft + 2 >= g_RegionStart && Deletions[C_S].BPLeft + 2 < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(g_RefCoverageRegion[Deletions[C_S].BPLeft + 2 - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(-1);
        }
    }
    if (Deletions[C_S].BPRight > g_RegionStart && Deletions[C_S].BPRight < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(g_RefCoverageRegion[Deletions[C_S].BPRight - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(-1);
        }
    }
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      DeletionOutf << "\t" << indexToSampleMap[i]
                   << " " << CoverageStart[i] << " " << CoverageEnd[i]
                   << " " << NumSupportPerTag[i].NumPlus
                   << " " << NumSupportPerTag[i].NumUPlus
                   << " " << NumSupportPerTag[i].NumMinus
                   << " " << NumSupportPerTag[i].NumUMinus;
   DeletionOutf << std::endl;
   LOG_DEBUG(*logStream << "d_7" << std::endl);

   DeletionOutf << TheInput.substr (Deletions[C_S].Left - g_reportLength + Deletions[C_S].BP + 1, g_reportLength);	// << endl;// g_reportLength
   if (Deletions[C_S].IndelSize >= 14) {
      DeletionOutf << Cap2Low (TheInput.substr (Deletions[C_S].Left + Deletions[C_S].BP + 1, 5)) << "<" << Deletions[C_S].IndelSize - 10 << ">" << 
										 Cap2Low (TheInput.substr (Deletions[C_S].Right - Deletions[C_S].getReadLength() + Deletions[C_S].BP - 3, 5));
   }
   else {
      DeletionOutf << Cap2Low (TheInput.substr (Deletions[C_S].Left + Deletions[C_S].BP + 1, GapSize));
	}
   DeletionOutf << TheInput.substr (Deletions[C_S].Left + Deletions[C_S].BP + 1 + Deletions[C_S].IndelSize, g_reportLength - GapSize) << std::endl;	// g_reportLength
   short SpaceBeforeReadSeq;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      SpaceBeforeReadSeq = g_reportLength - Deletions[GoodIndex].BP - 1;

      for (int i = 0; i < SpaceBeforeReadSeq; i++) {
         DeletionOutf << " ";
      }
      short SpaceBeforeD =
         g_reportLength + g_reportLength - SpaceBeforeReadSeq -
         Deletions[GoodIndex].getReadLength();
      if (Deletions[GoodIndex].MatchedD == Minus) {
         DeletionOutf << Deletions[GoodIndex].getUnmatchedSeq().substr (0, Deletions[GoodIndex].BP + 1);	// << endl;
         for (int i = 0; i < GapSize; i++) {
            DeletionOutf << " ";
         }
         DeletionOutf << Deletions[GoodIndex].getUnmatchedSeq().substr (Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].getReadLength() - Deletions[GoodIndex].BP);	// << endl;
      }
      else {
         DeletionOutf << Deletions[GoodIndex].getUnmatchedSeqRev().substr (0, Deletions[GoodIndex].BP + 1);	// << endl;
         for (int i = 0; i < GapSize; i++) {
            DeletionOutf << " ";
         }
         DeletionOutf << Deletions[GoodIndex].getUnmatchedSeqRev().substr (Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].getReadLength() - Deletions[GoodIndex].BP);	// << endl;
      }
      for (int i = 0; i < SpaceBeforeD; i++) {
         DeletionOutf << " ";
      }
      DeletionOutf << "\t" << Deletions[GoodIndex].MatchedD << "\t"
                   << Deletions[GoodIndex].MatchedRelPos
                   << "\t" << Deletions[GoodIndex].MS
                   << "\t" << Deletions[GoodIndex].Tag
                   << "\t" << Deletions[GoodIndex].Name << std::endl;
   }
}

void OutputInversions (std::vector < SPLIT_READ > &Inv,
                  const std::string & TheInput,
                  const unsigned int &C_S,
                  const unsigned int &C_E,
                  const unsigned int &RealStart,
                  const unsigned int &RealEnd, std::ofstream & InvOutf)
{
   //if (!(Inv[C_S].BPLeft + 2 >= g_RegionStart &&  Inv[C_S].BPRight + 2 < g_RegionEnd)) return;
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
   std::string LeftNT_str = "";
   std::string RightNT_str = "";
   if (LeftNT_index != -1) {
      LeftNT_size = Inv[LeftNT_index].NT_size;
      LeftNT_str = Inv[LeftNT_index].NT_str;
   }
   if (RightNT_index != -1) {
      RightNT_size = Inv[RightNT_index].NT_size;
      RightNT_str = Inv[RightNT_index].NT_str;
   }
   unsigned int NumberOfReads = 0;//C_E - C_S + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;

   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];

   for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
      NumSupportPerTag[i].NumPlus = 0;
      NumSupportPerTag[i].NumMinus = 0;
      NumSupportPerTag[i].NumUPlus = 0;
      NumSupportPerTag[i].NumUMinus = 0;
   }

   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( Inv, C_S, C_E, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( NumSupportPerTag, LeftS, LeftUNum, RightS, RightUNum );

   short NumberSupSamples = 0;
   short NumU_SupSamples = 0;
   int Num_U_Reads = 0;
   for (unsigned short i = 0; i < g_sampleNames.size (); ++i) {
      if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) {
         ++NumberSupSamples;
      }
      if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) {
         ++NumU_SupSamples;
      }
		NumberOfReads += NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus;
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }
   //*logStream << "+" << std::endl;
   unsigned int EasyScore = (LeftS + 1) * (RightS + 1 );
   CurrentChrMask[Inv[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[Inv[C_S].BPRight + g_SpacerBeforeAfter] = 'B';
   reportBreakDancerEvent(Inv[C_S].FragName, Inv[C_S].BPLeft, Inv[C_S].BPRight+2, Inv[C_S].IndelSize, "INV", g_numberOfInvInstances);
   InvOutf <<
           "####################################################################################################"
           << std::endl;
   InvOutf << g_numberOfInvInstances++ << "\tINV " << Inv[C_S].IndelSize	// << " bases "
           << "\tNT " << LeftNT_size << ":" << RightNT_size << " \"" << LeftNT_str << "\":\"" << RightNT_str << "\"" << "\tChrID " << Inv[C_S].FragName << "\tBP " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 << "\tBP_range " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS << "\t" << LeftUNum << "\t- " << RightS << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += Inv[i].MS;
   }
   InvOutf << "\tSUM_MS " << SUM_MS;
   InvOutf << "\t" << g_sampleNames.
           size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
           NumU_SupSamples;
   
    std::vector <int> CoverageStart, CoverageEnd; //if (!(TDs[C_S].BPLeft + 2 >= g_RegionStart &&  TDs[C_S].BPRight + 2 < g_RegionEnd)) return;
    if (Inv[C_S].BPLeft + 1 >= g_RegionStart && Inv[C_S].BPLeft + 1 < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(g_RefCoverageRegion[Inv[C_S].BPLeft + 1 - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(-1);
        }
    }
    if (Inv[C_S].BPRight + 1 > g_RegionStart && Inv[C_S].BPRight + 1 < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(g_RefCoverageRegion[Inv[C_S].BPRight + 1 - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(-1);
        }
    }
    
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      InvOutf << "\t" << indexToSampleMap[i]
              << " " << CoverageStart[i] << " " << CoverageEnd[i]
              << " " << NumSupportPerTag[i].NumPlus
              << " " << NumSupportPerTag[i].NumUPlus
              << " " << NumSupportPerTag[i].NumMinus
              << " " << NumSupportPerTag[i].NumUMinus;
   InvOutf << std::endl;

   short SpaceBeforeReadSeq;
   InvOutf << TheInput.substr (Inv[C_S].BPLeft + g_SpacerBeforeAfter - g_reportLength, g_reportLength);	//;// << endl;// g_reportLength
//   LOG_DEBUG(*logStream << Inv[C_S].NT_size << "\t" << Inv[C_S].NT_2size << std::endl);
   if (LeftNT_size) {
      for (int i = 0; i < LeftNT_size; i++) {
         InvOutf << " ";
      }
   }
   InvOutf <<
           Cap2Low (ReverseComplement
                    (TheInput.
                     substr (Inv[C_S].BPRight + 1 + g_SpacerBeforeAfter - g_reportLength, g_reportLength))) << std::endl;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      if (Inv[GoodIndex].MatchedD == Plus) {
         SpaceBeforeReadSeq = g_reportLength - Inv[GoodIndex].BP - 1;
         InvOutf << std::string( SpaceBeforeReadSeq, ' ' );
         if (Inv[GoodIndex].UP_Close[0].AbsLoc < Inv[GoodIndex].UP_Far[0].AbsLoc ) {

            InvOutf << Inv[GoodIndex].getUnmatchedSeqRev();
             
            InvOutf << std::string( Inv[GoodIndex].BP, ' ' );
         }
         else {
            InvOutf << Inv[GoodIndex].getUnmatchedSeq();
         }
         InvOutf	<< "\t" << Inv[GoodIndex].MatchedD << "\t"
                  << Inv[GoodIndex].MatchedRelPos
                  << "\t" << Inv[GoodIndex].MS
                  << "\t" << Inv[GoodIndex].Tag
                  << "\t" << Inv[GoodIndex].Name << std::endl;
      }
   }
   InvOutf <<
           "----------------------------------------------------------------------------------------------------"
           << std::endl;


   InvOutf <<
           Cap2Low (ReverseComplement
                    (TheInput.
                     substr (Inv[C_S].BPLeft + g_SpacerBeforeAfter, g_reportLength)));
   if (RightNT_size) {
      for (int i = 0; i < RightNT_size; i++) {
         InvOutf << " ";
      }
   }
   InvOutf << TheInput.substr (Inv[C_S].BPRight + 1 + g_SpacerBeforeAfter, g_reportLength) << std::endl;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      if (Inv[GoodIndex].MatchedD == Minus ) {
         SpaceBeforeReadSeq = g_reportLength - Inv[GoodIndex].BP - 1;
         InvOutf << std::string( SpaceBeforeReadSeq, ' ' );

         if (Inv[GoodIndex].UP_Close[0].AbsLoc > Inv[GoodIndex].UP_Far[0].AbsLoc) {
            InvOutf << Inv[GoodIndex].getUnmatchedSeq();
            InvOutf << std::string( Inv[GoodIndex].BP, ' ');
         }
         else {
            InvOutf << Inv[GoodIndex].getUnmatchedSeqRev();
         }
         InvOutf  << "\t" << Inv[GoodIndex].MatchedD << "\t"
                  << Inv[GoodIndex].MatchedRelPos
                  << "\t" << Inv[GoodIndex].MS
                  << "\t" << Inv[GoodIndex].Tag
                  << "\t" << Inv[GoodIndex].Name << std::endl;
      }
   }
}

void OutputSIs (std::vector < SPLIT_READ > &SIs,
           const std::string & TheInput,
           const unsigned int &C_S,
           const unsigned int &C_E,
           const unsigned int &RealStart,
           const unsigned int &RealEnd, std::ofstream & SIsOutf)
{
   //if (!(SIs[C_S].BPLeft + 2 >= g_RegionStart &&  SIs[C_S].BPRight + 2 < g_RegionEnd)) return;
   unsigned int NumberOfReads = 0;//C_E - C_S + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;

   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];

   for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
      NumSupportPerTag[i].NumPlus = 0;
      NumSupportPerTag[i].NumMinus = 0;
      NumSupportPerTag[i].NumUPlus = 0;
      NumSupportPerTag[i].NumUMinus = 0;
   }
   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( SIs, C_S, C_E, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( NumSupportPerTag, LeftS, LeftUNum, RightS, RightUNum );

   short NumberSupSamples = 0;
   short NumU_SupSamples = 0;
   int Num_U_Reads = 0;
   for (unsigned short i = 0; i < g_sampleNames.size (); ++i) {
      if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) {
         ++NumberSupSamples;
      }
      if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) {
         ++NumU_SupSamples;
      }
	NumberOfReads += NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus;
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = (LeftS + 1) * (RightS + 1 );
   std::string CurrentReadSeq;

   CurrentChrMask[SIs[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[SIs[C_S].BPRight + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[RealStart + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[RealEnd + g_SpacerBeforeAfter] = 'B';

   reportBreakDancerEvent(SIs[C_S].FragName, SIs[C_S].BPLeft+1, SIs[C_S].BPRight+1, SIs[C_S].IndelSize, "SI", NumberOfSIsInstances);

   SIsOutf <<
           "####################################################################################################"
           << std::endl;
   SIsOutf << NumberOfSIsInstances << "\tI " << SIs[C_S].IndelSize << "\tNT " << SIs[C_S].IndelSize << " \"" << GetConsensusInsertedStr(SIs, C_S, C_E) << "\"" << "\tChrID " << SIs[C_S].FragName << "\tBP " << SIs[C_S].BPLeft + 1 << "\t" << SIs[C_S].BPRight + 1 << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS << "\t" << LeftUNum << "\t- " << RightS << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += SIs[i].MS;
   }
   SIsOutf << "\tSUM_MS " << SUM_MS;

   SIsOutf << "\t" << g_sampleNames.
           size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
           NumU_SupSamples;
    
    std::vector <int> CoverageStart, CoverageEnd; //if (!(TDs[C_S].BPLeft + 2 >= g_RegionStart &&  TDs[C_S].BPRight + 2 < g_RegionEnd)) return;
    if (SIs[C_S].BPLeft + 2 >= g_RegionStart && SIs[C_S].BPLeft + 2 < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(g_RefCoverageRegion[SIs[C_S].BPLeft + 2 - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(-1);
        }
    }
    if (SIs[C_S].BPRight > g_RegionStart && SIs[C_S].BPRight < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(g_RefCoverageRegion[SIs[C_S].BPRight - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(-1);
        }
    }
    
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      SIsOutf << "\t" << indexToSampleMap[i]
              << " " << CoverageStart[i] << " " << CoverageEnd[i]
              << " " << NumSupportPerTag[i].NumPlus
              << " " << NumSupportPerTag[i].NumUPlus
              << " " << NumSupportPerTag[i].NumMinus
              << " " << NumSupportPerTag[i].NumUMinus;
   SIsOutf << std::endl;

   SIsOutf << TheInput.substr (SIs[C_S].Left - g_reportLength + SIs[C_S].BP + 1, g_reportLength);	// g_reportLength
   for (unsigned int i = 0; i < SIs[C_S].IndelSize; i++) {
      SIsOutf << " ";
   }
   SIsOutf << TheInput.substr (SIs[C_S].Left + SIs[C_S].BP + 1, g_reportLength) << std::endl;	// g_reportLength

   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      short SpaceBeforeReadSeq = g_reportLength - SIs[GoodIndex].BP - 1;
      for (short i = 0; i < SpaceBeforeReadSeq; i++) {
         SIsOutf << " ";
      }
      if (SIs[GoodIndex].MatchedD == Minus) {
         SIsOutf << SIs[GoodIndex].getUnmatchedSeq();
      }
      else {
         SIsOutf << SIs[GoodIndex].getUnmatchedSeqRev();
      }
      short SpaceBeforeD =
         g_reportLength + g_reportLength - SpaceBeforeReadSeq -
         SIs[GoodIndex].getReadLength();
      for (short i = 0; i < SpaceBeforeD; i++) {
         SIsOutf << " ";
      }
      SIsOutf << "\t" << SIs[GoodIndex].MatchedD
              << "\t" << SIs[GoodIndex].MatchedRelPos
              << "\t" << SIs[GoodIndex].MS
              << "\t" << SIs[GoodIndex].Tag << "\t" << SIs[GoodIndex].Name << std::endl;
   }
   NumberOfSIsInstances++;
}

void OutputDI (std::vector < SPLIT_READ > &DI,
          const std::string & TheInput,
          const unsigned int &C_S,
          const unsigned int &C_E,
          const unsigned int &RealStart,
          const unsigned int &RealEnd, std::ofstream & DeletionOutf)
{
   //if (!(DI[C_S].BPLeft + 2 >= g_RegionStart &&  DI[C_S].BPRight + 2 < g_RegionEnd)) return;
   unsigned int NumberOfReads = 0;//C_E - C_S + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;

   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];

   for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
      NumSupportPerTag[i].NumPlus = 0;
      NumSupportPerTag[i].NumMinus = 0;
      NumSupportPerTag[i].NumUPlus = 0;
      NumSupportPerTag[i].NumUMinus = 0;
   }

   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( DI, C_S, C_E, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( NumSupportPerTag, LeftS, LeftUNum, RightS, RightUNum );

   short NumberSupSamples = 0;
   short NumU_SupSamples = 0;
   int Num_U_Reads = 0;
   for (unsigned short i = 0; i < g_sampleNames.size (); ++i) {
      if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) {
         ++NumberSupSamples;
      }
      if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) {
         ++NumU_SupSamples;
      }
	NumberOfReads += NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus;
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = (LeftS + 1) * (RightS + 1 );
   CurrentChrMask[DI[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[DI[C_S].BPRight + g_SpacerBeforeAfter] = 'B';
   reportBreakDancerEvent(DI[C_S].FragName, DI[C_S].BPLeft+1, DI[C_S].BPRight+1, DI[C_S].IndelSize, "D", deletionFileData.getSvIndex());
   DeletionOutf <<
                "####################################################################################################"
                << std::endl;
   DeletionOutf << deletionFileData.getSvIndex() << "\tD " << DI[C_S].IndelSize << "\tNT " << DI[C_S].NT_size << " \"" << DI[C_S].NT_str
<< "\"" << "\tChrID " << DI[C_S].FragName << "\tBP " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 << "\tBP_range " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS << "\t" << LeftUNum << "\t- " << RightS << "\t" << RightUNum << "\tS1 " << EasyScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += DI[i].MS;
   }
   DeletionOutf << "\tSUM_MS " << SUM_MS;


   DeletionOutf << "\t" << g_sampleNames.
                size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
                NumU_SupSamples;
   
    std::vector <int> CoverageStart, CoverageEnd; //if (!(TDs[C_S].BPLeft + 2 >= g_RegionStart &&  TDs[C_S].BPRight + 2 < g_RegionEnd)) return;
    if (DI[C_S].BPLeft + 2 >= g_RegionStart && DI[C_S].BPLeft + 2 < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(g_RefCoverageRegion[DI[C_S].BPLeft + 2 - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(-1);
        }
    }
    if (DI[C_S].BPRight > g_RegionStart && DI[C_S].BPRight < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(g_RefCoverageRegion[DI[C_S].BPRight - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(-1);
        }
    }
    
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      DeletionOutf << "\t" << indexToSampleMap[i]
                   << " " << CoverageStart[i] << " " << CoverageEnd[i]
                   << " " << NumSupportPerTag[i].NumPlus
                   << " " << NumSupportPerTag[i].NumUPlus
                   << " " << NumSupportPerTag[i].NumMinus
                   << " " << NumSupportPerTag[i].NumUMinus;
   DeletionOutf << std::endl;
   DeletionOutf << TheInput.substr (DI[C_S].Left - g_reportLength + DI[C_S].BP + 1, g_reportLength);	// << endl;// g_reportLength

   for (short i = 0; i < DI[C_S].NT_size; i++) {
      DeletionOutf << " ";
   }
   DeletionOutf << TheInput.substr (DI[C_S].Left + DI[C_S].BP + 1 + DI[C_S].IndelSize, g_reportLength) << std::endl;	// g_reportLength
   short SpaceBeforeReadSeq;
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      SpaceBeforeReadSeq = g_reportLength - DI[GoodIndex].BP - 1;

      for (int i = 0; i < SpaceBeforeReadSeq; i++) {
         DeletionOutf << " ";
      }
      if (DI[GoodIndex].MatchedD == Minus) {
         DeletionOutf << DI[GoodIndex].getUnmatchedSeq() << "\t";
      }
      else {
         DeletionOutf << DI[GoodIndex].getUnmatchedSeqRev() << "\t";
      }
      DeletionOutf << "\t" << DI[GoodIndex].MatchedD << "\t" << DI[GoodIndex].
                   MatchedRelPos << "\t" << DI[GoodIndex].MS << "\t" << DI[GoodIndex].
                   Tag << "\t" << DI[GoodIndex].Name << std::endl;
   }
}

bool CompareFragName(const std::string & a, const std::string & b) {
    if (a.size() > b.size()) {
        for (unsigned pos = 0; pos < b.size(); pos++) {
            if ((short)a[pos] > (short)b[pos]) return true;
            else if ((short)a[pos] < (short)b[pos]) return false;
        }
    }
    else if (a.size() < b.size()) {
        for (unsigned pos = 0; pos < a.size(); pos++) {
            for (unsigned pos = 0; pos < b.size(); pos++) {
                if ((short)a[pos] > (short)b[pos]) return true;
                else if ((short)a[pos] < (short)b[pos]) return false;
            }
        }
    }
    else {
        for (unsigned pos = 0; pos < b.size(); pos++) {
            for (unsigned pos = 0; pos < b.size(); pos++) {
                if ((short)a[pos] > (short)b[pos]) return true;
                else if ((short)a[pos] < (short)b[pos]) return false;
            }
        }
    }
    return false;
}

bool smaller( const SPLIT_READ& firstRead, const SPLIT_READ& secondRead )
{
	if (firstRead.FragName != secondRead.FragName ) {
		return (CompareFragName(firstRead.FragName, secondRead.FragName));
	}
	if (firstRead.BPLeft != secondRead.BPLeft ) {
		return (firstRead.BPLeft < secondRead.BPLeft);
	}
	if (firstRead.BPRight != secondRead.BPRight ) {
		return (firstRead.BPRight < secondRead.BPRight);
	}
    if (firstRead.IndelSize != secondRead.IndelSize ) {
		return (firstRead.IndelSize < secondRead.IndelSize);
	}
    if (firstRead.NT_size != secondRead.NT_size ) {
		return (firstRead.NT_size < secondRead.NT_size);
	}
	if (firstRead.BP != secondRead.BP ) {
		return (firstRead.BP < secondRead.BP);
	}
	return false; // are equal, though we may want to do more stringent sorting
}


void bubblesortReads(const std::vector < SPLIT_READ >& Reads, std::vector < unsigned >& VariantIndices)
{
	unsigned int SIsNum = VariantIndices.size();
	for (unsigned int First = 0; First < SIsNum - 1; First++) {
		for (unsigned int Second = First + 1; Second < SIsNum; Second++) {
			if (!smaller(Reads[VariantIndices[First]], Reads[VariantIndices[Second]])) {
				std::swap( VariantIndices[First],VariantIndices[Second]); 
			}
		}
	}
}

// NOTES:						// since both lengths are the same, the second argument of OR is equal to the first)
					// firstRead.LeftMostPos + firstRead.getReadLength() == secondRead.LeftMostPos + secondRead.getReadLength()) {
void markDuplicates(std::vector < SPLIT_READ >& Reads, const std::vector < unsigned >& VariantIndices)
{
	// mark reads that are not unique 
	unsigned int Num = VariantIndices.size();

   for (unsigned int First = 0; First < Num - 1; First++) {
		SPLIT_READ& firstRead = Reads[VariantIndices[First]];
		if ( firstRead.UniqueRead == false ) { continue; }
      for (unsigned int Second = First + 1; Second < Num; Second++) {
			SPLIT_READ& secondRead = Reads[VariantIndices[Second]];		

         if ( firstRead.LeftMostPos == secondRead.LeftMostPos 
					&& firstRead.Tag == secondRead.Tag
					&& firstRead.MatchedD == secondRead.MatchedD
					&& firstRead.getReadLength() == secondRead.getReadLength() 
				 ) { // added EW210812

             secondRead.UniqueRead = false;
         }
      }
   }
}


void SortOutputSI (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
              std::vector < SPLIT_READ > &Reads, std::vector < unsigned >SIs[],
              std::ofstream & SIsOutf)
{
   LOG_INFO(*logStream << "Sorting and outputing short insertions ..." << std::endl);
   unsigned int SIsNum;
   std::vector < SPLIT_READ > GoodIndels;
   unsigned int GoodNum;
   std::vector < Indel4output > IndelEvents;
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

	// loop over all boxes
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
       //std::cout << "Box_index " << Box_index << std::endl;
		// does this box contain at least enough reads to warrant an event?
      if (SIs[Box_index].size () >= userSettings->NumRead2ReportCutOff) {
         //std::cout << "Bi 1" << std::endl;
         SIsNum = SIs[Box_index].size();
         LOG_DEBUG(*logStream << "SIsNum " << SIsNum << std::endl);
			//std::cout << "Number of reads in current box: " << SIsNum << "\n";
         //std::cout << "Bi 2" << std::endl;
         bubblesortReads( Reads, SIs[ Box_index ] );
         //std::cout << "Bi 3" << std::endl;
         markDuplicates( Reads, SIs[ Box_index ] );
         //std::cout << "Bi 4" << std::endl;
         GoodIndels.clear ();
         IndelEvents.clear ();
         LOG_DEBUG(*logStream << GoodIndels.size() << std::endl);

			// push all the reads of this box in the GoodIndels vector.
         for (unsigned int index = 0; index < SIsNum; index++) {
            GoodIndels.push_back (Reads[SIs[Box_index][index]]);
         }
         //std::cout << "Bi 5" << std::endl;
         GoodNum = GoodIndels.size();
         LOG_DEBUG(*logStream << Box_index << " " << GoodNum << std::endl);
         if (GoodNum == 0) {
            continue;
         }
         Indel4output OneIndelEvent;
         //std::cout << "Bi 6" << std::endl;
         OneIndelEvent.initialize( 0, GoodIndels[0] );
			//std::cout << "Number of Good indels " << GoodNum << "\n";
			//for (unsigned x=0; x<GoodNum; x++ ) { std::cout << GoodIndels[x].getUnmatchedSeq() << ":" << GoodIndels[x].BPLeft << "-" << GoodIndels[x].IndelSize << std::endl; }

			// here we're fusing different reads into events. Typical A while() [B A ] B loop; there was a trick for that...
         //std::cout << "Bi 7" << std::endl;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].IndelSize == OneIndelEvent.IndelSize) {

               OneIndelEvent.End = GoodIndex;
            }
            else {
               OneIndelEvent.complete();
               GetRealStart4Insertion (CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
               IndelEvents.push_back (OneIndelEvent);
               OneIndelEvent.initialize( GoodIndex, GoodIndels[GoodIndex]);        
            }
         }
         //std::cout << "Bi 8" << std::endl;
         OneIndelEvent.complete();
         //std::cout << "Bi 9" << std::endl;
         GetRealStart4Insertion (CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         //std::cout << "Bi 10" << std::endl;
         IndelEvents.push_back (OneIndelEvent);

			//std::cout << "Number of indel events' " << IndelEvents.size() << "\n";
			//for (unsigned x=0; x<IndelEvents.size(); x++ ) { std::cout << IndelEvents[x].Start << "-" << IndelEvents[x].End << std::endl; }

         if (IndelEvents.size ()) {

				// loop over all events; if two events are 'identical' (same Start, end, indelsize) the support is the biggest support available. [which is weird, why not sum?]
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size (); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
					//	const Indel4output& firstEvent = IndelEvents[EventIndex];
                  //unsigned int Max_Support = IndelEvents[EventIndex].Support;
                  //unsigned int Max_Support_Index = EventIndex;
                  /*
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size (); EventIndex_left++) {
							Indel4output& secondEvent = IndelEvents[EventIndex_left];
                     if ( ( secondEvent.WhetherReport == true ) 
								&& ( secondEvent.RealStart == firstEvent.RealStart )
								&& ( secondEvent.RealEnd == firstEvent.RealEnd )
								&& ( secondEvent.IndelSize == firstEvent.IndelSize )
								&& ( secondEvent.IndelStr == firstEvent.IndelStr ) ) {
 
                        secondEvent.WhetherReport = false;
								// EW: if the events are the same, should you not be adding their support together? And if the support is higher, why not report this one and make reporting the previous event false?
                         
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
									std::cout << "THIS SHOULD NOT HAPPEN@@@!\n";
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        }
                         
                     }
                  }
                   */
                  // report max one
                  LOG_DEBUG(*logStream << IndelEvents[EventIndex].Support << std::endl);
                  if (IndelEvents[EventIndex].Support >= userSettings->NumRead2ReportCutOff && IndelEvents[EventIndex].RealStart < IndelEvents[EventIndex].RealEnd) {
                      //std::cout << "Bi 11 " << IndelEvents[EventIndex].Start << " " << IndelEvents[EventIndex].End << " " << IndelEvents[EventIndex].RealStart << " " << IndelEvents[EventIndex].RealEnd << std::endl;
                      
                     OutputSIs (GoodIndels, CurrentChr, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End, IndelEvents[EventIndex].RealStart, IndelEvents[EventIndex].RealEnd, SIsOutf);
                      //std::cout << "Bi 12" << std::endl;

                  }
               }
            }
         }
      }												// if (!insertion[Box_index].empty())
   }
   LOG_INFO(*logStream << "Short insertions: " << NumberOfSIsInstances << std::endl << std::endl);
}

bool IsGoodTD(std::vector < SPLIT_READ > & GoodIndels, Indel4output & OneIndelEvent, unsigned RealStart, unsigned RealEnd, ControlState& currentState) {
    //std::cout << "Real Start and End: " << RealStart << " " << RealEnd
    //<< " " << RP_support_D(currentState, OneIndelEvent, RealStart, RealEnd)
    if (RealEnd < RealStart) return false;
    if (RealStart == 0) return false;
    //UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
    
    //if (userSettings->pindelConfigFileAsInput() ||userSettings->singlePindelFileAsInput()) return true;
    if (userSettings->pindelConfigFileAsInput() || userSettings->singlePindelFileAsInput() || userSettings->NormalSamples == false) return true;
    //return true;
    
    //<< std::endl;
    
    if (RealEnd - RealStart < (unsigned)(GoodIndels[0].getReadLength() * 2)) {
        //std::cout << "<1000 good" << std::endl;
        return true;
    }
    else {
        //std::cout << "entering > 2000" << std::endl;
        int chromosomeID = g_genome.getChrID(OneIndelEvent.ChrName);
        //std::cout << "chromosomeID " << chromosomeID << std::endl;
        if (chromosomeID == -1) {
            std::cout << "ID -1 " << OneIndelEvent.ChrName << " " << chromosomeID << std::endl;
            return false;
        }
        
        Genotyping OneDEL;
        OneDEL.Type = "TD";
        OneDEL.ChrA = OneIndelEvent.ChrName;
        OneDEL.ChrB = OneIndelEvent.ChrName;
        OneDEL.CI_A = 50;
        OneDEL.CI_B = 50;
        OneDEL.PosA = RealStart;
        OneDEL.PosB = RealEnd;
        std::vector <unsigned> SampleIDs;
        //std::cout << "Before SampleIDs.size() " << SampleIDs.size() << std::endl;
        UpdateSampleID(currentState, GoodIndels, OneIndelEvent, SampleIDs);
        //std::cout << "After SampleIDs.size() " << SampleIDs.size() << std::endl;
        //g_genome.getChr(chromosomeID);
        getRelativeCoverageForFiltering(chromosomeID, currentState, OneDEL, g_genome.getChr(chromosomeID), SampleIDs);
        //std::cout << ">2000, true " << OneIndelEvent.End - OneIndelEvent.Start  << std::endl;
        //std::cout << "After getRelativeCoverageForFiltering " << SampleIDs.size() << " " << OneDEL.RD_signals.size() << std::endl;
        unsigned CountGoodSamples = 0;
        for (unsigned RD_index = 0; RD_index < OneDEL.RD_signals.size(); RD_index++) {
            if (OneDEL.RD_signals[RD_index] >= 2.7) CountGoodSamples++;
            //std::cout << " " << std::fixed << OneDEL.RD_signals[RD_index] << " ";
        }
        if ((OneDEL.RD_signals.size() == 1 && CountGoodSamples == 1) || (OneDEL.RD_signals.size() > 1 && OneDEL.RD_signals.size() <= 4 && OneDEL.RD_signals.size() - CountGoodSamples <= 1) || (OneDEL.RD_signals.size() > 4 && (((float)CountGoodSamples / OneDEL.RD_signals.size()) > 0.66) ))
            return true;
        else return false;
    }
    return false;
}

void SortAndOutputTandemDuplications (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr, std::vector < SPLIT_READ > &AllReads, std::vector < unsigned >TDs[],
                                 std::ofstream & TDOutf, const bool nonTemplate)
{

   if (nonTemplate) {
      LOG_INFO(*logStream << "Sorting and outputing tandem duplications with non-template sequence ..." << std::endl);
   }
   else {
      LOG_INFO(*logStream << "Sorting and outputing tandem duplications ..." << std::endl);
   }
   unsigned int TDNum;
   //short CompareResult;
   //unsigned Temp4Exchange;
   int countTandemDuplications = 0;
   unsigned int GoodNum;
   std::vector < SPLIT_READ > GoodIndels;
   std::vector < Indel4output > IndelEvents;
	//UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (TDs[Box_index].size () >= userSettings->NumRead2ReportCutOff) {
         TDNum = TDs[Box_index].size ();

          bubblesortReads( AllReads, TDs[ Box_index ] );
          markDuplicates( AllReads, TDs[ Box_index ] );
          
         GoodIndels.clear ();
         IndelEvents.clear ();

         for (unsigned int First = 0; First < TDNum; First++) {
            GoodIndels.push_back (AllReads[TDs[Box_index][First]]);
         }

         GoodNum = GoodIndels.size ();
         if (GoodNum == 0) {
            continue;
         }
         Indel4output OneIndelEvent;
         OneIndelEvent.ChrName = GoodIndels[0].FragName;
         OneIndelEvent.Start = 0;
         OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight = GoodIndels[0].BPRight;
         OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                  && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight) {
               OneIndelEvent.End = GoodIndex;
            }
            else {
               OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
               OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
               OneIndelEvent.Support =
                  OneIndelEvent.End - OneIndelEvent.Start + 1;
               GetRealStart4Deletion (CurrentChr, OneIndelEvent.RealStart,
                                      OneIndelEvent.RealEnd);
               IndelEvents.push_back (OneIndelEvent);
               OneIndelEvent.Start = GoodIndex;
               OneIndelEvent.End = GoodIndex;
               OneIndelEvent.BPLeft = GoodIndels[GoodIndex].BPLeft;
               OneIndelEvent.BPRight = GoodIndels[GoodIndex].BPRight;
               OneIndelEvent.ChrName = GoodIndels[GoodIndex].FragName; 
            }
         }

         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion (CurrentChr, OneIndelEvent.RealStart,
                                OneIndelEvent.RealEnd);
         IndelEvents.push_back (OneIndelEvent);

         if (IndelEvents.size ()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
                  EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  unsigned int RealStart = IndelEvents[EventIndex].RealStart;
                  unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
                  //unsigned int Max_Support = IndelEvents[EventIndex].Support;
                  //unsigned int Max_Support_Index = EventIndex;
                   /*
                  for (unsigned EventIndex_left = 0;
                        EventIndex_left < IndelEvents.size ();
                        EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport ==
                           false) {
                        continue;
                     }
                     else if (IndelEvents[EventIndex_left].RealStart !=
                              RealStart) {
                        continue;
                     }
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) {
                        continue;
                     }
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        }
                     }
                  }
                   */
                  // report max one
                  if (IndelEvents[EventIndex].Support >= userSettings->NumRead2ReportCutOff && IsGoodTD(GoodIndels, IndelEvents[EventIndex], RealStart, RealEnd, currentState)) {
                     if (GoodIndels[ IndelEvents[EventIndex].Start].IndelSize < userSettings->BalanceCutoff) {
                        OutputTDs (GoodIndels, CurrentChr, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End, RealStart, RealEnd, TDOutf);
                        NumberOfTDInstances++;
                        countTandemDuplications++;
                     }
                     else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
                        OutputTDs (GoodIndels, CurrentChr, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End, RealStart, RealEnd, TDOutf);
                        NumberOfTDInstances++;
                        countTandemDuplications++;
                     }
                  }
               }
            }
         }
      } // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   if (nonTemplate) {
      LOG_INFO(*logStream << "Tandem duplications with non-template sequence (TD_NT): " <<
               countTandemDuplications << std::endl
               << std::endl);
   }
   else {
      LOG_INFO(*logStream << "Tandem duplications: " << NumberOfTDInstances << std::endl
               << std::endl);
   }
}

bool RP_support_D(ControlState& currentState, Indel4output & OneIndelEvent, unsigned start, unsigned end) {
    // RP_READ
    // Reads_RP_Discovery
    unsigned CountSupportingReadPairRead = 0;
    unsigned Cutoff = OneIndelEvent.Support * 2;
    if (Cutoff > 10) Cutoff = 10;
    for (unsigned index = 0; index < currentState.Reads_RP_Discovery.size(); index++) {// DA ～＝ DB ChrNameA == ChrNameB == CurrentChrName PosA < start PosB > end
        RP_READ & currentRead = currentState.Reads_RP_Discovery[index];
        if ((currentRead.PosA < currentRead.PosB && currentRead.PosA < start && currentRead.PosB > end) || (currentRead.PosA > currentRead.PosB && currentRead.PosB < start && currentRead.PosA > end)) {
            if (currentRead.DA != currentRead.DB) {
                if (currentRead.ChrNameA == OneIndelEvent.ChrName && currentRead.ChrNameB == OneIndelEvent.ChrName) {
                    CountSupportingReadPairRead++;
                }
            }
        }
        //std::cout << "CountSupportingReadPairRead " << CountSupportingReadPairRead << std::endl;
        if (CountSupportingReadPairRead >= Cutoff) return true;
    }
    return false;
}



bool IsGoodDeletion(std::vector < SPLIT_READ > & GoodIndels, Indel4output & OneIndelEvent, unsigned RealStart, unsigned RealEnd, ControlState& currentState) {
    //std::cout << "Real Start and End: " << RealStart << " " << RealEnd
    //<< " " << RP_support_D(currentState, OneIndelEvent, RealStart, RealEnd)
    //<< std::endl;
    
    //UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
    
    if (RealEnd < RealStart) return false;
    if (RealStart == 0) return false;
    
    if (userSettings->pindelConfigFileAsInput() || userSettings->singlePindelFileAsInput() || userSettings->NormalSamples == false) return true;
    //return true;
    //std::cout << "in IsGoodDeletion RealStart " << RealStart << " RealEnd " << RealEnd << std::endl;
    if (RealEnd < RealStart) return false;
    //std::cout << "should be returned already, you should not see RealStart " << RealStart << " RealEnd " << RealEnd << std::endl;
    if (RealEnd - RealStart < 1000) {
        //std::cout << "<1000 good" << std::endl;
        return true;
    }
    else if (RealEnd - RealStart < 2000) {
        
        if (RP_support_D(currentState, OneIndelEvent, RealStart, RealEnd)) {
            //std::cout << "<2000, true " << OneIndelEvent.End - OneIndelEvent.Start << std::endl;
            return true;   
        }
        else {
            //std::cout << "<2000, false " << OneIndelEvent.End - OneIndelEvent.Start << std::endl;
            return false;
        }
    }
    else if (RP_support_D(currentState, OneIndelEvent, RealStart, RealEnd)) {
        //std::cout << "entering > 2000" << std::endl;
        int chromosomeID = g_genome.getChrID(OneIndelEvent.ChrName);
        //std::cout << "chromosomeID " << chromosomeID << std::endl;
        if (chromosomeID == -1) {
            std::cout << "ID -1 " << OneIndelEvent.ChrName << " " << chromosomeID << std::endl;
            return false;
        }
                
        Genotyping OneDEL;
        OneDEL.Type = "DEL";
        OneDEL.ChrA = OneIndelEvent.ChrName;
        OneDEL.ChrB = OneIndelEvent.ChrName;
        OneDEL.CI_A = 50;
        OneDEL.CI_B = 50;
        OneDEL.PosA = RealStart;
        OneDEL.PosB = RealEnd;
        std::vector <unsigned> SampleIDs;
        //std::cout << "Before SampleIDs.size() " << SampleIDs.size() << std::endl;
        UpdateSampleID(currentState, GoodIndels, OneIndelEvent, SampleIDs);
        //std::cout << "After SampleIDs.size() " << SampleIDs.size() << std::endl;
        //g_genome.getChr(chromosomeID);
        getRelativeCoverageForFiltering(chromosomeID, currentState, OneDEL, g_genome.getChr(chromosomeID), SampleIDs);
        //std::cout << ">2000, true " << OneIndelEvent.End - OneIndelEvent.Start  << std::endl;
        //std::cout << "After getRelativeCoverageForFiltering " << SampleIDs.size() << " " << OneDEL.RD_signals.size() << std::endl;
        unsigned CountGoodSamples = 0;
        for (unsigned RD_index = 0; RD_index < OneDEL.RD_signals.size(); RD_index++) {
            if (OneDEL.RD_signals[RD_index] <= 1.3) CountGoodSamples++;
            //std::cout << " " << std::fixed << OneDEL.RD_signals[RD_index];
        }
        if ((OneDEL.RD_signals.size() == 1 && CountGoodSamples == 1) || (OneDEL.RD_signals.size() > 1 && OneDEL.RD_signals.size() <= 4 && OneDEL.RD_signals.size() - CountGoodSamples <= 1) || (OneDEL.RD_signals.size() > 4 && (((float)CountGoodSamples / OneDEL.RD_signals.size()) > 0.66) ))
            return true;
        else return false;
    }
    
    return false;
}

void SortOutputD (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
             std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Deletions[],
             std::ofstream & DeletionOutf)
{
   LOG_INFO(*logStream << "Sorting and outputing deletions ..." << std::endl);
   unsigned int DeletionsNum;
   //short CompareResult;
   //unsigned Temp4Exchange;
    

   unsigned int GoodNum;
   std::vector < SPLIT_READ > GoodIndels;
   std::vector < Indel4output > IndelEvents;
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
    //std::cout << "S1" << std::endl;
    //std::cout << "here" << std::endl;
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      LOG_DEBUG(*logStream << Box_index << "\t" << NumBoxes << "\t" << Deletions[Box_index].size() << std::endl);
      if (Deletions[Box_index].size () >= userSettings->NumRead2ReportCutOff) {
         //std::cout << "S2" << std::endl;
         DeletionsNum = Deletions[Box_index].size ();
          /*
          for (unsigned index = 0; index < DeletionsNum; index++) {
              std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
              std::cout << Reads[Deletions[ Box_index ][index]];
              std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
          }
          */
          //std::cout << "S3" << std::endl;
          bubblesortReads( Reads, Deletions[ Box_index ] );
          //std::cout << "S4" << std::endl;
          markDuplicates( Reads, Deletions[ Box_index ] );
          //std::cout << "S5" << std::endl;
         GoodIndels.clear ();
         IndelEvents.clear ();
         for (unsigned int First = 0; First < DeletionsNum; First++) {
            GoodIndels.push_back (Reads[Deletions[Box_index][First]]);
         }
         //std::cout << "S6" << std::endl;
         GoodNum = GoodIndels.size ();
         LOG_DEBUG(*logStream << Box_index << " box read size " << GoodNum << std::endl);
         if (GoodNum == 0) {
            continue;
         }
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
         OneIndelEvent.ChrName = GoodIndels[0].FragName;
         OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight = GoodIndels[0].BPRight;
         OneIndelEvent.WhetherReport = true;
         LOG_DEBUG(*logStream << "here" << std::endl);
         //std::cout << "S7" << std::endl;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                  && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight
                  && GoodIndels[GoodIndex].FragName == OneIndelEvent.ChrName
                  && GoodIndels[GoodIndex].FarFragName == OneIndelEvent.ChrName) {
               OneIndelEvent.End = GoodIndex;
            }
            else {
               OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
               OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
               OneIndelEvent.Support =
                  OneIndelEvent.End - OneIndelEvent.Start + 1;
               GetRealStart4Deletion (CurrentChr, OneIndelEvent.RealStart,
                                      OneIndelEvent.RealEnd);
               IndelEvents.push_back (OneIndelEvent);
            /*    std::cout << OneIndelEvent.RealStart << " "
                          << OneIndelEvent.Start << " "
                          << OneIndelEvent.End << " "
                          << OneIndelEvent.RealEnd << " "
                          << OneIndelEvent.BPLeft << " "
                          << OneIndelEvent.BPRight << " " << OneIndelEvent.Support << std::endl;*/
               OneIndelEvent.Start = GoodIndex;
               OneIndelEvent.End = GoodIndex;
               OneIndelEvent.BPLeft = GoodIndels[GoodIndex].BPLeft;
               OneIndelEvent.BPRight = GoodIndels[GoodIndex].BPRight;
               OneIndelEvent.ChrName = GoodIndels[GoodIndex].FragName;
            }
         }
         //std::cout << "S8" << std::endl;
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         // std::cout << OneIndelEvent.RealStart << " " << OneIndelEvent.RealEnd << std::endl;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion (CurrentChr, OneIndelEvent.RealStart,
                                OneIndelEvent.RealEnd);
         IndelEvents.push_back (OneIndelEvent);
         LOG_DEBUG(*logStream << "IndelEvent: " << IndelEvents.size() << std::endl);

         if (IndelEvents.size ()) {
            // std::cout << "IndelEvents.size ()" << IndelEvents.size () << std::endl;
            //std::cout << "S8a" << std::endl;
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
                  EventIndex++) {
               LOG_DEBUG(*logStream << EventIndex << " EventIndex" << std::endl);
               if (IndelEvents[EventIndex].WhetherReport) {
                  unsigned int RealStart = IndelEvents[EventIndex].RealStart;
                  unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
                  //unsigned int Max_Support = IndelEvents[EventIndex].Support;
                  //unsigned int Max_Support_Index = EventIndex;
                  /*
                  for (unsigned EventIndex_left = 0;
                        EventIndex_left < IndelEvents.size ();
                        EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport ==
                           false) {
                        continue;
                     }
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) {
                        continue;
                     }
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) {
                        continue;
                     }
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        }
                     }
                  }
                  */ 
                  // report max one
                  LOG_DEBUG(*logStream << "max" << std::endl);
                   //std::cout << "current D " << RealStart << " " << RealEnd << " " << RealEnd - RealStart << " " << IndelEvents[EventIndex].Support << " " << IsGoodDeletion(GoodIndels, IndelEvents[EventIndex], RealStart, RealEnd, currentState) << std::endl;
                  if (IndelEvents[EventIndex].Support >= userSettings->NumRead2ReportCutOff && IsGoodDeletion(GoodIndels, IndelEvents[EventIndex], RealStart, RealEnd, currentState))
                  {
                      //std::cout << "passed IsGoodDeletion" << std::endl;
                     LOG_DEBUG(*logStream << "aa" << std::endl);
                     // std::cout << "aa" << std::endl;
                     if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < userSettings->BalanceCutoff) {
                        LOG_DEBUG(*logStream << "ba" << std::endl);
                       //  std::cout << "ba" << std::endl;
                        OutputDeletions (GoodIndels, CurrentChr, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End, RealStart,
									RealEnd, DeletionOutf);
                        deletionFileData.increaseTemplateSvCounter();
                        LOG_DEBUG(*logStream << "bb" << std::endl);
                         //std::cout << "bb" << std::endl;
                     }
                     else if (ReportEvent( GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
                        LOG_DEBUG(*logStream << "ca" << std::endl);
                         //std::cout << "ca" << std::endl;
                        //std::cout << "there 1" << std::endl;
                        OutputDeletions (GoodIndels, CurrentChr, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End, RealStart, RealEnd,
                                         DeletionOutf);
                         //std::cout << "there 2" << std::endl;
                        deletionFileData.increaseTemplateSvCounter();
                        LOG_DEBUG(*logStream << "cb" << std::endl);
                         //std::cout << "cb" << std::endl;
                     }
                  }
               }
            }
            //std::cout << "S8b" << std::endl;
         }
         //std::cout << "S9" << std::endl;
      }												// if (!Deletions[Box_index].empty())
   }														// for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
    //std::cout << "there" << std::endl;
   LOG_INFO(*logStream << "Deletions: " << deletionFileData.getTemplateSvCounter() << std::endl << std::endl);
}

void SortOutputInv (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
               std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Inv[],
               std::ofstream & InvOutf)
{
   OutputSorter os(NumBoxes, CurrentChr, InvOutf);
   os.SortAndOutputInversions(currentState, Reads, Inv);
}

void SortOutputInv_NT (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
                  std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Inv[],
                  std::ofstream & InvOutf)
{
   OutputSorter os(NumBoxes, CurrentChr, InvOutf);
   os.SortAndOutputNonTemplateInversions(currentState, Reads, Inv);
}

void OutputShortInversion (std::vector < SPLIT_READ > &supportingReads,
                           const std::string &chromosome,
                           const unsigned int &indexOfFirstRead,
                           const unsigned int &indexOfLastRead,
                           const unsigned int &RealStart,
                           const unsigned int &RealEnd, std::ofstream & InversionOutF)
{
   //if (!(supportingReads[indexOfFirstRead].BPLeft + 2 >= g_RegionStart &&  supportingReads[indexOfFirstRead].BPRight + 2 < g_RegionEnd)) return;
   unsigned int NumberOfReads = 0;//indexOfLastRead - indexOfFirstRead + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;

   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];
   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( supportingReads, indexOfFirstRead, indexOfLastRead, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( NumSupportPerTag, LeftS, LeftUNum, RightS, RightUNum );

   short NumberSupSamples = 0;
   short NumU_SupSamples = 0;
   int Num_U_Reads = 0;
   for (unsigned short i = 0; i < g_sampleNames.size (); ++i) {
      if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) {
         ++NumberSupSamples;
      }
      if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) {
         ++NumU_SupSamples;
      }
		NumberOfReads += NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus;
      Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = (LeftS + 1) * (RightS + 1 );
   const SPLIT_READ& firstSupportingRead = supportingReads[indexOfFirstRead];
   CurrentChrMask[supportingReads[ indexOfFirstRead ].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[supportingReads[ indexOfFirstRead ].BPRight + g_SpacerBeforeAfter] = 'B';
   reportBreakDancerEvent( supportingReads[ indexOfFirstRead ].FragName, supportingReads[indexOfFirstRead].BPLeft+1,
                           supportingReads[indexOfFirstRead].BPRight+1, supportingReads[indexOfFirstRead].IndelSize, "INV", deletionFileData.getSvIndex());
   InversionOutF <<
                 "####################################################################################################"
                 << std::endl;
   InversionOutF << g_numberOfInvInstances++ << "\tINV " << supportingReads[indexOfFirstRead].IndelSize << "\tNT " << supportingReads[indexOfFirstRead].NT_size << " \"" <<
                 supportingReads[indexOfFirstRead].NT_str << "\"" << "\tChrID " << supportingReads[indexOfFirstRead].FragName << "\tBP " << supportingReads[indexOfFirstRead].BPLeft + 1 <<
                 "\t" << supportingReads[indexOfFirstRead].BPRight + 1 << "\tBP_range " << supportingReads[indexOfFirstRead].BPLeft + 1 << "\t" <<
                 supportingReads[indexOfFirstRead].BPRight + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS << "\t" << LeftUNum << "\t- " << RightS <<
                 "\t" << RightUNum << "\tS1 " << EasyScore;

   int SUM_MS = 0;
   for (unsigned int i = indexOfFirstRead; i <= indexOfLastRead; i++) {
      SUM_MS += supportingReads[i].MS;
   }
   InversionOutF << "\tSUM_MS " << SUM_MS;


   InversionOutF << "\t" << g_sampleNames.size() << "\tNumSupSamples " << NumberSupSamples << "\t" <<
                 NumU_SupSamples;
    
    std::vector <int> CoverageStart, CoverageEnd; //if (!(TDs[C_S].BPLeft + 2 >= g_RegionStart &&  TDs[C_S].BPRight + 2 < g_RegionEnd)) return;
    if (supportingReads[indexOfFirstRead].BPLeft + 2 >= g_RegionStart && supportingReads[indexOfFirstRead].BPLeft + 2 < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(g_RefCoverageRegion[supportingReads[indexOfFirstRead].BPLeft + 2 - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageStart.push_back(-1);
        }
    }
    if (supportingReads[indexOfFirstRead].BPRight > g_RegionStart && supportingReads[indexOfFirstRead].BPRight < g_RegionEnd) {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(g_RefCoverageRegion[supportingReads[indexOfFirstRead].BPRight - g_RegionStart].RefCoveragePerSample[i]);
        }
    }
    else {
        for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
            CoverageEnd.push_back(-1);
        }
    }
    
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      InversionOutF << "\t" << indexToSampleMap[i]
                    << " " << CoverageStart[i] << " " << CoverageEnd[i]
                    << " " << NumSupportPerTag[i].NumPlus
                    << " " << NumSupportPerTag[i].NumUPlus
                    << " " << NumSupportPerTag[i].NumMinus
                    << " " << NumSupportPerTag[i].NumUMinus;
   InversionOutF << std::endl;

   InversionOutF << chromosome.substr (supportingReads[indexOfFirstRead].Left - g_reportLength + supportingReads[indexOfFirstRead].BP + 1, g_reportLength);
   InversionOutF << Cap2Low( ReverseComplement( chromosome.substr( firstSupportingRead.Left + firstSupportingRead.BP + 1, firstSupportingRead.NT_size )));
   InversionOutF << chromosome.substr (firstSupportingRead.Left + firstSupportingRead.BP + 1 + firstSupportingRead.IndelSize, g_reportLength) << std::endl;	// ReportLength
   short SpaceBeforeReadSeq;
   for (unsigned int GoodIndex = indexOfFirstRead; GoodIndex <= indexOfLastRead; GoodIndex++) {
      SpaceBeforeReadSeq = g_reportLength - supportingReads[GoodIndex].BP - 1;

      for (int i = 0; i < SpaceBeforeReadSeq; i++) {
         InversionOutF << " ";
      }
      if (supportingReads[GoodIndex].MatchedD == Minus) {
         InversionOutF << supportingReads[GoodIndex].getUnmatchedSeq() << "\t";
      }
      else {
         InversionOutF << supportingReads[GoodIndex].getUnmatchedSeqRev() << "\t";
      }
      InversionOutF << "\t" << supportingReads[GoodIndex].MatchedD << "\t" << supportingReads[GoodIndex].MatchedRelPos << "\t" << supportingReads[GoodIndex].MS << "\t"
                    << supportingReads[GoodIndex].Tag << "\t" << supportingReads[GoodIndex].Name << std::endl;
   }
}


bool IsInversion( const SPLIT_READ& read, const std::string& chromosome )
{
   if (read.IndelSize == read.NT_size ) {
      std::string replacedString = chromosome.substr( g_SpacerBeforeAfter + 1 + read.BPLeft, read.NT_size );
      if ( ReverseComplement( replacedString ) == read.NT_str ) {
         return true;
      }
   }
   return false;
}

void SortOutputDI (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
                   std::vector < SPLIT_READ > &Reads, std::vector < unsigned >DI[],
                   std::ofstream & DIOutf, std::ofstream & InvOutf)
{
   LOG_INFO(*logStream << "Sorting and outputing deletions with non-template sequences ..." <<
            std::endl);
    std::cout << "Added: Sorting and outputing deletions with non-template sequences ..." << std::endl;
   unsigned int DINum;
   short CompareResult;
   unsigned Temp4Exchange;

   unsigned int GoodNum;
   std::vector < SPLIT_READ > GoodIndels;
   std::vector < Indel4output > IndelEvents;
	//UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      LOG_DEBUG(*logStream << "Box_index: "   << Box_index << std::endl);
      //if (DI[Box_index].size ())
       //std::cout << Box_index << " " << DI[Box_index].size () << std::endl;
      if (DI[Box_index].size () >= userSettings->NumRead2ReportCutOff) {
         DINum = DI[Box_index].size ();
         for (unsigned int First = 0; First < DINum - 1; First++) {
            {
               for (unsigned int Second = First + 1; Second < DINum; Second++) {
                  {  
                     CompareResult = 0;
                     if (Reads[DI[Box_index][First]].BPLeft < Reads[DI[Box_index][Second]].BPLeft) {
                        continue;
                     }
                     else if (Reads[DI[Box_index][First]].BPLeft > Reads[DI[Box_index][Second]].BPLeft) {
                        CompareResult = 1;
                     }
                     else if (Reads[DI[Box_index][First]].BPLeft == Reads[DI[Box_index][Second]].BPLeft) {
                        if (Reads[DI[Box_index][First]].BPRight < Reads[DI[Box_index][Second]].BPRight) {
                           continue;
                        }
                        else if (Reads[DI[Box_index][First]].BPRight > Reads[DI[Box_index][Second]].BPRight) {
                           CompareResult = 1;
                        }
                        else {
                           if (Reads[DI[Box_index][First]].NT_size < Reads[DI[Box_index][Second]].NT_size) {
                              continue;
                           }
                           else if (Reads[DI[Box_index][First]].NT_size >
                                    Reads[DI[Box_index][Second]].NT_size) {
                              CompareResult = 1;
                           }
                           else if (Reads[DI[Box_index][First]].BP > Reads[DI[Box_index][Second]].BP) CompareResult = 1; 
                        }
                     }
                     if (CompareResult == 1) {
                        Temp4Exchange = DI[Box_index][First];
                        DI[Box_index][First] = DI[Box_index][Second];
                        DI[Box_index][Second] = Temp4Exchange;
                     }
                  }
               }
            }
         }
          for (unsigned int First = 0; First < DINum - 1; First++) {
              for (unsigned int Second = First + 1; Second < DINum; Second++) {
                  if (Reads[DI[Box_index][First]].getReadLength() == Reads[DI[Box_index][Second]].getReadLength()) {
                      if (Reads[DI[Box_index][First]].LeftMostPos ==
                          Reads[DI[Box_index][Second]].LeftMostPos || Reads[DI[Box_index][First]].LeftMostPos + Reads[DI[Box_index][First]].getReadLength() ==
                          Reads[DI[Box_index][Second]].LeftMostPos + Reads[DI[Box_index][Second]].getReadLength()) {
                          if (Reads[DI[Box_index][First]].MatchedD == Reads[DI[Box_index][Second]].MatchedD)
                          Reads[DI[Box_index][Second]].UniqueRead = false;
                      }
                  }
              }
          }
         GoodIndels.clear ();
         IndelEvents.clear ();

         for (unsigned int First = 0; First < DINum; First++) {
            GoodIndels.push_back (Reads[DI[Box_index][First]]);
         }

         GoodNum = GoodIndels.size ();
         LOG_DEBUG(*logStream << Box_index << " " << GoodNum << std::endl);
         if (GoodNum == 0) {
            continue;
         }
         LOG_DEBUG(*logStream << "GoodNum: " << GoodNum << std::endl);
         Indel4output OneIndelEvent;
         OneIndelEvent.ChrName = GoodIndels[0].FragName;
         OneIndelEvent.Start = 0;
         OneIndelEvent.End = 0;
         OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.NT_size = GoodIndels[0].NT_size;
         OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight = GoodIndels[0].BPRight;
         OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                  && GoodIndels[GoodIndex].IndelSize ==
                  OneIndelEvent.IndelSize
                  && GoodIndels[GoodIndex].NT_size == OneIndelEvent.NT_size) {
               OneIndelEvent.End = GoodIndex;
            }
            else {
               IndelEvents.push_back (OneIndelEvent);
               OneIndelEvent.Start = GoodIndex;
               OneIndelEvent.End = GoodIndex;
               OneIndelEvent.BPLeft = GoodIndels[GoodIndex].BPLeft;
               OneIndelEvent.IndelSize = GoodIndels[GoodIndex].IndelSize;
               OneIndelEvent.NT_size = GoodIndels[GoodIndex].NT_size;
               OneIndelEvent.ChrName = GoodIndels[GoodIndex].FragName;
            }
         }

         IndelEvents.push_back (OneIndelEvent);
         unsigned int RealStart;
         unsigned int RealEnd;
         // std::cout << "DI IndelEvents " << IndelEvents.size() << std::endl;
         for (unsigned EventIndex = 0; EventIndex < IndelEvents.size (); EventIndex++) {
            if (IndelEvents[EventIndex].End - IndelEvents[EventIndex].Start + 1 >= userSettings->NumRead2ReportCutOff) {
               RealStart =	GoodIndels[IndelEvents[EventIndex].Start].BPLeft;
               RealEnd = GoodIndels[IndelEvents[EventIndex].Start].BPRight;
               if (  (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < userSettings->BalanceCutoff)
                     || (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) ) {

                  if (IsInversion(GoodIndels[IndelEvents[EventIndex].Start], CurrentChr )) {
                     OutputShortInversion(GoodIndels, CurrentChr,IndelEvents[EventIndex].Start,
                                          IndelEvents[EventIndex].End,RealStart, RealEnd, InvOutf);
                      //std::cout << "Small inv " << IndelEvents[EventIndex].End - IndelEvents[EventIndex].Start + 1 << std::endl;
                     // increase the number of inversion instances in the OutputShortInversion itself; it'd add
                     // needless complexity here.
                  }
                  else if (IsGoodDeletion(GoodIndels, IndelEvents[EventIndex], RealStart, RealEnd, currentState)) {
                     OutputDI (GoodIndels, CurrentChr,IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End,
                               RealStart, RealEnd, DIOutf);
                     deletionFileData.increaseNonTemplateSvCounter();
                     //std::cout << "Large del " << IndelEvents[EventIndex].End - IndelEvents[EventIndex].Start + 1 << std::endl;
                  }
               }
            }
         }
      }												// if (!Deletions[Box_index].empty())
   }														// for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   LOG_INFO(*logStream << "deletions with non-template sequences: " << deletionFileData.getNonTemplateSvCounter() <<
            std::endl << std::endl);
}



void SortOutputLI (ControlState& currentState, const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads, const SearchWindow& window, const std::string& filename)
{
   unsigned UP_Close_index;
   unsigned temp_AbsLoc;
   unsigned int LI_BORDER_BUFFER = 4 * g_maxInsertSize;
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
	time_t Time_LI_S, Time_LI_E;
   Time_LI_S = time(NULL);

	std::ofstream LargeInsertionOutf( filename.c_str(), std::ios::app);
   *logStream << "LI: maxInsertSize: " << g_maxInsertSize << std::endl;

   unsigned int absStartLIWindow = g_SpacerBeforeAfter + window.getStart();
   unsigned int absEndLIWindow = g_SpacerBeforeAfter + window.getEnd();

   if (absEndLIWindow > CurrentChr.size() - g_SpacerBeforeAfter ) {
      absEndLIWindow = CurrentChr.size() - g_SpacerBeforeAfter;
   }

   unsigned int absStartBuffered = absStartLIWindow - LI_BORDER_BUFFER;
   unsigned int absEndBuffered = absEndLIWindow + LI_BORDER_BUFFER;

   //*logStream << "SBA: " << g_SpacerBeforeAfter << " WS " << windowStart << " Total: " << g_SpacerBeforeAfter + windowStart - LI_BORDER_BUFFER << std::endl;
   ShiftedVector< uint8_t > plus_LI_Pos( absStartBuffered , absEndBuffered , 0 );
   ShiftedVector< uint8_t > minus_LI_Pos( absStartBuffered , absEndBuffered , 0 );
   ShiftedVector< int32_t > EventIndex_Pos( absStartBuffered , absEndBuffered , -1 );

   for (unsigned Index = 0; Index < Reads.size (); Index++) {
      if (Reads[Index].Used || !Reads[Index].UP_Far.empty ()) {
         continue;
      }
      temp_AbsLoc = Reads[Index].UP_Close[Reads[Index].UP_Close.size () - 1].AbsLoc;
      if (plus_LI_Pos[temp_AbsLoc] < Max_short) {
         if (Reads[Index].MatchedD == Plus) {
            plus_LI_Pos[temp_AbsLoc]++;
         }
      }
      if (minus_LI_Pos[temp_AbsLoc] < Max_short) {
         if (Reads[Index].MatchedD == Minus) {
            minus_LI_Pos[temp_AbsLoc]++;
         }
      }
   }
   std::vector < LI_Pos > LI_Positions;
   LI_Pos temp_LI_pos;
   bool SkipThis;
   int LI_Positions_Size = 0;
   bool SkipPlus;

   for (unsigned int Index_Minus = absStartBuffered; Index_Minus < absEndBuffered; Index_Minus++) {
      SkipPlus = false;
      for (unsigned int MaskedPosIndexMinus = Index_Minus + 10;
            MaskedPosIndexMinus >= Index_Minus - 10; MaskedPosIndexMinus--) {
         if (CurrentChrMask[MaskedPosIndexMinus] == 'B') {
            Index_Minus = MaskedPosIndexMinus + 10;
            SkipPlus = true;
            break;
         }
      }
      if (SkipPlus == false && minus_LI_Pos[Index_Minus] >= userSettings->NumRead2ReportCutOff) {
         for (unsigned int Index_Plus = Index_Minus - 1; Index_Plus <= Index_Minus + 30; Index_Plus++) {
            SkipThis = false;
            for (unsigned int MaskedPosIndexPlus = Index_Plus + 10; MaskedPosIndexPlus >= Index_Plus - 10; MaskedPosIndexPlus--) {
               if (CurrentChrMask[MaskedPosIndexPlus] == 'B') {
                  if ((MaskedPosIndexPlus + 10) > Index_Minus) {
                     Index_Minus = MaskedPosIndexPlus + 10;
                  }
                  SkipThis = true;
                  break;
               }
            }
            if (SkipThis == false && plus_LI_Pos[Index_Plus] >= userSettings->NumRead2ReportCutOff) {
               temp_LI_pos.Plus_Pos = Index_Plus;
               temp_LI_pos.Minus_Pos = Index_Minus;
               temp_LI_pos.WhetherReport = false;
               LI_Positions.push_back (temp_LI_pos);
               EventIndex_Pos[Index_Plus] = LI_Positions_Size;
               EventIndex_Pos[Index_Minus] = LI_Positions_Size;
               LI_Positions_Size++;
            }
         }
      }
   }
   LOG_DEBUG(*logStream << "LI: " << LI_Positions.size() << std::endl);
   static int Count_LI = 0;
   // find LI supporting reads


   for (unsigned Index = 0; Index < Reads.size (); Index++) {
      if (Reads[Index].Used || !Reads[Index].UP_Far.empty ()) {
         continue;
      }
      temp_AbsLoc =
         Reads[Index].UP_Close[Reads[Index].UP_Close.size () - 1].AbsLoc;
      if (EventIndex_Pos[temp_AbsLoc] == -1) {
         continue;
      }
      Reads[Index].Used = true;
      if (Reads[Index].MatchedD == Plus) {
         LI_Positions[EventIndex_Pos[temp_AbsLoc]].Plus_Reads.
         push_back (Index);
      }
      else {
         LI_Positions[EventIndex_Pos[temp_AbsLoc]].Minus_Reads.
         push_back (Index);
      }
   }

   std::vector < SPLIT_READ > temp_Plus_Reads, temp_Minus_Reads;

   bool temp_BalancedPlus_Plus, temp_BalancedPlus_Minus,
        temp_BalancedMinus_Plus, temp_BalancedMinus_Minus;
   short temp_LengthStr;

   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;

   for (unsigned LI_index = 0; LI_index < LI_Positions.size (); LI_index++) {
      if (LI_Positions[LI_index].Minus_Reads.empty ()
            || LI_Positions[LI_index].Plus_Reads.empty ()) {
         continue;
      }
      temp_BalancedPlus_Plus = false;
      temp_BalancedPlus_Minus = false;
      temp_BalancedMinus_Plus = false;
      temp_BalancedMinus_Minus = false;
      temp_Plus_Reads.clear ();
      temp_Minus_Reads.clear ();
      LOG_DEBUG(*logStream << "Here: " << LI_index <<
                "\t" << LI_Positions[LI_index].Minus_Reads.size() <<
                "\t" << LI_Positions[LI_index].Plus_Reads.size() << std::endl);
      for (unsigned int i = 0; i < LI_Positions[LI_index].Minus_Reads.size (); i++) {
         temp_Minus_Reads.
         push_back (Reads[LI_Positions[LI_index].Minus_Reads[i]]);
      }
      for (unsigned int i = 0; i < LI_Positions[LI_index].Minus_Reads.size (); i++) {
         UP_Close_index = temp_Minus_Reads[i].UP_Close.size () - 1;
         temp_LengthStr =
            temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
         if ((float) temp_LengthStr > temp_Minus_Reads[i].getReadLength() * 0.5) {
            temp_BalancedMinus_Plus = true;
         }
         else if ((float) temp_LengthStr <
                  temp_Minus_Reads[i].getReadLength() * 0.5) {
            temp_BalancedMinus_Minus = true;
         }
      }
      for (unsigned int i = 0; i < LI_Positions[LI_index].Plus_Reads.size (); i++) {
         temp_Plus_Reads.
         push_back (Reads[LI_Positions[LI_index].Plus_Reads[i]]);
      }
      for (unsigned int i = 0; i < LI_Positions[LI_index].Plus_Reads.size (); i++) {
         UP_Close_index = temp_Plus_Reads[i].UP_Close.size () - 1;
         temp_LengthStr =
            temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
         if ((float) temp_LengthStr > temp_Plus_Reads[i].getReadLength() * 0.5) {
            temp_BalancedPlus_Plus = true;
         }
         else if ((float) temp_LengthStr <
                  temp_Plus_Reads[i].getReadLength() * 0.5) {
            temp_BalancedPlus_Minus = true;
         }
      }

      unsigned int NumSupportPerTagPlus[g_sampleNames.size ()];
      unsigned int NumSupportPerTagMinus[g_sampleNames.size ()];
      for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
         NumSupportPerTagPlus[i] = 0;
         NumSupportPerTagMinus[i] = 0;
      }
      for (unsigned int i = 0; i < temp_Minus_Reads.size (); i++) {
         std::string currentTag = temp_Minus_Reads[i].Tag;
         int tagIndex = sampleToIndexMap[ currentTag ];
         NumSupportPerTagMinus[tagIndex]++;
      }
      for (unsigned int i = 0; i < temp_Plus_Reads.size (); i++) {
         std::string currentTag = temp_Plus_Reads[i].Tag;
         int tagIndex = sampleToIndexMap[ currentTag ];
         NumSupportPerTagPlus[tagIndex]++;
      }

      bool SupportedByOneSample = false;
      for (unsigned short j = 0; j < g_sampleNames.size (); j++) {
         LOG_DEBUG(*logStream << NumSupportPerTagPlus[j] << "\t" << NumSupportPerTagMinus[j] << std::endl);
         if (NumSupportPerTagPlus[j] > 0 && NumSupportPerTagMinus[j] > 0) {
            SupportedByOneSample = true;
            break;
         }
      }

      short PositiveBool = 0;
      if (temp_BalancedPlus_Plus) {
         PositiveBool++;
      }
      if (temp_BalancedPlus_Minus) {
         PositiveBool++;
      }
      if (temp_BalancedMinus_Plus) {
         PositiveBool++;
      }
      if (temp_BalancedMinus_Minus) {
         PositiveBool++;
      }

      if (SupportedByOneSample && PositiveBool >= 3) {
         {
            reportBreakDancerEvent(temp_Plus_Reads[0].FragName, LI_Positions[LI_index].
                                   Plus_Pos - g_SpacerBeforeAfter + 1, LI_Positions[LI_index].Minus_Pos -
                                   g_SpacerBeforeAfter +
                                   1 , -1, "LI", Count_LI);
            LargeInsertionOutf <<
                               "########################################################" <<
                               std::endl;
            LargeInsertionOutf << Count_LI++ << "\tLI\tChrID " <<
                               temp_Plus_Reads[0].FragName << "\t" << LI_Positions[LI_index].
                               Plus_Pos - g_SpacerBeforeAfter +
                               1 << "\t+ " << temp_Plus_Reads.
                               size () << "\t" << LI_Positions[LI_index].Minus_Pos -
                               g_SpacerBeforeAfter +
                               1 << "\t- " << temp_Minus_Reads.size ();
             
            for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
                LargeInsertionOutf << "\t" << indexToSampleMap[i] 
                                   << " + " << NumSupportPerTagPlus[i]
                                   << " - " << NumSupportPerTagMinus[i];
            }
            LargeInsertionOutf << std::endl;

            LargeInsertionOutf << (CurrentChr.substr (LI_Positions[LI_index].Plus_Pos - g_reportLength + 1, g_reportLength)) 
                               << Cap2Low (CurrentChr.substr (LI_Positions[LI_index].Plus_Pos + 1, g_reportLength)) << std::endl;
            for (unsigned int i = 0; i < temp_Plus_Reads.size (); i++) {
               UP_Close_index = temp_Plus_Reads[i].UP_Close.size () - 1;
               temp_LengthStr =
                  temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
               for (int j = 0; j < g_reportLength - temp_LengthStr; j++) {
                  LargeInsertionOutf << " ";
               }
               LargeInsertionOutf << temp_Plus_Reads[i].getUnmatchedSeqRev() << "\t"
                                  << temp_Plus_Reads[i].MatchedD << "\t"
                                  << temp_Plus_Reads[i].MatchedRelPos << "\t"
                                  << temp_Plus_Reads[i].MS << "\t"
                                  << temp_Plus_Reads[i].Tag << "\t" 
                                  << temp_Plus_Reads[i].Name << std::endl;
            }

            LargeInsertionOutf <<
                               "--------------------------------------------------------" <<
                               std::endl;
            LargeInsertionOutf << Cap2Low (CurrentChr.
                                           substr (LI_Positions[LI_index].
                                                   Minus_Pos - g_reportLength,
                                                   g_reportLength)) <<
                               (CurrentChr.
                                substr (LI_Positions[LI_index].Minus_Pos,
                                        g_reportLength)) << std::endl;
            for (unsigned int i = 0; i < temp_Minus_Reads.size (); i++) {
               UP_Close_index = temp_Minus_Reads[i].UP_Close.size () - 1;
               temp_LengthStr =
                  temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
               for (int j = 0;
                     j <
                     g_reportLength + temp_LengthStr -
                     temp_Minus_Reads[i].getReadLength(); j++) {
                  LargeInsertionOutf << " ";
               }
               LargeInsertionOutf << (temp_Minus_Reads[i].getUnmatchedSeq()) 
                                  << temp_Minus_Reads[i].MatchedD << "\t"
                                  << temp_Minus_Reads[i].MatchedRelPos << "\t"
                                  << temp_Minus_Reads[i].MS << "\t"
                                  << temp_Minus_Reads[i].Tag << "\t" 
                                  << temp_Minus_Reads[i].Name << std::endl;                
            }
         }

      }
   }
   LOG_INFO(*logStream << "Breakpoints for large insertions (LI): " << Count_LI << std::endl << std::endl);
	LargeInsertionOutf.close();
   Time_LI_E = time(NULL);
	*logStream << "Mining, Sorting and output LI results: " << (unsigned int) difftime(Time_LI_E, Time_LI_S) << " seconds." << std::endl << std::endl;
}

/** writes (appends) breakpoints to the file with the name "filename" */
void SortOutputRest (ControlState& currentState, const std::string & CurrentChr, 
                std::vector < SPLIT_READ > &Reads,
               // std::vector < SPLIT_READ > &BP_Reads, 
					const SearchWindow& window,
					 const std::string& filename)
{
   SPLIT_READ one_BP_read;
   unsigned UP_Close_index;
   unsigned temp_AbsLoc;
   LOG_DEBUG(*logStream << "1" << std::endl);
   // find LI combinations
	time_t Time_BP_S, Time_BP_E;
   Time_BP_S = time(NULL);


   const int BP_BORDER_BUFFER = 4 * g_maxInsertSize;
   std::ofstream Outf_Rest(filename.c_str(), std::ios::app);

   unsigned int absStartBPWindow = g_SpacerBeforeAfter + window.getStart() ;
   unsigned int absEndBPWindow = g_SpacerBeforeAfter + window.getEnd() ;
   if (absEndBPWindow > CurrentChr.size()- g_SpacerBeforeAfter ) {
      absEndBPWindow = CurrentChr.size()- g_SpacerBeforeAfter;
   }

   unsigned int absStartBuffered = absStartBPWindow - BP_BORDER_BUFFER;
   unsigned int absEndBuffered = absEndBPWindow + BP_BORDER_BUFFER;

   //*logStream << "SBA: " << g_SpacerBeforeAfter << " WS " << windowStart << " Total: " << g_SpacerBeforeAfter + windowStart - BP_BORDER_BUFFER << std::endl;
   ShiftedVector< uint8_t > plus_LI_Pos( absStartBuffered , absEndBuffered , 0 );
   ShiftedVector< uint8_t > minus_LI_Pos( absStartBuffered , absEndBuffered , 0 );

   for (unsigned Index = 0; Index < Reads.size (); Index++) {
      if (Reads[Index].Used || !Reads[Index].UP_Far.empty ()) {
         continue;
      }
      UP_Close_index = Reads[Index].UP_Close.size () - 1;
      temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
      if (plus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP) {
         if (Reads[Index].MatchedD == Plus) {
            plus_LI_Pos[temp_AbsLoc]++;
         }
      }
      if (minus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP) {
         if (Reads[Index].MatchedD == Minus) {
            minus_LI_Pos[temp_AbsLoc]++;
         }
      }
   }
   std::vector < Rest_Pos > Rest_Positions;
   Rest_Pos temp_Rest_pos;

   LOG_DEBUG(*logStream << "2" << std::endl);

   bool SkipThisPos;
   for (unsigned int Index = absStartBuffered; Index < absEndBuffered; Index++) {
         SkipThisPos = false;
         if (SkipThisPos == true) {
            continue;
         }
         if (plus_LI_Pos[Index] >= NumRead2ReportCutOff_BP) {
            temp_Rest_pos.Strand = Plus;
            temp_Rest_pos.Pos = Index;
            Rest_Positions.push_back (temp_Rest_pos);
         }
         if (minus_LI_Pos[Index] >= NumRead2ReportCutOff_BP) {
            temp_Rest_pos.Strand = Minus;
            temp_Rest_pos.Pos = Index;
            Rest_Positions.push_back (temp_Rest_pos);
         }
   }

      // find supporting reads
      for (unsigned Index = 0; Index < Reads.size (); Index++) {
         if (Reads[Index].Used || !Reads[Index].UP_Far.empty ()) {
            continue;
         }
         UP_Close_index = Reads[Index].UP_Close.size () - 1;
         temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
         for (unsigned Pos_index = 0; Pos_index < Rest_Positions.size ();
            Pos_index++) {
            if (Reads[Index].MatchedD == Rest_Positions[Pos_index].Strand) {
               if (temp_AbsLoc == Rest_Positions[Pos_index].Pos) {
                  Reads[Index].Used = true;
                  Rest_Positions[Pos_index].Pos_Reads.push_back (Index);	// copy index to save memory
               }
            }
         }
      }

   LOG_DEBUG(*logStream << "Other unassigned breakpoints (BP): " << Rest_Positions.size() << std::endl << std::endl);
   int Count_BP = 0;
   bool temp_BalancedPlus, temp_BalancedMinus;
   short temp_LengthStr;

    std::map<std::string,int> sampleToIndexMap;
    std::map<int,std::string> indexToSampleMap;
    createMaps( sampleToIndexMap, indexToSampleMap) ;
   for (unsigned LI_index = 0; LI_index < Rest_Positions.size (); LI_index++) {
      temp_BalancedPlus = false;
      temp_BalancedMinus = false;
      for (unsigned int i = 0; i < Rest_Positions[LI_index].Pos_Reads.size (); i++) {
         SPLIT_READ & CurrentSupportingRead = Reads[Rest_Positions[LI_index].Pos_Reads[i]];  
         UP_Close_index = CurrentSupportingRead.UP_Close.size () - 1;
         temp_LengthStr = CurrentSupportingRead.UP_Close[UP_Close_index].LengthStr;
         if ((float) temp_LengthStr > CurrentSupportingRead.getReadLength() * 0.5) {
            temp_BalancedPlus = true;
         }
         else if ((float) temp_LengthStr < CurrentSupportingRead.getReadLength() * 0.5) {
            temp_BalancedMinus = true;
         }
      }
      if (temp_BalancedPlus && temp_BalancedMinus) {
         Count_BP++;
          SupportPerSample NumSupportPerTag[g_sampleNames.size ()];
          for (unsigned short i = 0; i < g_sampleNames.size (); i++) {
              NumSupportPerTag[i].NumPlus = 0;
              NumSupportPerTag[i].NumMinus = 0;
          }

          for (unsigned int readIndex = 0; readIndex < Rest_Positions[LI_index].Pos_Reads.size(); readIndex++) {
             SPLIT_READ & CurrentSupportingRead = Reads[Rest_Positions[LI_index].Pos_Reads[readIndex]];
             std::string currentTag = CurrentSupportingRead.Tag;
             int tagIndex = sampleToIndexMap[ currentTag ];
             if (CurrentSupportingRead.MatchedD == Plus)	{
                NumSupportPerTag[tagIndex].NumPlus++;
             }
             else {
                NumSupportPerTag[tagIndex].NumMinus++;
             }
          }

          
          
          
         if (Rest_Positions[LI_index].Strand == Plus) {
            reportBreakDancerEvent(Reads[Rest_Positions[LI_index].Pos_Reads[0]].FragName,  0, 0, -1, "BP", -1);
            Outf_Rest <<
                      "########################################################" <<
                      std::endl;
            Outf_Rest << "ChrID " << Reads[Rest_Positions[LI_index].Pos_Reads[0]].FragName 
                      << "\t" << Rest_Positions[LI_index].Pos - g_SpacerBeforeAfter + 1 
                      << "\t+ " << Rest_Positions[LI_index].Pos_Reads.size();
            for (unsigned short SampleIndex = 0; SampleIndex < g_sampleNames.size (); ++SampleIndex) {
                Outf_Rest << "\t" << indexToSampleMap[SampleIndex] << " " << NumSupportPerTag[SampleIndex].NumPlus;
            }
            Outf_Rest << std::endl;

            Outf_Rest << (CurrentChr.substr (Rest_Positions[LI_index].Pos - g_reportLength + 1, g_reportLength)) 
                      << Cap2Low (CurrentChr.substr(Rest_Positions[LI_index].Pos + 1, g_reportLength))
                      << std::endl;
            for (unsigned int i = 0; i < Rest_Positions[LI_index].Pos_Reads.size(); i++) {
               SPLIT_READ & CurrentSupportingRead = Reads[Rest_Positions[LI_index].Pos_Reads[i]]; 
               UP_Close_index = CurrentSupportingRead.UP_Close.size () - 1;
               temp_LengthStr = CurrentSupportingRead.UP_Close[UP_Close_index].LengthStr;
               for (int j = 0; j < g_reportLength - temp_LengthStr; j++) {
                  Outf_Rest << " ";
               }
               Outf_Rest << CurrentSupportingRead.getUnmatchedSeqRev() << "\t" <<
                         CurrentSupportingRead.MatchedD << "\t" << CurrentSupportingRead.
                         MatchedRelPos << "\t" << CurrentSupportingRead.
                         MS << "\t" << CurrentSupportingRead.
                         Tag << "\t" << CurrentSupportingRead.Name << std::endl;
            }

         }
         else {
            reportBreakDancerEvent(Reads[Rest_Positions[LI_index].Pos_Reads[0]].FragName,  0,  0, -1, "BP", -1);
            Outf_Rest <<
                      "########################################################" <<
                      std::endl;
            Outf_Rest << "ChrID " << Reads[Rest_Positions[LI_index].Pos_Reads[0]].
                      FragName << "\t" << Rest_Positions[LI_index].Pos -
                      g_SpacerBeforeAfter +
                      1 << "\t- " << Rest_Positions[LI_index].Pos_Reads.
                      size ();// << "\t" << std::endl;
            for (unsigned short SampleIndex = 0; SampleIndex < g_sampleNames.size (); ++SampleIndex) {
                Outf_Rest << "\t" << indexToSampleMap[SampleIndex] << " " << NumSupportPerTag[SampleIndex].NumMinus;
            } 
            Outf_Rest << std::endl;
            Outf_Rest << Cap2Low (CurrentChr.
                                  substr (Rest_Positions[LI_index].Pos -
                                          g_reportLength,
                                          g_reportLength)) << (CurrentChr.
                                                substr
                                                (Rest_Positions
                                                 [LI_index].
                                                 Pos,
                                                 g_reportLength))
                      << std::endl;
            for (unsigned int i = 0; i < Rest_Positions[LI_index].Pos_Reads.size (); i++) {
                SPLIT_READ & CurrentSupportingRead = Reads[Rest_Positions[LI_index].Pos_Reads[i]];
               UP_Close_index = CurrentSupportingRead.UP_Close.size () - 1;
               temp_LengthStr =
                  CurrentSupportingRead.UP_Close[UP_Close_index].LengthStr;
               for (int j = 0; j < g_reportLength + temp_LengthStr - CurrentSupportingRead.getReadLength(); j++) {
                  Outf_Rest << " ";
               }
               Outf_Rest << (CurrentSupportingRead.getUnmatchedSeq())
                         << "\t" << CurrentSupportingRead.MatchedD
                         << "\t" << CurrentSupportingRead.MatchedRelPos
                         << "\t" << CurrentSupportingRead.MS
                         << "\t" << CurrentSupportingRead.Tag
                         << "\t" << CurrentSupportingRead.Name << std::endl;
            }
         }
      }
   }
   LOG_INFO(*logStream << "Other unassigned breakpoints (BP): " << Count_BP << std::endl << std::endl);
	Outf_Rest.close();
   Time_BP_E = time(NULL);
   *logStream << "Mining, Sorting and output BP results: " << (unsigned int) difftime(Time_BP_E, Time_BP_S) << " seconds." << std::endl << std::endl;
}

std::string GetConsensusInsertedStr(const std::vector <SPLIT_READ> & Reads, const int & StartIndex, const int & EndIndex) {
    // InsertedStr
    std::map<std::string,int> NT_str_2_count;
    std::map<std::string,int>::iterator it;
    for (int i = StartIndex; i <= EndIndex; i++) {
        it = NT_str_2_count.find(Reads[i].NT_str);
        if (it == NT_str_2_count.end()) {
            NT_str_2_count[Reads[i].NT_str] = 1;
        }
        else it->second++;
    }
    int Max = 0;
    std::string OutputStr = "";
    for (std::map<std::string,int>::iterator it=NT_str_2_count.begin(); it!=NT_str_2_count.end(); it++ ) {
        if (it->second > Max) {
            Max = it->second;
            OutputStr = it->first;
        }
    }
    return OutputStr;
}

std::string IntToString ( int number )
{
    std::ostringstream oss;
    
    // Works just like cout
    oss<< number;
    
    // Return the underlying string
    return oss.str();
}

void UpdateInterChromosomeCallAndSupport(std::map<std::string,int> & CallAndSupport, std::string & tempResult) {
    std::map<std::string,int> ::iterator it;
    it = CallAndSupport.find(tempResult);
    if (it == CallAndSupport.end()) { // no such element
        CallAndSupport.insert ( std::pair<std::string,int>(tempResult,1) );
    }
    else {
        (*it).second++;
    }
}

std::string OtherStrand(const char input) {
    if (input == '+') return "-";
    else if (input == '-') return "+";
    else return "";
}

std::string SameStrand(const char input) {
    if (input == '+') return "+";
    else if (input == '-') return "-";
    else return "";
}

void SortAndReportInterChromosomalEvents(ControlState& current_state, Genome& genome, UserDefinedSettings* user_settings)
{
    
    std::ofstream INToutputfile(user_settings->getINTOutputFilename().c_str(), std::ios::app);

	bool GoodRead;

	std::string InsertedSequence;

    //std::map Result;
    std::map<std::string,int> CallAndSupport;
    std::string tempResult;
    std::set<std::string> ChrNames;
	std::set<std::string> ReadNames;
	 std::set<std::string>::iterator it;
    std::set<std::string>::iterator first, second;
    std::cout << "reporting interchromosome variants" << std::endl;
    for (unsigned index = 0; index < current_state.InterChromosome_SR.size(); index++) {
	
        //std::cout << current_state.InterChromosome_SR[index];//.FragName << " " << current_state.InterChromosome_SR[index].FarFragName << std::endl;
        ChrNames.insert(current_state.InterChromosome_SR[index].FragName);
        ChrNames.insert(current_state.InterChromosome_SR[index].FarFragName);
    }
    //std::cout << "ChrNames.size() " << ChrNames.size() << std::endl;
    for (first = ChrNames.begin(); first != ChrNames.end(); first++) {
        //second = first;
        //second++;
        //if (second == ChrNames.end()) break;
        for (second = first; second != ChrNames.end(); second++) {
            if (first == second) {
                continue;
            }
            for (unsigned index = 0; index < current_state.InterChromosome_SR.size(); index++) {
		GoodRead = false;

                SPLIT_READ & currentRead = current_state.InterChromosome_SR[index];
		it = ReadNames.find(currentRead.Name);
		if (it != ReadNames.end()) {
			// this reads is already seen. skip it
			//std::cout << "skip one interchr read" << std::endl;
			continue;
		}
		else ReadNames.insert(currentRead.Name);
                if (currentRead.FragName == *first && currentRead.FarFragName == *second) {
                    //std::cout << *first << " " << *second << std::endl;
                    if (currentRead.MatchedD == '+') { // close smaller
                        for (unsigned CloseIndex = 0; CloseIndex < currentRead.UP_Close.size(); CloseIndex++) {
				if (currentRead.Used) break;
                            for (int FarIndex = currentRead.UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
				if (currentRead.Used) break;
                                if (currentRead.UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr == currentRead.getReadLength()) {
                                    currentRead.Used = true;
                                    currentRead.Left = currentRead. UP_Close[CloseIndex].AbsLoc - currentRead. UP_Close[CloseIndex].LengthStr + 1;
                                    if (currentRead.MatchedFarD == '+') 
                                        currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc - currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    else currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc + currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    currentRead.BP = currentRead.UP_Close[CloseIndex].LengthStr - 1;
                                    currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                                    currentRead.BPRight = currentRead.UP_Far[FarIndex]. AbsLoc - g_SpacerBeforeAfter;//0;//
					GoodRead = true;
					InsertedSequence = "\"\"";
					if (currentRead.BPLeft == 0 || currentRead.BPRight == 0) {
						std::cout << currentRead;
                    				//std::cout << tempResult << " " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
						//	<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;
					}
                                }
                            }
                        }
                    }
                    else { //close bigger
                        for (int CloseIndex = currentRead.UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
				if (currentRead.Used) break;
                            for (unsigned FarIndex = 0; FarIndex < currentRead.UP_Far.size(); FarIndex++) {
				if (currentRead.Used) break;
                                if (currentRead.UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr == currentRead.getReadLength()) {
                                    currentRead.Used = true;
                                    currentRead.Left = currentRead. UP_Close[CloseIndex].AbsLoc + currentRead. UP_Close[CloseIndex].LengthStr - 1;
                                    if (currentRead.MatchedFarD == '+')
                                        currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc - currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    else currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc + currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    currentRead.BP = currentRead.UP_Close[CloseIndex].LengthStr - 1;
                                    currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                                    currentRead.BPRight = currentRead.UP_Far[FarIndex]. AbsLoc - g_SpacerBeforeAfter;//0;//
					GoodRead = true;
					InsertedSequence = "\"\"";
					if (currentRead.BPLeft == 0 || currentRead.BPRight == 0) {
						std::cout << currentRead;
                    				//std::cout << tempResult << " " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
						//	<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;
					}
                                }
                            }
                        }
                    }
			if (GoodRead == false) {
				unsigned EffectiveLength = currentRead.UP_Close[currentRead.UP_Close.size() - 1].LengthStr + currentRead.UP_Far[currentRead.UP_Far.size() - 1].LengthStr;
				if (EffectiveLength >= 30 && 
					currentRead.UP_Close[currentRead.UP_Close.size() - 1].LengthStr >= 10 && currentRead.UP_Far[currentRead.UP_Far.size() - 1].LengthStr >= 10) {
					//if (currentRead.MatchedD == '+') {
					InsertedSequence = currentRead.UnmatchedSeq.substr(currentRead.UP_Far[currentRead.UP_Far.size() - 1].LengthStr, currentRead.getReadLength() - EffectiveLength);
					//}
					InsertedSequence = "\"" + InsertedSequence + "\"";
                              	      currentRead.BPLeft = currentRead.UP_Close[currentRead.UP_Close.size() - 1].AbsLoc - g_SpacerBeforeAfter;
                              	      currentRead.BPRight = currentRead.UP_Far[currentRead.UP_Far.size() - 1]. AbsLoc - g_SpacerBeforeAfter;//0;//
					GoodRead = true;
				}
			}

			if (GoodRead == true) {
                    		tempResult = "Anchor " + SameStrand(currentRead.MatchedD) + " " +
                                 currentRead.FragName + " " + IntToString(currentRead.BPLeft) + " " + OtherStrand(currentRead.MatchedD) + " " +
                                 currentRead.FarFragName + " " + IntToString(currentRead.BPRight)  + " " + SameStrand(currentRead.MatchedFarD) + " " + InsertedSequence;
				if (currentRead.BPLeft == 0 || currentRead.BPRight == 0) {
					//std::cout << currentRead;
                    			std::cout << tempResult << " " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
						<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;
				}
                    		UpdateInterChromosomeCallAndSupport(CallAndSupport, tempResult);
			}
			//else std::cout << "non-template in read " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
			//			<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;

                }
                else if (currentRead.FragName == *second && currentRead.FarFragName == *first) {
                    //std::cout << *second << " " << *first << std::endl;
                    if (currentRead.MatchedFarD == '-') { // close smaller
                        for (unsigned CloseIndex = 0; CloseIndex < currentRead.UP_Close.size(); CloseIndex++) {
				if (currentRead.Used) break;
                            for (int FarIndex = currentRead.UP_Far.size() - 1; FarIndex >= 0; FarIndex--) {
				if (currentRead.Used) break;
                                if (currentRead.UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr == currentRead.getReadLength()) {
                                    currentRead.Used = true;
                                    currentRead.Left = currentRead. UP_Far[FarIndex].AbsLoc - currentRead. UP_Close[CloseIndex].LengthStr + 1;
                                    if (currentRead.MatchedD == '+')
                                        currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc - currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    else currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc + currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    currentRead.BP = currentRead.UP_Close[CloseIndex].LengthStr - 1;
                                    currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                                    currentRead.BPRight = currentRead.UP_Far[FarIndex]. AbsLoc - g_SpacerBeforeAfter;//0;//
					GoodRead = true;
					InsertedSequence = "\"\"";
					if (currentRead.BPLeft == 0 || currentRead.BPRight == 0) {
						std::cout << currentRead;
                    				//std::cout << tempResult << " " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
						//	<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;
					}
                                }
                            }
                        }
                    }
                    else { //close bigger
                        for (int CloseIndex = currentRead.UP_Close.size() - 1; CloseIndex >= 0; CloseIndex--) {
				if (currentRead.Used) break;
                            for (unsigned FarIndex = 0; FarIndex < currentRead.UP_Far.size(); FarIndex++) {
				if (currentRead.Used) break;
                                if (currentRead.UP_Close[CloseIndex].LengthStr + currentRead.UP_Far[FarIndex].LengthStr == currentRead.getReadLength()) {
                                    currentRead.Used = true;
                                    currentRead.Left = currentRead. UP_Close[CloseIndex].AbsLoc + currentRead. UP_Close[CloseIndex].LengthStr - 1;
                                    if (currentRead.MatchedFarD == '+')
                                        currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc - currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    else currentRead.Right = currentRead.UP_Far[FarIndex]. AbsLoc + currentRead.UP_Far[FarIndex]. LengthStr - 1;
                                    currentRead.BP = currentRead.UP_Close[CloseIndex].LengthStr - 1;
                                    currentRead.BPLeft = currentRead.UP_Close[CloseIndex].AbsLoc - g_SpacerBeforeAfter;
                                    currentRead.BPRight = currentRead.UP_Far[FarIndex]. AbsLoc - g_SpacerBeforeAfter;//0;//
					GoodRead = true;
					InsertedSequence = "\"\"";
					if (currentRead.BPLeft == 0 || currentRead.BPRight == 0) {
						std::cout << currentRead;
                    				//std::cout << tempResult << " " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
						//	<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;
					}
                                }
                            }
                        }
                    }
			if (GoodRead == false) {
				unsigned EffectiveLength = currentRead.UP_Close[currentRead.UP_Close.size() - 1].LengthStr + currentRead.UP_Far[currentRead.UP_Far.size() - 1].LengthStr;
				if (EffectiveLength >= 30 && 
					currentRead.UP_Close[currentRead.UP_Close.size() - 1].LengthStr >= 10 && currentRead.UP_Far[currentRead.UP_Far.size() - 1].LengthStr >= 10) {
					//if (currentRead.MatchedD == '+') {
					InsertedSequence = currentRead.UnmatchedSeq.substr(currentRead.UP_Far[currentRead.UP_Far.size() - 1].LengthStr, currentRead.getReadLength() - EffectiveLength);
					//}
					InsertedSequence = "\"" + InsertedSequence + "\"";
                              	      currentRead.BPLeft = currentRead.UP_Close[currentRead.UP_Close.size() - 1].AbsLoc - g_SpacerBeforeAfter;
                              	      currentRead.BPRight = currentRead.UP_Far[currentRead.UP_Far.size() - 1]. AbsLoc - g_SpacerBeforeAfter;//0;//
					GoodRead = true;
				}
			}
			if (GoodRead == true) {
                    		tempResult = "Anchor " + SameStrand(currentRead.MatchedD) + " " +
                                 currentRead.FragName + " " + IntToString(currentRead.BPLeft) + " " + OtherStrand(currentRead.MatchedD) + " " +
                                 currentRead.FarFragName + " " + IntToString(currentRead.BPRight)  + " " + SameStrand(currentRead.MatchedFarD) + " " + InsertedSequence;
				if (currentRead.BPLeft == 0 || currentRead.BPRight == 0) {
					//std::cout << currentRead;
                    			std::cout << tempResult << " " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
						<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;
				}
                    		UpdateInterChromosomeCallAndSupport(CallAndSupport, tempResult);
			}
			//else std::cout << "non-template in read " << currentRead.FragName << " " << currentRead.MatchedRelPos << " " << currentRead.BPLeft << " " << currentRead.MatchedD << " " 
			//			<< currentRead.FarFragName << " " << currentRead.BPRight << " " << currentRead.MatchedFarD << " " << currentRead.Left << " " << currentRead.Right << std::endl;
                }
            }
        }
    }
    
    for ( std::map<std::string,int> ::iterator it = CallAndSupport.begin() ; it != CallAndSupport.end(); it++ ) {
        if ((*it).second >= 2) {
            INToutputfile << (*it).first << "\tsupport: " << (*it).second << std::endl;
		//std::cout << (*it).first << "\tsupport: " << (*it).second << std::endl;
	}
    }
        
    
}

void GetVariants(std::ifstream & input_file, std::vector <Variant> & variants) {
    Variant TempOne;
    TempOne.Report = true;
    unsigned NumberOfSamples;
    std::string tempStr;
    int tempInt;
    char tempChar;
    unsigned RefStartSupport, RefEndSupport, AlleleStartSupport, AlleleEndSupport;
    while (input_file >> tempStr) {
        if (tempStr.size() == 0) return;
        if (tempStr.substr(0, 2) == "##") {
            input_file >> tempInt >> TempOne.VariantType >> TempOne.Length
                       >> tempStr >> TempOne.NT_length >> TempOne.NT_str
                       >> tempStr >> TempOne.ChrName
                       >> tempStr >> TempOne.Start >> TempOne.End // BP
                       >> tempStr >> TempOne.Start_Range >> TempOne.End_Range // BP_range
                       >> tempStr >> TempOne.AlleleSupport >> tempInt // support
                       >> tempChar >> tempInt >> tempInt // +
                       >> tempChar >> tempInt >> tempInt // -
                       >> tempStr >> tempInt // S
                       >> tempStr >> tempInt // SUM_MS
                       >> NumberOfSamples // total number of samples
                       >> tempStr >> tempInt >> tempInt
                       >> tempStr >> RefStartSupport >> RefEndSupport
                       >> AlleleStartSupport >> tempInt >> AlleleEndSupport >> tempInt;
            getline(input_file, tempStr);
            if (NumberOfSamples != 1) {
                std::cout << "NumberOfSamples != 1: " << NumberOfSamples << std::endl;
                return;
            }
            for (unsigned i = 0; i <= TempOne.AlleleSupport; i++) {
                getline(input_file, tempStr);
            }
            if (RefStartSupport >= RefEndSupport) TempOne.RefSupport = RefStartSupport;
            else TempOne.RefSupport = RefEndSupport;
            TempOne.AlleleSupport = AlleleStartSupport + AlleleEndSupport;
            if (TempOne.AlleleSupport > 10 && TempOne.AlleleSupport > TempOne.RefSupport * 1.5) {
                
                variants.push_back(TempOne);
                //std::cout << TempOne.VariantType << " " << TempOne.Length << " " << TempOne.ChrName << " " << TempOne.Start << " " << TempOne.End << " " << TempOne.RefSupport << " " << TempOne.AlleleSupport << std::endl;
            }
        }
        else {
            std::cout << "someting is wrong here: " << tempStr << std::endl;
        }
    }
    //std::cout << "There are " << variants.size() << " variants reported." << std::endl;
}

bool CompareTwoVariants(const Variant & first, const Variant & second) {
    if (first.ChrName.size() != second.ChrName.size()) { // name length !=
        if (first.ChrName.size() < second.ChrName.size()) return true;
        else return false;
    }
    else if (first.ChrName != second.ChrName) { // name length equal
        for (unsigned CharIndex = 0; CharIndex < first.ChrName.size(); CharIndex++) {
            if ((int)first.ChrName[CharIndex] < (int)second.ChrName[CharIndex]) return true;
            else if ((int)first.ChrName[CharIndex] > (int)second.ChrName[CharIndex]) return false;
        }
    }
    else {
        if (first.Start_Range < second.Start_Range) return true;
        else if (first.Start_Range > second.Start_Range) return false;
        else if (first.End_Range < second.End_Range) return true;
        else if (first.End_Range > second.End_Range) return false;
        else {
            std::cout << "we shall never be here: WhetherExchange() in reporter.cpp" << std::endl;
            return true;
        }
    }
    return true;
}

void ModifyRefAndOutput(std::map<std::string, short> ChrName2Index, std::vector <Variant> & VariantsForReport, std::vector <VariantsPerChr> & WholeGenomeVariants, std::ofstream & Ploidy_Contig_Output_File) {
    //std::cout << "entering ModifyRefAndOutput" << std::endl;
    std::string NT_String, Left, Right;
    for (unsigned index = 0; index < VariantsForReport.size(); index++) {
        
        //std::cout << "Adding " << index << std::endl;
        WholeGenomeVariants[ChrName2Index.find(VariantsForReport[index].ChrName)->second].IndexOfVariants.push_back(index);
    }
    for (unsigned ChrIndex = 0; ChrIndex < WholeGenomeVariants.size(); ChrIndex++) {
        std::string & Current_ChrName = WholeGenomeVariants[ChrIndex].ChrName;
        std::string & Current_ChrSeq = WholeGenomeVariants[ChrIndex].ChrSeq;
        Current_ChrSeq = Current_ChrSeq.substr(g_SpacerBeforeAfter, Current_ChrSeq.size() - 2 * g_SpacerBeforeAfter);
        //std::cout << "Before: " << Current_ChrName << " " << Current_ChrSeq.size() << std::endl;
        if (WholeGenomeVariants[ChrIndex].IndexOfVariants.size()) {
            for (int VariantIndex = (int)WholeGenomeVariants[ChrIndex].IndexOfVariants.size() - 1; VariantIndex >= 0; VariantIndex--) {
                Variant & Current_Variant = VariantsForReport[WholeGenomeVariants[ChrIndex].IndexOfVariants[VariantIndex]];
                //std::cout << Current_Variant.VariantType << " " << Current_Variant.Length << " " << Current_Variant.NT_length << " " << Current_Variant.NT_str << " " << Current_Variant.ChrName << " " << Current_Variant.Start << " " << Current_Variant.End << std::endl;
                if (Current_Variant.NT_length == 0) NT_String = "";
                else NT_String = Current_Variant.NT_str.substr(1, Current_Variant.NT_length);
                //std::cout << "NT_String = |" << NT_String << "|" << std::endl;
                Left = Current_ChrSeq.substr(0, Current_Variant.Start);
                //std::cout << "Left = |" << Left.size() << "|" << std::endl;
                Right = Current_ChrSeq.substr(Current_Variant.End - 1, Current_ChrSeq.size() - Current_Variant.End + 1);
                //std::cout << "Right = |" << Right.size() << "|" << std::endl;
                Current_ChrSeq = Left + NT_String + Right;
            }
        }
        //std::cout << "After: " << Current_ChrName << " " << Current_ChrSeq.size() << std::endl;
        Ploidy_Contig_Output_File << ">" << Current_ChrName << "\n" << Current_ChrSeq << std::endl;
    }
    //std::cout << "leaving ModifyRefAndOutput" << std::endl;
}

void GetConsensusBasedOnPloidy(ControlState& current_state, Genome& genome, UserDefinedSettings* user_settings) {
    //for ( std::map<std::string, unsigned>::iterator it = g_ChrName2Ploidy.begin(); it != g_ChrName2Ploidy.end(); it++ ) {
    //    if ((*it).second != 1) {
    //        std::cout << "only accept ploidy 1 for testing at this stage" << std::endl;
    //        return;
    //    }
    //}
    std::ifstream D_Input_file(user_settings->getDOutputFilename().c_str());
    std::ifstream SI_input_file(user_settings->getSIOutputFilename().c_str());
    std::ofstream Ploidy_Variant_Output_File( (user_settings->getIndelConsensusOutputFilename()).c_str() );
    std::ofstream Ploidy_Contig_Output_File( (user_settings->getContigOutputFilename()).c_str() );
    std::vector <Variant> Variants;
    std::cout << "1" << std::endl;
    GetVariants(D_Input_file, Variants);
    GetVariants(SI_input_file, Variants);
    std::cout << "2" << std::endl;
    if (Variants.size() == 0) return;
    for (unsigned first = 0; first < Variants.size() - 1; first++) {
        if (Variants[first].Report == false) continue;
        for (unsigned second = first + 1; second < Variants.size(); second++) {
            if (Variants[second].Report == false) continue;
            if (Variants[first].ChrName == Variants[second].ChrName) {
                if (Variants[first].Start >= Variants[second].End || Variants[second].Start >= Variants[first].End) { // no overlap
                    continue;
                }
                else {
                    if (Variants[first].AlleleSupport >= Variants[second].AlleleSupport) {
                        Variants[second].Report = false;
                    }
                    else {
                        Variants[first].Report = false;
                        break;
                    }
                }
            }
            else continue;
        }
    }
    std::cout << "3" << std::endl;
    std::vector <Variant> VariantsForReport;
    for (unsigned index = 0; index < Variants.size(); index++) {
        if (Variants[index].Report) VariantsForReport.push_back(Variants[index]);
    }
    std::cout << "4" << std::endl;
    sort(VariantsForReport.begin(), VariantsForReport.end(), CompareTwoVariants);
    std::cout << "5" << std::endl;
    for (unsigned index = 0; index < VariantsForReport.size(); index++) {
            Ploidy_Variant_Output_File << "Type " << VariantsForReport[index].VariantType << "\t"
                                       << "Size " << VariantsForReport[index].Length << "\t"
                                       << "NT_Length " << VariantsForReport[index].NT_length << "\t"
                                       << "NT_String " << VariantsForReport[index].NT_str << "\t"
                                       << "ChrName " << VariantsForReport[index].ChrName << "\t"
                                       << "Start " << VariantsForReport[index].Start << "\t"
                                       << "End " << VariantsForReport[index].End << "\t"
                                       << "RefSupport " << VariantsForReport[index].RefSupport << "\t"
                                       << "VariantSupport " << VariantsForReport[index].AlleleSupport << std::endl;
    }
    std::cout << "6" << std::endl;
    std::vector <VariantsPerChr> WholeGenomeVariants;
    g_genome.reset();
    std::map<std::string, short> ChrName2Index;
    short ChrIndex = 0;
    while (const Chromosome* currentChromosome = g_genome.getNextChromosome()) {
        //std::cout << "ChrName: " << currentChromosome->getName() << "\tSize: " << (currentChromosome->getSeq()).size() - 2 * g_SpacerBeforeAfter << std::endl;
        VariantsPerChr tempOne;
        tempOne.ChrName = currentChromosome->getName();
        tempOne.ChrSeq = currentChromosome->getSeq();
        ChrName2Index.insert(std::pair<std::string, short>(tempOne.ChrName, ChrIndex) );
        ChrIndex++;
        WholeGenomeVariants.push_back(tempOne);
    }
    std::cout << "7" << std::endl;
    ModifyRefAndOutput(ChrName2Index, VariantsForReport, WholeGenomeVariants, Ploidy_Contig_Output_File);
    std::cout << "8" << std::endl;
}
