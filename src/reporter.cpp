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

// Pindel header files
#include "logstream.h"
#include "pindel.h"
#include "logdef.h"
#include "output_file_data.h"
#include "output_sorter.h"
#include "reporter.h"
#include "shifted_vector.h"

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
void calculateSupportPerTag( const std::vector< SPLIT_READ >& reads, const unsigned int firstReadIndex, const unsigned int lastReadIndex,
                             std::map<std::string,int>& sampleToIndexMap, SupportPerSample* NumSupportPerTag )
{
   for (unsigned int readIndex = firstReadIndex; readIndex <= lastReadIndex; readIndex++) {
      std::string currentTag = reads[readIndex].Tag;
      int tagIndex = sampleToIndexMap[ currentTag ];

      if (reads[readIndex].MatchedD == Plus)	{
         NumSupportPerTag[tagIndex].NumPlus++;
         if (reads[readIndex].UniqueRead) {
            NumSupportPerTag[tagIndex].NumUPlus++;
         }
      }
      else {
         NumSupportPerTag[tagIndex].NumMinus++;
         if (reads[readIndex].UniqueRead) {
            NumSupportPerTag[tagIndex].NumUMinus++;
         }
      }
   }
}

/* 'calculateSupportPerStrand (EWL, Aug31st2011) calculates the number of reads per strand (sample) supporting the event. */
void calculateSupportPerStrand( const std::vector< SPLIT_READ >& reads, const unsigned int firstReadIndex, const unsigned int lastReadIndex,
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

void OutputTDs (const std::vector < SPLIT_READ > &TDs,
           const std::string & TheInput,
           const unsigned int &C_S,
           const unsigned int &C_E,
           const unsigned int &RealStart,
           const unsigned int &RealEnd, std::ofstream & TDOutf)
{

   unsigned int NumberOfReads = C_E - C_S + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;

   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];
   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( TDs, C_S, C_E, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( TDs, C_S, C_E, LeftS, LeftUNum, RightS, RightUNum );

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
      Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = LeftS * RightS;

   CurrentChrMask[TDs[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[TDs[C_S].BPRight + g_SpacerBeforeAfter] = 'B';

   // reports BreakDancer event if applicable (-Q-option set and both ends of the SV referring to the same event)
   reportBreakDancerEvent(TDs[C_S].FragName, TDs[C_S].BPLeft, TDs[C_S].BPRight, TDs[C_S].IndelSize, "TD", NumberOfTDInstances);
   TDOutf <<
          "####################################################################################################"
          << std::endl;
   TDOutf << NumberOfTDInstances << "\tTD " << TDs[C_S].IndelSize	// << " bases "
          << "\tNT " << TDs[C_S].NT_size << " \"" << TDs[C_S].NT_str << "\"" << "\tChrID " << TDs[C_S].FragName << "\tBP " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2 << "\tBP_range " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111 << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += TDs[i].MS;
   }
   TDOutf << "\tSUM_MS " << SUM_MS;

   TDOutf << "\t" << g_sampleNames.size() << "\tNumSupSamples " << NumberSupSamples << "\t" <<
          NumU_SupSamples;
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      TDOutf << "\t" << indexToSampleMap[i] << " " << NumSupportPerTag[i].NumPlus
             << " " << NumSupportPerTag[i].NumUPlus
             << " " << NumSupportPerTag[i].NumMinus
             << " " << NumSupportPerTag[i].NumUMinus;
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
         TDOutf << ReverseComplement (TDs[GoodIndex].getUnmatchedSeq()) << std::endl;
      }
      TDOutf << "\t" << TDs[GoodIndex].MatchedD << "\t"
             << TDs[GoodIndex].MatchedRelPos
             << "\t" << TDs[GoodIndex].MS
             << "\t" << TDs[GoodIndex].Tag << "\t" << TDs[GoodIndex].Name << std::endl;
   }
}

void
OutputDeletions (const std::vector < SPLIT_READ > &Deletions,
                 const std::string & TheInput,
                 const unsigned int &C_S,
                 const unsigned int &C_E,
                 const unsigned int &RealStart,
                 const unsigned int &RealEnd, std::ofstream & DeletionOutf)
{
   LOG_DEBUG(*logStream << "d_1" << std::endl);
   unsigned int NumberOfReads = C_E - C_S + 1;
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
   calculateSupportPerStrand( Deletions, C_S, C_E, LeftS, LeftUNum, RightS, RightUNum );


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
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }
   LOG_DEBUG(*logStream << "d_5" << std::endl);
   unsigned int EasyScore = LeftS * RightS;
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
                << "\tNT " << Deletions[C_S].NT_size << " \"" << Deletions[C_S].NT_str << "\"" << "\tChrID " << Deletions[C_S].FragName << "\tBP " << Deletions[C_S].BPLeft + 1 << "\t" << Deletions[C_S].BPRight + 1 << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore;
   LOG_DEBUG(*logStream << "d_6" << std::endl);
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += Deletions[i].MS;
   }
   DeletionOutf << "\tSUM_MS " << SUM_MS;

   DeletionOutf << "\t" << g_sampleNames.
                size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
                NumU_SupSamples;
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      DeletionOutf << "\t" << indexToSampleMap[i] << " " << NumSupportPerTag[i].NumPlus
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
         DeletionOutf << ReverseComplement (Deletions[GoodIndex].getUnmatchedSeq()).substr (0, Deletions[GoodIndex].BP + 1);	// << endl;
         for (int i = 0; i < GapSize; i++) {
            DeletionOutf << " ";
         }
         DeletionOutf << ReverseComplement (Deletions[GoodIndex].getUnmatchedSeq()).substr (Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].getReadLength() - Deletions[GoodIndex].BP);	// << endl;
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

void OutputInversions (const std::vector < SPLIT_READ > &Inv,
                  const std::string & TheInput,
                  const unsigned int &C_S,
                  const unsigned int &C_E,
                  const unsigned int &RealStart,
                  const unsigned int &RealEnd, std::ofstream & InvOutf)
{
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
   unsigned int NumberOfReads = C_E - C_S + 1;
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
   calculateSupportPerStrand( Inv, C_S, C_E, LeftS, LeftUNum, RightS, RightUNum );

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
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }
   //*logStream << "+" << std::endl;
   unsigned int EasyScore = LeftS * RightS;
   CurrentChrMask[Inv[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[Inv[C_S].BPRight + g_SpacerBeforeAfter] = 'B';
   reportBreakDancerEvent(Inv[C_S].FragName, Inv[C_S].BPLeft, Inv[C_S].BPRight+2, Inv[C_S].IndelSize, "INV", g_numberOfInvInstances);
   InvOutf <<
           "####################################################################################################"
           << std::endl;
   InvOutf << g_numberOfInvInstances++ << "\tINV " << Inv[C_S].IndelSize	// << " bases "
           << "\tNT " << LeftNT_size << ":" << RightNT_size << " \"" << LeftNT_str << "\":\"" << RightNT_str << "\"" << "\tChrID " << Inv[C_S].FragName << "\tBP " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 << "\tBP_range " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += Inv[i].MS;
   }
   InvOutf << "\tSUM_MS " << SUM_MS;
   InvOutf << "\t" << g_sampleNames.
           size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
           NumU_SupSamples;
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      InvOutf << "\t" << indexToSampleMap[i] << " " << NumSupportPerTag[i].NumPlus
              << " " << NumSupportPerTag[i].NumUPlus
              << " " << NumSupportPerTag[i].NumMinus
              << " " << NumSupportPerTag[i].NumUMinus;
   InvOutf << std::endl;

   short SpaceBeforeReadSeq;
   InvOutf << TheInput.substr (Inv[C_S].BPLeft + g_SpacerBeforeAfter - g_reportLength, g_reportLength);	//;// << endl;// g_reportLength
   LOG_DEBUG(*logStream << Inv[C_S].NT_size << "\t" << Inv[C_S].NT_2size << std::endl);
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

            InvOutf << ReverseComplement (Inv[GoodIndex].getUnmatchedSeq());
             
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
            InvOutf << ReverseComplement (Inv[GoodIndex].getUnmatchedSeq());
         }
         InvOutf  << "\t" << Inv[GoodIndex].MatchedD << "\t"
                  << Inv[GoodIndex].MatchedRelPos
                  << "\t" << Inv[GoodIndex].MS
                  << "\t" << Inv[GoodIndex].Tag
                  << "\t" << Inv[GoodIndex].Name << std::endl;
      }
   }
}

void OutputSIs (const std::vector < SPLIT_READ > &SIs,
           const std::string & TheInput,
           const unsigned int &C_S,
           const unsigned int &C_E,
           const unsigned int &RealStart,
           const unsigned int &RealEnd, std::ofstream & SIsOutf)
{
   unsigned int NumberOfReads = C_E - C_S + 1;
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
   calculateSupportPerStrand( SIs, C_S, C_E, LeftS, LeftUNum, RightS, RightUNum );

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
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = LeftS * RightS;
   std::string CurrentReadSeq;

   CurrentChrMask[SIs[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[SIs[C_S].BPRight + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[RealStart + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[RealEnd + g_SpacerBeforeAfter] = 'B';

   reportBreakDancerEvent(SIs[C_S].FragName, SIs[C_S].BPLeft+1, SIs[C_S].BPRight+1, SIs[C_S].IndelSize, "SI", NumberOfSIsInstances);

   SIsOutf <<
           "####################################################################################################"
           << std::endl;
   SIsOutf << NumberOfSIsInstances << "\tI " << SIs[C_S].IndelSize << "\tNT " << SIs[C_S].IndelSize << " \"" << GetConsensusInsertedStr(SIs, C_S, C_E) << "\"" << "\tChrID " << SIs[C_S].FragName << "\tBP " << SIs[C_S].BPLeft + 1 << "\t" << SIs[C_S].BPRight + 1 << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += SIs[i].MS;
   }
   SIsOutf << "\tSUM_MS " << SUM_MS;

   SIsOutf << "\t" << g_sampleNames.
           size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
           NumU_SupSamples;
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      SIsOutf << "\t" << indexToSampleMap[i] << " " << NumSupportPerTag[i].NumPlus
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
         SIsOutf << ReverseComplement (SIs[GoodIndex].getUnmatchedSeq());
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
}

void
OutputDI (const std::vector < SPLIT_READ > &DI,
          const std::string & TheInput,
          const unsigned int &C_S,
          const unsigned int &C_E,
          const unsigned int &RealStart,
          const unsigned int &RealEnd, std::ofstream & DeletionOutf)
{
   unsigned int NumberOfReads = C_E - C_S + 1;
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
   calculateSupportPerStrand( DI, C_S, C_E, LeftS, LeftUNum, RightS, RightUNum );

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
      Num_U_Reads +=
         NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = LeftS * RightS;
   CurrentChrMask[DI[C_S].BPLeft + g_SpacerBeforeAfter] = 'B';
   CurrentChrMask[DI[C_S].BPRight + g_SpacerBeforeAfter] = 'B';
   reportBreakDancerEvent(DI[C_S].FragName, DI[C_S].BPLeft+1, DI[C_S].BPRight+1, DI[C_S].IndelSize, "D", deletionFileData.getSvIndex());
   DeletionOutf <<
                "####################################################################################################"
                << std::endl;
   DeletionOutf << deletionFileData.getSvIndex() << "\tD " << DI[C_S].IndelSize << "\tNT " << DI[C_S].NT_size << " \"" << DI[C_S].NT_str
<< "\"" << "\tChrID " << DI[C_S].FragName << "\tBP " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 << "\tBP_range " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;

   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) {
      SUM_MS += DI[i].MS;
   }
   DeletionOutf << "\tSUM_MS " << SUM_MS;


   DeletionOutf << "\t" << g_sampleNames.
                size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
                NumU_SupSamples;
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      DeletionOutf << "\t" << indexToSampleMap[i] << " " << NumSupportPerTag[i].NumPlus
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
         DeletionOutf << ReverseComplement (DI[GoodIndex].getUnmatchedSeq()) << "\t";
      }
      DeletionOutf << "\t" << DI[GoodIndex].MatchedD << "\t" << DI[GoodIndex].
                   MatchedRelPos << "\t" << DI[GoodIndex].MS << "\t" << DI[GoodIndex].
                   Tag << "\t" << DI[GoodIndex].Name << std::endl;
   }
}

void
SortOutputSI (const unsigned &NumBoxes, const std::string & CurrentChr,
              std::vector < SPLIT_READ > &Reads, std::vector < unsigned >SIs[],
              std::ofstream & SIsOutf)
{
   LOG_INFO(*logStream << "Sorting and outputing short insertions ..." << std::endl);
   unsigned int SIsNum;
   short CompareResult;
   unsigned Temp4Exchange;
   std::vector < SPLIT_READ > GoodIndels;
   unsigned int GoodNum;
   std::vector < Indel4output > IndelEvents;
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (SIs[Box_index].size () >= userSettings->NumRead2ReportCutOff) {
         SIsNum = SIs[Box_index].size ();
         LOG_DEBUG(*logStream << "SIsNum " << SIsNum << std::endl);
         for (unsigned int First = 0; First < SIsNum - 1; First++) {
            {
               for (unsigned int Second = First + 1; Second < SIsNum; Second++) {
                  {  
                      CompareResult = 0;
                     if (Reads[SIs[Box_index][First]].BPLeft < Reads[SIs[Box_index][Second]].BPLeft) {
                        continue;
                     }
                     else if (Reads[SIs[Box_index][First]].BPLeft > Reads[SIs[Box_index][Second]].BPLeft) {
                        CompareResult = 1;
                     }
                     else {
                        if (Reads[SIs[Box_index][First]].IndelSize < Reads[SIs[Box_index][Second]].IndelSize) {
                           continue;
                        }
                        else if (Reads[SIs[Box_index][First]].IndelSize > Reads[SIs[Box_index][Second]].IndelSize) {
                           CompareResult = 1;
                        }
                        else if (Reads[SIs[Box_index][First]].BP > Reads[SIs[Box_index][Second]].BP) {
                            CompareResult = 1; 
                        }
                     }
                     if (CompareResult == 1) {
                        Temp4Exchange = SIs[Box_index][First];
                        SIs[Box_index][First] = SIs[Box_index][Second];
                        SIs[Box_index][Second] = Temp4Exchange;
                     }
                  }
               }
            }
         }
          
         for (unsigned int First = 0; First < SIsNum - 1; First++) {
             for (unsigned int Second = First + 1; Second < SIsNum; Second++) {
                 if (Reads[SIs[Box_index][First]].getReadLength() == Reads[SIs[Box_index][Second]].getReadLength()) {
                     if (Reads[SIs[Box_index][First]].LeftMostPos ==
                         Reads[SIs[Box_index][Second]].LeftMostPos || Reads[SIs[Box_index][First]].LeftMostPos + Reads[SIs[Box_index][First]].getReadLength() ==
                         Reads[SIs[Box_index][Second]].LeftMostPos + Reads[SIs[Box_index][Second]].getReadLength()) {
                         if (Reads[SIs[Box_index][First]].MatchedD == Reads[SIs[Box_index][Second]].MatchedD) 
                         Reads[SIs[Box_index][Second]].UniqueRead = false;
                     }
                 }
             }
         }

          
         GoodIndels.clear ();
         IndelEvents.clear ();
         LOG_DEBUG(*logStream << GoodIndels.size() << std::endl);

         for (unsigned int First = 0; First < SIsNum; First++) {
            GoodIndels.push_back (Reads[SIs[Box_index][First]]);
         }

         GoodNum = GoodIndels.size ();
         LOG_DEBUG(*logStream << Box_index << " " << GoodNum << std::endl);
         if (GoodNum == 0) {
            continue;
         }
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
         OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.IndelStr = GoodIndels[0].NT_str;
         OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight = GoodIndels[0].BPRight;
         OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                  && GoodIndels[GoodIndex].IndelSize ==
                  OneIndelEvent.IndelSize)

            {
               OneIndelEvent.End = GoodIndex;
            }
            else {
               OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
               OneIndelEvent.RealEnd = OneIndelEvent.BPRight;

               OneIndelEvent.Support =
                  OneIndelEvent.End - OneIndelEvent.Start + 1;
               GetRealStart4Insertion (CurrentChr, OneIndelEvent.IndelStr,
                                       OneIndelEvent.RealStart,
                                       OneIndelEvent.RealEnd);
               IndelEvents.push_back (OneIndelEvent);
               OneIndelEvent.Start = GoodIndex;
               OneIndelEvent.End = GoodIndex;
               OneIndelEvent.BPLeft = GoodIndels[GoodIndex].BPLeft;
               OneIndelEvent.BPRight = GoodIndels[GoodIndex].BPRight;
               OneIndelEvent.IndelSize = GoodIndels[GoodIndex].IndelSize;
               OneIndelEvent.IndelStr = GoodIndels[GoodIndex].NT_str;
            }
         }
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;

         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Insertion (CurrentChr, OneIndelEvent.IndelStr,
                                 OneIndelEvent.RealStart,
                                 OneIndelEvent.RealEnd);
         IndelEvents.push_back (OneIndelEvent);

         if (IndelEvents.size ()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
                  EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  unsigned int RealStart = IndelEvents[EventIndex].RealStart;
                  unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
                  unsigned int IndelSize = IndelEvents[EventIndex].IndelSize;
                  unsigned int Max_Support = IndelEvents[EventIndex].Support;
                  unsigned int Max_Support_Index = EventIndex;

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
                     else if (IndelEvents[EventIndex_left].IndelSize != IndelSize) {
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
                  // report max one
                  LOG_DEBUG(*logStream << Max_Support << std::endl);
                  if (Max_Support >= userSettings->NumRead2ReportCutOff) {
                     OutputSIs (GoodIndels, CurrentChr, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End, RealStart, RealEnd, SIsOutf);
                     NumberOfSIsInstances++;
                  }
               }
            }
         }
      }												// if (!insertion[Box_index].empty())
   }
   LOG_INFO(*logStream << "Short insertions: " << NumberOfSIsInstances << std::endl << std::endl);
}

void SortAndOutputTandemDuplications (const unsigned &NumBoxes, const std::string & CurrentChr, std::vector < SPLIT_READ > &AllReads, std::vector < unsigned >TDs[],
                                 std::ofstream & TDOutf, const bool nonTemplate)
{

   if (nonTemplate) {
      LOG_INFO(*logStream << "Sorting and outputing tandem duplications with non-template sequence ..." << std::endl);
   }
   else {
      LOG_INFO(*logStream << "Sorting and outputing tandem duplications ..." << std::endl);
   }
   unsigned int TDNum;
   short CompareResult;
   unsigned Temp4Exchange;
   int countTandemDuplications = 0;
   unsigned int GoodNum;
   std::vector < SPLIT_READ > GoodIndels;
   std::vector < Indel4output > IndelEvents;
	UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (TDs[Box_index].size () >= userSettings->NumRead2ReportCutOff) {
         TDNum = TDs[Box_index].size ();

         for (unsigned int First = 0; First < TDNum - 1; First++) {
            {
               for (unsigned int Second = First + 1; Second < TDNum;  Second++) {
                  {  
                      CompareResult = 0;
                     if (AllReads[TDs[Box_index][First]].BPLeft < AllReads[TDs[Box_index][Second]].BPLeft) {
                        continue;
                     }
                     else if (AllReads[TDs[Box_index][First]].BPLeft > AllReads[TDs[Box_index][Second]].BPLeft) {
                        CompareResult = 1;
                     }
                     else if (AllReads[TDs[Box_index][First]].BPLeft == AllReads[TDs[Box_index][Second]].BPLeft) {
                        if (AllReads[TDs[Box_index][First]].BPRight < AllReads[TDs[Box_index][Second]].BPRight) {
                           continue;
                        }
                        else if (AllReads[TDs[Box_index][First]].BPRight > AllReads[TDs[Box_index][Second]].BPRight) {
                           CompareResult = 1;
                        }
                        else if (nonTemplate) {
                           if (AllReads[TDs[Box_index][First]].NT_size <
                                 AllReads[TDs[Box_index][Second]].NT_size) {
                              continue;
                           }
                           else if (AllReads[TDs[Box_index][First]].
                                    NT_size >
                                    AllReads[TDs[Box_index][Second]].
                                    NT_size) {
                              CompareResult = 1;
                           }
                           else if (AllReads[TDs[Box_index][First]].BP > AllReads[TDs[Box_index][Second]].BP) CompareResult = 1; 
                        }
                        else if (AllReads[TDs[Box_index][First]].BP > AllReads[TDs[Box_index][Second]].BP) CompareResult = 1; 
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
          
          for (unsigned int First = 0; First < TDNum - 1; First++) {
              for (unsigned int Second = First + 1; Second < TDNum; Second++) {
                  if (AllReads[TDs[Box_index][First]].getReadLength() == AllReads[TDs[Box_index][Second]].getReadLength()) {
                      if (AllReads[TDs[Box_index][First]].LeftMostPos ==
                          AllReads[TDs[Box_index][Second]].LeftMostPos || AllReads[TDs[Box_index][First]].LeftMostPos + AllReads[TDs[Box_index][First]].getReadLength() ==
                          AllReads[TDs[Box_index][Second]].LeftMostPos + AllReads[TDs[Box_index][Second]].getReadLength()) {
                          if (AllReads[TDs[Box_index][First]].MatchedD == AllReads[TDs[Box_index][Second]].MatchedD) 
                          AllReads[TDs[Box_index][Second]].UniqueRead = false;
                      }
                  }
              }
          }
          
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
                  unsigned int Max_Support = IndelEvents[EventIndex].Support;
                  unsigned int Max_Support_Index = EventIndex;
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
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= userSettings->NumRead2ReportCutOff) {
                     if (GoodIndels[ IndelEvents[Max_Support_Index].Start].IndelSize < userSettings->BalanceCutoff) {
                        OutputTDs (GoodIndels, CurrentChr, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End, RealStart, RealEnd, TDOutf);
                        NumberOfTDInstances++;
                        countTandemDuplications++;
                     }
                     else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
                        OutputTDs (GoodIndels, CurrentChr, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End, RealStart, RealEnd, TDOutf);
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


void SortOutputD (const unsigned &NumBoxes, const std::string & CurrentChr,
             std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Deletions[],
             std::ofstream & DeletionOutf)
{
   LOG_INFO(*logStream << "Sorting and outputing deletions ..." << std::endl);
   unsigned int DeletionsNum;
   short CompareResult;
   unsigned Temp4Exchange;

   unsigned int GoodNum;
   std::vector < SPLIT_READ > GoodIndels;
   std::vector < Indel4output > IndelEvents;
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      LOG_DEBUG(*logStream << Box_index << "\t" << NumBoxes << "\t" << Deletions[Box_index].size() << std::endl);
      if (Deletions[Box_index].size () >= userSettings->NumRead2ReportCutOff) {
         DeletionsNum = Deletions[Box_index].size ();

         for (unsigned int First = 0; First < DeletionsNum - 1; First++) {
            {
               for (unsigned int Second = First + 1; Second < DeletionsNum; Second++) {
                  { 
                      CompareResult = 0;
                     if (Reads[Deletions[Box_index][First]].BPLeft < Reads[Deletions[Box_index][Second]].BPLeft) {
                        continue;
                     }
                     else if (Reads[Deletions[Box_index][First]].BPLeft > Reads[Deletions[Box_index][Second]].BPLeft) {
                        CompareResult = 1;
                     }
                     else if (Reads[Deletions[Box_index][First]].BPLeft == Reads[Deletions[Box_index][Second]].BPLeft) {
                        if (Reads[Deletions[Box_index][First]].BPRight < Reads[Deletions[Box_index][Second]].BPRight) {
                           continue;
                        }
                        else if (Reads[Deletions[Box_index][First]].BPRight >Reads[Deletions[Box_index][Second]].BPRight) {
                           CompareResult = 1;
                        }
                        else if (Reads[Deletions[Box_index][First]].BP > Reads[Deletions[Box_index][Second]].BP) CompareResult = 1;
                     }
                     if (CompareResult == 1) {
                        Temp4Exchange = Deletions[Box_index][First];
                        Deletions[Box_index][First] =
                           Deletions[Box_index][Second];
                        Deletions[Box_index][Second] = Temp4Exchange;
                     }
                  }
               }
            }
         }
          for (unsigned int First = 0; First < DeletionsNum - 1; First++) {
              for (unsigned int Second = First + 1; Second < DeletionsNum; Second++) {
                  if (Reads[Deletions[Box_index][First]].getReadLength() == Reads[Deletions[Box_index][Second]].getReadLength()) {
                      if (Reads[Deletions[Box_index][First]].LeftMostPos ==
                          Reads[Deletions[Box_index][Second]].LeftMostPos || Reads[Deletions[Box_index][First]].LeftMostPos + Reads[Deletions[Box_index][First]].getReadLength() ==
                          Reads[Deletions[Box_index][Second]].LeftMostPos + Reads[Deletions[Box_index][Second]].getReadLength()) {
                          if (Reads[Deletions[Box_index][First]].MatchedD == Reads[Deletions[Box_index][Second]].MatchedD)
                          Reads[Deletions[Box_index][Second]].UniqueRead = false;
                      }
                  }
              }
          }
         GoodIndels.clear ();
         IndelEvents.clear ();
         for (unsigned int First = 0; First < DeletionsNum; First++) {
            GoodIndels.push_back (Reads[Deletions[Box_index][First]]);
         }

         GoodNum = GoodIndels.size ();
         LOG_DEBUG(*logStream << Box_index << " box read size " << GoodNum << std::endl);
         if (GoodNum == 0) {
            continue;
         }
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
         OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight = GoodIndels[0].BPRight;
         OneIndelEvent.WhetherReport = true;
         LOG_DEBUG(*logStream << "here" << std::endl);
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
            }
         }

         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion (CurrentChr, OneIndelEvent.RealStart,
                                OneIndelEvent.RealEnd);
         IndelEvents.push_back (OneIndelEvent);
         LOG_DEBUG(*logStream << "IndelEvent: " << IndelEvents.size() << std::endl);

         if (IndelEvents.size ()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
                  EventIndex++) {
               LOG_DEBUG(*logStream << EventIndex << " EventIndex" << std::endl);
               if (IndelEvents[EventIndex].WhetherReport) {
                  unsigned int RealStart = IndelEvents[EventIndex].RealStart;
                  unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
                  unsigned int Max_Support = IndelEvents[EventIndex].Support;
                  unsigned int Max_Support_Index = EventIndex;
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
                  // report max one
                  LOG_DEBUG(*logStream << "max" << std::endl);
                  if (IndelEvents[Max_Support_Index].Support >= userSettings->NumRead2ReportCutOff) {
                     LOG_DEBUG(*logStream << "aa" << std::endl);
                     if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < userSettings->BalanceCutoff) {
                        LOG_DEBUG(*logStream << "ba" << std::endl);
                        OutputDeletions (GoodIndels, CurrentChr, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End, RealStart,
									RealEnd, DeletionOutf);
                        deletionFileData.increaseTemplateSvCounter();
                        LOG_DEBUG(*logStream << "bb" << std::endl);
                     }
                     else if (ReportEvent( GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
                        LOG_DEBUG(*logStream << "ca" << std::endl);
                        OutputDeletions (GoodIndels, CurrentChr, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End, RealStart, RealEnd,
                                         DeletionOutf);
                        deletionFileData.increaseTemplateSvCounter();
                        LOG_DEBUG(*logStream << "cb" << std::endl);
                     }
                  }
               }
            }
         }
      }												// if (!Deletions[Box_index].empty())
   }														// for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   LOG_INFO(*logStream << "Deletions: " << deletionFileData.getTemplateSvCounter() << std::endl << std::endl);
}

void SortOutputInv (const unsigned &NumBoxes, const std::string & CurrentChr,
               std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Inv[],
               std::ofstream & InvOutf)
{
   OutputSorter os(NumBoxes, CurrentChr, InvOutf);
   os.SortAndOutputInversions(Reads, Inv);
}

void SortOutputInv_NT (const unsigned &NumBoxes, const std::string & CurrentChr,
                  std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Inv[],
                  std::ofstream & InvOutf)
{
   OutputSorter os(NumBoxes, CurrentChr, InvOutf);
   os.SortAndOutputNonTemplateInversions(Reads, Inv);
}

void OutputShortInversion (const std::vector < SPLIT_READ > &supportingReads,
                           const std::string &chromosome,
                           const unsigned int &indexOfFirstRead,
                           const unsigned int &indexOfLastRead,
                           const unsigned int &RealStart,
                           const unsigned int &RealEnd, std::ofstream & InversionOutF)
{
   unsigned int NumberOfReads = indexOfLastRead - indexOfFirstRead + 1;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
   unsigned int LeftUNum = 0;
   unsigned int RightUNum = 0;

   SupportPerSample NumSupportPerTag[g_sampleNames.size ()];
   std::map<std::string,int> sampleToIndexMap;
   std::map<int,std::string> indexToSampleMap;
   createMaps( sampleToIndexMap, indexToSampleMap) ;
   calculateSupportPerTag( supportingReads, indexOfFirstRead, indexOfLastRead, sampleToIndexMap, NumSupportPerTag );
   calculateSupportPerStrand( supportingReads, indexOfFirstRead, indexOfLastRead, LeftS, LeftUNum, RightS, RightUNum );

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
      Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
   }

   unsigned int EasyScore = LeftS * RightS;
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
                 supportingReads[indexOfFirstRead].BPRight + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 <<
                 "\t" << RightUNum << "\tS1 " << EasyScore;

   int SUM_MS = 0;
   for (unsigned int i = indexOfFirstRead; i <= indexOfLastRead; i++) {
      SUM_MS += supportingReads[i].MS;
   }
   InversionOutF << "\tSUM_MS " << SUM_MS;


   InversionOutF << "\t" << g_sampleNames.size() << "\tNumSupSamples " << NumberSupSamples << "\t" <<
                 NumU_SupSamples;
   for (unsigned short i = 0; i < g_sampleNames.size (); i++)
      InversionOutF << "\t" << indexToSampleMap[i] << " " << NumSupportPerTag[i].NumPlus
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
         InversionOutF << ReverseComplement (supportingReads[GoodIndex].getUnmatchedSeq()) << "\t";
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

void SortOutputDI (const unsigned &NumBoxes, const std::string & CurrentChr,
                   std::vector < SPLIT_READ > &Reads, std::vector < unsigned >DI[],
                   std::ofstream & DIOutf, std::ofstream & InvOutf)
{
   LOG_INFO(*logStream << "Sorting and outputing deletions with non-template sequences ..." <<
            std::endl);
   unsigned int DINum;
   short CompareResult;
   unsigned Temp4Exchange;

   unsigned int GoodNum;
   std::vector < SPLIT_READ > GoodIndels;
   std::vector < Indel4output > IndelEvents;
	UserDefinedSettings *userSettings = UserDefinedSettings::Instance();

   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      LOG_DEBUG(*logStream << "Box_index: "   << Box_index << std::endl);
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
            }
         }

         IndelEvents.push_back (OneIndelEvent);
         unsigned int RealStart;
         unsigned int RealEnd;
         for (unsigned EventIndex = 0; EventIndex < IndelEvents.size (); EventIndex++) {
            if (IndelEvents[EventIndex].End - IndelEvents[EventIndex].Start + 1 >= userSettings->NumRead2ReportCutOff) {
               RealStart =	GoodIndels[IndelEvents[EventIndex].Start].BPLeft;
               RealEnd = GoodIndels[IndelEvents[EventIndex].Start].BPRight;
               if (  (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < userSettings->BalanceCutoff)
                     || (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) ) {

                  if (IsInversion(GoodIndels[IndelEvents[EventIndex].Start], CurrentChr )) {
                     OutputShortInversion(GoodIndels, CurrentChr,IndelEvents[EventIndex].Start,
                                          IndelEvents[EventIndex].End,RealStart, RealEnd, InvOutf);
                     // increase the number of inversion instances in the OutputShortInversion itself; it'd add
                     // needless complexity here.
                  }
                  else {
                     OutputDI (GoodIndels, CurrentChr,IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End,
                               RealStart, RealEnd, DIOutf);
                     deletionFileData.increaseNonTemplateSvCounter();
                  }
               }
            }
         }
      }												// if (!Deletions[Box_index].empty())
   }														// for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   LOG_INFO(*logStream << "deletions with non-template sequences: " << deletionFileData.getNonTemplateSvCounter() <<
            std::endl << std::endl);
}



void SortOutputLI (const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads, const SearchWindow& window, const std::string& filename)
{
   unsigned UP_Close_index;
   unsigned temp_AbsLoc;
   unsigned int LI_BORDER_BUFFER = 4 * g_maxInsertSize;
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
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
               LargeInsertionOutf << ReverseComplement (temp_Plus_Reads[i].getUnmatchedSeq()) << "\t"
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
void SortOutputRest (const std::string & CurrentChr, 
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
               Outf_Rest << ReverseComplement (CurrentSupportingRead.
                                               getUnmatchedSeq()) << "\t" <<
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
