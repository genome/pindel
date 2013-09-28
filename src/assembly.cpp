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

#include <iostream>
#include <vector>
#include "assembly.h"
#include "reader.h"
#include "pindel.h"
#include "farend_searcher.h"
#include "control_state.h"
#include "user_defined_settings.h"
#include <map>
#include <string>
#include <utility>

void doAssembly (ControlState & CurrentState, UserDefinedSettings* userSettings){// std::ifstream& FastaFile ) {
        /*
        // step 1: define output file
       // std:string AssemblyOutputFilename = inf_AssemblyInputFilename + "_ASM";
        std::ofstream AssemblyOutput(AssemblyOutputFilename.c_str()); 
        
        // step 1: get whole genome as vector of chr
        
        // step 2: get a list of breakpoints 
         */
    
    
    //std::ifstream FastaFile((userSettings->getRefFilename()).c_str());
    std::cout << "Entering assembly mode ..." << std::endl; 
    std::map<std::string,int> ChrName2Index;
    
    std::cout << "Get whole genome sequence..." << std::endl;
    // step 1. get the whole genome sequence

    g_genome.loadAll( userSettings->getRefFilename() );


    //std::vector <AChromosome> AllChromosomes;
    //getWholeGenome(FastaFile, AllChromosomes);
    unsigned CountChr = 0;
    do {
        const Chromosome* currentChromosome = g_genome.getNextChromosome();
        if (currentChromosome == NULL) {
			break;
		}
        std::cout << "ChrName " << currentChromosome->getName() << "\tChrSeqSize " << currentChromosome->getBiolSize() << std::endl;
        ChrName2Index[currentChromosome->getName()] = CountChr;
        CountChr++;
        
    } while (true);

    // step 2. get all SVs
    //CurrentState.inf_AssemblyInput.open(par.inf_AssemblyInputFilename.c_str());
    std::cout << "\nGet all SVs to assemble..." << std::endl;
    std::vector <Assembly> AllSV4Assembly;
    Assembly OneSV;
    unsigned SV_Count = 0;
    while (CurrentState.inf_AssemblyInput >> OneSV.Type
                                          >> OneSV.ChrA
                                          >> OneSV.PosA
                                          >> OneSV.CI_A
                                          >> OneSV.ChrB
                                          >> OneSV.PosB
                                          >> OneSV.CI_B) {
        OneSV.Index = SV_Count++;
        if (OneSV.ChrA == OneSV.ChrB) {
            if (OneSV.PosA > OneSV.PosB) {
                unsigned Exchange = OneSV.PosA;
                OneSV.PosA = OneSV.PosB;
                OneSV.PosB = Exchange;
            }
        }
        AllSV4Assembly.push_back(OneSV);
    }
        
    
    // step 3. per SV, find BP, sort and output
    // step 3.1 define output
    //std::string ASM_OutputFileName = userSettings->getAOutputFolder + "_ASM";
    std::ofstream ASM_Output(userSettings->getASMOutputFilename().c_str());
    // step 3.2 per SV, collect reads around BP, sort, combine and output.
    for (unsigned SV_index = 0; SV_index < AllSV4Assembly.size(); SV_index++) {
        std::cout << AllSV4Assembly[SV_index].Type << "\t" << AllSV4Assembly[SV_index].ChrA << "\t" << AllSV4Assembly[SV_index].PosA << "\t" << AllSV4Assembly[SV_index].CI_A << "\t" << AllSV4Assembly[SV_index].ChrB << "\t" << AllSV4Assembly[SV_index].PosB << "\t" << AllSV4Assembly[SV_index].CI_B << std::endl;
        AssembleOneSV(ChrName2Index, CurrentState, AllSV4Assembly[SV_index], ASM_Output);
    }
    std::cout << "Leaving assembly mode and terminating this run." << std::endl;
     
    return;
}


short AssembleOneSV(std::map<std::string,int> &ChrName2Index, ControlState & CurrentState, const Assembly & OneSV, std::ofstream & ASM_Output) {
    //std::cout << "AssembleOneSV 1" << std::endl;
    //short Max_NT_Size = 30;
    bool WhetherFirstBP = true;
    std::vector <SPLIT_READ> First, Second;
    unsigned SearchCenter;
    unsigned SearchRange;
    std::cout << "Current SV: " << OneSV.Index << " " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " << OneSV.CI_A 
              << "\t" << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV.CI_B << std::endl;
    // get first BP
    CurrentState.Reads_SR.clear();
    CurrentState.OneEndMappedReads.clear();
    //std::cout << "AssembleOneSV 2" << std::endl;
    const Chromosome* currentChromosomeFirst = g_genome.getChr(OneSV.ChrA);
    const Chromosome* currentChromosomeSecond = g_genome.getChr(OneSV.ChrB);
    std::string FirstChrSeq = currentChromosomeFirst->getSeq();
    std::string SecondChrSeq = currentChromosomeSecond->getSeq();
    
    if (CurrentState.CurrentChrName != OneSV.ChrA) {
        CurrentState.CurrentChrName = OneSV.ChrA;
        CurrentState.CurrentChrSeq = FirstChrSeq;
        //CurrentState.CurrentChrSeq = AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq; // change later, copying one chrseq for each SV is expensive.
    }

    //const unsigned CONS_Chr_Size = currentChromosome->getBiolSize(); // #################
    //std::cout << "CONS_Chr_Size " << CONS_Chr_Size << std::endl;
    g_maxPos = 0; // #################
    unsigned Left, Right;
	unsigned int lowerBinBorder = 1;
    if (OneSV.PosA > OneSV.CI_A + 1000)  
        lowerBinBorder = OneSV.PosA - OneSV.CI_A - 1000; //CurrentState.
    unsigned int upperBinBorder = OneSV.PosA + OneSV.CI_A + 1000;
	SearchWindow window( currentChromosomeFirst, lowerBinBorder, upperBinBorder );
    Left = OneSV.PosA + g_SpacerBeforeAfter - OneSV.CI_A;
    Right = OneSV.PosA + g_SpacerBeforeAfter + OneSV.CI_A;

    std::cout << "\nFirst BP\tChrName " << CurrentState.CurrentChrName << "\tRange " <<  lowerBinBorder << " " << upperBinBorder << std::endl;
    get_SR_Reads(CurrentState, window );

    //std::cout << "First size: " << CurrentState.Reads_SR.size() << std::endl;
    CombineAndSort(ChrName2Index, CurrentState, OneSV, First, lowerBinBorder, upperBinBorder, WhetherFirstBP);

    CleanUpCloseEnd(First, Left, Right); // vector of reads

    std::cout << "\nFirst size " << First.size() << std::endl;
    SearchRange = OneSV.CI_B + 1000;
    SearchCenter = OneSV.PosB + g_SpacerBeforeAfter;
	std::vector< SearchWindow > searchCluster;
	searchCluster.push_back( SearchWindow( currentChromosomeSecond, SearchCenter-SearchRange, SearchCenter+SearchRange ) );
    Left = OneSV.PosB + g_SpacerBeforeAfter - OneSV.CI_B;
    Right = OneSV.PosB + g_SpacerBeforeAfter + OneSV.CI_B;
    
	for (unsigned ReadIndex = 0; ReadIndex < First.size(); ReadIndex++) {
      First[ReadIndex].FarFragName = OneSV.ChrB;
      //SearchFarEndAtPos(AllChromosomes[ChrName2Index.find(OneSV.ChrB)->second].ChrSeq, First[ReadIndex], SearchCenter, SearchRange);
		SearchFarEndAtPos(First[ReadIndex], searchCluster);
    }
    //std::cout << "AssembleOneSV 7" << std::endl;

    CleanUpFarEnd(First, Left, Right);
    //std::cout << "AssembleOneSV 8" << std::endl;


	for (unsigned ReadIndex = 0; ReadIndex < First.size(); ReadIndex++) {
        if (First[ReadIndex].UP_Close.size()) {
			if (First[ReadIndex].UP_Far.size()) {
				OutputCurrentRead(ChrName2Index, CurrentState, OneSV, First[ReadIndex], ASM_Output);
            }
        }
    }
    
    
    // get second BP
    CurrentState.Reads_SR.clear();
    WhetherFirstBP = false;
    CurrentState.CurrentChrName = OneSV.ChrB;
    
    CurrentState.CurrentChrSeq = SecondChrSeq;
    //CONS_Chr_Size = CurrentState.CurrentChrSeq.size() - 2 * g_SpacerBeforeAfter; // #################
    g_maxPos = 0; // #################
	lowerBinBorder = 1;
    if (OneSV.PosB > OneSV.CI_B + 1000)  
        lowerBinBorder = OneSV.PosB - OneSV.CI_B - 1000;
    upperBinBorder = OneSV.PosB + OneSV.CI_B + 1000;
	SearchWindow window2( currentChromosomeSecond, lowerBinBorder, upperBinBorder );
    Left = OneSV.PosB + g_SpacerBeforeAfter - OneSV.CI_B;
    Right = OneSV.PosB + g_SpacerBeforeAfter + OneSV.CI_B;

    std::cout << "\nSecond BP\tChrName " << CurrentState.CurrentChrName << "\tRange " << lowerBinBorder << " " << upperBinBorder << std::endl;    
    get_SR_Reads(CurrentState, window2 );

    CombineAndSort(ChrName2Index, CurrentState, OneSV, Second, lowerBinBorder, upperBinBorder, WhetherFirstBP);

    CleanUpCloseEnd(Second, Left, Right);

    std::cout << "\nSecond size " << Second.size() << std::endl;
    SearchRange = OneSV.CI_A + 1000;
    SearchCenter = OneSV.PosA + g_SpacerBeforeAfter;
    Left = OneSV.PosA + g_SpacerBeforeAfter - OneSV.CI_A;
    Right = OneSV.PosA + g_SpacerBeforeAfter + OneSV.CI_A;
	searchCluster.clear();
	searchCluster.push_back( SearchWindow( currentChromosomeFirst, SearchCenter-SearchRange, SearchCenter+SearchRange ) );
  
   for (unsigned ReadIndex = 0; ReadIndex < Second.size(); ReadIndex++) {
      Second[ReadIndex].FarFragName = OneSV.ChrA;
      //SearchFarEndAtPos(AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq, Second[ReadIndex], SearchCenter, SearchRange);
		SearchFarEndAtPos(Second[ReadIndex], searchCluster);
   }

    CleanUpFarEnd(Second, Left, Right);

   for (unsigned ReadIndex = 0; ReadIndex < Second.size(); ReadIndex++) {
      if (Second[ReadIndex].UP_Close.size()) {
         if (Second[ReadIndex].UP_Far.size()) {
            OutputCurrentRead(ChrName2Index, CurrentState, OneSV, Second[ReadIndex], ASM_Output);
         }
      }
   }

    unsigned SumSize = 0;
    for (unsigned ReadIndex = 0; ReadIndex < First.size(); ReadIndex++) {
        SumSize += First[ReadIndex].UP_Far.size(); 
    }
    for (unsigned ReadIndex = 0; ReadIndex < Second.size(); ReadIndex++) {
        SumSize += Second[ReadIndex].UP_Far.size(); 
    }
    if (SumSize == 0 && OneSV.ChrA == OneSV.ChrB) {
        TryLI(ChrName2Index, CurrentState, OneSV, First, Second, ASM_Output);
    }
    return 0;
}

void CombineAndSort(std::map<std::string,int> &ChrName2Index, ControlState & CurrentState, const Assembly & OneSV, std::vector <SPLIT_READ> & output_reads, const unsigned & lowerBinBorder, const unsigned & upperBinBorder, const bool & First) {
    
    const unsigned AssemblyCutOff = 3;
    const char Plus = '+';
    const char Minus = '-';
    unsigned WindowSize = (upperBinBorder - lowerBinBorder) * 3;
    unsigned OffSet;
    if (lowerBinBorder * 2 > upperBinBorder) OffSet = lowerBinBorder * 2 - upperBinBorder;
    else OffSet = 0;
    std::vector<unsigned int> Index_Plus[WindowSize], Index_Minus[WindowSize];
    std::cout << "\nin CombineAndSort..." << std::endl;
    // get reads anchored at plus and minus
    std::vector <int> Plus_Read_Index, Minus_Read_Index;
    for (unsigned i = 0; i < CurrentState.Reads_SR.size(); i++) {
        if (First) {
            if (CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc + OneSV.CI_A + CurrentState.Reads_SR[i].getReadLength() > g_SpacerBeforeAfter + OneSV.PosA && CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc < g_SpacerBeforeAfter + OneSV.PosA + OneSV.CI_A + CurrentState.Reads_SR[i].getReadLength()) {
                if (CurrentState.Reads_SR[i].MatchedD == Plus) Index_Plus[CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc - OffSet - g_SpacerBeforeAfter].push_back(i);
                else if (CurrentState.Reads_SR[i].MatchedD == Minus) Index_Minus[CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc - OffSet - g_SpacerBeforeAfter].push_back(i);
            }
        }
        else { // second BP
            if (CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc + OneSV.CI_B + CurrentState.Reads_SR[i].getReadLength() > g_SpacerBeforeAfter + OneSV.PosB && CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc < g_SpacerBeforeAfter+ OneSV.PosB + OneSV.CI_B + CurrentState.Reads_SR[i].getReadLength()) {
                if (CurrentState.Reads_SR[i].MatchedD == Plus) Index_Plus[CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc - OffSet - g_SpacerBeforeAfter].push_back(i);
                else if (CurrentState.Reads_SR[i].MatchedD == Minus) Index_Minus[CurrentState.Reads_SR[i].UP_Close[CurrentState.Reads_SR[i].UP_Close.size() - 1].AbsLoc - OffSet - g_SpacerBeforeAfter].push_back(i);
            }
        }
    }
    for (unsigned i = 0; i < WindowSize; i++) if (Index_Plus[i].size() >= AssemblyCutOff || Index_Minus[i].size() >= AssemblyCutOff) {
       std::cout << "Candidate: " << i << " " << i + OffSet << "\t+ " << Index_Plus[i].size() << "\t-" << Index_Minus[i].size() << std::endl;    
       if (Index_Plus[i].size() >= AssemblyCutOff) 
          CombineReads(CurrentState.CurrentChrSeq, Plus, CurrentState.Reads_SR, Index_Plus[i], output_reads);
       if (Index_Minus[i].size() >= AssemblyCutOff) CombineReads(CurrentState.CurrentChrSeq, Minus, CurrentState.Reads_SR, Index_Minus[i], output_reads); 
    }
    //std::cout << Plus_Read_Index.size() << " " << Minus_Read_Index.size() << std::endl;
}

void CombineReads(const std::string & CurrentChrSeq, const char & Strand, const std::vector <SPLIT_READ> & input_reads, const std::vector <unsigned int> & Index_Of_Useful_Reads, std::vector <SPLIT_READ> & output_reads) {
    //std::cout << "start of CombineReads" << std::endl;
    std::string Spacer = "";
    unsigned Max_ReadLength = 0;
    unsigned Max_AssembledLength = 0;
    unsigned Min_LeftMostPos = input_reads[Index_Of_Useful_Reads[0]].LeftMostPos;
    SPLIT_READ output_one_read;// = input_reads[Index_Of_Useful_Reads[0]];
    unsigned Index2Read4Copy = 0;
    //unsigned Min_Close_Size = 10000;
    for (unsigned ReadIndex = 0; ReadIndex < Index_Of_Useful_Reads.size(); ReadIndex++) {
        //if (input_reads[Index_Of_Useful_Reads[ReadIndex]].UP_Close.size() < Min_Close_Size) {
        //    Min_Close_Size = input_reads[Index_Of_Useful_Reads[ReadIndex]].UP_Close.size();
        //    Index2Read4Copy = ReadIndex; // input_reads[Index_Of_Useful_Reads[Index2Read4Copy]]
        //}
        //std::cout << Strand << " " << input_reads[Index_Of_Useful_Reads[ReadIndex]].UP_Close.size() << std::endl;
        //std::cout << input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos << " " << input_reads[Index_Of_Useful_Reads[ReadIndex]].UnmatchedSeq << std::endl;
        if (input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos < (int)Min_LeftMostPos)
            Min_LeftMostPos = input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos;
        if (input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength() > (short)Max_ReadLength)
            Max_ReadLength = input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength();
        if (input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos + (unsigned)input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength() > Max_AssembledLength) 
            Max_AssembledLength = input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos + input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength();
    }
    Max_AssembledLength = Max_AssembledLength - Min_LeftMostPos;
    if ((float)Max_AssembledLength < Max_ReadLength * 1.3) return;
    std::cout << "Max_AssembledLength " << Max_AssembledLength << std::endl;
    for (unsigned ReadIndex = 0; ReadIndex < Index_Of_Useful_Reads.size(); ReadIndex++) {
        Spacer.clear();
        
        if (Strand == '+') {
            for (unsigned SpacerIndex = 0; SpacerIndex < Max_AssembledLength + Min_LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength(); SpacerIndex++) Spacer += " ";
            //std::cout << Spacer << (input_reads[Index_Of_Useful_Reads[ReadIndex]].UnmatchedSeq) << std::endl;
        }
        else {
            for (unsigned SpacerIndex = 0; SpacerIndex < input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - Min_LeftMostPos; SpacerIndex++) 
                Spacer += " ";
            //std::cout << Spacer << (input_reads[Index_Of_Useful_Reads[ReadIndex]].UnmatchedSeq) << std::endl;
        }
    }
    if (Strand == '+') { // UnmatchedSeq
        std::cout << "+ Max_ReadLength " << Max_ReadLength << "\t" << "Min_LeftMostPos " << Min_LeftMostPos 
        << "\nref: \n" << ReverseComplement(CurrentChrSeq.substr(Min_LeftMostPos, Max_AssembledLength)) << std::endl;
    }
    else if (Strand == '-') {
        std::cout << "- Max_ReadLength " << Max_ReadLength << "\t" << "Min_LeftMostPos " << Min_LeftMostPos 
        << "\nref: \n" << CurrentChrSeq.substr(Min_LeftMostPos, Max_AssembledLength) << std::endl;
    }
    unsigned Count[5][Max_AssembledLength];
    float Ratio[5][Max_AssembledLength];
    for (short i = 0; i < 5; i++) {
        for (unsigned j = 0; j < Max_AssembledLength; j++) {
            Count[i][j] = 0;
            Ratio[i][j] = 0.0;
        }
    }
    if (Strand == '+') {
        for (unsigned ReadIndex = 0; ReadIndex < Index_Of_Useful_Reads.size(); ReadIndex++) {
            //std::cout << std::endl;
            for (short BaseIndex = 0; BaseIndex < input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength(); BaseIndex++) {
                //std::cout << input_reads[Index_Of_Useful_Reads[ReadIndex]].UnmatchedSeq[BaseIndex];
                switch (input_reads[Index_Of_Useful_Reads[ReadIndex]].getUnmatchedSeq()[BaseIndex]) {
                    case 'A':
                        Count[0][Max_AssembledLength + Min_LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength() + BaseIndex]++;
                        break;	// 00000000
                    case 'C':
                        Count[1][Max_AssembledLength + Min_LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength() + BaseIndex]++;
                        break;	// 00010000
                    case 'G':
                        Count[2][Max_AssembledLength + Min_LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength() + BaseIndex]++;
                        break;	// 00100000
                    case 'T':
                        Count[3][Max_AssembledLength + Min_LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength() + BaseIndex]++;
                        break;	// 00110000
                    default:
                        Count[4][Max_AssembledLength + Min_LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength() + BaseIndex]++;
                        // 01000000
                }
            }
        }
    }
    else if (Strand == '-') {
        for (unsigned ReadIndex = 0; ReadIndex < Index_Of_Useful_Reads.size(); ReadIndex++) {
            //std::cout << std::endl;
            for (short BaseIndex = 0; BaseIndex < input_reads[Index_Of_Useful_Reads[ReadIndex]].getReadLength(); BaseIndex++) {
                //std::cout << input_reads[Index_Of_Useful_Reads[ReadIndex]].UnmatchedSeq[BaseIndex];
                switch (input_reads[Index_Of_Useful_Reads[ReadIndex]].getUnmatchedSeq()[BaseIndex]) {
                    case 'A':
                        Count[0][input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - Min_LeftMostPos + BaseIndex]++;
                        break;	// 00000000
                    case 'C':
                        Count[1][input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - Min_LeftMostPos + BaseIndex]++;
                        break;	// 00010000
                    case 'G':
                        Count[2][input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - Min_LeftMostPos + BaseIndex]++;
                        break;	// 00100000
                    case 'T':
                        Count[3][input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - Min_LeftMostPos + BaseIndex]++;
                        break;	// 00110000
                    default:
                        Count[4][input_reads[Index_Of_Useful_Reads[ReadIndex]].LeftMostPos - Min_LeftMostPos + BaseIndex]++;
                        // 01000000
                }
            }
        }
    }

    float Sum;
    for (unsigned PosIndex = 0; PosIndex < Max_AssembledLength; PosIndex++) {
        Sum = Count[0][PosIndex] + Count[1][PosIndex] + Count[2][PosIndex] + Count[3][PosIndex] + Count[4][PosIndex];
        //std::cout << Count[0][PosIndex] << " " << Count[1][PosIndex] << " " << Count[2][PosIndex] << " " << Count[3][PosIndex] << " " << Count[4][PosIndex] << std::endl;
        for (unsigned BaseIndex = 0; BaseIndex < 5; BaseIndex++) {
            Ratio[BaseIndex][PosIndex] = Count[BaseIndex][PosIndex] / Sum;
        }
    }
    std::string OutputOneStr = "";
    const float RatioCutoff = 0.66;
    unsigned Max_Base_Count = 0;
    short Max_Base_Count_Index = -1;
    for (unsigned PosIndex = 0; PosIndex < Max_AssembledLength; PosIndex++) {
        if (Ratio[0][PosIndex] > RatioCutoff) {
            OutputOneStr += "A";
            continue;
        }
        if (Ratio[1][PosIndex] > RatioCutoff) {
            OutputOneStr += "C";
            continue;
        }
        if (Ratio[2][PosIndex] > RatioCutoff) {
            OutputOneStr += "G";
            continue;
        }
        if (Ratio[3][PosIndex] > RatioCutoff) {
            OutputOneStr += "T";
            continue;
        }
        for (short BaseIndex = 0; BaseIndex < 4; BaseIndex++) {
            if (Count[BaseIndex][PosIndex] > Max_Base_Count && Count[BaseIndex][PosIndex] >= 3)
                Max_Base_Count_Index = BaseIndex;
        }
        //if (Max_Base_Count_Index != -1) {
            switch (Max_Base_Count_Index) {
                case 0:
                    OutputOneStr += "A";
                    break;
                case 1:
                    OutputOneStr += "C";
                    break;
                case 2:
                    OutputOneStr += "G";
                    break;
                case 3:
                    OutputOneStr += "T";
                    break;
                case -1:
                    OutputOneStr += "N";
                    break;
                default:
                    break;
            }
        //}
        //else OutputOneStr += "N";
    }
    if (Strand == '+') {
        std::cout << "Final merged string +: original\n" << (OutputOneStr) << std::endl;
        std::cout << "Final merged string +: convert to ref\n" << ReverseComplement(OutputOneStr) << std::endl;
    }
    else {
        std::cout << "Final merged string -: original\n" << (OutputOneStr) << std::endl;
        std::cout << "Final merged string : convert to ref\n" << (OutputOneStr) << std::endl;
    }//std::cout << "Final merged string: -\n" << (OutputOneStr) << std::endl;
    //std::cout << "here1" << std::endl;
    unsigned Count_N = 0;
    for (unsigned pos_index = 0;  pos_index < OutputOneStr.size(); pos_index++) {
        if (OutputOneStr[pos_index] == 'N') Count_N++;
    }
    if ((float)Count_N < OutputOneStr.size() * 0.05) {
        
        unsigned Min_Close_Size = 10000;
        Index2Read4Copy = 0; // if the best one cannot be found due to N or whatever reasons, use the first read as the template for copy.
        //std::cout << "Original Index2Read4Copy: " << Index2Read4Copy << std::endl;
        for (unsigned ReadIndex = 0; ReadIndex < Index_Of_Useful_Reads.size(); ReadIndex++) {
            if (input_reads[Index_Of_Useful_Reads[ReadIndex]].UP_Close.size() < Min_Close_Size && OutputOneStr.find(input_reads[Index_Of_Useful_Reads[ReadIndex]].getUnmatchedSeq()) !=std::string::npos) { // quick fix here: need more work
                Min_Close_Size = input_reads[Index_Of_Useful_Reads[ReadIndex]].UP_Close.size();
                Index2Read4Copy = ReadIndex; // input_reads[Index_Of_Useful_Reads[Index2Read4Copy]]
                //std::cout << "Changed Index2Read4Copy: " << Index2Read4Copy << std::endl;
            }
        }
        //std::cout << "here2" << std::endl;
        output_one_read = input_reads[Index_Of_Useful_Reads[Index2Read4Copy]];
        //std::cout << "here2a" << std::endl;
        output_one_read.setUnmatchedSeq( OutputOneStr );
        //update std::map <std::string, int> ReadCountPerSample;
        GetReadCountPerSample(input_reads, Index_Of_Useful_Reads, output_one_read);
        //std::cout << "here2b" << std::endl;
        //std::cout << "Before: " << output_one_read.UP_Close.size() << std::endl;
        output_one_read.UP_Close.clear();
        //std::cout << "here3" << std::endl;
        output_one_read.Thickness = Index_Of_Useful_Reads.size();
        //std::cout << "Thickness " << output_one_read.Thickness << std::endl;
        GetCloseEnd(CurrentChrSeq, output_one_read);
        //std::cout << "After: " << output_one_read.UP_Close.size() << std::endl;
        output_reads.push_back(output_one_read);
        //std::cout << "here4" << std::endl;
        //std::cout << "end of CombineReads" << std::endl;
    }
}

void OutputCurrentRead(std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Assembly & OneSV, SPLIT_READ & OneRead, std::ofstream & ASM_Output) {

    //std::cout << "start of OutputCurrentRead" << std::endl;
    //if (OneRead.UP_Far_backup.size())
    //std::cout << "OutputCurrentRead start " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].LengthStr << std::endl;
    if (OneRead.UP_Far.size()) {
        //std::cout << "UP_Far start" << std::endl;
        ASM_Output << OneSV.Index + 1 << " " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " << OneSV.CI_A << "\t" << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV .CI_B 
        << "\tA " << OneRead.MatchedD << " " << OneRead.MatchedRelPos << " " << OneRead.Thickness << "\t" << OneRead.FragName  
        << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].Strand 
        << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].Direction 
        << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].LengthStr 
        << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].AbsLoc - g_SpacerBeforeAfter + 1 
        << " | " << OneRead.FarFragName // << " " << OneRead.MatchedD 
        << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].Strand 
        << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].Direction 
        << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].LengthStr 
        << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].AbsLoc - g_SpacerBeforeAfter + 1;
        int CountSupportingSamples = 0;
        for (std::map<std::string,int>::iterator it = OneRead.ReadCountPerSample.begin(); it != OneRead.ReadCountPerSample.end(); ++it )
            CountSupportingSamples++;
        ASM_Output << "\t#samples " << CountSupportingSamples << " "; 
        for (std::map<std::string,int>::iterator it = OneRead.ReadCountPerSample.begin(); it != OneRead.ReadCountPerSample.end(); ++it ) {
            ASM_Output << "\t" << (*it).first << " " << (*it).second;
        }
        ASM_Output << "\t NT_Size: 0\tNT_Str: \"\"" << std::endl;
        //std::cout << "UP_Far end" << std::endl;
    }   
}

void TryLI(std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Assembly & OneSV,  std::vector <SPLIT_READ> & First, std::vector <SPLIT_READ> & Second, std::ofstream & ASM_Output) {
    short MinimumOverlap = 10;
    short MaximumOverlap;// = min();
    short MaxMismatch = 3;
    short CountMismatch;
    short FirstLength, SecondLength;
    std::string FirstOne, SecondOne, MergedString;
    std::cout << "TryLI" << std::endl;
    for (unsigned ReadIndex = 0; ReadIndex < First.size(); ReadIndex++) {
        std::cout << First[ReadIndex].MatchedD << " " << First[ReadIndex].MatchedRelPos << std::endl;
    }
    for (unsigned ReadIndex = 0; ReadIndex < Second.size(); ReadIndex++) {
        std::cout << Second[ReadIndex].MatchedD << " " << Second[ReadIndex].MatchedRelPos << std::endl;
    }
    for (unsigned ReadIndex_Plus = 0; ReadIndex_Plus < First.size(); ReadIndex_Plus++) {
        if (First[ReadIndex_Plus].MatchedD == '-') continue;
        for (unsigned ReadIndex_Minus = 0; ReadIndex_Minus < Second.size(); ReadIndex_Minus++) {
            if (Second[ReadIndex_Minus].MatchedD == '+') continue;
            MaximumOverlap = std::min(First[ReadIndex_Plus].getReadLength(), Second[ReadIndex_Minus].getReadLength());
            std::cout << MaximumOverlap << std::endl;
            FirstOne = First[ReadIndex_Plus].getUnmatchedSeqRev();
            SecondOne = Second[ReadIndex_Minus].getUnmatchedSeq();
            FirstLength = FirstOne.size();
            SecondLength = SecondOne.size();
            std::cout << FirstOne << "\n" << SecondOne << "\n";
            for (short OverlapCount = MinimumOverlap; OverlapCount < MaximumOverlap; OverlapCount++) {
                CountMismatch = 0;
                for (short pos_index = 0; pos_index < OverlapCount; pos_index++) {
                    if (FirstOne[FirstLength - OverlapCount + pos_index] != SecondOne[pos_index]) {
                        ++CountMismatch;
                    }
                    if (CountMismatch > MaxMismatch) {
                        break;    
                    }
                } 
                if (CountMismatch <= MaxMismatch) {
                    MergedString = FirstOne.substr(0, FirstLength - OverlapCount) + SecondOne;
                    std::cout << "MergedString: " << OverlapCount << " " << MergedString << std::endl;
                    //ReportLI(AllChromosomes, ChrName2Index, CurrentState, OneSV, First[ReadIndex_Plus], Second[ReadIndex_Minus], MergedString, OverlapCount, ASM_Output);
                }
            }
        }
    }
}

/*void ReportLI(const std::vector <AChromosome> & AllChromosomes, std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Assembly & OneSV,  SPLIT_READ & First, SPLIT_READ & Second, const std::string & MergedString, short & OverlapCount, std::ofstream & ASM_Output){
    
    unsigned CountSupportingSamples = 0;
    SPLIT_READ OneRead = Second;
    OneRead.setUnmatchedSeq( MergedString );
    
    OneRead.UP_Close.clear();
    OneRead.UP_Far.clear();
    //OneRead.UP_Far_backup.clear();
    GetCloseEnd(AllChromosomes[ChrName2Index.find(OneSV.ChrB)->second].ChrSeq, OneRead); //        
    unsigned SearchRange = OneSV.CI_A + 1000;
    unsigned SearchCenter = OneSV.PosA + g_SpacerBeforeAfter;
    SearchFarEndAtPos(AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq, OneRead, SearchCenter, SearchRange);
    std::cout << OneRead.UP_Close.size() << " " << OneRead.UP_Far.size() << " " << "\n"; // OneRead.UP_Far_backup.size() << std::endl;
    if (OneRead.UP_Far_backup.size()) {
        if (OneRead.MatchedD == '+') {
            std::cout << "something is wrong: assembly::ReportLI, line 567" << std::endl;
            return;
        }
        else {
            short NT_Size = OneRead.getReadLength() - OneRead.UP_Close[OneRead.UP_Close.size() - 1].LengthStr - OneRead.UP_Far_backup[OneRead.UP_Far_backup.size() - 1].LengthStr;
            std::string NT_Str = OneRead.getUnmatchedSeq().substr(OneRead.UP_Far_backup[OneRead.UP_Far_backup.size() - 1].LengthStr, NT_Size);
            
            ASM_Output << OneSV.Index + 1 << " " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " << OneSV.CI_A << "\t" << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV .CI_B 
            << "\tLIA " << OneRead.MatchedD << " " << OneRead.MatchedRelPos << " " << OneRead.Thickness << "\t" << OneRead.FragName  
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].Strand 
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].Direction 
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].LengthStr 
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].AbsLoc - g_SpacerBeforeAfter + 1
            << " | " << OneRead.FarFragName // << " " << OneRead.MatchedD 
            << " " << OneRead.UP_Far_backup[OneRead.UP_Far_backup.size() - 1].Strand 
            << " " << OneRead.UP_Far_backup[OneRead.UP_Far_backup.size() - 1].Direction 
            << " " << OneRead.UP_Far_backup[OneRead.UP_Far_backup.size() - 1].LengthStr 
            << " " << OneRead.UP_Far_backup[OneRead.UP_Far_backup.size() - 1].AbsLoc - g_SpacerBeforeAfter + 1;
            CountSupportingSamples = 0;
            for (std::map<std::string,int>::iterator it = OneRead.ReadCountPerSample.begin(); it != OneRead.ReadCountPerSample.end(); it++ )
                CountSupportingSamples++;
            ASM_Output << "\t#samples " << CountSupportingSamples << " "; 
            for (std::map<std::string,int>::iterator it = OneRead.ReadCountPerSample.begin(); it != OneRead.ReadCountPerSample.end(); it++ ) {
                ASM_Output << "\t" << (*it).first << " " << (*it).second;
            }
            ASM_Output << "\t" << "NT_Size: " << NT_Size << "\tNT_Str: " << NT_Str
            << std::endl;
        }
    }
}*/

/*
void ReportLI(ControlState & CurrentState, const Assembly & OneSV,  SPLIT_READ & First, SPLIT_READ & Second, const std::string & MergedString, short & OverlapCount, std::ofstream & ASM_Output){
    
    unsigned CountSupportingSamples = 0;
    SPLIT_READ OneRead = Second;
    OneRead.setUnmatchedSeq( MergedString );
    
    OneRead.UP_Close.clear();
    OneRead.UP_Far.clear();
    //OneRead.UP_Far_backup.clear();
    GetCloseEnd((g_genome.getChr(OneSV.ChrB))->getSeq(), OneRead); //
    unsigned SearchRange = OneSV.CI_A + 1000;
    unsigned SearchCenter = OneSV.PosA + g_SpacerBeforeAfter;
	std::vector< SearchWindow > searchCluster;
	searchCluster.push_back( SearchWindow( AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrName, SearchCenter-SearchRange, SearchCenter+SearchRange ) );

   //SearchFarEndAtPos(AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq, OneRead, SearchCenter, SearchRange);
	SearchFarEndAtPos(AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq, OneRead, searchCluster);
    std::cout << OneRead.UP_Close.size() << " " << OneRead.UP_Far.size() << " " << "\n"; // OneRead.UP_Far.size() << std::endl;
    if (OneRead.UP_Far.size()) {
        if (OneRead.MatchedD == '+') {
            std::cout << "something is wrong: assembly::ReportLI, line 567" << std::endl;
            return;
        }
        else {
            short NT_Size = OneRead.getReadLength() - OneRead.UP_Close[OneRead.UP_Close.size() - 1].LengthStr - OneRead.UP_Far[OneRead.UP_Far.size() - 1].LengthStr;
            std::string NT_Str = OneRead.getUnmatchedSeq().substr(OneRead.UP_Far[OneRead.UP_Far.size() - 1].LengthStr, NT_Size);
            
            ASM_Output << OneSV.Index + 1 << " " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " << OneSV.CI_A << "\t" << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV .CI_B 
            << "\tLIA " << OneRead.MatchedD << " " << OneRead.MatchedRelPos << " " << OneRead.Thickness << "\t" << OneRead.FragName  
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].Strand 
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].Direction 
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].LengthStr 
            << " " << OneRead.UP_Close[OneRead.UP_Close.size() - 1].AbsLoc - g_SpacerBeforeAfter + 1
            << " | " << OneRead.FarFragName // << " " << OneRead.MatchedD 
            << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].Strand 
            << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].Direction 
            << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].LengthStr 
            << " " << OneRead.UP_Far[OneRead.UP_Far.size() - 1].AbsLoc - g_SpacerBeforeAfter + 1;
            CountSupportingSamples = 0;
            for (std::map<std::string,int>::iterator it = OneRead.ReadCountPerSample.begin(); it != OneRead.ReadCountPerSample.end(); it++ )
                CountSupportingSamples++;
            ASM_Output << "\t#samples " << CountSupportingSamples << " "; 
            for (std::map<std::string,int>::iterator it = OneRead.ReadCountPerSample.begin(); it != OneRead.ReadCountPerSample.end(); it++ ) {
                ASM_Output << "\t" << (*it).first << " " << (*it).second;
            }
            ASM_Output << "\t" << "NT_Size: " << NT_Size << "\tNT_Str: " << NT_Str
            << std::endl;
        }
    }
}
 */

void GetReadCountPerSample(const std::vector <SPLIT_READ> & input_reads, const std::vector <unsigned int> & Index_Of_Useful_Reads, SPLIT_READ & output_one_read) {
    //std::map <std::string, int> ReadCountPerSample;
    std::map<std::string,int>::iterator it;
    for (unsigned ReadIndex = 0; ReadIndex < Index_Of_Useful_Reads.size(); ReadIndex++) {
        it = output_one_read.ReadCountPerSample.find(input_reads[Index_Of_Useful_Reads[ReadIndex]].Tag);
        if (it == output_one_read.ReadCountPerSample.end()) {
            output_one_read.ReadCountPerSample.insert ( std::pair<std::string,int>(input_reads[Index_Of_Useful_Reads[ReadIndex]].Tag, 1) );
        }
        else {
            (*it).second++; 
        }
    }
}

void CleanUpCloseEnd(std::vector <SPLIT_READ> & input, const unsigned & Left, const unsigned & Right) {
    std::vector <SPLIT_READ> output;
     

    for (unsigned ReadIndex = 0; ReadIndex < input.size(); ReadIndex++) {
        if (input[ReadIndex].UP_Close.size()) {
            if ((input[ReadIndex].UP_Close[input[ReadIndex].UP_Close.size() - 1].AbsLoc >= Left && input[ReadIndex].UP_Close[input[ReadIndex].UP_Close.size() - 1].AbsLoc <= Right) || (input[ReadIndex].UP_Close[0].AbsLoc >= Left && input[ReadIndex].UP_Close[0].AbsLoc <= Right)) {
                output.push_back(input[ReadIndex]);
            }
        }
    }
    input.swap(output);
}

void CleanUpFarEnd(std::vector <SPLIT_READ> & input, const unsigned & Left, const unsigned & Right) {
    std::vector <SPLIT_READ> output;
    for (unsigned ReadIndex = 0; ReadIndex < input.size(); ReadIndex++) {
        if (input[ReadIndex].UP_Far.size()) {
           if ((input[ReadIndex].UP_Far[input[ReadIndex].UP_Far.size() - 1].AbsLoc >= Left && input[ReadIndex].UP_Far[input[ReadIndex].UP_Far.size() - 1].AbsLoc <= Right) || (input[ReadIndex].UP_Far[0].AbsLoc >= Left && input[ReadIndex].UP_Far[0].AbsLoc <= Right)) {
              output.push_back(input[ReadIndex]);
           }
        }
    }
    input.swap(output);
}
