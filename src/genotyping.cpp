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
#include "genotyping.h"
#include "bam2depth.h"
#include "assembly.h"
#include "reader.h"
#include "pindel.h"
#include "farend_searcher.h"
#include <map>
#include <set>
#include <string>
#include <utility>
#include <algorithm>
#include <math.h>


void doGenotyping (ControlState & CurrentState, ParCollection & par) {
    // step 1 load whole genome sequences into memory
    
    std::map<std::string,int> ChrName2Index;
    
    std::cout << "Get whole genome sequence..." << std::endl;
    // step 1. get the whole genome sequence
    
    std::vector <Chromosome> AllChromosomes;
    getWholeGenome(CurrentState, AllChromosomes);
    
    for (unsigned i = 0; i < AllChromosomes.size(); i++) {
        //std::cout << "ChrName " << AllChromosomes[i].ChrName << "\tChrSeqSize " << AllChromosomes[i].ChrSeq.size() << std::endl;
        ChrName2Index[AllChromosomes[i].ChrName] = i;
    }
    
    
    // step 2 load all variants into memory
    
    // step 2. get all SVs
    //CurrentState.inf_AssemblyInput.open(par.inf_AssemblyInputFilename.c_str());
    std::cout << "\nGet all SVs to genotype ..." << std::endl;
    std::vector <Genotyping> AllSV4Genotyping;
    Genotyping OneSV;
    unsigned SV_Count = 0;
    while (CurrentState.inf_GenotypingInput >> OneSV.Type
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
        //std::cout << "getting OneSV " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " 
        //          << OneSV.CI_A << " " << OneSV.ChrB << " " << OneSV.PosB << " " 
        //          << OneSV.CI_B << std::endl;
        AllSV4Genotyping.push_back(OneSV);
    }
    std::cout << "\nAllSV4Genotyping size " << AllSV4Genotyping.size() << "\n" << std::endl;
    // step 3 define output
    std::string GT_OutputFileName = CurrentState.OutputFolder + "_GT";
    std::ofstream GT_Output(GT_OutputFileName.c_str());

    // step 4 for each variant, do genotyping
    for (unsigned SV_index =0; SV_index < AllSV4Genotyping.size(); SV_index++) {
        // step 4.1 if type == DEL, GenotypeDel

        if (AllSV4Genotyping[SV_index].Type == "DEL") GenotypingOneDEL(AllChromosomes, ChrName2Index, CurrentState, par, AllSV4Genotyping[SV_index], GT_Output);

        // step 4.2 if type == DUP, GenotypeDup
        
        // step 4.3 if type == INV, GenotypeINV
        
        // step 4.4 if type == ITX, GenotypeINV
        
        // step 4.5 if type == CTX, GenotypeINV
    }
}

short GenotypingOneDEL(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> &ChrName2Index, ControlState & CurrentState, ParCollection & par, Genotyping & OneSV, std::ofstream & GT_Output) {
    std::cout << "Genotyping " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " 
    << OneSV.CI_A << " " << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV.CI_B << std::endl;
    // get RD signals
    const std::string & CurrentChrSeq = AllChromosomes[ ChrName2Index[ OneSV.ChrA ]].ChrSeq;
    CurrentState.CurrentChrName = OneSV.ChrA;
    getRelativeCoverage(CurrentChrSeq, ChrName2Index[OneSV.ChrA], CurrentState, OneSV);
    GetRP4OnDEL(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, GT_Output);
    //short AssembleOneSV(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, ParCollection & par, const Assembly & OneSV, std::ofstream & ASM_Output);
    //getRP_counts4DEL(CurrentChrSeq, ChrName2Index[OneSV.ChrA], CurrentState, OneSV);
    //std::cout << "after getRelativeCoverage " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " 
    //<< OneSV.CI_A << " " << OneSV.ChrB << " " << OneSV.PosB << " " 
    //<< OneSV.CI_B << std::endl;
    //CountRP();
    
    return 0;
}

void getAverageAndSTDE(const std::vector <unsigned> & Distances, unsigned & Average, unsigned & STDE) {
    float float_average = 0;
    unsigned Sum = 0;
    for (unsigned i = 0; i < Distances.size(); i++) Sum += Distances[i];
    float_average = (float)Sum / Distances.size();
    Average = (unsigned) float_average;
    float Diff = 0;
    for (unsigned i = 0; i < Distances.size(); i++) Diff += pow(Distances[i] - float_average, 2);
    STDE = (unsigned)(sqrt(Diff / Distances.size()));
    std::cout << Average << " " << STDE << std::endl;
}

void getMAD(const std::vector <unsigned> & Distances, const unsigned & Median, unsigned & MAD) {
    std::vector <unsigned> Diff;
    unsigned TempDiff;
    for (unsigned i = 0; i < Distances.size(); i++) {
        if (Distances[i] > Median) TempDiff = Distances[i] - Median;
        else TempDiff = Median - Distances[i];
        Diff.push_back(TempDiff);
    }
    sort(Diff.begin(), Diff.end());
    MAD = Diff[Diff.size() / 2];
    std::cout << "MAD: " << Median << " " << MAD << std::endl; 
}

void CountREF_RP(const std::vector <RP_READ> & Reads_RP, const std::vector <unsigned> & RP_READ_Index, unsigned lower, unsigned upper, unsigned Cutoff, unsigned & CountREF) {
    CountREF = 0;
    //std::cout << "lower and upper: " << lower << " " << upper << std::endl;
    for (unsigned i = 0; i < RP_READ_Index.size(); i++) {
        if ((unsigned)Reads_RP[RP_READ_Index[i]].Distance <= Cutoff) {
            //std::cout << "here CountREF_RP " << Reads_RP[RP_READ_Index[i]].Distance << " " << Reads_RP[RP_READ_Index[i]].PosA << " " << Reads_RP[RP_READ_Index[i]].PosB << std::endl;
            if (Reads_RP[RP_READ_Index[i]].PosA < Reads_RP[RP_READ_Index[i]].PosB) {
                if (Reads_RP[RP_READ_Index[i]].PosA <= lower && Reads_RP[RP_READ_Index[i]].PosB >= upper) {
                    CountREF++;
                }
            }
            else {
                if (Reads_RP[RP_READ_Index[i]].PosB <= lower && Reads_RP[RP_READ_Index[i]].PosA >= upper) {
                    CountREF++;
                }
            }
        }
    }
}

void CountALT_RP(const std::vector <RP_READ> & Reads_RP, const std::vector <unsigned> & RP_READ_Index, unsigned lower, unsigned upper, unsigned Cutoff, unsigned & CountALT) {
    CountALT = 0;
    for (unsigned i = 0; i < RP_READ_Index.size(); i++) {
        if ((unsigned)Reads_RP[RP_READ_Index[i]].Distance > Cutoff) {
            //std::cout << "here CountALT_RP " << Reads_RP[RP_READ_Index[i]].Distance << " " << Reads_RP[RP_READ_Index[i]].PosA << " " << Reads_RP[RP_READ_Index[i]].PosB << std::endl;
            if (Reads_RP[RP_READ_Index[i]].PosA < Reads_RP[RP_READ_Index[i]].PosB) {
                if (Reads_RP[RP_READ_Index[i]].PosA <= lower && Reads_RP[RP_READ_Index[i]].PosB >= upper) {
                    CountALT++;
                }
            }
            else {
                if (Reads_RP[RP_READ_Index[i]].PosB <= lower && Reads_RP[RP_READ_Index[i]].PosA >= upper) {
                    CountALT++;
                }
            }
        }
    }
}

void CountRPSupport4DEL(const std::vector <RP_READ> & Reads_RP, const std::vector <unsigned> RP_READ_Index, const Genotyping & OneSV, const unsigned Median, const unsigned MAD, unsigned & CountREF_A, unsigned & CountREF_B, unsigned & CountALT) {
    unsigned Cutoff = Median + 5 * MAD;
    std::cout << "Cutoff " << Cutoff << std::endl;
    CountREF_RP(Reads_RP, RP_READ_Index, OneSV.PosA - OneSV.CI_A, OneSV.PosA + OneSV.CI_A, Cutoff, CountREF_A);
    CountREF_RP(Reads_RP, RP_READ_Index, OneSV.PosB - OneSV.CI_B, OneSV.PosB + OneSV.CI_B, Cutoff, CountREF_B);
    CountALT_RP(Reads_RP, RP_READ_Index, OneSV.PosA - OneSV.CI_A, OneSV.PosB + OneSV.CI_B, Cutoff, CountALT);
    std::cout << "REF A: " << CountREF_A << "\tREF B: " << CountREF_B << "\t ALT: " << CountALT << std::endl;
    if (CountREF_A + CountREF_B + CountALT)
        std::cout << "Genotype from RP information " << (float)(CountREF_A + CountREF_B) * 2 / (CountREF_A + CountREF_B + CountALT * 2) << std::endl;
    else std::cout << "Genotype undefined from RP" << std::endl;
}

short GetRP4OnDEL(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, ParCollection & par, const Genotyping & OneSV, std::ofstream & GT_Output) {
    //std::vector <RP_READ> ALL_RP_reads;
    short Min_MQ = 20;
    std::set<std::string> ReadNames;
    
    if (CurrentState.CurrentChrName != OneSV.ChrA) {
        CurrentState.CurrentChrName = OneSV.ChrA;
        CurrentState.CurrentChrSeq = AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq; // change later, copying one chrseq for each SV is expensive. 
    }
    
    //unsigned SearchCenter;
    //unsigned SearchRange;
    
    CurrentState.Reads_RP.clear();
    unsigned Overhead = 1000;
    if (OneSV.PosA > OneSV.CI_A + Overhead)  
        CurrentState.lowerBinBorder = OneSV.PosA - OneSV.CI_A - Overhead; //CurrentState.
    else CurrentState.lowerBinBorder = 1;
    CurrentState.upperBinBorder = OneSV.PosB + OneSV.CI_B + Overhead;
    
    get_RP_Reads(CurrentState, par);
    //std::cout << "Reads around BP 1 " << CurrentState.Reads_RP.size() << std::endl;
    
    /*
    if (CurrentState.CurrentChrName != OneSV.ChrB) {
        CurrentState.CurrentChrName = OneSV.ChrB;
        CurrentState.CurrentChrSeq = AllChromosomes[ChrName2Index.find(OneSV.ChrB)->second].ChrSeq; // change later, copying one chrseq for each SV is expensive. 
    }
    
    if (OneSV.PosB > OneSV.CI_B + Overhead)  
        CurrentState.lowerBinBorder = OneSV.PosB - OneSV.CI_B - Overhead; //CurrentState.
    else CurrentState.lowerBinBorder = 1;
    CurrentState.upperBinBorder = OneSV.PosB + OneSV.CI_B + Overhead;
    
    get_RP_Reads(CurrentState, par);
    */
    std::cout << "Reads around both breakpoints " << CurrentState.Reads_RP.size() << std::endl;
    std::vector <unsigned> Distances, RP_READ_Index;
    unsigned TempDistance;
    for (unsigned ReadIndex = 0; ReadIndex < CurrentState.Reads_RP.size(); ReadIndex++) {
        if (CurrentState.Reads_RP[ReadIndex].ChrNameA == CurrentState.Reads_RP[ReadIndex].ChrNameB 
            && CurrentState.Reads_RP[ReadIndex].ChrNameA == OneSV.ChrA) {
            if (CurrentState.Reads_RP[ReadIndex].PosA == CurrentState.Reads_RP[ReadIndex].PosB) continue;
            if (CurrentState.Reads_RP[ReadIndex].MQA >= Min_MQ && CurrentState.Reads_RP[ReadIndex].MQB >= Min_MQ) {
                RP_READ_Index.push_back(ReadIndex);
                
                if (CurrentState.Reads_RP[ReadIndex].PosA > CurrentState.Reads_RP[ReadIndex].PosB) TempDistance = CurrentState.Reads_RP[ReadIndex].PosA - CurrentState.Reads_RP[ReadIndex].PosB;
                else TempDistance = CurrentState.Reads_RP[ReadIndex].PosB - CurrentState.Reads_RP[ReadIndex].PosA;
                CurrentState.Reads_RP[ReadIndex].Distance = TempDistance;
                //if (TempDistance > 10000)
                //std::cout << CurrentState.Reads_RP[ReadIndex].PosA << " " << CurrentState.Reads_RP[ReadIndex].PosB << " " << TempDistance << std::endl;
                Distances.push_back(TempDistance);
                
            }
        }
    }
    sort(Distances.begin(), Distances.end());
    //for (unsigned ReadIndex = 0; ReadIndex < Distances.size(); ReadIndex++)
    //    std::cout << Distances[ReadIndex] << " ";
    //std::cout << std::endl;
    unsigned Median = Distances[Distances.size() / 2]; 
    std::cout << "Insert size = " << Median << std::endl;
    unsigned Average, STDE;
    getAverageAndSTDE(Distances, Average, STDE);
    unsigned MAD;
    getMAD(Distances, Median, MAD);
    unsigned CountREF_A, CountREF_B, CountALT;
    CountRPSupport4DEL(CurrentState.Reads_RP, RP_READ_Index, OneSV, Median, MAD, CountREF_A, CountREF_B, CountALT);
    /*
    std::vector <RP_READ> ALL_RP_reads;
    unsigned SearchCenter;
    unsigned SearchRange;
    std::cout << "Current SV: " << OneSV.Index << " " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " << OneSV.CI_A 
    << "\t" << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV.CI_B << std::endl;
    // get first BP
    CurrentState.Reads.clear();
    //std::cout << "AssembleOneSV 2" << std::endl;
    
    CurrentState.CurrentChrName = OneSV.ChrA;
    CurrentState.CurrentChrSeq = AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq;
    CONS_Chr_Size = CurrentState.CurrentChrSeq.size() - 2 * g_SpacerBeforeAfter; // #################
    //std::cout << "CONS_Chr_Size " << CONS_Chr_Size << std::endl;
    g_maxPos = 0; // #################
    g_NumReadInWindow = 0; // #################
    g_InWinPlus = 0; // #################
    g_InWinMinus = 0; // #################
    g_CloseMappedPlus = 0; // #################
    g_CloseMappedMinus = 0; // #################
    unsigned Left, Right;
    if (OneSV.PosA > OneSV.CI_A + 1000)  
        CurrentState.lowerBinBorder = OneSV.PosA - OneSV.CI_A - 1000; //CurrentState.
    else CurrentState.lowerBinBorder = 1;
    CurrentState.upperBinBorder = OneSV.PosA + OneSV.CI_A + 1000;
    Left = OneSV.PosA + g_SpacerBeforeAfter - OneSV.CI_A;
    Right = OneSV.PosA + g_SpacerBeforeAfter + OneSV.CI_A;
    
    std::cout << "\nFirst BP\tChrName " << CurrentState.CurrentChrName << "\tRange " << CurrentState.lowerBinBorder << " " << CurrentState.upperBinBorder << std::endl;
    get_SR_Reads(CurrentState, par);
    
    //std::cout << "First size: " << CurrentState.Reads.size() << std::endl;
    CombineAndSort(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, First, CurrentState.lowerBinBorder, CurrentState.upperBinBorder, WhetherFirstBP);
    
    CleanUpCloseEnd(First, Left, Right); // vector of reads
    
    std::cout << "\nFirst size " << First.size() << std::endl;
    SearchRange = OneSV.CI_B + 1000;
    SearchCenter = OneSV.PosB + g_SpacerBeforeAfter;
    Left = OneSV.PosB + g_SpacerBeforeAfter - OneSV.CI_B;
    Right = OneSV.PosB + g_SpacerBeforeAfter + OneSV.CI_B;
    
    for (unsigned ReadIndex = 0; ReadIndex < First.size(); ReadIndex++) {
        First[ReadIndex].FarFragName = OneSV.ChrB;
        SearchFarEndAtPos(AllChromosomes[ChrName2Index.find(OneSV.ChrB)->second].ChrSeq, First[ReadIndex], SearchCenter, SearchRange);
    }
    //std::cout << "AssembleOneSV 7" << std::endl;
    
    CleanUpFarEnd(First, Left, Right);
    //std::cout << "AssembleOneSV 8" << std::endl;
    
    
    for (unsigned ReadIndex = 0; ReadIndex < First.size(); ReadIndex++) {
        if (First[ReadIndex].UP_Close.size()) {
            if (First[ReadIndex].UP_Far.size()) {
                //if (First[ReadIndex].UP_Far[0].LengthStr < 0) continue;
                //std::cout << "First UP_Far: ";
                //std::cout << First[ReadIndex].UP_Far.size() << std::endl;
                //std::cout << "First[ReadIndex].UP_Far.size() " << First[ReadIndex].UP_Close[First[ReadIndex].UP_Close.size() - 1].LengthStr << std::endl;
                if (First[ReadIndex].UP_Far[First[ReadIndex].UP_Far.size() - 1].LengthStr + First[ReadIndex].UP_Close[First[ReadIndex].UP_Close.size() - 1].LengthStr + Max_NT_Size >= First[ReadIndex].ReadLength) OutputCurrentRead(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, First[ReadIndex], ASM_Output);
            }
            else if (First[ReadIndex].UP_Far_backup.size()) {
                //std::cout << "First UP_Far_backup ";
                //std::cout << First[ReadIndex].UP_Far_backup.size() << std::endl; //First[ReadIndex].UP_Far_backup[First[ReadIndex].UP_Far_backup.size() - 1].LengthStr << " " << First[ReadIndex].UP_Close[First[ReadIndex].UP_Close.size() - 1].LengthStr << " " << Max_NT_Size << " " << First[ReadIndex].ReadLength << std::endl;
                //if (First[ReadIndex].UP_Far_backup[First[ReadIndex].UP_Far_backup.size() - 1].LengthStr + First[ReadIndex].UP_Close[First[ReadIndex].UP_Close.size() - 1].LengthStr + Max_NT_Size >= First[ReadIndex].ReadLength) 
                //std::cout << "First[ReadIndex].UP_Far_backup.size() # " << First[ReadIndex].UP_Close[First[ReadIndex].UP_Close.size() - 1].LengthStr << std::endl;
                OutputCurrentRead(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, First[ReadIndex], ASM_Output);
            }
        }
    }
    
    
    // get second BP
    CurrentState.Reads.clear();
    WhetherFirstBP = false;
    CurrentState.CurrentChrName = OneSV.ChrB;
    CurrentState.CurrentChrSeq = AllChromosomes[ChrName2Index.find(OneSV.ChrB)->second].ChrSeq;
    CONS_Chr_Size = CurrentState.CurrentChrSeq.size() - 2 * g_SpacerBeforeAfter; // #################
    g_maxPos = 0; // #################
    g_NumReadInWindow = 0; // #################
    g_InWinPlus = 0; // #################
    g_InWinMinus = 0; // #################
    g_CloseMappedPlus = 0; // #################
    g_CloseMappedMinus = 0; // #################
    if (OneSV.PosB > OneSV.CI_B + 1000)  
        CurrentState.lowerBinBorder = OneSV.PosB - OneSV.CI_B - 1000;
    else CurrentState.lowerBinBorder = 1;
    CurrentState.upperBinBorder = OneSV.PosB + OneSV.CI_B + 1000;
    Left = OneSV.PosB + g_SpacerBeforeAfter - OneSV.CI_B;
    Right = OneSV.PosB + g_SpacerBeforeAfter + OneSV.CI_B;
    
    std::cout << "\nSecond BP\tChrName " << CurrentState.CurrentChrName << "\tRange " << CurrentState.lowerBinBorder << " " << CurrentState.upperBinBorder << std::endl;    
    get_SR_Reads(CurrentState, par);
    
    CombineAndSort(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, Second, CurrentState.lowerBinBorder, CurrentState.upperBinBorder, WhetherFirstBP);
    
    CleanUpCloseEnd(Second, Left, Right);
    
    std::cout << "\nSecond size " << Second.size() << std::endl;
    SearchRange = OneSV.CI_A + 1000;
    SearchCenter = OneSV.PosA + g_SpacerBeforeAfter;
    Left = OneSV.PosA + g_SpacerBeforeAfter - OneSV.CI_A;
    Right = OneSV.PosA + g_SpacerBeforeAfter + OneSV.CI_A;
    
    for (unsigned ReadIndex = 0; ReadIndex < Second.size(); ReadIndex++) {
        Second[ReadIndex].FarFragName = OneSV.ChrA;
        SearchFarEndAtPos(AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq, Second[ReadIndex], SearchCenter, SearchRange);
    }
    
    CleanUpFarEnd(Second, Left, Right);
    
    for (unsigned ReadIndex = 0; ReadIndex < Second.size(); ReadIndex++) {
        if (Second[ReadIndex].UP_Close.size()) {
            if (Second[ReadIndex].UP_Far.size()) {
                //std::cout << "Second UP_Far: " << Second[ReadIndex].UP_Far.size() << std::endl;
                //std::cout << "Second[ReadIndex].UP_Far.size() " << Second[ReadIndex].UP_Close[Second[ReadIndex].UP_Close.size() - 1].LengthStr << std::endl;
                if (Second[ReadIndex].UP_Far[Second[ReadIndex].UP_Far.size() - 1].LengthStr + Second[ReadIndex].UP_Close[Second[ReadIndex].UP_Close.size() - 1].LengthStr + Max_NT_Size >= Second[ReadIndex].ReadLength) OutputCurrentRead(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, Second[ReadIndex], ASM_Output);
            }
            else if (Second[ReadIndex].UP_Far_backup.size()) {
                //std::cout << "Second UP_Far_backup " << Second[ReadIndex].UP_Far_backup.size() << std::endl; //<< Second[ReadIndex].UP_Far_backup[Second[ReadIndex].UP_Far_backup.size() - 1].LengthStr << " " << Second[ReadIndex].UP_Close[Second[ReadIndex].UP_Close.size() - 1].LengthStr << " " << Max_NT_Size << " " << Second[ReadIndex].ReadLength << std::endl;
                //if (Second[ReadIndex].UP_Far_backup[Second[ReadIndex].UP_Far_backup.size() - 1].LengthStr + Second[ReadIndex].UP_Close[Second[ReadIndex].UP_Close.size() - 1].LengthStr + Max_NT_Size >= Second[ReadIndex].ReadLength) 
                //std::cout << "Second[ReadIndex].UP_Far_backup.size() " << Second[ReadIndex].UP_Close[Second[ReadIndex].UP_Close.size() - 1].LengthStr << std::endl;
                OutputCurrentRead(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, Second[ReadIndex], ASM_Output);
            }
        }
    }
    //std::cout << "AssembleOneSV 15" << std::endl;
    
    unsigned SumSize = 0;
    for (unsigned ReadIndex = 0; ReadIndex < First.size(); ReadIndex++) {
        SumSize += First[ReadIndex].UP_Far.size() + First[ReadIndex].UP_Far_backup.size();
    }
    for (unsigned ReadIndex = 0; ReadIndex < Second.size(); ReadIndex++) {
        SumSize += Second[ReadIndex].UP_Far.size() + Second[ReadIndex].UP_Far_backup.size();
    }
    if (SumSize == 0 && OneSV.ChrA == OneSV.ChrB) {
        TryLI(AllChromosomes, ChrName2Index, CurrentState, par, OneSV, First, Second, ASM_Output);
    }

    */
    return 0;
}

