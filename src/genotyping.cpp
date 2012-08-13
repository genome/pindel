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
//#include <pair.h>
#include <string>
#include <utility>
#include <algorithm>
#include <math.h>


void doGenotyping (ControlState & CurrentState, std::ifstream& FastaFile ) {
    const int SV_Genotype_Cutoff = 1000;
    
    // step 1 load whole genome sequences into memory
    
    std::map<std::string,int> ChrName2Index;
    
    std::cout << "Get whole genome sequence..." << std::endl;
    // step 1. get the whole genome sequence
    
    std::vector <Chromosome> AllChromosomes;
    getWholeGenome(FastaFile, AllChromosomes);
    
    for (unsigned i = 0; i < AllChromosomes.size(); i++) {
        //std::cout << "ChrName " << AllChromosomes[i].ChrName << "\tChrSeqSize " << AllChromosomes[i].ChrSeq.size() << std::endl;
        ChrName2Index[AllChromosomes[i].ChrName] = i;
    }
    
    std::set<std::string> SampleNameAsSet;
    std::map<std::string, unsigned> SampleName2IndexAsMap;
    std::vector<std::string> SampleNameAsVector;
    for (unsigned BamIndex = 0; BamIndex < CurrentState.bams_to_parse.size(); BamIndex++) {
        if (SampleNameAsSet.find(CurrentState.bams_to_parse[BamIndex].Tag) == SampleNameAsSet.end()) { // not in the set
            SampleName2IndexAsMap.insert ( std::pair<std::string,int>(CurrentState.bams_to_parse[BamIndex].Tag, SampleNameAsSet.size() ) );
            SampleNameAsSet.insert(CurrentState.bams_to_parse[BamIndex].Tag);
            SampleNameAsVector.push_back(CurrentState.bams_to_parse[BamIndex].Tag);
        }
        else {
            std::cout << "Two BAM files with the same sample name.\n";
            exit(EXIT_FAILURE);
        }
    }
    std::cout << "There are " << SampleNameAsVector.size() << " samples.\n";
    std::cout << "Samples:";
    for (unsigned SampleIndex = 0; SampleIndex < SampleNameAsVector.size(); SampleIndex++) {
        std::cout << " " << SampleNameAsVector[SampleIndex];
    }
    std::cout << std::endl;
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
    
    std::cout << "Samples:";
    for (unsigned SampleIndex = 0; SampleIndex < SampleNameAsVector.size(); SampleIndex++) {
        std::cout << " " << SampleNameAsVector[SampleIndex];
    }
    std::cout << "\n\n";
    
    // step 3 define output
   // std::string GT_OutputFileName = CurrentState.OutputFolder + "_GT";
    std::ofstream GT_Output(UserDefinedSettings::Instance()->getGTOutputFilename().c_str());

    // step 4 for each variant, do genotyping
    for (unsigned SV_index =0; SV_index < AllSV4Genotyping.size(); SV_index++) {
        // step 4.1 if type == DEL, GenotypeDel
        if (AllSV4Genotyping[SV_index].ChrA == AllSV4Genotyping[SV_index].ChrB && abs(AllSV4Genotyping[SV_index].PosA - AllSV4Genotyping[SV_index].PosB) < SV_Genotype_Cutoff) {
            std::cout << "Skip One SV " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " 
                      << OneSV.CI_A << " " << OneSV.ChrB << " " << OneSV.PosB << " " 
                      << OneSV.CI_B << std::endl;
        }
        if (AllSV4Genotyping[SV_index].Type == "DEL") GenotypingOneDEL(AllChromosomes, ChrName2Index, CurrentState, AllSV4Genotyping[SV_index], SampleName2IndexAsMap, GT_Output);

        // step 4.2 if type == DUP, GenotypeDup
        if (AllSV4Genotyping[SV_index].Type == "DUP" || AllSV4Genotyping[SV_index].Type == "TD" || AllSV4Genotyping[SV_index].Type == "GT") GenotypingOneDUP(AllChromosomes, ChrName2Index, CurrentState, AllSV4Genotyping[SV_index], SampleName2IndexAsMap, GT_Output);
        // step 4.3 if type == INV, GenotypeINV
        
        // step 4.4 if type == ITX, GenotypeINV
        
        // step 4.5 if type == CTX, GenotypeINV
    }
}

short GenotypingOneDEL(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> &ChrName2Index, ControlState & CurrentState, Genotyping & OneSV, std::map<std::string, unsigned> & SampleName2IndexAsMap, std::ofstream & GT_Output) {
    std::cout << "\nGenotyping " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " 
    << OneSV.CI_A << " " << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV.CI_B << std::endl;
    // get RD signals
    const std::string & CurrentChrSeq = AllChromosomes[ ChrName2Index[ OneSV.ChrA ]].ChrSeq;
    CurrentState.CurrentChrName = OneSV.ChrA;
    //std::cout << "1" << std::endl;
    getRelativeCoverage(CurrentChrSeq, ChrName2Index[OneSV.ChrA], CurrentState, OneSV);
    //std::cout << "2" << std::endl;
    GetRP4OnDEL(AllChromosomes, ChrName2Index, CurrentState, OneSV, SampleName2IndexAsMap, GT_Output);
    //std::cout << "3" << std::endl;
    //short AssembleOneSV(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, ParCollection & par, const Assembly & OneSV, std::ofstream & ASM_Output);
    //getRP_counts4DEL(CurrentChrSeq, ChrName2Index[OneSV.ChrA], CurrentState, OneSV);
    //std::cout << "after getRelativeCoverage " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " " 
    //<< OneSV.CI_A << " " << OneSV.ChrB << " " << OneSV.PosB << " " 
    //<< OneSV.CI_B << std::endl;
    //CountRP();
    
    return 0;
}

short GenotypingOneDUP(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> &ChrName2Index, ControlState & CurrentState, Genotyping & OneSV, std::map<std::string, unsigned> & SampleName2IndexAsMap, std::ofstream & GT_Output) {
    std::cout << "\nGenotyping " << OneSV.Type << " " << OneSV.ChrA << " " << OneSV.PosA << " "
    << OneSV.CI_A << " " << OneSV.ChrB << " " << OneSV.PosB << " " << OneSV.CI_B << std::endl;
    const std::string & CurrentChrSeq = AllChromosomes[ ChrName2Index[ OneSV.ChrA ]].ChrSeq;
    CurrentState.CurrentChrName = OneSV.ChrA;
    getRelativeCoverage(CurrentChrSeq, ChrName2Index[OneSV.ChrA], CurrentState, OneSV);
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
    //std::cout << Average << " " << STDE << std::endl;
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
    //std::cout << "MAD: " << Median << " " << MAD << std::endl; 
}

void CountREF_RP_DEL(const std::vector <RPVector> & Reads_RP, const std::vector <std::vector <unsigned> > & RP_READ_Index, unsigned lower, unsigned upper, unsigned * Cutoff, unsigned * CountREF, std::map<std::string, unsigned> & SampleName2IndexAsMap) {
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++)
       CountREF[SampleIndex] = 0;
    //std::cout << "lower and upper: " << lower << " " << upper << std::endl;
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++) {
        if (Cutoff[SampleIndex] == 0) continue;
        for (unsigned i = 0; i < RP_READ_Index[SampleIndex].size(); i++) {
            if ((unsigned)Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].Distance <= Cutoff[SampleIndex]) {
                //std::cout << "here CountREF_RP " << Reads_RP[RP_READ_Index[i]].Distance << " " << Reads_RP[RP_READ_Index[i]].PosA << " " << Reads_RP[RP_READ_Index[i]].PosB << std::endl;
                if (Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosA < Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosB) {
                    if (Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosA <= lower && Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosB >= upper) {
                        CountREF[SampleIndex]++;
                    }
                }
                else {
                    if (Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosB <= lower && Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosA >= upper) {
                        CountREF[SampleIndex]++;
                    }
                }
            }
        }
    }
}

void CountALT_RP_DEL(const std::vector <RPVector> & Reads_RP, const std::vector <std::vector <unsigned> > & RP_READ_Index, unsigned lower, unsigned upper, unsigned * Cutoff, unsigned * CountALT, std::map<std::string, unsigned> & SampleName2IndexAsMap) {
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++)
        CountALT[SampleIndex] = 0;
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++) {
        if (Cutoff[SampleIndex] == 0) continue;
        for (unsigned i = 0; i < RP_READ_Index[SampleIndex].size(); i++) {
            if ((unsigned)Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].Distance > Cutoff[SampleIndex]) {
                //std::cout << "here CountREF_RP " << Reads_RP[RP_READ_Index[i]].Distance << " " << Reads_RP[RP_READ_Index[i]].PosA << " " << Reads_RP[RP_READ_Index[i]].PosB << std::endl;
                if (Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosA < Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosB) {
                    if (Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosA <= lower && Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosB >= upper) {
                        CountALT[SampleIndex]++;
                    }
                }
                else {
                    if (Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosB <= lower && Reads_RP[SampleIndex][RP_READ_Index[SampleIndex][i]].PosA >= upper) {
                        CountALT[SampleIndex]++;
                    }
                }
            }
        }
    }
}

void CountRPSupport4DEL(const std::vector <RPVector> & Reads_RP, const std::vector <std::vector <unsigned> >  RP_READ_Index, const Genotyping & OneSV, const unsigned * Median, const unsigned * MAD, std::map<std::string, unsigned> & SampleName2IndexAsMap) {
    //unsigned Cutoff = Median + 5 * MAD;
    //std::cout << "entering CountRPSupport4DEL ..." << std::endl;
    unsigned cutoff[SampleName2IndexAsMap.size()];
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++) {
        cutoff[SampleIndex] = Median[SampleIndex] + 5 * MAD[SampleIndex];
    }
    unsigned CountREF_A[SampleName2IndexAsMap.size()], CountREF_B[SampleName2IndexAsMap.size()], CountALT[SampleName2IndexAsMap.size()];
    CountREF_RP_DEL(Reads_RP, RP_READ_Index, OneSV.PosA - OneSV.CI_A, OneSV.PosA + OneSV.CI_A, cutoff, CountREF_A, SampleName2IndexAsMap);
    CountREF_RP_DEL(Reads_RP, RP_READ_Index, OneSV.PosB - OneSV.CI_B, OneSV.PosB + OneSV.CI_B, cutoff, CountREF_B, SampleName2IndexAsMap);
    CountALT_RP_DEL(Reads_RP, RP_READ_Index, OneSV.PosA - OneSV.CI_A, OneSV.PosB + OneSV.CI_B, cutoff, CountALT, SampleName2IndexAsMap);
    //std::cout << "REF A: " << CountREF_A << "\tREF B: " << CountREF_B << "\t ALT: " << CountALT << std::endl;
    std::cout << "Genotype_Based_On_RP:";
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++) {
        if (CountREF_A[SampleIndex] + CountREF_B[SampleIndex] + CountALT[SampleIndex])
            std::cout << " " << (float)(CountREF_A[SampleIndex] + CountREF_B[SampleIndex]) * 2 / (CountREF_A[SampleIndex] + CountREF_B[SampleIndex] + CountALT[SampleIndex] * 2);
        else std::cout << " -1";// << std::endl;
    }
    std::cout << std::endl;
    //std::cout << "leaving CountRPSupport4DEL ..." << std::endl;
}

short GetRP4OnDEL(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Genotyping & OneSV, std::map<std::string, unsigned> & SampleName2IndexAsMap, std::ofstream & GT_Output) {
    //std::vector <RP_READ> ALL_RP_reads;
    //std::cout << "GetRP4OnDEL 1" << std::endl;
    short Min_MQ = 20;
    std::set<std::string> ReadNames;
    
    if (CurrentState.CurrentChrName != OneSV.ChrA) {
        CurrentState.CurrentChrName = OneSV.ChrA;
        CurrentState.CurrentChrSeq = AllChromosomes[ChrName2Index.find(OneSV.ChrA)->second].ChrSeq; // change later, copying one chrseq for each SV is expensive. 
    }
    //std::cout << "GetRP4OnDEL 2" << std::endl;
    //unsigned SearchCenter;
    //unsigned SearchRange;
    
    CurrentState.Reads_RP.clear();
    unsigned Overhead = 1000;
	unsigned int lowerBinBorder = 1;
    if (OneSV.PosA > OneSV.CI_A + Overhead)  
        lowerBinBorder = OneSV.PosA - OneSV.CI_A - Overhead; //CurrentState.
   unsigned int upperBinBorder = OneSV.PosB + OneSV.CI_B + Overhead;
	SearchWindow window(CurrentState.CurrentChrName,lowerBinBorder, upperBinBorder); 
    //std::cout << "GetRP4OnDEL 3" << std::endl;
    get_RP_Reads(CurrentState, window );
    //std::cout << "Reads around BP 1 " << CurrentState.Reads_RP.size() << std::endl;
    //std::cout << "GetRP4OnDEL 4" << std::endl;
    //std::cout << "Reads around both breakpoints " << CurrentState.Reads_RP.size() << std::endl;
    typedef std::vector <unsigned> VectorOfUnsigned;
    std::vector <VectorOfUnsigned> Distances, RP_READ_Index;
    VectorOfUnsigned TempDistances, TempRP_READ_Index;
    //std::cout << "GetRP4OnDEL 5" << std::endl;
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++) {
        Distances.push_back(TempDistances);
        RP_READ_Index.push_back(TempRP_READ_Index);
    }
    unsigned TempDistance;
    //std::cout << "GetRP4OnDEL 6" << std::endl;
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++) {
        for (unsigned ReadIndex = 0; ReadIndex < CurrentState.Reads_RP[SampleIndex].size(); ReadIndex++) {
            if (CurrentState.Reads_RP[SampleIndex][ReadIndex].ChrNameA == CurrentState.Reads_RP[SampleIndex][ReadIndex].ChrNameB 
                && CurrentState.Reads_RP[SampleIndex][ReadIndex].ChrNameA == OneSV.ChrA) {
                if (CurrentState.Reads_RP[SampleIndex][ReadIndex].PosA == CurrentState.Reads_RP[SampleIndex][ReadIndex].PosB) continue;
                if (CurrentState.Reads_RP[SampleIndex][ReadIndex].MQA >= Min_MQ && CurrentState.Reads_RP[SampleIndex][ReadIndex].MQB >= Min_MQ) {
                    RP_READ_Index[SampleIndex].push_back(ReadIndex);
                    
                    if (CurrentState.Reads_RP[SampleIndex][ReadIndex].PosA > CurrentState.Reads_RP[SampleIndex][ReadIndex].PosB) TempDistance = CurrentState.Reads_RP[SampleIndex][ReadIndex].PosA - CurrentState.Reads_RP[SampleIndex][ReadIndex].PosB;
                    else TempDistance = CurrentState.Reads_RP[SampleIndex][ReadIndex].PosB - CurrentState.Reads_RP[SampleIndex][ReadIndex].PosA;
                    CurrentState.Reads_RP[SampleIndex][ReadIndex].Distance = TempDistance;
                    //if (TempDistance > 10000)
                    //std::cout << CurrentState.Reads_RP[ReadIndex].PosA << " " << CurrentState.Reads_RP[ReadIndex].PosB << " " << TempDistance << std::endl;
                    Distances[SampleIndex].push_back(TempDistance);
                    
                }
            }
        }
    }
    //std::cout << "GetRP4OnDEL 7" << std::endl;

    unsigned Median[SampleName2IndexAsMap.size()], Average[SampleName2IndexAsMap.size()], STDE[SampleName2IndexAsMap.size()], MAD[SampleName2IndexAsMap.size()];
    for (unsigned SampleIndex = 0; SampleIndex < SampleName2IndexAsMap.size(); SampleIndex++) {
        if (Distances[SampleIndex].size()) {
            sort(Distances[SampleIndex].begin(), Distances[SampleIndex].end());
            Median[SampleIndex] = Distances[SampleIndex][Distances[SampleIndex].size() / 2]; 
            getAverageAndSTDE(Distances[SampleIndex], Average[SampleIndex], STDE[SampleIndex]);
            getMAD(Distances[SampleIndex], Median[SampleIndex], MAD[SampleIndex]);
        }
        else {
            Median[SampleIndex] = 0;
            Average[SampleIndex] = 0;
            STDE[SampleIndex] = 0;
            MAD[SampleIndex] = 0;
        }
    }
    //std::cout << "GetRP4OnDEL 8" << std::endl;
    CountRPSupport4DEL(CurrentState.Reads_RP, RP_READ_Index, OneSV, Median, MAD, SampleName2IndexAsMap);
    //std::cout << "GetRP4OnDEL 9" << std::endl;
    return 0;
}

