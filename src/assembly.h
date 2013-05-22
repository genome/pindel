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

#include "control_state.h"
#include "user_defined_settings.h"
#include "pindel.h"
#include <vector>
#include <map>
#include <string>

#ifndef ASSEMBLY_H
#define ASSEMBLY_H

/*
struct AChromosome {
    AChromosome() {
        ChrName = "";
        ChrSeq = ""; 
    }
    std::string ChrName;
    std::string ChrSeq;
};
 */

struct Assembly {
    Assembly() {
        Index = 0;
        Type = "";
        ChrA = "";
        PosA = 0;
        CI_A = 0;
        ChrB = "";
        PosB = 0;
        CI_B = 0;        
    }
    unsigned Index;
    std::string Type;
    std::string ChrA;
    unsigned PosA;
    unsigned CI_A;
    std::string ChrB;
    unsigned PosB;
    unsigned CI_B;
};

void doAssembly (ControlState & CurrentState, UserDefinedSettings* userSettings);//std::ifstream& FastaFile );
//short getWholeGenome(std::ifstream& FastaFile, std::vector <AChromosome> & AllChromosomes);
short AssembleOneSV(std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Assembly & OneSV, std::ofstream & ASM_Output);
void CombineAndSort(std::map<std::string,int> &ChrName2Index, ControlState & CurrentState, const Assembly & OneSV, std::vector <SPLIT_READ> & output_reads, const unsigned & lowerBinBorder, const unsigned & upperBinBorder, const bool & First);
void CombineReads(const std::string & CurrentChrSeq, const char & Strand, const std::vector <SPLIT_READ> & input_reads, const std::vector <unsigned int> & Index_Of_Useful_Reads, std::vector <SPLIT_READ> & output_reads);
void OutputCurrentRead(std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Assembly & OneSV, SPLIT_READ & OneRead, std::ofstream & ASM_Output);
void TryLI(std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Assembly & OneSV,  std::vector <SPLIT_READ> & First, std::vector <SPLIT_READ> & Second, std::ofstream & ASM_Output);
void ReportLI(std::map<std::string,int> & ChrName2Index, ControlState & CurrentState, const Assembly & OneSV,  SPLIT_READ & First, SPLIT_READ & Second, const std::string & MergedString, short & OverlapCount, std::ofstream & ASM_Output);
void GetReadCountPerSample(const std::vector <SPLIT_READ> & input_reads, const std::vector <unsigned int> & Index_Of_Useful_Reads, SPLIT_READ & output_one_read);
void CleanUpCloseEnd(std::vector <SPLIT_READ> & input, const unsigned & Left, const unsigned & Right);
void CleanUpFarEnd(std::vector <SPLIT_READ> & input, const unsigned & Left, const unsigned & Right);

#endif /* ASSEMBLY_H */
