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
#include <map>
#include <string>

void doAssembly (ControlState & CurrentState, ParCollection & par) {
    
    std::map<std::string,int> ChrName2Index;
    
    
    std::cout << "In assembly mode ..." << std::endl;
    // step 1. get the whole genome sequence
    //currentState.inf_Seq.open(par.referenceFileName.c_str());
    std::vector <Chromosome> AllChromosomes;
    getWholeGenome(CurrentState, AllChromosomes);
    CurrentState.inf_Seq.close();
    for (unsigned i = 0; i < AllChromosomes.size(); i++) {
        std::cout << AllChromosomes[i].ChrName << "\t" << AllChromosomes[i].ChrSeq.size() << std::endl;
        ChrName2Index[AllChromosomes[i].ChrName] = i;
    }


    // step 2. get all SVs
    //CurrentState.inf_AssemblyInput.open(par.inf_AssemblyInputFilename.c_str());
    std::vector <Assembly> AllSV4Assembly;
    Assembly OneSV;
    while (CurrentState.inf_AssemblyInput >> OneSV.Type
                                          >> OneSV.ChrA
                                          >> OneSV.PosA
                                          >> OneSV.CI_A
                                          >> OneSV.ChrB
                                          >> OneSV.PosB
                                          >> OneSV.CI_B)
        AllSV4Assembly.push_back(OneSV);
    
    // step 3. per SV, find BP, sort and output
    //CurrentState.OutputFolder = par.outputFileName;
    // step 3.1 define output
    std::string ASM_OutputFileName = CurrentState.OutputFolder + "_ASM";
    std::ofstream ASM_Output(ASM_OutputFileName.c_str());
    // step 3.2 per SV, collect reads around BP, sort, combine and output.
    for (unsigned SV_index = 0; SV_index < AllSV4Assembly.size(); SV_index++) {
        std::cout << AllSV4Assembly[SV_index].Type << "\t" << AllSV4Assembly[SV_index].ChrA << "\t" << AllSV4Assembly[SV_index].PosA << "\t" << AllSV4Assembly[SV_index].CI_A << "\t" << AllSV4Assembly[SV_index].ChrB << "\t" << AllSV4Assembly[SV_index].PosB << "\t" << AllSV4Assembly[SV_index].CI_B << std::endl;
        AssembleOneSV(AllChromosomes, ChrName2Index, CurrentState, par, AllSV4Assembly[SV_index], ASM_Output);
    }
    
    return;
}

short getWholeGenome(ControlState & CurrentState, std::vector <Chromosome> & AllChromosomes) {
    short Diff2UpperCase = 'A' - 'a';
    Chromosome OneChr;
    std::string TempLine;
    char TempChar;
    while (CurrentState.inf_Seq >> OneChr.ChrName) {
        getline(CurrentState.inf_Seq, TempLine);
        while (CurrentState.inf_Seq >> TempChar) {
            if (TempChar != '\n' && TempChar != '\r') {
                if (TempChar == '>') {
                    AllChromosomes.push_back(OneChr);
                    OneChr.ChrSeq = "";
                    OneChr.ChrName = "";
                    break;
                }
                else {
                    if ('a' <= TempChar && TempChar <= 'z') {
                        TempChar = TempChar + Diff2UpperCase;
                    }
                    switch (TempChar) {
                        case 'A':
                            OneChr.ChrSeq += 'A';
                            break;	// 00000000
                        case 'C':
                            OneChr.ChrSeq += 'C';
                            break;	// 00010000
                        case 'G':
                            OneChr.ChrSeq += 'G';
                            break;	// 00100000
                        case 'T':
                            OneChr.ChrSeq += 'T';
                            break;	// 00110000
                        default:
                            OneChr.ChrSeq += 'N';
                            // 01000000
                    }
                }						
            }
        }
    }
    AllChromosomes.push_back(OneChr);
    return 0;
}

short AssembleOneSV(const std::vector <Chromosome> & AllChromosomes, std::map<std::string,int> &ChrName2Index, ControlState & CurrentState, ParCollection & par, const Assembly & OneSV, std::ofstream & ASM_Output) {
    // get first BP
    CurrentState.Reads.clear();
    /*
    currentState.CurrentChrName = OneSV.ChrA;
    currentState.CurrentChrSeq = OneSV.ChrA;
    currentState.lowerBinBorder = ;
    currentState.upperBinBorder;
     */
    getReads(CurrentState, par);
    
    // get second BP
    CurrentState.Reads.clear();
    /*
    currentState.lowerBinBorder;
    currentState.upperBinBorder;
     */
    getReads(CurrentState, par);
    return 0;
}

