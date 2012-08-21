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
#include "do_rd_rp.h"
#include "reader.h"
#include "pindel.h"
#include "farend_searcher.h"
#include <map>
#include <string>
#include <utility>

/*
 
 struct BAM_Path_IS {
 BAM_Path_IS() {
 BamFile = "";
 InsertSize = 0;
 }
 std::string BamFile;
 int InsertSize;
 };
 
 struct SampleAndBAMFiles {
 std::string SampleName;
 std::vector <BAM_Path_IS> BAMs;
 };
 
*/


void do_rd_rp (ControlState & CurrentState, ParCollection & par) {
    // bin size and other parameters
    unsigned BinSize = 10000000; // 10M
    unsigned Skip = 100;
    unsigned ExtraRegionSize = 10000; // 10k
    unsigned CoverageCount[BinSize / Skip + 1];
    unsigned binStart, binEnd, binStart_Extra, binEnd_Extra;
    //unsigned ReadCount[BinSize];
    std::vector <DiscordantRP> RP_oneChr; 
    short RD_MQ = 20;
    short RP_MQ = 20;
    
    // get sample names and their associated BAM files
    std::map<std::string, int> SampleName2ID;
    std::map<std::string, int>::iterator it;
    std::vector <SampleAndBAMFiles> AllSamples;
    BAM_Path_IS TempOneBam;
    for (unsigned int BAMindex = 0; BAMindex < CurrentState.bams_to_parse.size(); BAMindex++) {
        it = SampleName2ID.find(CurrentState.bams_to_parse[BAMindex].Tag);
        if (it != SampleName2ID.end()) { // already in the list, add it
           TempOneBam.BamFile = CurrentState.bams_to_parse[BAMindex].BamFile;
           TempOneBam.InsertSize = CurrentState.bams_to_parse[BAMindex].InsertSize;
           AllSamples[(*it).second].BAMs.push_back(TempOneBam); 
        }
        else { // not in the list add one new entry // mymap.insert ( pair<char,int>('a',100) );
            SampleName2ID.insert(std::pair <std::string, int> (CurrentState.bams_to_parse[BAMindex].Tag, SampleName2ID.size())); 
            SampleAndBAMFiles OneSample;
            OneSample.SampleName = CurrentState.bams_to_parse[BAMindex].Tag;
            TempOneBam.BamFile = CurrentState.bams_to_parse[BAMindex].BamFile;
            TempOneBam.InsertSize = CurrentState.bams_to_parse[BAMindex].InsertSize;
            OneSample.BAMs.push_back(TempOneBam);
            AllSamples.push_back(OneSample);
        }
    }

    
    
    // get chromsome sizes using .fai
    std::vector <ChrNameAndSize> ChrNameSizes;
    getWholeGenomeSize(CurrentState, ChrNameSizes);
    // for each sample, get rd signal for every n bases and all discordant reads (the right read must be in the window); write to disk for each sample.

    for (unsigned SampleIndex = 0; SampleIndex < AllSamples.size(); SampleIndex++) {
        std::cout << SampleIndex << " " << AllSamples[SampleIndex].SampleName << " " << AllSamples[SampleIndex].BAMs.size() << std::endl;
        
        for (unsigned ChrIndex = 0; ChrIndex < ChrNameSizes.size(); ChrIndex++) {
            RP_oneChr.clear();
            binStart = 1;
            binEnd = BinSize;
            binStart_Extra = 1;
            binEnd_Extra = binEnd + ExtraRegionSize;
            unsigned NumberOfScan = ChrNameSizes[ChrIndex].ChrSize / BinSize + 1;
            for (unsigned ScanIndex = 0; ScanIndex < NumberOfScan; ScanIndex++) {
                // get all reads for this sample in this window & update RD RP
                for (unsigned posIndex = 0; posIndex < BinSize / Skip + 1; posIndex++) CoverageCount[posIndex] = 0;
                for (unsigned BAMindex = 0; BAMindex < AllSamples[SampleIndex].BAMs.size(); BAMindex++) {
                    std::cout << "\t" << AllSamples[SampleIndex].BAMs[BAMindex].BamFile << " " << AllSamples[SampleIndex].BAMs[BAMindex].InsertSize << std::endl;
                    UpdateRDRP(AllSamples[SampleIndex].BAMs[BAMindex].BamFile, AllSamples[SampleIndex].BAMs[BAMindex].InsertSize, ChrNameSizes[ChrIndex].ChrName, binStart, binEnd, binStart_Extra, binEnd_Extra, RD_MQ, RP_MQ, CoverageCount, RP_oneChr);
                }
                
                // report RD
                                    
                // adjust windows
                binStart += BinSize;
                binEnd += BinSize;
                if (binStart > ExtraRegionSize) binStart_Extra = binStart - ExtraRegionSize;
                else binStart_Extra = 1;
                if (binEnd + ExtraRegionSize < ChrNameSizes[ChrIndex].ChrSize) binEnd_Extra = binEnd + ExtraRegionSize;
                else binEnd_Extra = ChrNameSizes[ChrIndex].ChrSize;
            }
            // report RP
        }
    }
    // for each sample {
    //    clear data 
    //    for each chromosome {
    //       for each 10MB get reads {
    //          // count reads
    //          // get discordant reads
    //       }
    //    }
    //    write the data to hard drive: *_RD and *_RP
    // }

    return;
}

short getWholeGenomeSize(ControlState & CurrentState, std::vector <ChrNameAndSize> & ChrNameSizes) {
    ChrNameAndSize TempOneChr;
    std::string TempStr;
    while (CurrentState.inf_Seq_Fai >> TempOneChr.ChrName >> TempOneChr.ChrSize) {
       safeGetline(CurrentState.inf_Seq_Fai, TempStr); 
       std::cout << "ChrName: " << TempOneChr.ChrName << "\tChrSize: " << TempOneChr.ChrSize << std::endl;
       ChrNameSizes.push_back(TempOneChr); 
    }

    CurrentState.inf_Seq_Fai.close();
    return 0;
}
