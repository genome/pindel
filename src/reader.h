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

#ifndef READER_H
#define	READER_H

#include "read_buffer.h"
#include "control_state.h"

class PindelReadReader;

void GetOneChrSeq (std::ifstream & inf_Seq, std::string & CurrentChr, bool WhetherBuildUp);
bool ReadInBamReads_SR (const char *bam_path,
                        const std::string & FragName,
						const std::string * CurrentChr,
                        std::vector < SPLIT_READ > &LeftReads,
                        std::vector < SPLIT_READ > & OneEndMappedReads,
                        std::vector <REF_READ> &RefSupportingReads,
						int InsertSize,
                        std::string Tag,
                        const SearchWindow& currentWindow,
                        ReadBuffer& readBuffer,
                        bool verbose=true);
/*bool ReadInBamReads_RP (const char *bam_path, const std::string & FragName,
                        std::string * CurrentChr, std::vector <RP_READ> &LeftReads, 
                        int InsertSize, std::string Tag, const SearchWindow& currentWindow, ReadBuffer& readBuffer);
//short ReadInRead (PindelReadReader & inf_ReadSeq,
									const std::string & CurrentFragName,
									const std::string & CurrentFrag,
									std::vector < SPLIT_READ > &Reads, const unsigned int lowerBinBorder,
									const unsigned int upperBinBorder);*/

short get_SR_Reads(ControlState& currentState, const SearchWindow& currentWindow );
short get_RP_Reads(ControlState& currentState, const SearchWindow& currentWindow );
short get_RP_Reads_Discovery(ControlState& currentState, const SearchWindow& currentWindow );
//short get_RP_Reads(ControlState& currentState, const SearchWindow& currentWindow );
//void readInPindelReads(PindelReadReader &reader, const std::string& pindelFilename, ControlState& currentState, const SearchWindow& currentWindow  );

// Get the read group ID as a string for an alignment.
void get_read_group(const bam1_t* alignment, std::string& read_group);

#endif /* READER_H */
