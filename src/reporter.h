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

#ifndef REPORTER_H
#define	REPORTER_H

#include <fstream>
#include <string>
#include <vector>
#include "pindel.h"
#include "shifted_vector.h"
#include "control_state.h"
#include "bddata.h"

struct Variant {
    std::string VariantType;
    unsigned Length;
    unsigned NT_length;
    std::string NT_str;
    std::string ChrName;
    unsigned Start;
    unsigned End;
    unsigned Start_Range;
    unsigned End_Range;
    unsigned RefSupport;
    unsigned AlleleSupport;
    bool Report;
};

struct VariantsPerChr {
    
    std::string ChrName;
    std::string ChrSeq;
    std::vector <unsigned> IndexOfVariants;
};


void SortOutputD (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
									std::vector < SPLIT_READ > &AllReads,
									std::vector < unsigned >Deletions[],
									std::ofstream & DeletionOutf);
void SortOutputSI (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
									 std::vector < SPLIT_READ > &AllReads,
									 std::vector < unsigned >SIs[], std::ofstream & SIsOutf);
void SortAndOutputTandemDuplications (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
									  std::vector < SPLIT_READ > &AllReads, std::vector < unsigned >TDs[],
									  std::ofstream & TDOutf, const bool nonTemplate);
void SortOutputInv (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
										std::vector < SPLIT_READ > &AllReads,
										std::vector < unsigned >Inv[], std::ofstream & InvOutf);
void SortOutputInv_NT (ControlState& currentState, const unsigned &NumBoxes,
											 const std::string & CurrentChr,
											 std::vector < SPLIT_READ > &AllReads,
											 std::vector < unsigned >Inv[],
											 std::ofstream & InvOutf);
void SortOutputDI (ControlState& currentState, const unsigned &NumBoxes, const std::string & CurrentChr,
							std::vector < SPLIT_READ > &Reads, std::vector < unsigned >DI[],
							std::ofstream & DIOutf, std::ofstream & InvOutf);
void SortOutputLI (ControlState& currentState, const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads, const SearchWindow& currentWindow, const std::string& filename);
void SortOutputRest (ControlState& currentState, const std::string & CurrentChr,  std::vector < SPLIT_READ > &Reads,  const SearchWindow& currentWindow, const std::string& filename);

void SortAndReportInterChromosomalEvents(ControlState& current_state, Genome& genome, UserDefinedSettings* user_settings);

void GetConsensusBasedOnPloidy(ControlState& current_state, Genome& genome, UserDefinedSettings* user_settings);

#endif /* REPORTER_H */
