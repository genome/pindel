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

void SortOutputD (const unsigned &NumBoxes, const std::string & CurrentChr,
									std::vector < SPLIT_READ > &AllReads,
									std::vector < unsigned >Deletions[],
									std::ofstream & DeletionOutf);
void SortOutputSI (const unsigned &NumBoxes, const std::string & CurrentChr,
									 std::vector < SPLIT_READ > &AllReads,
									 std::vector < unsigned >SIs[], std::ofstream & SIsOutf);
void SortOutputTD (const unsigned &NumBoxes, const std::string & CurrentChr,
									 std::vector < SPLIT_READ > &AllReads,
									 std::vector < unsigned >TDs[], std::ofstream & TDOutf);
void SortOutputTD_NT (const unsigned &NumBoxes,
											const std::string & CurrentChr,
											std::vector < SPLIT_READ > &AllReads,
											std::vector < unsigned >TDs[], std::ofstream & TDOutf);
void SortOutputInv (const unsigned &NumBoxes, const std::string & CurrentChr,
										std::vector < SPLIT_READ > &AllReads,
										std::vector < unsigned >Inv[], std::ofstream & InvOutf);
void SortOutputInv_NT (const unsigned &NumBoxes,
											 const std::string & CurrentChr,
											 std::vector < SPLIT_READ > &AllReads,
											 std::vector < unsigned >Inv[],
											 std::ofstream & InvOutf);
void SortOutputDI (const unsigned &NumBoxes, const std::string & CurrentChr,
									 std::vector < SPLIT_READ > &AllReads,
									 std::vector < unsigned >DI[], std::ofstream & DIOutf);
void SortOutputLI (const std::string & CurrentChr,
									 std::vector < SPLIT_READ > &AllReads,
									 std::ofstream & Outf_LI);
void SortOutputRest (const std::string & CurrentChr,
										 std::vector < SPLIT_READ > &AllReads,
										 std::vector < SPLIT_READ > &BP_Reads,
										 std::ofstream & Outf_Rest);


#endif /* REPORTER_H */
