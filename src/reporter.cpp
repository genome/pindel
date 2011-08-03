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

// Pindel header files
#include "pindel.h"
#include "logdef.h"
#include "output_file_data.h"
#include "output_sorter.h"
#include "reporter.h"

OutputFileData deletionFileData;

void
OutputTDs (const std::vector < SPLIT_READ > &TDs,
					 const std::string & TheInput,
					 const unsigned int &C_S,
					 const unsigned int &C_E,
					 const unsigned int &RealStart,
					 const unsigned int &RealEnd, std::ofstream & TDOutf)
{

	//short ReadLength = Deletions[C_S].ReadLength;
	//short ReadLengthMinus = ReadLength - 1;
	unsigned int NumberOfReads = C_E - C_S + 1;
	//float LeftScore = 0;
	//float RightScore = 0;
	unsigned int LeftS = 1;
	unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;

	SupportPerSample NumSupportPerTag[VectorTag.size ()];

	for (unsigned short i = 0; i < VectorTag.size (); i++)
		{
			NumSupportPerTag[i].NumPlus = 0;
			NumSupportPerTag[i].NumMinus = 0;
			NumSupportPerTag[i].NumUPlus = 0;
			NumSupportPerTag[i].NumUMinus = 0;
		}

	for (unsigned int i = C_S; i <= C_E; i++)
		{
			for (unsigned short j = 0; j < VectorTag.size (); j++)
				{
					if (TDs[i].Tag == VectorTag[j])
						{
							if (TDs[i].MatchedD == Plus)
								{
									NumSupportPerTag[j].NumPlus++;
									if (TDs[i].Unique)
										NumSupportPerTag[j].NumUPlus++;
								}
							else
								{
									NumSupportPerTag[j].NumMinus++;
									if (TDs[i].Unique)
										NumSupportPerTag[j].NumUMinus++;
								}
						}
				}
			if (TDs[i].MatchedD == Plus)
				{
				  //LeftScore += TDs[i].score;
					LeftS++;
					if (TDs[i].Unique)
						LeftUNum++;
				}
			else
				{
				  //RightScore += TDs[i].score;
					RightS++;
					if (TDs[i].Unique)
						RightUNum++;
				}
		}

	short NumberSupSamples = 0;
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
	for (unsigned short i = 0; i < VectorTag.size (); ++i)
		{
			if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus)
				++NumberSupSamples;
			if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus)
				++NumU_SupSamples;
			Num_U_Reads +=
				NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
		}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples

	unsigned int EasyScore = LeftS * RightS;
	//double PreciseScore = (LeftScore + RightScore) * (-1);
	//short GapSize = 0;
	//if (TDs[C_S].IndelSize < 14)
	//	GapSize = TDs[C_S].IndelSize;
	//else
	//	GapSize = 13 + (int) log10 (TDs[C_S].IndelSize - 10);
	CurrentChrMask[TDs[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[TDs[C_S].BPRight + SpacerBeforeAfter] = 'B';

	TDOutf <<
		"####################################################################################################"
		<< std::endl;
	TDOutf << NumberOfTDInstances << "\tTD " << TDs[C_S].IndelSize	// << " bases " 
		<< "\tNT " << TDs[C_S].NT_size << " \"" << TDs[C_S].NT_str << "\"" << "\tChrID " << TDs[C_S].FragName << "\tBP " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2 << "\tBP_range " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111 << "\tS2 " << PreciseScore; 

	int SUM_MS = 0;
	for (unsigned int i = C_S; i <= C_E; i++)
		SUM_MS += TDs[i].MS;
	TDOutf << "\tSUM_MS " << SUM_MS;

	TDOutf << "\t" << VectorTag.
		size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
		NumU_SupSamples;
	for (unsigned short i = 0; i < VectorTag.size (); i++)
		TDOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
			<< " " << NumSupportPerTag[i].NumUPlus
			<< " " << NumSupportPerTag[i].NumMinus
			<< " " << NumSupportPerTag[i].NumUMinus;
	TDOutf << std::endl;

	//DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
	TDOutf << TheInput.substr (TDs[C_S].BPRight + SpacerBeforeAfter - ReportLength + 1, ReportLength);	// << endl;// 
	if (TDs[C_S].NT_size)
		{
			for (short i = 0; i < TDs[C_S].NT_size; i++)
				TDOutf << " ";
		}
	//ReportLength 
	TDOutf << Cap2Low (TheInput.
										 substr (TDs[C_S].BPLeft + SpacerBeforeAfter,
														 ReportLength)) << std::endl;

	short SpaceBeforeReadSeq;
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++)
		{
			SpaceBeforeReadSeq = ReportLength - TDs[GoodIndex].BP - 1;

			for (int i = 0; i < SpaceBeforeReadSeq; i++)
				TDOutf << " ";
			//short SpaceBeforeD =
			//	ReportLength + ReportLength - SpacerBeforeAfter -
			//	TDs[GoodIndex].ReadLength;
			if (TDs[GoodIndex].MatchedD == Minus)
				{
					TDOutf << TDs[GoodIndex].UnmatchedSeq << std::endl;
					//for (int i = 0; i < GapSize; i++) TDOutf << " ";    
					//TDOutf << TDs[GoodIndex].UnmatchedSeq.substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
				}
			else
				{
					TDOutf << ReverseComplement (TDs[GoodIndex].UnmatchedSeq) << std::endl;
					//for (int i = 0; i < GapSize; i++) TDOutf << " ";  
					//TDOutf << ReverseComplement(TDs[GoodIndex].UnmatchedSeq).substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
				}
			//for (int i = 0; i < SpaceBeforeD; i++) TDOutf << " ";
			TDOutf << "\t" << TDs[GoodIndex].MatchedD << "\t"
				<< TDs[GoodIndex].MatchedRelPos
				<< "\t" << TDs[GoodIndex].MS
				<< "\t" << TDs[GoodIndex].Tag << "\t" << TDs[GoodIndex].Name << std::endl;
			//<< "\t" << TDs[C_S].BPLeft
			//<< "\t" << TDs[C_S].BPRight << endl; 
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
	//short ReadLength = Deletions[C_S].ReadLength;
	//short ReadLengthMinus = ReadLength - 1;
	LOG_DEBUG(std::cout << "d_1" << std::endl);
	unsigned int NumberOfReads = C_E - C_S + 1;
	//float LeftScore = 0;
	//float RightScore = 0;
	unsigned int LeftS = 1;
	unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	LOG_DEBUG(std::cout << "d_2" << std::endl);
	SupportPerSample NumSupportPerTag[VectorTag.size ()];

	DeletionOutf << "EWL BREAKING JENKINS! (see if it still works)\n";

	for (unsigned short i = 0; i < VectorTag.size (); i++)
		{
			NumSupportPerTag[i].NumPlus = 0;
			NumSupportPerTag[i].NumMinus = 0;
			NumSupportPerTag[i].NumUPlus = 0;
			NumSupportPerTag[i].NumUMinus = 0;
		}
	LOG_DEBUG(std::cout << "d_3" << std::endl);
	for (unsigned int i = C_S; i <= C_E; i++)
		{
			for (unsigned short j = 0; j < VectorTag.size (); j++)
				{
					if (Deletions[i].Tag == VectorTag[j])
						{
							if (Deletions[i].MatchedD == Plus)
								{
									NumSupportPerTag[j].NumPlus++;
									if (Deletions[i].Unique)
										NumSupportPerTag[j].NumUPlus++;
								}
							else
								{
									NumSupportPerTag[j].NumMinus++;
									if (Deletions[i].Unique)
										NumSupportPerTag[j].NumUMinus++;
								}
						}
				}
			if (Deletions[i].MatchedD == Plus)
				{
				  //LeftScore += Deletions[i].score;
					LeftS++;
					if (Deletions[i].Unique)
						LeftUNum++;
				}
			else
				{
				  //RightScore += Deletions[i].score;
					RightS++;
					if (Deletions[i].Unique)
						RightUNum++;
				}
		}
	LOG_DEBUG(std::cout << "d_4" << std::endl);
	short NumberSupSamples = 0;
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
	for (unsigned short i = 0; i < VectorTag.size (); ++i)
		{
			if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus)
				++NumberSupSamples;
			if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus)
				++NumU_SupSamples;
			Num_U_Reads +=
				NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
		}
	//  << "\t" << Num_U_Reads
	LOG_DEBUG(std::cout << "d_5" << std::endl);
	unsigned int EasyScore = LeftS * RightS;
	//double PreciseScore = (LeftScore + RightScore) * (-1);
	short GapSize = 0;
	if (Deletions[C_S].IndelSize < 14)
		GapSize = Deletions[C_S].IndelSize;
	else
		GapSize = 13 + (int) log10 (Deletions[C_S].IndelSize - 10);
	LOG_DEBUG(std::cout << "d_5a" << std::endl);
	CurrentChrMask[Deletions[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	LOG_DEBUG(std::cout << "d_5b" << std::endl);
	CurrentChrMask[Deletions[C_S].BPRight + SpacerBeforeAfter] = 'B';
	LOG_DEBUG(std::cout << "d_5c" << std::endl);
	CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	LOG_DEBUG(std::cout << "d_5d" << std::endl);
	CurrentChrMask[RealEnd + SpacerBeforeAfter] = 'B';
	LOG_DEBUG(std::cout << "d_5e" << std::endl);
	DeletionOutf <<
		"####################################################################################################"
		<< std::endl;
	DeletionOutf << deletionFileData.getSvIndex() << "\tD " << Deletions[C_S].IndelSize	// << " bases " 
		<< "\tNT " << Deletions[C_S].NT_size << " \"" << Deletions[C_S].NT_str << "\"" << "\tChrID " << Deletions[C_S].FragName << "\tBP " << Deletions[C_S].BPLeft + 1 << "\t" << Deletions[C_S].BPRight + 1 << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore; 
	LOG_DEBUG(std::cout << "d_6" << std::endl);
	int SUM_MS = 0;
	for (unsigned int i = C_S; i <= C_E; i++)
		SUM_MS += Deletions[i].MS;
	DeletionOutf << "\tSUM_MS " << SUM_MS;

	DeletionOutf << "\t" << VectorTag.
		size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
		NumU_SupSamples;
	for (unsigned short i = 0; i < VectorTag.size (); i++)
		DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
			<< " " << NumSupportPerTag[i].NumUPlus
			<< " " << NumSupportPerTag[i].NumMinus
			<< " " << NumSupportPerTag[i].NumUMinus;
	DeletionOutf << std::endl;
	LOG_DEBUG(std::cout << "d_7" << std::endl);
	//DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
	DeletionOutf << TheInput.substr (Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength);	// << endl;// ReportLength    
	if (Deletions[C_S].IndelSize >= 14)
		{
			DeletionOutf << Cap2Low (TheInput.
															 substr (Deletions[C_S].Left +
																			 Deletions[C_S].BP + 1,
																			 5)) << "<" << Deletions[C_S].
				IndelSize -
				10 << ">" << Cap2Low (TheInput.
															substr (Deletions[C_S].Right -
																			Deletions[C_S].ReadLength +
																			Deletions[C_S].BP - 3, 5));
		}
	else
		DeletionOutf << Cap2Low (TheInput.
														 substr (Deletions[C_S].Left + Deletions[C_S].BP +
																		 1, GapSize));
	DeletionOutf << TheInput.substr (Deletions[C_S].Left + Deletions[C_S].BP + 1 + Deletions[C_S].IndelSize, ReportLength - GapSize) << std::endl;	// ReportLength
	short SpaceBeforeReadSeq;
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++)
		{
			SpaceBeforeReadSeq = ReportLength - Deletions[GoodIndex].BP - 1;

			for (int i = 0; i < SpaceBeforeReadSeq; i++)
				DeletionOutf << " ";
			short SpaceBeforeD =
				ReportLength + ReportLength - SpaceBeforeReadSeq -
				Deletions[GoodIndex].ReadLength;
			if (Deletions[GoodIndex].MatchedD == Minus)
				{
					DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr (0, Deletions[GoodIndex].BP + 1);	// << endl;
					for (int i = 0; i < GapSize; i++)
						DeletionOutf << " ";
					DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr (Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);	// << endl;
				}
			else
				{
					DeletionOutf << ReverseComplement (Deletions[GoodIndex].UnmatchedSeq).substr (0, Deletions[GoodIndex].BP + 1);	// << endl; 
					for (int i = 0; i < GapSize; i++)
						DeletionOutf << " ";
					DeletionOutf << ReverseComplement (Deletions[GoodIndex].UnmatchedSeq).substr (Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);	// << endl;
				}
			for (int i = 0; i < SpaceBeforeD; i++)
				DeletionOutf << " ";
			DeletionOutf << "\t" << Deletions[GoodIndex].MatchedD << "\t"
				<< Deletions[GoodIndex].MatchedRelPos
				<< "\t" << Deletions[GoodIndex].MS
				<< "\t" << Deletions[GoodIndex].Tag
				<< "\t" << Deletions[GoodIndex].Name << std::endl;
			//<< "\t" << Deletions[C_S].BPLeft
			//<< "\t" << Deletions[C_S].BPRight << endl; 
		}
}

void
OutputInversions (const std::vector < SPLIT_READ > &Inv,
									const std::string & TheInput,
									const unsigned int &C_S,
									const unsigned int &C_E,
									const unsigned int &RealStart,
									const unsigned int &RealEnd, std::ofstream & InvOutf)
{
	int LeftNT_index = -1;
	int RightNT_index = -1;
	for (unsigned Index = C_S; Index <= C_E; Index++)
		{
			if (Inv[Index].MatchedD == Plus)
				{
					LeftNT_index = Index;
					break;
				}
		}
	for (unsigned Index = C_S; Index <= C_E; Index++)
		{
			if (Inv[Index].MatchedD == Minus)
				{
					RightNT_index = Index;
					break;
				}
		}
	short LeftNT_size = 0;
	short RightNT_size = 0;
	std::string LeftNT_str = "";
	std::string RightNT_str = "";
	if (LeftNT_index != -1)
		{
			LeftNT_size = Inv[LeftNT_index].NT_size;
			LeftNT_str = Inv[LeftNT_index].NT_str;
		}
	if (RightNT_index != -1)
		{
			RightNT_size = Inv[RightNT_index].NT_size;
			RightNT_str = Inv[RightNT_index].NT_str;
		}
	unsigned int NumberOfReads = C_E - C_S + 1;
	//float LeftScore = 0;
	//float RightScore = 0;
	unsigned int LeftS = 1;
	unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;

	SupportPerSample NumSupportPerTag[VectorTag.size ()];

	for (unsigned short i = 0; i < VectorTag.size (); i++)
		{
			NumSupportPerTag[i].NumPlus = 0;
			NumSupportPerTag[i].NumMinus = 0;
			NumSupportPerTag[i].NumUPlus = 0;
			NumSupportPerTag[i].NumUMinus = 0;
		}

	for (unsigned int i = C_S; i <= C_E; i++)
		{
			for (unsigned short j = 0; j < VectorTag.size (); j++)
				{
					if (Inv[i].Tag == VectorTag[j])
						{
							if (Inv[i].MatchedD == Plus)
								{
									NumSupportPerTag[j].NumPlus++;
									if (Inv[i].Unique)
										NumSupportPerTag[j].NumUPlus++;
								}
							else
								{
									NumSupportPerTag[j].NumMinus++;
									if (Inv[i].Unique)
										NumSupportPerTag[j].NumUMinus++;
								}
						}
				}
			if (Inv[i].MatchedD == Plus)
				{
				  //LeftScore += Inv[i].score;
					LeftS++;
					if (Inv[i].Unique)
						LeftUNum++;
				}
			else
				{
				  //RightScore += Inv[i].score;
					RightS++;
					if (Inv[i].Unique)
						RightUNum++;
				}
		}

	short NumberSupSamples = 0;
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
	for (unsigned short i = 0; i < VectorTag.size (); ++i)
		{
			if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus)
				++NumberSupSamples;
			if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus)
				++NumU_SupSamples;
			Num_U_Reads +=
				NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
		}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples

	unsigned int EasyScore = LeftS * RightS;
	//double PreciseScore = (LeftScore + RightScore) * (-1);
	//short GapSize = 0;
	//if (Inv[C_S].IndelSize < 14)
	//	GapSize = Inv[C_S].IndelSize;
	//else
	//	GapSize = 13 + (int) log10 (Inv[C_S].IndelSize - 10);
	CurrentChrMask[Inv[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[Inv[C_S].BPRight + SpacerBeforeAfter] = 'B';
	InvOutf <<
		"####################################################################################################"
		<< std::endl;
	InvOutf << NumberOfInvInstances << "\tINV " << Inv[C_S].IndelSize	// << " bases " 
		<< "\tNT " << LeftNT_size << ":" << RightNT_size << " \"" << LeftNT_str << "\":\"" << RightNT_str << "\"" << "\tChrID " << Inv[C_S].FragName << "\tBP " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 << "\tBP_range " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore; 

	int SUM_MS = 0;
	for (unsigned int i = C_S; i <= C_E; i++)
		SUM_MS += Inv[i].MS;
	InvOutf << "\tSUM_MS " << SUM_MS;

	InvOutf << "\t" << VectorTag.
		size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
		NumU_SupSamples;
	for (unsigned short i = 0; i < VectorTag.size (); i++)
		InvOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
			<< " " << NumSupportPerTag[i].NumUPlus
			<< " " << NumSupportPerTag[i].NumMinus
			<< " " << NumSupportPerTag[i].NumUMinus;
	InvOutf << std::endl;

	short SpaceBeforeReadSeq;
	//DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
	InvOutf << TheInput.substr (Inv[C_S].BPLeft + SpacerBeforeAfter - ReportLength, ReportLength);	//;// << endl;// ReportLength  
	LOG_DEBUG(std::cout << Inv[C_S].NT_size << "\t" << Inv[C_S].NT_2size << std::endl);
	if (LeftNT_size)
		{
			for (int i = 0; i < LeftNT_size; i++)
				{
					InvOutf << " ";
				}
		}
	InvOutf <<
		Cap2Low (ReverseComplement
						 (TheInput.
							substr (Inv[C_S].BPRight + 1 + SpacerBeforeAfter - ReportLength,
											ReportLength))) << std::endl;
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++)
		{
			if (Inv[GoodIndex].MatchedD == Plus
					&& Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPLeft)
				{
					SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
					for (int i = 0; i < SpaceBeforeReadSeq; i++)
						InvOutf << " ";
					InvOutf << ReverseComplement (Inv[GoodIndex].UnmatchedSeq);
					for (int i = 0; i < Inv[GoodIndex].BP; i++)
						InvOutf << " ";
					InvOutf								//<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
						<< "\t" << Inv[GoodIndex].MatchedD << "\t"
						<< Inv[GoodIndex].MatchedRelPos
						<< "\t" << Inv[GoodIndex].MS
						<< "\t" << Inv[GoodIndex].Tag
						<< "\t" << Inv[GoodIndex].Name << std::endl;
					//<< "\t" << Deletions[C_S].BPLeft
					//<< "\t" << Deletions[C_S].BPRight << endl; 
				}
			else if (Inv[GoodIndex].MatchedD == Plus
							 && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPLeft)
				{
					SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
					for (int i = 0; i < SpaceBeforeReadSeq; i++)
						InvOutf << " ";
					InvOutf << Inv[GoodIndex].UnmatchedSeq;
					InvOutf								//<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
						<< "\t" << Inv[GoodIndex].MatchedD << "\t"
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
							substr (Inv[C_S].BPLeft + SpacerBeforeAfter, ReportLength)));
	if (RightNT_size)
		{
			for (int i = 0; i < RightNT_size; i++)
				{
					InvOutf << " ";
				}
		}
	InvOutf << TheInput.substr (Inv[C_S].BPRight + 1 + SpacerBeforeAfter,
															ReportLength) << std::endl;
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++)
		{
			if (Inv[GoodIndex].MatchedD == Minus
					&& Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPRight)
				{
					SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
					for (int i = 0; i < SpaceBeforeReadSeq; i++)
						InvOutf << " ";
					InvOutf << ReverseComplement (Inv[GoodIndex].UnmatchedSeq);
					for (int i = 0; i < Inv[GoodIndex].BP; i++)
						InvOutf << " ";
					InvOutf								//<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str 
						<< "\t" << Inv[GoodIndex].MatchedD << "\t"
						<< Inv[GoodIndex].MatchedRelPos
						<< "\t" << Inv[GoodIndex].MS
						<< "\t" << Inv[GoodIndex].Tag
						<< "\t" << Inv[GoodIndex].Name << std::endl;
					//<< "\t" << Deletions[C_S].BPLeft
					//<< "\t" << Deletions[C_S].BPRight << endl; 
				}
			else if (Inv[GoodIndex].MatchedD == Minus
							 && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPRight)
				{
					SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
					for (int i = 0; i < SpaceBeforeReadSeq; i++)
						InvOutf << " ";
					InvOutf << Inv[GoodIndex].UnmatchedSeq;
					InvOutf								//<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
						<< "\t" << Inv[GoodIndex].MatchedD << "\t"
						<< Inv[GoodIndex].MatchedRelPos
						<< "\t" << Inv[GoodIndex].MS
						<< "\t" << Inv[GoodIndex].Tag
						<< "\t" << Inv[GoodIndex].Name << std::endl;
				}
		}


}

void
OutputSIs (const std::vector < SPLIT_READ > &SIs,
					 const std::string & TheInput,
					 const unsigned int &C_S,
					 const unsigned int &C_E,
					 const unsigned int &RealStart,
					 const unsigned int &RealEnd, std::ofstream & SIsOutf)
{
	//short ReadLength = SIs[C_S].ReadLength;
	//short ReadLengthMinus = ReadLength - 1;               
	unsigned int NumberOfReads = C_E - C_S + 1;
	//float LeftScore = 0;
	//float RightScore = 0;
	unsigned int LeftS = 1;
	unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;

	SupportPerSample NumSupportPerTag[VectorTag.size ()];

	for (unsigned short i = 0; i < VectorTag.size (); i++)
		{
			NumSupportPerTag[i].NumPlus = 0;
			NumSupportPerTag[i].NumMinus = 0;
			NumSupportPerTag[i].NumUPlus = 0;
			NumSupportPerTag[i].NumUMinus = 0;
		}

	for (unsigned int i = C_S; i <= C_E; i++)
		{
			for (unsigned short j = 0; j < VectorTag.size (); j++)
				{
					if (SIs[i].Tag == VectorTag[j])
						{
							if (SIs[i].MatchedD == Plus)
								{
									NumSupportPerTag[j].NumPlus++;
									if (SIs[i].Unique)
										NumSupportPerTag[j].NumUPlus++;
								}
							else
								{
									NumSupportPerTag[j].NumMinus++;
									if (SIs[i].Unique)
										NumSupportPerTag[j].NumUMinus++;
								}
						}
				}
			if (SIs[i].MatchedD == Plus)
				{
				  //LeftScore += SIs[i].score;
					LeftS++;
					if (SIs[i].Unique)
						LeftUNum++;
				}
			else
				{
				  //RightScore += SIs[i].score;
					RightS++;
					if (SIs[i].Unique)
						RightUNum++;
				}
		}


	short NumberSupSamples = 0;
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
	for (unsigned short i = 0; i < VectorTag.size (); ++i)
		{
			if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus)
				++NumberSupSamples;
			if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus)
				++NumU_SupSamples;
			Num_U_Reads +=
				NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
		}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples

	unsigned int EasyScore = LeftS * RightS;
	//double PreciseScore = (LeftScore + RightScore) * (-1);
	std::string CurrentReadSeq;

	CurrentChrMask[SIs[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[SIs[C_S].BPRight + SpacerBeforeAfter] = 'B';
	CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	CurrentChrMask[RealEnd + SpacerBeforeAfter] = 'B';

	SIsOutf <<
		"####################################################################################################"
		<< std::endl;
	SIsOutf << NumberOfSIsInstances << "\tI " << SIs[C_S].IndelSize << "\tNT " << SIs[C_S].IndelSize << " \"" << SIs[C_S].InsertedStr << "\"" << "\tChrID " << SIs[C_S].FragName << "\tBP " << SIs[C_S].BPLeft + 1 << "\t" << SIs[C_S].BPRight + 1 << "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	//EWL070111  << "\tS2 " << PreciseScore; 

	int SUM_MS = 0;
	for (unsigned int i = C_S; i <= C_E; i++)
		SUM_MS += SIs[i].MS;
	SIsOutf << "\tSUM_MS " << SUM_MS;

	SIsOutf << "\t" << VectorTag.
		size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
		NumU_SupSamples;
	for (unsigned short i = 0; i < VectorTag.size (); i++)
		SIsOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
			<< " " << NumSupportPerTag[i].NumUPlus
			<< " " << NumSupportPerTag[i].NumMinus
			<< " " << NumSupportPerTag[i].NumUMinus;
	SIsOutf << std::endl;

	SIsOutf << TheInput.substr (SIs[C_S].Left - ReportLength + SIs[C_S].BP + 1, ReportLength);	// ReportLength
	for (unsigned int i = 0; i < SIs[C_S].IndelSize; i++)
		SIsOutf << " ";
	SIsOutf << TheInput.substr (SIs[C_S].Left + SIs[C_S].BP + 1, ReportLength) << std::endl;	// ReportLength

	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++)
		{
			short SpaceBeforeReadSeq = ReportLength - SIs[GoodIndex].BP - 1;
			for (short i = 0; i < SpaceBeforeReadSeq; i++)
				SIsOutf << " ";
			if (SIs[GoodIndex].MatchedD == Minus)
				SIsOutf << SIs[GoodIndex].UnmatchedSeq;
			else
				SIsOutf << ReverseComplement (SIs[GoodIndex].UnmatchedSeq);
			short SpaceBeforeD =
				ReportLength + ReportLength - SpaceBeforeReadSeq -
				SIs[GoodIndex].ReadLength;
			for (short i = 0; i < SpaceBeforeD; i++)
				SIsOutf << " ";
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
	//short ReadLength = DI[C_S].ReadLength;
	//short ReadLengthMinus = ReadLength - 1;
	unsigned int NumberOfReads = C_E - C_S + 1;
	//float LeftScore = 0;
	//float RightScore = 0;
	unsigned int LeftS = 1;
	unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;

	SupportPerSample NumSupportPerTag[VectorTag.size ()];

	for (unsigned short i = 0; i < VectorTag.size (); i++)
		{
			NumSupportPerTag[i].NumPlus = 0;
			NumSupportPerTag[i].NumMinus = 0;
			NumSupportPerTag[i].NumUPlus = 0;
			NumSupportPerTag[i].NumUMinus = 0;
		}

	for (unsigned int i = C_S; i <= C_E; i++)
		{
			for (unsigned short j = 0; j < VectorTag.size (); j++)
				{
					if (DI[i].Tag == VectorTag[j])
						{
							if (DI[i].MatchedD == Plus)
								{
									NumSupportPerTag[j].NumPlus++;
									if (DI[i].Unique)
										NumSupportPerTag[j].NumUPlus++;
								}
							else
								{
									NumSupportPerTag[j].NumMinus++;
									if (DI[i].Unique)
										NumSupportPerTag[j].NumUMinus++;
								}
						}
				}
			if (DI[i].MatchedD == Plus)
				{
				  //LeftScore += DI[i].score;
					LeftS++;
					if (DI[i].Unique)
						LeftUNum++;
				}
			else
				{
				  //RightScore += DI[i].score;
					RightS++;
					if (DI[i].Unique)
						RightUNum++;
				}
		}

	short NumberSupSamples = 0;
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
	for (unsigned short i = 0; i < VectorTag.size (); ++i)
		{
			if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus)
				++NumberSupSamples;
			if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus)
				++NumU_SupSamples;
			Num_U_Reads +=
				NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
		}
	//  << "\t" << Num_U_Reads

	unsigned int EasyScore = LeftS * RightS;
	//double PreciseScore = (LeftScore + RightScore) * (-1);
	//short GapSize =0;
	//if (DI[C_S].IndelSize < 14) GapSize = DI[C_S].IndelSize;
	//else GapSize = 13 + (int)log10(Deletions[C_S].IndelSize - 10);
	CurrentChrMask[DI[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[DI[C_S].BPRight + SpacerBeforeAfter] = 'B';
	//CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	//CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B'; 
	DeletionOutf <<
		"####################################################################################################"
		<< std::endl;
	DeletionOutf << deletionFileData.getSvIndex() << "\tD " << DI[C_S].IndelSize << "\tNT " << DI[C_S].NT_size << " \"" << DI[C_S].NT_str << "\"" << "\tChrID " << DI[C_S].FragName << "\tBP " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 << "\tBP_range " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 << "\tSupports " << NumberOfReads << "\t" << Num_U_Reads << "\t+ " << LeftS - 1 << "\t" << LeftUNum << "\t- " << RightS - 1 << "\t" << RightUNum << "\tS1 " << EasyScore;	// << "\tS2 0.0";// << PreciseScore << "\t"; 
	//EWL070111 << "\tS2 " << PreciseScore; 

	int SUM_MS = 0;
	for (unsigned int i = C_S; i <= C_E; i++)
		SUM_MS += DI[i].MS;
	DeletionOutf << "\tSUM_MS " << SUM_MS;


	DeletionOutf << "\t" << VectorTag.
		size () << "\tNumSupSamples " << NumberSupSamples << "\t" <<
		NumU_SupSamples;
	for (unsigned short i = 0; i < VectorTag.size (); i++)
		DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
			<< " " << NumSupportPerTag[i].NumUPlus
			<< " " << NumSupportPerTag[i].NumMinus
			<< " " << NumSupportPerTag[i].NumUMinus;
	DeletionOutf << std::endl;

	//DeletionOutf << TheInput.substr(DI[C_S].Left - ReportLength + DI[C_S].BP + 1, 2 * ReportLength) << endl;       
	DeletionOutf << TheInput.substr (DI[C_S].Left - ReportLength + DI[C_S].BP + 1, ReportLength);	// << endl;// ReportLength    

	for (short i = 0; i < DI[C_S].NT_size; i++)
		DeletionOutf << " ";
	DeletionOutf << TheInput.substr (DI[C_S].Left + DI[C_S].BP + 1 + DI[C_S].IndelSize, ReportLength) << std::endl;	// ReportLength
	short SpaceBeforeReadSeq;
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++)
		{
			SpaceBeforeReadSeq = ReportLength - DI[GoodIndex].BP - 1;

			for (int i = 0; i < SpaceBeforeReadSeq; i++)
				DeletionOutf << " ";
			//short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - Deletions[GoodIndex].ReadLength;
			if (DI[GoodIndex].MatchedD == Minus)
				{
					DeletionOutf << DI[GoodIndex].UnmatchedSeq << "\t";
					//for (int i = 0; i < GapSize; i++) DeletionOutf << " ";    
					//DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
				}
			else
				{
					DeletionOutf << ReverseComplement (DI[GoodIndex].
																						 UnmatchedSeq) << "\t";
					//for (int i = 0; i < GapSize; i++) DeletionOutf << " ";  
					//DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
				}
			//for (int i = 0; i < SpaceBeforeD; i++) DeletionOutf << " ";
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
	LOG_INFO(std::cout << "Sorting and outputing short insertions ..." << std::endl);
	unsigned int SIsNum;
	short CompareResult;
	unsigned Temp4Exchange;

	//vector <SPLIT_READ> InputIndels; 
	std::vector < SPLIT_READ > GoodIndels;
	unsigned int GoodNum;
	std::vector < Indel4output > IndelEvents;

	for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
		{
			if (SIs[Box_index].size () >= NumRead2ReportCutOff)
				{
					//InputIndels.clear();
					SIsNum = SIs[Box_index].size ();
					LOG_DEBUG(std::cout << "SIsNum " << SIsNum << std::endl);
					//for (int i = 0; i < SIsNum; i++) {
					//  InputIndels.push_back(Reads[SIs[Box_index][i]]);
					//}
					for (unsigned int First = 0; First < SIsNum - 1; First++)
						{
							//if (InputIndels[First].Unique) 
							{
								for (unsigned int Second = First + 1; Second < SIsNum;
										 Second++)
									{
										//if (InputIndels[Second].Unique) 
										{
											if (Reads[SIs[Box_index][First]].ReadLength ==
													Reads[SIs[Box_index][Second]].ReadLength)
												{
													if (Reads[SIs[Box_index][First]].LeftMostPos ==
															Reads[SIs[Box_index][Second]].LeftMostPos)
														Reads[SIs[Box_index][Second]].Unique = false;
												}

											if (Reads[SIs[Box_index][First]].BPLeft <
													Reads[SIs[Box_index][Second]].BPLeft)
												continue;
											else if (Reads[SIs[Box_index][First]].BPLeft >
															 Reads[SIs[Box_index][Second]].BPLeft)
												{
													CompareResult = 1;
												}
											else
												{
													if (Reads[SIs[Box_index][First]].IndelSize <
															Reads[SIs[Box_index][Second]].IndelSize)
														continue;
													else if (Reads[SIs[Box_index][First]].IndelSize >
																	 Reads[SIs[Box_index][Second]].IndelSize)
														{
															CompareResult = 1;
														}
													//else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
													//    short Compare2Str = CompareTwoString(InputIndels[First].InsertedStr, InputIndels[Second].InsertedStr );
													//if (Compare2Str > 0) CompareResult = 1;
													//else if (Compare2Str == 0) CompareResult = 2;
													//else continue;

													//}
													//else if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
													//  if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
													//    InputIndels[Second].Unique = false;
													//  }
													//}
												}
											//CompareResult = CompareTwoReads(InputIndels[First], InputIndels[Second]);
											if (CompareResult == 1)
												{
													Temp4Exchange = SIs[Box_index][First];
													SIs[Box_index][First] = SIs[Box_index][Second];
													SIs[Box_index][Second] = Temp4Exchange;
												}
											//else if (CompareResult == 2) InputIndels[Second].Unique = false;
										}
									}
							}
						}
					GoodIndels.clear ();
					IndelEvents.clear ();
					LOG_DEBUG(std::cout << GoodIndels.size() << std::endl);

					for (unsigned int First = 0; First < SIsNum; First++)
						{
							//if (InputIndels[First].Unique) 
							GoodIndels.push_back (Reads[SIs[Box_index][First]]);
						}

					GoodNum = GoodIndels.size ();
					LOG_DEBUG(std::cout << Box_index << " " << GoodNum << std::endl);
					if (GoodNum == 0)
						continue;
					Indel4output OneIndelEvent;
					OneIndelEvent.Start = 0;
					OneIndelEvent.End = 0;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
					OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
					OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
					OneIndelEvent.BPRight = GoodIndels[0].BPRight;
					//OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
					OneIndelEvent.WhetherReport = true;
					for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++)
						{
							if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
									&& GoodIndels[GoodIndex].IndelSize ==
									OneIndelEvent.IndelSize)
								//&& OneIndelEvent.IndelStr == GoodIndels[GoodIndex].InsertedStr )

								OneIndelEvent.End = GoodIndex;
							else
								{
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
									OneIndelEvent.IndelStr = GoodIndels[GoodIndex].InsertedStr;
									//OneIndelEvent.IndelStr = GoodIndels[GoodIndex].InsertedStr;
								}
						}

					//if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff) 
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;

					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Insertion (CurrentChr, OneIndelEvent.IndelStr,
																	OneIndelEvent.RealStart,
																	OneIndelEvent.RealEnd);
					IndelEvents.push_back (OneIndelEvent);

					if (IndelEvents.size ())
						{
							for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
									 EventIndex++)
								{
									if (IndelEvents[EventIndex].WhetherReport)
										{
											unsigned int RealStart = IndelEvents[EventIndex].RealStart;
											unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
											unsigned int IndelSize = IndelEvents[EventIndex].IndelSize;
											unsigned int Max_Support = IndelEvents[EventIndex].Support;
											unsigned int Max_Support_Index = EventIndex;

											for (unsigned EventIndex_left = 0;
													 EventIndex_left < IndelEvents.size ();
													 EventIndex_left++)
												{
													if (IndelEvents[EventIndex_left].WhetherReport ==
															false)
														continue;
													else if (IndelEvents[EventIndex_left].RealStart !=
																	 RealStart)
														continue;
													else if (IndelEvents[EventIndex_left].RealEnd !=
																	 RealEnd)
														continue;
													else if (IndelEvents[EventIndex_left].IndelSize !=
																	 IndelSize)
														continue;
													else
														{
															IndelEvents[EventIndex_left].WhetherReport =
																false;
															if (IndelEvents[EventIndex_left].Support >
																	Max_Support)
																{
																	Max_Support =
																		IndelEvents[EventIndex_left].Support;
																	Max_Support_Index = EventIndex_left;
																}
														}
												}
											// report max one
											LOG_DEBUG(std::cout << Max_Support << std::endl);
											if (Max_Support >= NumRead2ReportCutOff)
												{
													OutputSIs (GoodIndels, CurrentChr,
																		 IndelEvents[Max_Support_Index].Start,
																		 IndelEvents[Max_Support_Index].End,
																		 RealStart, RealEnd, SIsOutf);
													NumberOfSIsInstances++;
												}
										}
								}
						}
				}												// if (!insertion[Box_index].empty())
		}
	LOG_INFO(std::cout << "Short insertions: " << NumberOfSIsInstances << std::endl << std::endl);
}

void
SortAndOutputTandemDuplications (const unsigned &NumBoxes, const std::string & CurrentChr,
								 std::vector < SPLIT_READ > &AllReads, std::vector < unsigned >TDs[],
								 std::ofstream & TDOutf, const bool nonTemplate)
{

  if (nonTemplate)
	{
	  LOG_INFO(std::cout << "Sorting and outputing tandem duplications with non-template sequence ..." << std::endl);
	}
  else
	{
	  LOG_INFO(std::cout << "Sorting and outputing tandem duplications ..." << std::endl);
	}
  unsigned int TDNum;
  short CompareResult;
  unsigned Temp4Exchange;
  int countTandemDuplications = 0;
  unsigned int GoodNum;
  std::vector < SPLIT_READ > GoodIndels;
  std::vector < Indel4output > IndelEvents;

  for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
	{
	  if (TDs[Box_index].size () >= NumRead2ReportCutOff)
		{
		  TDNum = TDs[Box_index].size ();

		  for (unsigned int First = 0; First < TDNum - 1; First++)
			{
			  {
				for (unsigned int Second = First + 1; Second < TDNum;
					Second++)
				  {
					{
					  if (AllReads[TDs[Box_index][First]].ReadLength ==
						  AllReads[TDs[Box_index][Second]].ReadLength)
						{
						  if (AllReads[TDs[Box_index][First]].LeftMostPos ==
							  AllReads[TDs[Box_index][Second]].LeftMostPos)
							AllReads[TDs[Box_index][Second]].Unique = false;
						}
					  if (AllReads[TDs[Box_index][First]].BPLeft <
						  AllReads[TDs[Box_index][Second]].BPLeft)
						continue;
					  else if (AllReads[TDs[Box_index][First]].BPLeft >
							   AllReads[TDs[Box_index][Second]].BPLeft)
						{
						  CompareResult = 1;
						}
					  else if (AllReads[TDs[Box_index][First]].BPLeft ==
		 AllReads[TDs[Box_index][Second]].BPLeft)
						{
						  if (AllReads[TDs[Box_index][First]].BPRight <
							  AllReads[TDs[Box_index][Second]].BPRight)
							continue;
						  else if (AllReads[TDs[Box_index][First]].BPRight >
								   AllReads[TDs[Box_index][Second]].BPRight)
							{
							  CompareResult = 1;
							}
						  else if (nonTemplate)
							{ // InputIndels[First].BPRight == InputIndels[Second].BPRight
							  if (AllReads[TDs[Box_index][First]].NT_size <
								  AllReads[TDs[Box_index][Second]].NT_size)
								continue;
							  else if (AllReads[TDs[Box_index][First]].
									   NT_size >
									   AllReads[TDs[Box_index][Second]].
									   NT_size)
								CompareResult = 1;
							}
						}
					  if (CompareResult == 1)
						{
						  Temp4Exchange = TDs[Box_index][First];
						  TDs[Box_index][First] = TDs[Box_index][Second];
						  TDs[Box_index][Second] = Temp4Exchange;
						}
					}
				  }
			  }
			}
		  GoodIndels.clear ();
		  IndelEvents.clear ();

		  for (unsigned int First = 0; First < TDNum; First++)
			{
			  GoodIndels.push_back (AllReads[TDs[Box_index][First]]);
			}

		  GoodNum = GoodIndels.size ();
		  if (GoodNum == 0)
			continue;
		  Indel4output OneIndelEvent;
		  OneIndelEvent.Start = 0;
		  OneIndelEvent.End = 0;
		  OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
		  OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
		  OneIndelEvent.BPRight = GoodIndels[0].BPRight;
		  OneIndelEvent.WhetherReport = true;
		  for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++)
			{
			  if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
				  && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
				OneIndelEvent.End = GoodIndex;
			  else
				{
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

		  if (IndelEvents.size ())
			{
			  for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
				  EventIndex++)
				{
				  if (IndelEvents[EventIndex].WhetherReport)
					{
					  unsigned int RealStart = IndelEvents[EventIndex].RealStart;
					  unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
					  unsigned int Max_Support = IndelEvents[EventIndex].Support;
					  unsigned int Max_Support_Index = EventIndex;
					  for (unsigned EventIndex_left = 0;
						  EventIndex_left < IndelEvents.size ();
						  EventIndex_left++)
						{
						  if (IndelEvents[EventIndex_left].WhetherReport ==
							  false)
							continue;
						  else if (IndelEvents[EventIndex_left].RealStart !=
								   RealStart)
							continue;
						  else if (IndelEvents[EventIndex_left].RealEnd !=
								   RealEnd)
							continue;
						  else
							{
							  IndelEvents[EventIndex_left].WhetherReport =
								  false;
							  if (IndelEvents[EventIndex_left].Support >
								  Max_Support)
								{
								  Max_Support =
									  IndelEvents[EventIndex_left].Support;
								  Max_Support_Index = EventIndex_left;
								}
							}
						}
					  // report max one
					  if (IndelEvents[Max_Support_Index].Support >=
						  NumRead2ReportCutOff)
						{
						  if (GoodIndels
							  [IndelEvents[Max_Support_Index].Start].
							  IndelSize < BalanceCutoff)
							{
							  OutputTDs (GoodIndels, CurrentChr,
										 IndelEvents[Max_Support_Index].Start,
										 IndelEvents[Max_Support_Index].End,
										 RealStart, RealEnd, TDOutf);
							  NumberOfTDInstances++;
							  countTandemDuplications++;
							}
						  else
							if (ReportEvent
								(GoodIndels,
								IndelEvents[Max_Support_Index].Start,
								IndelEvents[Max_Support_Index].End))
							{
							  OutputTDs (GoodIndels, CurrentChr,
										 IndelEvents[Max_Support_Index].Start,
										 IndelEvents[Max_Support_Index].End,
										 RealStart, RealEnd, TDOutf);
							  NumberOfTDInstances++;
							  countTandemDuplications++;
							}
						}
					}
				}
			}
		} // if (!Deletions[Box_index].empty())
	} // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
  if (nonTemplate)
	{
	  LOG_INFO(std::cout << "Tandem duplications with non-template sequence (TD_NT): " <<
						 countTandemDuplications << std::endl
						 << std::endl);
	}
  else
	{
	  LOG_INFO(std::cout << "Tandem duplications: " << NumberOfTDInstances << std::endl
						 << std::endl);
	}
}


void
SortOutputD (const unsigned &NumBoxes, const std::string & CurrentChr,
						 std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Deletions[],
						 std::ofstream & DeletionOutf)
{
	LOG_INFO(std::cout << "Sorting and outputing deletions ..." << std::endl);
	unsigned int DeletionsNum;
	short CompareResult;
	unsigned Temp4Exchange;

	unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
	std::vector < SPLIT_READ > GoodIndels;
	std::vector < Indel4output > IndelEvents;

	for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
		{
			LOG_DEBUG(std::cout << Box_index << "\t" << NumBoxes << "\t" << Deletions[Box_index].size() << std::endl);
			//if (Deletions[Box_index].size() >= NumRead2ReportCutOff)
			//   cout << Box_index << "\t" << Deletions[Box_index].size() << endl;
			if (Deletions[Box_index].size () >= NumRead2ReportCutOff)
				{
					DeletionsNum = Deletions[Box_index].size ();

					for (unsigned int First = 0; First < DeletionsNum - 1; First++)
						{
							//if (InputIndels[First].Unique) 
							//cout << DeletionsNum << " First " << First << endl;
							{
								for (unsigned int Second = First + 1; Second < DeletionsNum;
										 Second++)
									{
										//if (InputIndels[Second].Unique) 
										//cout << DeletionsNum << " Second " << Second << endl;
										{
											if (Reads[Deletions[Box_index][First]].ReadLength ==
													Reads[Deletions[Box_index][Second]].ReadLength)
												{
													if (Reads[Deletions[Box_index][First]].
															LeftMostPos ==
															Reads[Deletions[Box_index][Second]].LeftMostPos)
														Reads[Deletions[Box_index][Second]].Unique =
															false;
												}
											if (Reads[Deletions[Box_index][First]].BPLeft <
													Reads[Deletions[Box_index][Second]].BPLeft)
												continue;
											else if (Reads[Deletions[Box_index][First]].BPLeft >
															 Reads[Deletions[Box_index][Second]].BPLeft)
												{
													CompareResult = 1;
												}
											else if (Reads[Deletions[Box_index][First]].BPLeft ==
															 Reads[Deletions[Box_index][Second]].BPLeft)
												{
													if (Reads[Deletions[Box_index][First]].BPRight <
															Reads[Deletions[Box_index][Second]].BPRight)
														continue;
													else if (Reads[Deletions[Box_index][First]].
																	 BPRight >
																	 Reads[Deletions[Box_index][Second]].
																	 BPRight)
														{
															CompareResult = 1;
														}
													//else CompareResult = 2;
													//else {
													//  if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
													//    if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
													//      InputIndels[Second].Unique = false;
													//    }
													//      
													//  }
													//}
												}
											if (CompareResult == 1)
												{
													Temp4Exchange = Deletions[Box_index][First];
													Deletions[Box_index][First] =
														Deletions[Box_index][Second];
													Deletions[Box_index][Second] = Temp4Exchange;
												}
											//else if (CompareResult == 2) {
											//  Temp4Exchange = InputIndels[First + 1];
											//   InputIndels[First + 1] = InputIndels[Second];
											//   InputIndels[Second] = Temp4Exchange;
											//}
										}
									}
							}
						}
					GoodIndels.clear ();
					IndelEvents.clear ();
					//for (int i = 0 ; i < DeletionsNum; i ++) {
					//  InputIndels.push_back(Reads[Deletions[Box_index][i]]);
					//}Reads[Deletions[Box_index][First]]
					//cout << "GoodIndels" << endl;
					for (unsigned int First = 0; First < DeletionsNum; First++)
						{
							//if (InputIndels[First].Unique) 
							GoodIndels.push_back (Reads[Deletions[Box_index][First]]);
						}

					GoodNum = GoodIndels.size ();
					LOG_DEBUG(std::cout << Box_index << " box read size " << GoodNum << std::endl);
					if (GoodNum == 0) continue;
					//    cout << GoodNum << endl;
					Indel4output OneIndelEvent;
					OneIndelEvent.Start = 0;
					OneIndelEvent.End = 0;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
					OneIndelEvent.BPRight = GoodIndels[0].BPRight;
					OneIndelEvent.WhetherReport = true;
					LOG_DEBUG(std::cout << "here" << std::endl);
					for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++)
						{
							if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
									&& GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
								OneIndelEvent.End = GoodIndex;
							else
								{
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
					LOG_DEBUG(std::cout << "IndelEvent: " << IndelEvents.size() << std::endl);

					if (IndelEvents.size ())
						{
							for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
									 EventIndex++)
								{
									LOG_DEBUG(std::cout << EventIndex << " EventIndex" << std::endl);
									if (IndelEvents[EventIndex].WhetherReport)
										{
											unsigned int RealStart = IndelEvents[EventIndex].RealStart;
											unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
											unsigned int Max_Support = IndelEvents[EventIndex].Support;
											unsigned int Max_Support_Index = EventIndex;
											for (unsigned EventIndex_left = 0;
													 EventIndex_left < IndelEvents.size ();
													 EventIndex_left++)
												{
													if (IndelEvents[EventIndex_left].WhetherReport ==
															false)
														continue;
													else if (IndelEvents[EventIndex_left].RealStart !=
																	 RealStart)
														continue;
													else if (IndelEvents[EventIndex_left].RealEnd !=
																	 RealEnd)
														continue;
													else
														{
															IndelEvents[EventIndex_left].WhetherReport =
																false;
															if (IndelEvents[EventIndex_left].Support >
																	Max_Support)
																{
																	Max_Support =
																		IndelEvents[EventIndex_left].Support;
																	Max_Support_Index = EventIndex_left;
																}
														}
												}
											// report max one
											LOG_DEBUG(std::cout << "max" << std::endl);
											if (IndelEvents[Max_Support_Index].Support >=
													NumRead2ReportCutOff)
												{
													LOG_DEBUG(std::cout << "aa" << std::endl);
													if (GoodIndels
															[IndelEvents[Max_Support_Index].Start].
															IndelSize < BalanceCutoff)
														{
															LOG_DEBUG(std::cout << "ba" << std::endl);
															OutputDeletions (GoodIndels, CurrentChr,
																							 IndelEvents[Max_Support_Index].
																							 Start,
																							 IndelEvents[Max_Support_Index].
																							 End, RealStart, RealEnd,
																							 DeletionOutf);
															deletionFileData.increaseTemplateSvCounter();
															LOG_DEBUG(std::cout << "bb" << std::endl);
														}
													else
														if (ReportEvent
																(GoodIndels,
																 IndelEvents[Max_Support_Index].Start,
																 IndelEvents[Max_Support_Index].End))
														{
															LOG_DEBUG(std::cout << "ca" << std::endl);
															OutputDeletions (GoodIndels, CurrentChr,
																							 IndelEvents[Max_Support_Index].
																							 Start,
																							 IndelEvents[Max_Support_Index].
																							 End, RealStart, RealEnd,
																							 DeletionOutf);
															deletionFileData.increaseTemplateSvCounter();
															LOG_DEBUG(std::cout << "cb" << std::endl);
														}
												}
										}
								}
						}
				}												// if (!Deletions[Box_index].empty())
		}														// for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
	LOG_INFO(std::cout << "Deletions: " << deletionFileData.getTemplateSvCounter() << std::endl << std::endl);
}

void
SortOutputInv (const unsigned &NumBoxes, const std::string & CurrentChr,
							 std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Inv[],
							 std::ofstream & InvOutf)
{
  OutputSorter os(NumBoxes, CurrentChr, InvOutf);
  os.SortAndOutputInversions(Reads, Inv);
}

void
SortOutputInv_NT (const unsigned &NumBoxes, const std::string & CurrentChr,
									std::vector < SPLIT_READ > &Reads, std::vector < unsigned >Inv[],
									std::ofstream & InvOutf)
{
  OutputSorter os(NumBoxes, CurrentChr, InvOutf);
  os.SortAndOutputNonTemplateInversions(Reads, Inv);
}

void
SortOutputDI (const unsigned &NumBoxes, const std::string & CurrentChr,
							std::vector < SPLIT_READ > &Reads, std::vector < unsigned >DI[],
							std::ofstream & DIOutf)
{
	LOG_INFO(std::cout << "Sorting and outputing deletions with non-template sequences ..." <<
					 std::endl);
	unsigned int DINum;
	short CompareResult;
	unsigned Temp4Exchange;
	/*
	   unsigned int C_S = 0;
	   unsigned int C_E = 0;
	   unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	   unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	   unsigned int C_Indelsize;
	 */
	unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
	std::vector < SPLIT_READ > GoodIndels;
	std::vector < Indel4output > IndelEvents;

	for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
		{
			LOG_DEBUG(std::cout << "Box_index: "   << Box_index << std::endl);
			if (DI[Box_index].size () >= NumRead2ReportCutOff)
				{
					DINum = DI[Box_index].size ();
					//InputIndels.clear();
					//for (int i = 0; i < DINum; i++) InputIndels.push_back(Reads[DI[Box_index][i]]);
					for (unsigned int First = 0; First < DINum - 1; First++)
						{
							//if (InputIndels[First].Unique) 
							{
								for (unsigned int Second = First + 1; Second < DINum;
										 Second++)
									{
										//if (InputIndels[Second].Unique) 
										{
											if (Reads[DI[Box_index][First]].ReadLength ==
													Reads[DI[Box_index][Second]].ReadLength)
												{
													if (Reads[DI[Box_index][First]].LeftMostPos ==
															Reads[DI[Box_index][Second]].LeftMostPos)
														Reads[DI[Box_index][Second]].Unique = false;
												}
											if (Reads[DI[Box_index][First]].BPLeft <
													Reads[DI[Box_index][Second]].BPLeft)
												continue;
											else if (Reads[DI[Box_index][First]].BPLeft >
															 Reads[DI[Box_index][Second]].BPLeft)
												{
													CompareResult = 1;
												}
											else if (Reads[DI[Box_index][First]].BPLeft ==
															 Reads[DI[Box_index][Second]].BPLeft)
												{
													if (Reads[DI[Box_index][First]].BPRight
															< Reads[DI[Box_index][Second]].BPRight)
														continue;
													else if (Reads[DI[Box_index][First]].BPRight >
																	 Reads[DI[Box_index][Second]].BPRight)
														{
															CompareResult = 1;
														}
													else
														{
															if (Reads[DI[Box_index][First]].NT_size <
																	Reads[DI[Box_index][Second]].NT_size)
																continue;
															else if (Reads[DI[Box_index][First]].NT_size >
																			 Reads[DI[Box_index][Second]].NT_size)
																CompareResult = 1;
															//else CompareResult = 2;
															//else {
															//  if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
															//    if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
															//      InputIndels[Second].Unique = false;
															//    }
															//      
															//  }
															//}
														}
												}
											if (CompareResult == 1)
												{
													Temp4Exchange = DI[Box_index][First];
													DI[Box_index][First] = DI[Box_index][Second];
													DI[Box_index][Second] = Temp4Exchange;
												}
											//else if (CompareResult == 2) {
											//  Temp4Exchange = InputIndels[First + 1];
											//   InputIndels[First + 1] = InputIndels[Second];
											//   InputIndels[Second] = Temp4Exchange;
											//}
										}
									}
							}
						}
					GoodIndels.clear ();
					IndelEvents.clear ();

					for (unsigned int First = 0; First < DINum; First++)
						{
							//if (InputIndels[First].Unique) 
							GoodIndels.push_back (Reads[DI[Box_index][First]]);
						}

					GoodNum = GoodIndels.size ();
					LOG_DEBUG(std::cout << Box_index << " " << GoodNum << std::endl);
					if (GoodNum == 0)
						continue;
					LOG_DEBUG(std::cout << "GoodNum: " << GoodNum << std::endl);   
					//string InsertedStr;   
					//string NT_str;  
					//short NT_size;
					Indel4output OneIndelEvent;
					OneIndelEvent.Start = 0;
					OneIndelEvent.End = 0;
					OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
					OneIndelEvent.NT_size = GoodIndels[0].NT_size;
					OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
					OneIndelEvent.BPRight = GoodIndels[0].BPRight;
					//OneIndelEvent.IndelStr = GoodIndels[0].NT_str;
					OneIndelEvent.WhetherReport = true;
					for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++)
						{
							if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
									&& GoodIndels[GoodIndex].IndelSize ==
									OneIndelEvent.IndelSize
									&& GoodIndels[GoodIndex].NT_size == OneIndelEvent.NT_size)
								//&& OneIndelEvent.IndelStr == GoodIndels[GoodIndex].NT_str)
								OneIndelEvent.End = GoodIndex;
							else
								{
									IndelEvents.push_back (OneIndelEvent);
									OneIndelEvent.Start = GoodIndex;
									OneIndelEvent.End = GoodIndex;
									OneIndelEvent.BPLeft = GoodIndels[GoodIndex].BPLeft;
									OneIndelEvent.IndelSize = GoodIndels[GoodIndex].IndelSize;
									OneIndelEvent.NT_size = GoodIndels[GoodIndex].NT_size;
									//OneIndelEvent.IndelStr = GoodIndels[GoodIndex].NT_str;
								}
						}

					//if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff) 
					IndelEvents.push_back (OneIndelEvent);
					unsigned int RealStart;
					unsigned int RealEnd;
					for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
							 EventIndex++)
						{
							if (IndelEvents[EventIndex].End -
									IndelEvents[EventIndex].Start + 1 >= NumRead2ReportCutOff)
								{
									RealStart =
										GoodIndels[IndelEvents[EventIndex].Start].BPLeft;
									RealEnd = GoodIndels[IndelEvents[EventIndex].Start].BPRight;
									//if (IndelEvents[EventIndex].IndelSize < 100) {}
									//if (IndelSize < BalanceCutoff) {
									//        OutputDI(GoodIndels, CurrentChr, 
									//                 IndelEvents[EventIndex].Start, 
									//                IndelEvents[EventIndex].End, 
									//    RealStart, RealEnd, DIOutf);
									//        NumberOfDIInstances++;
									//}
									if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize <
											BalanceCutoff)
										{
											OutputDI (GoodIndels, CurrentChr,
																IndelEvents[EventIndex].Start,
																IndelEvents[EventIndex].End,
																RealStart, RealEnd, DIOutf);
											deletionFileData.increaseNonTemplateSvCounter();
										}
									else
										if (ReportEvent
												(GoodIndels, IndelEvents[EventIndex].Start,
												 IndelEvents[EventIndex].End))
										{
											OutputDI (GoodIndels, CurrentChr,
																IndelEvents[EventIndex].Start,
																IndelEvents[EventIndex].End,
																RealStart, RealEnd, DIOutf);
											deletionFileData.increaseNonTemplateSvCounter();
										}
								}
						}
				}												// if (!Deletions[Box_index].empty())
		}														// for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
	LOG_INFO(std::cout << "deletions with non-template sequences: " << deletionFileData.getNonTemplateSvCounter() <<
					 std::endl << std::endl);
}

void
SortOutputLI (const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads,
							std::ofstream & LargeInsertionOutf)
{
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	// find LI combinations
	uint8_t *plus_LI_Pos = new uint8_t[CurrentChr.size () + 1];
	uint8_t *minus_LI_Pos = new uint8_t[CurrentChr.size () + 1];
	int32_t *EventIndex_Pos = new int32_t[CurrentChr.size () + 1];
	for (unsigned i = 0; i < CurrentChr.size () + 1; i++)
		{
			plus_LI_Pos[i] = 0;
			minus_LI_Pos[i] = 0;
			EventIndex_Pos[i] = -1;
		}
	for (unsigned Index = 0; Index < Reads.size (); Index++)
		{
			//UP_Close_index = ;
			if (Reads[Index].Found || Reads[Index].Used
					|| !Reads[Index].UP_Far.empty ())
				continue;
			temp_AbsLoc =
				Reads[Index].UP_Close[Reads[Index].UP_Close.size () - 1].AbsLoc;
			if (plus_LI_Pos[temp_AbsLoc] < Max_short)
				{
					if (Reads[Index].MatchedD == Plus)
						plus_LI_Pos[temp_AbsLoc]++;
				}
			if (minus_LI_Pos[temp_AbsLoc] < Max_short)
				{
					if (Reads[Index].MatchedD == Minus)
						minus_LI_Pos[temp_AbsLoc]++;
				}
		}
	std::vector < LI_Pos > LI_Positions;
	LI_Pos temp_LI_pos;
	bool SkipThis;
	int LI_Positions_Size = 0;
	bool SkipPlus;
	for (unsigned int Index_Minus = SpacerBeforeAfter;
			 Index_Minus < CurrentChr.size () - SpacerBeforeAfter; Index_Minus++)
		{
			SkipPlus = false;
			for (unsigned int MaskedPosIndexMinus = Index_Minus + 10;
					 MaskedPosIndexMinus >= Index_Minus - 10; MaskedPosIndexMinus--)
				{
					if (CurrentChrMask[MaskedPosIndexMinus] == 'B')
						{
							Index_Minus = MaskedPosIndexMinus + 10;
							SkipPlus = true;
							break;
						}
				}
			if (SkipPlus == false
					&& minus_LI_Pos[Index_Minus] >= NumRead2ReportCutOff)
				{
					for (unsigned int Index_Plus = Index_Minus - 1;
							 Index_Plus <= Index_Minus + 30; Index_Plus++)
						{
							SkipThis = false;
							for (unsigned int MaskedPosIndexPlus = Index_Plus + 10;
									 MaskedPosIndexPlus >= Index_Plus - 10;
									 MaskedPosIndexPlus--)
								{
									if (CurrentChrMask[MaskedPosIndexPlus] == 'B')
										{
											if ((MaskedPosIndexPlus + 10) > Index_Minus)
												{
													Index_Minus = MaskedPosIndexPlus + 10;
												}
											SkipThis = true;
											break;
										}
								}
							if (SkipThis == false
									&& plus_LI_Pos[Index_Plus] >= NumRead2ReportCutOff)
								{
									temp_LI_pos.Plus_Pos = Index_Plus;
									temp_LI_pos.Minus_Pos = Index_Minus;
									// The following two were changed to assignents from comparisons
									// IS THAT CORRECT? (RWWH 20110523)
									// TODO Kai check if the CurrentChrMask really should be set
									// to 'B' here.
									// CurrentChrMask[Index_Plus] = 'B';
									// CurrentChrMask[Index_Minus] = 'B';
									//cout << Index_Plus << "\t" << (short)plus_LI_Pos[Index_Plus] << "\t" 
									//     << Index_Minus << "\t" << (short)minus_LI_Pos[Index_Minus] << endl;
									temp_LI_pos.WhetherReport = false;
									LI_Positions.push_back (temp_LI_pos);
									EventIndex_Pos[Index_Plus] = LI_Positions_Size;
									EventIndex_Pos[Index_Minus] = LI_Positions_Size;
									LI_Positions_Size++;
									//Index_Minus += 30;
								}
						}
				}
		}
	LOG_DEBUG(std::cout << "LI: " << LI_Positions.size() << std::endl);
	static int Count_LI = 0;
	// find LI supporting reads


	for (unsigned Index = 0; Index < Reads.size (); Index++)
		{
			//UP_Close_index = Reads[Index].UP_Close.size() - 1;
			if (Reads[Index].Used || !Reads[Index].UP_Far.empty ())
				continue;
			temp_AbsLoc =
				Reads[Index].UP_Close[Reads[Index].UP_Close.size () - 1].AbsLoc;
			if (EventIndex_Pos[temp_AbsLoc] == -1)
				continue;
			Reads[Index].Used = true;
			if (Reads[Index].MatchedD == Plus)
				{
					LI_Positions[EventIndex_Pos[temp_AbsLoc]].Plus_Reads.
						push_back (Index);
				}
			else
				{
					LI_Positions[EventIndex_Pos[temp_AbsLoc]].Minus_Reads.
						push_back (Index);
				}
		}

	std::vector < SPLIT_READ > temp_Plus_Reads, temp_Minus_Reads;

	bool temp_BalancedPlus_Plus, temp_BalancedPlus_Minus,
		temp_BalancedMinus_Plus, temp_BalancedMinus_Minus;
	short temp_LengthStr;
	for (unsigned LI_index = 0; LI_index < LI_Positions.size (); LI_index++)
		{
			if (LI_Positions[LI_index].Minus_Reads.empty ()
					|| LI_Positions[LI_index].Plus_Reads.empty ())
				continue;
			temp_BalancedPlus_Plus = false;
			temp_BalancedPlus_Minus = false;
			temp_BalancedMinus_Plus = false;
			temp_BalancedMinus_Minus = false;
			temp_Plus_Reads.clear ();
			temp_Minus_Reads.clear ();
			LOG_DEBUG(std::cout << "Here: " << LI_index << 
								"\t" << LI_Positions[LI_index].Minus_Reads.size() << 
								"\t" << LI_Positions[LI_index].Plus_Reads.size() << std::endl);
			for (unsigned int i = 0; i < LI_Positions[LI_index].Minus_Reads.size (); i++)
				{
					temp_Minus_Reads.
						push_back (Reads[LI_Positions[LI_index].Minus_Reads[i]]);
				}
			//CheckConsistancy(temp_Minus_Reads);
			for (unsigned int i = 0; i < LI_Positions[LI_index].Minus_Reads.size (); i++)
				{
					UP_Close_index = temp_Minus_Reads[i].UP_Close.size () - 1;
					temp_LengthStr =
						temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
					if ((float) temp_LengthStr > temp_Minus_Reads[i].ReadLength * 0.5)
						temp_BalancedMinus_Plus = true;
					else if ((float) temp_LengthStr <
									 temp_Minus_Reads[i].ReadLength * 0.5)
						temp_BalancedMinus_Minus = true;
				}
			for (unsigned int i = 0; i < LI_Positions[LI_index].Plus_Reads.size (); i++)
				{
					temp_Plus_Reads.
						push_back (Reads[LI_Positions[LI_index].Plus_Reads[i]]);
				}
			//CheckConsistancy(temp_Plus_Reads);
			for (unsigned int i = 0; i < LI_Positions[LI_index].Plus_Reads.size (); i++)
				{
					UP_Close_index = temp_Plus_Reads[i].UP_Close.size () - 1;
					temp_LengthStr =
						temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
					if ((float) temp_LengthStr > temp_Plus_Reads[i].ReadLength * 0.5)
						temp_BalancedPlus_Plus = true;
					else if ((float) temp_LengthStr <
									 temp_Plus_Reads[i].ReadLength * 0.5)
						temp_BalancedPlus_Minus = true;
				}

			unsigned int NumSupportPerTagPlus[VectorTag.size ()];
			unsigned int NumSupportPerTagMinus[VectorTag.size ()];
			for (unsigned short i = 0; i < VectorTag.size (); i++)
				{
					NumSupportPerTagPlus[i] = 0;
					NumSupportPerTagMinus[i] = 0;
				}
			for (unsigned int i = 0; i < temp_Minus_Reads.size (); i++)
				{
					for (unsigned short j = 0; j < VectorTag.size (); j++)
						{
							if (temp_Minus_Reads[i].Tag == VectorTag[j])
								{
									NumSupportPerTagMinus[j]++;
									break;
								}
						}
				}
			for (unsigned int i = 0; i < temp_Plus_Reads.size (); i++)
				{
					for (unsigned short j = 0; j < VectorTag.size (); j++)
						{
							if (temp_Plus_Reads[i].Tag == VectorTag[j])
								{
									NumSupportPerTagPlus[j]++;
									break;
								}
						}
				}
			bool SupportedByOneSample = false;
			for (unsigned short j = 0; j < VectorTag.size (); j++)
				{
					LOG_DEBUG(std::cout << NumSupportPerTagPlus[j] << "\t" << NumSupportPerTagMinus[j] << std::endl);
					if (NumSupportPerTagPlus[j] > 0 && NumSupportPerTagMinus[j] > 0)
						{
							SupportedByOneSample = true;
							break;
						}
				}

			short PositiveBool = 0;
			if (temp_BalancedPlus_Plus)
				PositiveBool++;
			if (temp_BalancedPlus_Minus)
				PositiveBool++;
			if (temp_BalancedMinus_Plus)
				PositiveBool++;
			if (temp_BalancedMinus_Minus)
				PositiveBool++;

			if (SupportedByOneSample && PositiveBool >= 3)
				{
					//if (LI_Positions[LI_index].WhetherReport) 
					{
						LargeInsertionOutf <<
							"########################################################" <<
							std::endl;
						LargeInsertionOutf << Count_LI++ << "\tLI\tChrID " <<
							temp_Plus_Reads[0].FragName << "\t" << LI_Positions[LI_index].
							Plus_Pos - SpacerBeforeAfter +
							1 << "\t" << temp_Plus_Reads.
							size () << "\t" << LI_Positions[LI_index].Minus_Pos -
							SpacerBeforeAfter +
							1 << "\t" << temp_Minus_Reads.size () << std::endl;

						LargeInsertionOutf << (CurrentChr.
																	 substr (LI_Positions[LI_index].Plus_Pos -
																					 ReportLength + 1,
																					 ReportLength)) <<
							Cap2Low (CurrentChr.
											 substr (LI_Positions[LI_index].Plus_Pos + 1,
															 ReportLength)) << std::endl;
						for (unsigned int i = 0; i < temp_Plus_Reads.size (); i++)
							{
								UP_Close_index = temp_Plus_Reads[i].UP_Close.size () - 1;
								temp_LengthStr =
									temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
								for (int j = 0; j < ReportLength - temp_LengthStr; j++)
									{
										LargeInsertionOutf << " ";
									}
								LargeInsertionOutf << ReverseComplement (temp_Plus_Reads[i].
																												 UnmatchedSeq) <<
									std::endl;
							}

						LargeInsertionOutf <<
							"--------------------------------------------------------" <<
							std::endl;
						//LargeInsertionOutf << "-\t" << minus_LI_Pos[Index_Minus].NumReads << endl;
						LargeInsertionOutf << Cap2Low (CurrentChr.
																					 substr (LI_Positions[LI_index].
																									 Minus_Pos - ReportLength,
																									 ReportLength)) <<
							(CurrentChr.
							 substr (LI_Positions[LI_index].Minus_Pos,
											 ReportLength)) << std::endl;
						for (unsigned int i = 0; i < temp_Minus_Reads.size (); i++)
							{
								UP_Close_index = temp_Minus_Reads[i].UP_Close.size () - 1;
								temp_LengthStr =
									temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
								for (int j = 0;
										 j <
										 ReportLength + temp_LengthStr -
										 temp_Minus_Reads[i].ReadLength; j++)
									{
										LargeInsertionOutf << " ";
									}
								LargeInsertionOutf << (temp_Minus_Reads[i].
																			 UnmatchedSeq) << std::endl;
							}
					}

				}

			//LI_Positions[LI_index].WhetherReport = true; 
		}
	// output
	delete[]plus_LI_Pos;
	delete[]minus_LI_Pos;
	delete[]EventIndex_Pos;

	LOG_INFO(std::cout << "Breakpoints for large insertions (LI): " << Count_LI << std::endl << std::endl);
}


void
SortOutputRest (const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads,
								std::vector < SPLIT_READ > &BP_Reads, std::ofstream & Outf_Rest)
{
	SPLIT_READ one_BP_read;
	std::string HalfMapped, HalfUnmapped;
	int HalfMappedIndex, HalfUnmappedIndex;
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	LOG_DEBUG(std::cout << "1" << std::endl);
	// find LI combinations
	uint8_t *plus_LI_Pos = new uint8_t[CurrentChr.size () + 1];
	uint8_t *minus_LI_Pos = new uint8_t[CurrentChr.size () + 1];
	for (unsigned i = 0; i < CurrentChr.size () + 1; i++)
		{
			plus_LI_Pos[i] = 0;
			minus_LI_Pos[i] = 0;
		}
	for (unsigned Index = 0; Index < Reads.size (); Index++)
		{
			if (Reads[Index].Found || Reads[Index].Used
					|| !Reads[Index].UP_Far.empty ())
				continue;
			UP_Close_index = Reads[Index].UP_Close.size () - 1;
			temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
			if (plus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP)
				{
					if (Reads[Index].MatchedD == Plus)
						plus_LI_Pos[temp_AbsLoc]++;
				}
			if (minus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP)
				{
					if (Reads[Index].MatchedD == Minus)
						minus_LI_Pos[temp_AbsLoc]++;
				}
		}
	std::vector < Rest_Pos > Rest_Positions;
	Rest_Pos temp_Rest_pos;

	LOG_DEBUG(std::cout << "2" << std::endl);

	bool SkipThisPos;
	for (unsigned int Index = SpacerBeforeAfter;
			 Index < CurrentChr.size () - SpacerBeforeAfter; Index++)
		{
			SkipThisPos = false;
			//for (int MaskIndex = Index + 10; MaskIndex >= Index - 10; MaskIndex--) {
			//  if (CurrentChrMask[MaskIndex] == 'B') {
			//    Index = MaskIndex + 10;
			//    SkipThisPos = true;
			//    break;
			//  }
			//}
			if (SkipThisPos == true)
				continue;
			if (plus_LI_Pos[Index] >= NumRead2ReportCutOff_BP)
				{
					temp_Rest_pos.Strand = Plus;
					temp_Rest_pos.Pos = Index;
					Rest_Positions.push_back (temp_Rest_pos);
				}
			if (minus_LI_Pos[Index] >= NumRead2ReportCutOff_BP)
				{
					temp_Rest_pos.Strand = Minus;
					temp_Rest_pos.Pos = Index;
					Rest_Positions.push_back (temp_Rest_pos);
				}
		}

	// find supporting reads
	for (unsigned Index = 0; Index < Reads.size (); Index++)
		{
			if (Reads[Index].Used || !Reads[Index].UP_Far.empty ())
				continue;
			UP_Close_index = Reads[Index].UP_Close.size () - 1;
			temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
			for (unsigned Pos_index = 0; Pos_index < Rest_Positions.size ();
					 Pos_index++)
				{
					if (Reads[Index].MatchedD == Rest_Positions[Pos_index].Strand)
						{
							if (temp_AbsLoc == Rest_Positions[Pos_index].Pos)
								{
									Reads[Index].Used = true;
									Rest_Positions[Pos_index].Pos_Reads.push_back (Index);	// copy index to save memory
								}
						}
				}
		}

	LOG_DEBUG(std::cout << "Other unassigned breakpoints (BP): " << Rest_Positions.size() << std::endl << std::endl);
	int Count_BP = 0;
	bool temp_BalancedPlus, temp_BalancedMinus;
	short temp_LengthStr;
	std::vector < SPLIT_READ > temp_Pos_Reads;
	for (unsigned LI_index = 0; LI_index < Rest_Positions.size (); LI_index++)
		{
			temp_Pos_Reads.clear ();
			temp_BalancedPlus = false;
			temp_BalancedMinus = false;
			for (unsigned int i = 0; i < Rest_Positions[LI_index].Pos_Reads.size (); i++)
				{
					temp_Pos_Reads.
						push_back (Reads[Rest_Positions[LI_index].Pos_Reads[i]]);
					UP_Close_index =
						Reads[Rest_Positions[LI_index].Pos_Reads[i]].UP_Close.size () - 1;
					temp_LengthStr =
						Reads[Rest_Positions[LI_index].Pos_Reads[i]].
						UP_Close[UP_Close_index].LengthStr;
					if ((float) temp_LengthStr >
							Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength * 0.5)
						temp_BalancedPlus = true;
					else if ((float) temp_LengthStr <
									 Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength *
									 0.5)
						temp_BalancedMinus = true;
				}
			if (temp_BalancedPlus && temp_BalancedMinus)
				{
					Count_BP++;
					if (Rest_Positions[LI_index].Strand == Plus)
						{
							Outf_Rest <<
								"########################################################" <<
								std::endl;
							Outf_Rest << "ChrID " << temp_Pos_Reads[0].
								FragName << "\t" << Rest_Positions[LI_index].Pos -
								SpacerBeforeAfter +
								1 << "\t" << Rest_Positions[LI_index].Pos_Reads.
								size () << "\t+" << std::endl;

							Outf_Rest << (CurrentChr.
														substr (Rest_Positions[LI_index].Pos -
																		ReportLength + 1,
																		ReportLength)) << Cap2Low (CurrentChr.
																															 substr
																															 (Rest_Positions
																																[LI_index].
																																Pos + 1,
																																ReportLength))
								<< std::endl;
							HalfMappedIndex = 0;
							HalfUnmappedIndex = 0;
							for (unsigned int i = 0; i < temp_Pos_Reads.size (); i++)
								{
									UP_Close_index = temp_Pos_Reads[i].UP_Close.size () - 1;
									temp_LengthStr =
										temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
									if (temp_LengthStr >
											temp_Pos_Reads[HalfMappedIndex].
											UP_Close[temp_Pos_Reads[HalfMappedIndex].UP_Close.
															 size () - 1].LengthStr)
										HalfMappedIndex = i;
									if (temp_Pos_Reads[i].ReadLength - temp_LengthStr >
											temp_Pos_Reads[HalfUnmappedIndex].ReadLength -
											temp_Pos_Reads[HalfUnmappedIndex].
											UP_Close[temp_Pos_Reads[HalfUnmappedIndex].UP_Close.
															 size () - 1].LengthStr)
										HalfUnmappedIndex = i;
									for (int j = 0; j < ReportLength - temp_LengthStr; j++)
										{
											Outf_Rest << " ";
										}
									Outf_Rest << ReverseComplement (temp_Pos_Reads[i].
																									UnmatchedSeq) << "\t" <<
										temp_Pos_Reads[i].MatchedD << "\t" << temp_Pos_Reads[i].
										MatchedRelPos << "\t" << temp_Pos_Reads[i].
										MS << "\t" << temp_Pos_Reads[i].
										Tag << "\t" << temp_Pos_Reads[i].Name << std::endl;
									//temp_Pos_Reads[0].
									//BP_Reads.push_back();
								}

						}
					else
						{
							Outf_Rest <<
								"########################################################" <<
								std::endl;
							Outf_Rest << "ChrID " << temp_Pos_Reads[0].
								FragName << "\t" << Rest_Positions[LI_index].Pos -
								SpacerBeforeAfter +
								1 << "\t" << Rest_Positions[LI_index].Pos_Reads.
								size () << "\t-" << std::endl;
							Outf_Rest << Cap2Low (CurrentChr.
																		substr (Rest_Positions[LI_index].Pos -
																						ReportLength,
																						ReportLength)) << (CurrentChr.
																															 substr
																															 (Rest_Positions
																																[LI_index].
																																Pos,
																																ReportLength))
								<< std::endl;
							for (unsigned int i = 0; i < temp_Pos_Reads.size (); i++)
								{
									UP_Close_index = temp_Pos_Reads[i].UP_Close.size () - 1;
									temp_LengthStr =
										temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
									for (int j = 0;
											 j <
											 ReportLength + temp_LengthStr -
											 temp_Pos_Reads[i].ReadLength; j++)
										{
											Outf_Rest << " ";
										}
									Outf_Rest << (temp_Pos_Reads[i].UnmatchedSeq)
										<< "\t" << temp_Pos_Reads[i].MatchedD
										<< "\t" << temp_Pos_Reads[i].MatchedRelPos
										<< "\t" << temp_Pos_Reads[i].MS
										<< "\t" << temp_Pos_Reads[i].Tag
										<< "\t" << temp_Pos_Reads[i].Name << std::endl;
								}
						}
				}
		}
	delete[]plus_LI_Pos;
	delete[]minus_LI_Pos;
	LOG_INFO(std::cout << "Other unassigned breakpoints (BP): " << Count_BP << std::endl << std::endl);
}
