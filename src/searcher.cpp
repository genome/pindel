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

// Pindel header files
#include "pindel.h"
#include "searcher.h"


void
CheckLeft_Close (const SPLIT_READ & OneRead,
								 const std::string & TheInput,
								 const std::string & CurrentReadSeq,
								 const std::vector < unsigned int >Left_PD[],
								 const short &BP_Left_Start,
								 const short &BP_Left_End,
								 const short &CurrentLength, std::vector < UniquePoint > &LeftUP)
{
	int Sum;
	if (CurrentLength >= BP_Left_Start && CurrentLength <= BP_Left_End)
		{
			// put it to LeftUP if unique
			for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++)
				{
					if (Left_PD[i].size () == 1 && CurrentLength >= BP_Left_Start + i)
						{
							Sum = 0;
							if (ADDITIONAL_MISMATCH)
								for (short j = 0; j <= i + ADDITIONAL_MISMATCH; j++)
									Sum += Left_PD[j].size ();

							if (Sum == 1
									&& i <= (short) (CurrentLength * Seq_Error_Rate + 1))
								{
									UniquePoint TempOne;
									TempOne.LengthStr = CurrentLength;
									TempOne.AbsLoc = Left_PD[i][0];
									TempOne.Direction = FORWARD;
									TempOne.Strand = ANTISENSE;
									TempOne.Mismatches = i; 
                                    if (CheckMismatches(TheInput, OneRead.UnmatchedSeq, TempOne)) {
                                        LeftUP.push_back (TempOne);
                                        break;
                                    }
								}
						}
				}
		}
	if (CurrentLength < BP_Left_End)
		{
			std::vector < unsigned int >Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
			for (int CheckedIndex = 0;
					 CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++)
				{
					Left_PD_Output[CheckedIndex].reserve (Left_PD[CheckedIndex].
																								size ());
				}
			const char CurrentChar = CurrentReadSeq[CurrentLength];
			//const int SizeOfCurrent = Left_PD.size();
			unsigned int pos;
			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++)
				{
					int SizeOfCurrent = Left_PD[i].size ();
					if (CurrentChar == 'N')
						{
							//#pragma omp parallel default(shared)
							{
								// #pragma omp for
								for (int j = 0; j < SizeOfCurrent; j++)
									{
										pos = Left_PD[i][j] + 1;
										if (Match2N[(short) TheInput[pos]] == 'N')
											Left_PD_Output[i].push_back (pos);
									}
							}
						}
					else
						{
							//#pragma omp parallel default(shared)
							{
								//#pragma omp for
								for (int j = 0; j < SizeOfCurrent; j++)
									{
										pos = Left_PD[i][j] + 1;
										if (TheInput[pos] == CurrentChar)
											Left_PD_Output[i].push_back (pos);
										else
											Left_PD_Output[i + 1].push_back (pos);
									}
							}
						}
				}

			int SizeOfCurrent =
				Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size ();
			if (CurrentChar == 'N')
				{
					//#pragma omp parallel default(shared)
					{
						// #pragma omp for
						for (int j = 0; j < SizeOfCurrent; j++)
							{
								pos = Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
								if (Match2N[(short) TheInput[pos]] == 'N')
									Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].
										push_back (pos);
							}
					}
				}
			else
				{
					//#pragma omp parallel default(shared)
					{
						// #pragma omp for
						for (int j = 0; j < SizeOfCurrent; j++)
							{
								pos = Left_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] + 1;
								if (TheInput[pos] == CurrentChar)
									Left_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].
										push_back (pos);
								//else Left_PD_Output[i + 1].push_back(pos); 
							}
					}
				}
			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++)
				{
					Sum += Left_PD_Output[i].size ();
				}
			if (Sum)
				{
					const short CurrentLengthOutput = CurrentLength + 1;
					CheckLeft_Close (OneRead, TheInput, CurrentReadSeq, Left_PD_Output,
													 BP_Left_Start, BP_Left_End,
													 CurrentLengthOutput, LeftUP);
				}
			else
				return;
		}
	else
		return;
}

void
CheckRight_Close (const SPLIT_READ & OneRead,
									const std::string & TheInput,
									const std::string & CurrentReadSeq,
									const std::vector < unsigned int >Right_PD[],
									const short &BP_Right_Start,
									const short &BP_Right_End,
									const short &CurrentLength, std::vector < UniquePoint > &RightUP)
{
	//cout << CurrentLength << "\t" << RightUP.size() << "\t" << Right_PD[0].size() << "\t" << Right_PD[1].size() << endl;
	short ReadLengthMinus = CurrentReadSeq.size () - 1;
	int Sum;
	if (CurrentLength >= BP_Right_Start && CurrentLength <= BP_Right_End)
		{
			for (short i = 0; i <= OneRead.MAX_SNP_ERROR; i++)
				{
					if (Right_PD[i].size () == 1 && CurrentLength >= BP_Right_Start + i)
						{
							Sum = 0;
							if (ADDITIONAL_MISMATCH)
								for (short j = 0; j <= i+ ADDITIONAL_MISMATCH; j++)
									Sum += Right_PD[j].size ();

							if (Sum == 1
									&& i <= (short) (CurrentLength * Seq_Error_Rate + 1))
								{
									UniquePoint TempOne;
									TempOne.LengthStr = CurrentLength;
									TempOne.AbsLoc = Right_PD[i][0];
									TempOne.Direction = BACKWARD;
									TempOne.Strand = SENSE;
									TempOne.Mismatches = i;
                                    if (CheckMismatches(TheInput, OneRead.UnmatchedSeq, TempOne)) {
                                        RightUP.push_back (TempOne);
                                        break;
                                    } // ###################################
								}
						}
				}
		}

	if (CurrentLength < BP_Right_End)
		{
			std::vector < unsigned int >Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED];
			for (int CheckedIndex = 0;
					 CheckedIndex < OneRead.TOTAL_SNP_ERROR_CHECKED; CheckedIndex++)
				{
					Right_PD_Output[CheckedIndex].reserve (Right_PD[CheckedIndex].
																								 size ());
				}
			const char CurrentChar =
				CurrentReadSeq[ReadLengthMinus - CurrentLength];
			unsigned int pos;

			for (int i = 0; i < OneRead.TOTAL_SNP_ERROR_CHECKED_Minus; i++)
				{
					int SizeOfCurrent = Right_PD[i].size ();
					if (CurrentChar == 'N')
						{
							//#pragma omp parallel default(shared)
							{
								// #pragma omp for
								for (int j = 0; j < SizeOfCurrent; j++)
									{
										pos = Right_PD[i][j] - 1;
										if (Match2N[(short) TheInput[pos]] == 'N')
											Right_PD_Output[i].push_back (pos);
									}
							}
						}
					else
						{
							//#pragma omp parallel default(shared)
							{
								// #pragma omp for
								for (int j = 0; j < SizeOfCurrent; j++)
									{
										pos = Right_PD[i][j] - 1;
										if (TheInput[pos] == CurrentChar)
											Right_PD_Output[i].push_back (pos);
										else
											Right_PD_Output[i + 1].push_back (pos);
									}
							}
						}
				}

			int SizeOfCurrent =
				Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].size ();
			if (CurrentChar == 'N')
				{
					for (int j = 0; j < SizeOfCurrent; j++)
						{
							pos = Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
							if (Match2N[(short) TheInput[pos]] == 'N')
								Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].
									push_back (pos);
						}
				}
			else
				{
					for (int j = 0; j < SizeOfCurrent; j++)
						{
							pos = Right_PD[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus][j] - 1;
							if (TheInput[pos] == CurrentChar)
								Right_PD_Output[OneRead.TOTAL_SNP_ERROR_CHECKED_Minus].
									push_back (pos);
							//else Left_PD_Output[i + 1].push_back(pos); 
						}
				}

			Sum = 0;
			for (int i = 0; i <= OneRead.MAX_SNP_ERROR; i++)
				{
					Sum += Right_PD_Output[i].size ();
				}
			if (Sum)
				{
					short CurrentLength_output = CurrentLength + 1;
					CheckRight_Close (OneRead, TheInput, CurrentReadSeq,
														Right_PD_Output, BP_Right_Start, BP_Right_End,
														CurrentLength_output, RightUP);
				}
			else
				return;
		}
	else
		return;
}

bool
CheckMismatches (const std::string & TheInput, const std::string & InputReadSeq,
								 //const unsigned int & Start,
								 const UniquePoint & UP)
{
	//return true;  short LengthStr;
	//unsigned int AbsLoc; 
	//cout << "CheckMismatches1" << endl;
	std::string CurrentReadSeq;
	if (UP.Strand == SENSE)
		CurrentReadSeq = InputReadSeq;
	else
		CurrentReadSeq = ReverseComplement (InputReadSeq);
	short CurrentReadLength = CurrentReadSeq.size ();
	unsigned int Start = 0;
	//cout << "CheckMismatches2" << endl;
	std::string BP_On_Read, BP_On_Ref;
	if (UP.Direction == FORWARD)
		{
			//cout << "+s" << endl;
			//cout << "UP.AbsLoc: " << UP.AbsLoc << "\t" << "UP.LengthStr: " << UP.LengthStr << endl;

			Start = UP.AbsLoc - UP.LengthStr + 1;
			//cout << "Start: " << Start << endl;
			//cout << "Min_Perfect_Match_Around_BP: " << Min_Perfect_Match_Around_BP << endl;
			if (UP.LengthStr <= Min_Perfect_Match_Around_BP)
				return false;
			BP_On_Read =
				CurrentReadSeq.substr (UP.LengthStr - Min_Perfect_Match_Around_BP,
															 Min_Perfect_Match_Around_BP);
			//cout << "BP_On_Read: " << BP_On_Read << endl;
			//if (UP.AbsLoc < Min_Perfect_Match_Around_BP) return false;
			BP_On_Ref =
				TheInput.substr (UP.AbsLoc - Min_Perfect_Match_Around_BP + 1,
												 Min_Perfect_Match_Around_BP);
			//cout << "BP_On_Ref: " << BP_On_Ref << endl;
			if (BP_On_Read != BP_On_Ref)
				return false;
			//cout << "+e" << endl;
		}
	else if (UP.Direction == BACKWARD)
		{
			//cout << "-s" << endl;
			Start = UP.AbsLoc + UP.LengthStr - CurrentReadLength;
			if (CurrentReadLength < UP.LengthStr)
				return false;
			BP_On_Read =
				CurrentReadSeq.substr (CurrentReadLength - UP.LengthStr,
															 Min_Perfect_Match_Around_BP);
			BP_On_Ref = TheInput.substr (UP.AbsLoc, Min_Perfect_Match_Around_BP);
			if (BP_On_Read != BP_On_Ref)
				return false;
			//cout << "-e" << endl;
		}
	//cout << "CheckMismatches3" << endl;
	short MAX_ALLOWED_MISMATCHES = (short) (CurrentReadSeq.size () * MaximumAllowedMismatchRate + 1);	// 

	short NumMismatches = 0;			// Match2N[(short)'A'] = 'N';    

	for (short i = 0; i < CurrentReadLength; i++)
		{
			if (CurrentReadSeq[i] == N_char)
				{
					if (Match2N[(short) TheInput[Start + i]] != N_char)
						NumMismatches++;
				}
			else
				{
					if (TheInput[Start + i] != CurrentReadSeq[i])
						NumMismatches++;
				}
		}
	//cout << "CheckMismatches4" << endl;
	if (NumMismatches > MAX_ALLOWED_MISMATCHES)
		return true;
	else
		return false;
}
