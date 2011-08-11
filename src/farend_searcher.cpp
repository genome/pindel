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

#include <math.h>

#include "pindel.h"
#include "searcher.h"

#include "farend_searcher.h"

FarEndSearcher::FarEndSearcher(const std::string & CurrentChr_in, SPLIT_READ & Temp_One_Read_in) {
	CurrentChr = &CurrentChr_in;
	Temp_One_Read = &Temp_One_Read_in;
}

void FarEndSearcher::GetFarEnd(const int &in_start, const int &in_end) {
	GetFarEnd_General(in_start, in_end, false, -1);
}

void FarEndSearcher::GetFarEnd_OtherStrand(const short &RangeIndex) {
	GetFarEnd_General(-1, -1, true, RangeIndex);
}

void FarEndSearcher::GetFarEnd_General(const int &in_start, const int &in_end, const bool &UseRangeIndex,
				const short &RangeIndex) {
	short BP_Start = Temp_One_Read->MinClose + (UseRangeIndex ? RangeIndex : 0);
	short BP_End = Temp_One_Read->ReadLengthMinus;
	std::vector<UniquePoint> UP;
	std::vector<unsigned int> PD_Plus[Temp_One_Read->TOTAL_SNP_ERROR_CHECKED];
	std::vector<unsigned int> PD_Minus[Temp_One_Read->TOTAL_SNP_ERROR_CHECKED];

	int Start;
	int End;

	if (UseRangeIndex) {
		if (Temp_One_Read->MatchedD == Plus) {
			Start = Temp_One_Read->MatchedRelPos + SpacerBeforeAfter - Temp_One_Read->ReadLength - 2
							* Temp_One_Read->InsertSize - DSizeArray[RangeIndex];
			End = Temp_One_Read->MatchedRelPos + SpacerBeforeAfter - Temp_One_Read->ReadLength + 3
							* Temp_One_Read->InsertSize + DSizeArray[RangeIndex];
		} else {
			Start = Temp_One_Read->MatchedRelPos + SpacerBeforeAfter - 3 * Temp_One_Read->InsertSize
							+ Temp_One_Read->ReadLength - DSizeArray[RangeIndex];
			End = Temp_One_Read->MatchedRelPos + SpacerBeforeAfter + 2 * Temp_One_Read->InsertSize
							+ Temp_One_Read->ReadLength + DSizeArray[RangeIndex];
		}
	} else {
		Start = in_start;
		End = in_end;
	}

	for (int CheckIndex = 0; CheckIndex < Temp_One_Read->TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(End - Start + 1);
		PD_Minus[CheckIndex].reserve(End - Start + 1);
	}

	char CurrentBase = Temp_One_Read->UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];

	if (UseRangeIndex) {
		if (Temp_One_Read->TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (Temp_One_Read->MatchedD == Plus) {
				if (CurrentBase != 'N') {
					for (int pos = Start; pos < End; pos++) {
						if (CurrentChr->at(pos) == CurrentBase)
							PD_Plus[0].push_back(pos);
						else
							PD_Plus[1].push_back(pos);
					}
				} else {
					for (int pos = Start; pos < End; pos++) {
						if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
							PD_Plus[0].push_back(pos);
						} else {
							PD_Plus[1].push_back(pos);
						}
					}
				}
			} else {
				if (CurrentBase != 'N') {
					for (int pos = Start; pos < End; pos++) {
						if (CurrentChr->at(pos) == CurrentBaseRC)
							PD_Minus[0].push_back(pos);
						else
							PD_Minus[1].push_back(pos);
					}
				} else {
					for (int pos = Start; pos < End; pos++) {
						if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
							PD_Minus[0].push_back(pos);
						} else {
							PD_Minus[1].push_back(pos);
						}
					}
				}
			}
		} else {
			if (Temp_One_Read->MatchedD == Plus) {
				if (CurrentBase != 'N') {
					for (int pos = Start; pos < End; pos++) {
						if (CurrentChr->at(pos) == CurrentBase)
							PD_Plus[0].push_back(pos);
					}
				} else {
					for (int pos = Start; pos < End; pos++) {
						if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
							PD_Plus[0].push_back(pos);
						}
					}
				}
			} else {
				if (CurrentBase != 'N') {
					for (int pos = Start; pos < End; pos++) {
						if (CurrentChr->at(pos) == CurrentBaseRC)
							PD_Minus[0].push_back(pos);
					}
				} else {
					for (int pos = Start; pos < End; pos++) {
						if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
							PD_Minus[0].push_back(pos);
						}
					}
				}
			}
		}
	} else {
		if (Temp_One_Read->TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == CurrentBase)
						PD_Plus[0].push_back(pos);
					else
						PD_Plus[1].push_back(pos);
					if (CurrentChr->at(pos) == CurrentBaseRC)
						PD_Minus[0].push_back(pos);
					else
						PD_Minus[1].push_back(pos);
				}
			} else {
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					} else {
						PD_Plus[1].push_back(pos);
						PD_Minus[1].push_back(pos);
					}
				}
			}
		} else {
			if (CurrentBase != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == CurrentBase)
						PD_Plus[0].push_back(pos);
					if (CurrentChr->at(pos) == CurrentBaseRC)
						PD_Minus[0].push_back(pos);
				}
			} else {
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD_Plus[0].push_back(pos);
						PD_Minus[0].push_back(pos);
					}
				}
			}
		}
	}

	CheckBoth(*Temp_One_Read, *CurrentChr, Temp_One_Read->UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase,
					UP);

	if (!UP.empty() && ((UP[UP.size() - 1].LengthStr + Temp_One_Read->MinClose) < Temp_One_Read->ReadLength)) {
		for (unsigned UP_index = 0; UP_index < UP.size(); UP_index++) {
			if (CheckMismatches(*CurrentChr, Temp_One_Read->UnmatchedSeq, UP[UP_index]))
				Temp_One_Read->UP_Far.push_back(UP[UP_index]);
		}
	}

	UP.clear();
}

void
FarEndSearcher::GetFarEnd_SingleStrandUpStream(const short &RangeIndex) {
	Temp_One_Read->ReadLength = Temp_One_Read->UnmatchedSeq.size();
	Temp_One_Read->ReadLengthMinus = Temp_One_Read->ReadLength - 1;
	std::string CurrentReadSeq;
	std::vector<unsigned int> PD[Temp_One_Read->TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read->TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read->InsertSize * 2 + DSizeArray[RangeIndex]);
	}
	std::vector<UniquePoint> UP;
	int Start, End;
	short BP_Start = Temp_One_Read->MinClose + RangeIndex;
	short BP_End = Temp_One_Read->ReadLengthMinus;

	if (Temp_One_Read->MatchedD == Minus) {
		CurrentReadSeq = Temp_One_Read->UnmatchedSeq;

		Start = Temp_One_Read->UP_Close[0].AbsLoc + Temp_One_Read->UP_Close[0].LengthStr;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read->InsertSize * 2;
		if (End > (int)CurrentChr->size() - (int)SpacerBeforeAfter)
			End = (int)CurrentChr->size() - (int)SpacerBeforeAfter;

		char LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read->TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == LeftChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else {
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		} else {
			if (LeftChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == LeftChar) {
						PD[0].push_back(pos);
					}
				}
			} else {
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		}

		CheckLeft_Far(*Temp_One_Read, *CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);

		if (!UP.empty() && (UP[UP.size() - 1].LengthStr + Temp_One_Read->MinClose < Temp_One_Read->ReadLength)) {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				if (CheckMismatches(*CurrentChr, Temp_One_Read->UnmatchedSeq, UP[LeftUP_index])) {
					Temp_One_Read->Used = false;
					Temp_One_Read->UP_Far.push_back(UP[LeftUP_index]);
				}
			}
		}

		UP.clear();
	} else if (Temp_One_Read->MatchedD == Plus) {
		char RightChar;
		CurrentReadSeq = ReverseComplement(Temp_One_Read->UnmatchedSeq);

		End = Temp_One_Read->UP_Close[0].AbsLoc - Temp_One_Read->UP_Close[0].LengthStr;

		if (End > (int)DSizeArray[RangeIndex] + (int)Temp_One_Read->InsertSize * 2 + (int)SpacerBeforeAfter)
			Start = End - DSizeArray[RangeIndex] - Temp_One_Read->InsertSize * 2;
		else
			Start = SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read->ReadLengthMinus];
		if (Temp_One_Read->TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == RightChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else {
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		} else {
			if (RightChar != 'N') {
				for (int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == RightChar) {
						PD[0].push_back(pos);
					}
				}
			} else {
				for (int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		}

		CheckRight_Far(*Temp_One_Read, *CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);

		if (!UP.empty() && (UP[UP.size() - 1].LengthStr + Temp_One_Read->MinClose < Temp_One_Read->ReadLength)) {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
				if (CheckMismatches(*CurrentChr, Temp_One_Read->UnmatchedSeq, UP[RightUP_index])) {
					Temp_One_Read->Used = false;
					Temp_One_Read->UP_Far.push_back(UP[RightUP_index]);
				}
			}
		}

		UP.clear();
	}
}

void
FarEndSearcher::GetFarEnd_SingleStrandDownStream(const short &RangeIndex) {
	std::string CurrentReadSeq;
	std::vector<unsigned int> PD[Temp_One_Read->TOTAL_SNP_ERROR_CHECKED];
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read->TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD[CheckIndex].reserve(Temp_One_Read->InsertSize * 2 + DSizeArray[RangeIndex]);
	}
	std::vector<UniquePoint> UP;
	unsigned int Start, End;
	short BP_Start = Temp_One_Read->MinClose + RangeIndex;
	short BP_End = Temp_One_Read->ReadLengthMinus;

	if (Temp_One_Read->MatchedD == Minus) {
		char LeftChar;
		CurrentReadSeq = Temp_One_Read->UnmatchedSeq;

		End = Temp_One_Read->UP_Close[0].AbsLoc
				+ Temp_One_Read->UP_Close[0].LengthStr
				- Temp_One_Read->ReadLength;
		if (End > SpacerBeforeAfter + Temp_One_Read->InsertSize * 2 + DSizeArray[RangeIndex])
			Start = End - DSizeArray[RangeIndex] - Temp_One_Read->InsertSize * 2;
		else
			Start = SpacerBeforeAfter;

		if (End > CurrentChr->size() - SpacerBeforeAfter)
			End = CurrentChr->size() - SpacerBeforeAfter;

		LeftChar = CurrentReadSeq[0];
		if (Temp_One_Read->TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (LeftChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == LeftChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		} else {
			if (LeftChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == LeftChar) {
						PD[0].push_back(pos);
					}
				}
			} else {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		}

		CheckLeft_Far(*Temp_One_Read, *CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);

		if (!UP.empty() && (UP[UP.size() - 1].LengthStr + Temp_One_Read->MinClose < Temp_One_Read->ReadLength)) {
			for (unsigned LeftUP_index = 0; LeftUP_index < UP.size(); LeftUP_index++) {
				if (CheckMismatches(*CurrentChr, Temp_One_Read->UnmatchedSeq, UP[LeftUP_index])) {
					Temp_One_Read->Used = false;
					Temp_One_Read->UP_Far.push_back(UP[LeftUP_index]);
				}
			}
		}

		UP.clear();
	} else if (Temp_One_Read->MatchedD == Plus) {
		char RightChar;
		CurrentReadSeq = ReverseComplement(Temp_One_Read->UnmatchedSeq);

		Start = Temp_One_Read->UP_Close[0].AbsLoc
				- Temp_One_Read->UP_Close[0].LengthStr
				+ Temp_One_Read->ReadLength;
		End = Start + DSizeArray[RangeIndex] + Temp_One_Read->InsertSize * 2;
		if (End > CurrentChr->size() - SpacerBeforeAfter)
			End = CurrentChr->size() - SpacerBeforeAfter;

		RightChar = CurrentReadSeq[Temp_One_Read->ReadLengthMinus];
		if (Temp_One_Read->TOTAL_SNP_ERROR_CHECKED_Minus) {
			if (RightChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == RightChar) {
						PD[0].push_back(pos);
					} else
						PD[1].push_back(pos);
				}
			} else {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		} else {
			if (RightChar != 'N') {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (CurrentChr->at(pos) == RightChar) {
						PD[0].push_back(pos);
					}
				}
			} else {
				for (unsigned int pos = Start; pos < End; pos++) {
					if (Match2N[(short) CurrentChr->at(pos)] == 'N') {
						PD[0].push_back(pos);
					}
				}
			}
		}

		CheckRight_Far(*Temp_One_Read, *CurrentChr, CurrentReadSeq, PD, BP_Start, BP_End, FirstBase, UP);

		if (!UP.empty() && (UP[UP.size() - 1].LengthStr + Temp_One_Read->MinClose < Temp_One_Read->ReadLength)) {
			for (unsigned RightUP_index = 0; RightUP_index < UP.size(); RightUP_index++) {
				if (CheckMismatches(*CurrentChr, Temp_One_Read->UnmatchedSeq, UP[RightUP_index])) {
					Temp_One_Read->Used = false;
					Temp_One_Read->UP_Far.push_back(UP[RightUP_index]);
				}
			}
		}

		UP.clear();
	}
}

void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const unsigned int SearchCenter, const unsigned int Range) {

    //short(log((double)(Temp_One_Read.InsertSize * 2 + Range))/log(4.0) + 0.8) + 3;
    // short (log ((double) (Temp_One_Read.InsertSize * 3)) / log (4.0) + 0.8) + 3;

if (Temp_One_Read.UP_Far.size()>0 ) { std::cout << "KAI1108 UP_Far.size() == " << Temp_One_Read.UP_Far.size() << std::endl; }

   short BP_Start = short(log((double)(Range))/log(4.0) + 0.8) + 3;  // required minimum length of match
	//Temp_One_Read->MinClose + (UseRangeIndex ? RangeIndex : 0);
	short BP_End = Temp_One_Read.ReadLengthMinus; // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
	std::vector<UniquePoint> UP; // temporary container for unique far ends
	std::vector<unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
	std::vector<unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
    
	int Start = SearchCenter - Range;
	int End = SearchCenter + Range;
    
	for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
		PD_Plus[CheckIndex].reserve(End - Start + 1);
		PD_Minus[CheckIndex].reserve(End - Start + 1);
	}
    
	char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
	char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];

	//if (Temp_One_Read.MatchedD == Plus)  {

	if (CurrentBase != 'N') {
	 	for (int pos = Start; pos < End; pos++) {
		   if (chromosome.at(pos) == CurrentBase) {
				PD_Plus[0].push_back(pos); // else 
		  	}	
		   else if (chromosome.at(pos) == CurrentBaseRC) {
			   PD_Minus[0].push_back(pos);
			}
      } 
	}

   if (PD_Minus[0].size() + PD_Plus[0].size() > 0) { // skip all reads starting with 'N'
	   CheckBoth(Temp_One_Read, chromosome, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
	}
    
	if (UP.empty()) {}
    //else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.MinClose >= Temp_One_Read.ReadLength) {} // match too long
	else if (UP[UP.size() - 1].LengthStr + Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size() - 1].LengthStr < Temp_One_Read.ReadLength) { // should put into UP_Far_backup
        if (Temp_One_Read.UP_Far_backup.empty()) { // UP_Far_backup is empty, put it straightforwards 
            Temp_One_Read.UP_Far_backup.swap(UP);
        }
        else if (UP[UP.size() - 1].LengthStr > Temp_One_Read.UP_Far_backup[Temp_One_Read.UP_Far_backup.size() - 1].LengthStr) { // check whether the new result is better
            Temp_One_Read.UP_Far_backup.clear();  
            Temp_One_Read.UP_Far_backup.swap(UP);
        } // if UP[UP.size()
    } // if read is incompletely mapped
	else { // should put into UP_Far
        if (Temp_One_Read.UP_Far.empty()) {
            Temp_One_Read.UP_Far.swap(UP);
        }
        else {
            std::cout << "We shouldn't get here: farend_searcher.cpp at line ~516" << std::endl;
        }
	}
	UP.clear();
}

FarEndSearcher::~FarEndSearcher() {
	// Empty.
}
