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

#include "pindel.h"
#include "searcher.h"

#include "FarEndSearcher.h"

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

FarEndSearcher::~FarEndSearcher() {
	// Empty.
}
