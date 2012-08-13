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

void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const unsigned int SearchCenter, const unsigned int Range)
{
   std::vector<unsigned int> PD_Plus[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
   std::vector<unsigned int> PD_Minus[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];

   int Start = SearchCenter - Range - Temp_One_Read.getReadLength();
   int End = SearchCenter + Range + Temp_One_Read.getReadLength();

   for (int CheckIndex = 0; CheckIndex < Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(); CheckIndex++) {
      PD_Plus[CheckIndex].reserve(End - Start + 1);
      PD_Minus[CheckIndex].reserve(End - Start + 1);
   }

   char CurrentBase = Temp_One_Read.getUnmatchedSeq()[0];
   char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];
	
	if (CurrentBase == 'N' || Temp_One_Read.MaxLenCloseEnd() == 0) return; // ask Kai: correct?

   for (int pos = Start; pos < End; pos++) {
      if (chromosome.at(pos) == CurrentBase) {
         PD_Plus[0].push_back(pos); // else
      }
      else if (chromosome.at(pos) == CurrentBaseRC) {
         PD_Minus[0].push_back(pos);
      }
   }

   if (PD_Minus[0].size() + PD_Plus[0].size() > 0) { 
	 	short BP_Start = 10; // perhaps use global constant like "g_MinimumLengthToReportMatch"
   	short BP_End = Temp_One_Read.getReadLengthMinus(); // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
	   SortedUniquePoints UP; // temporary container for unique far ends
      CheckBoth(Temp_One_Read, chromosome, Temp_One_Read.getUnmatchedSeq(), PD_Plus, PD_Minus, BP_Start, BP_End, 1, UP);
   
		if ( UP.MaxLen() > Temp_One_Read.MaxLenFarEnd() ) {
			Temp_One_Read.UP_Far.swap(UP);
		}
	   UP.clear(); // may not be necessary as this is deleted from the stack anyway
	}
}

void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions )
{
    std::vector<unsigned int> PD_Plus[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
    std::vector<unsigned int> PD_Minus[Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED()];
    
    int TotalSize = 0;
    for (unsigned int RegionIndex = 0; RegionIndex < Regions.size(); RegionIndex++) {
        TotalSize += Regions[RegionIndex].getEnd() - Regions[RegionIndex].getStart() + 1;        
    }
    
    for (int CheckIndex = 0; CheckIndex < Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(); CheckIndex++) {
        PD_Plus[CheckIndex].reserve(TotalSize);
        PD_Minus[CheckIndex].reserve(TotalSize);
    }
    
    char CurrentBase = Temp_One_Read.getUnmatchedSeq()[0];
    char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];
	
	if (CurrentBase == 'N' || Temp_One_Read.MaxLenCloseEnd() == 0) return; 
    int CurrentReadLength = Temp_One_Read.getReadLength();
    for (unsigned RegionIndex = 0; RegionIndex < Regions.size(); RegionIndex++) {
        int Start = Regions[RegionIndex].getStart() - CurrentReadLength;
        int End = Regions[RegionIndex].getEnd() + CurrentReadLength;
        for (int pos = Start; pos < End; pos++) {
            if (chromosome.at(pos) == CurrentBase) {
                PD_Plus[0].push_back(pos); // else
            }
            else if (chromosome.at(pos) == CurrentBaseRC) {
                PD_Minus[0].push_back(pos);
            }
        }
    }

    if (PD_Minus[0].size() + PD_Plus[0].size() > 0) {
	 	short BP_Start = 10; // perhaps use global constant like "g_MinimumLengthToReportMatch"
        short BP_End = Temp_One_Read.getReadLengthMinus(); // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
        SortedUniquePoints UP; // temporary container for unique far ends
        CheckBoth(Temp_One_Read, chromosome, Temp_One_Read.getUnmatchedSeq(), PD_Plus, PD_Minus, BP_Start, BP_End, 1, UP);
        
		/*if ((unsigned)(MaxEndSize(UP) + Temp_One_Read.MaxLenCloseEnd()) >= (unsigned)Temp_One_Read.getReadLength() && Temp_One_Read.UP_Far.empty()) {
			Temp_One_Read.UP_Far.swap(UP);
		}
		else */ if ( UP.MaxLen() > Temp_One_Read.MaxLenFarEnd() ) {
			Temp_One_Read.UP_Far.swap(UP);
		}
        UP.clear(); // may not be necessary as this is deleted from the stack anyway
	}
}

FarEndSearcher::~FarEndSearcher()
{
}
