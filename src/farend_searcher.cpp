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
   //if (Temp_One_Read.UP_Far.size()>0 ) {
   //   std::cout << "KAI1108 UP_Far.size() == " << Temp_One_Read.UP_Far.size() << std::endl;
   //}
    //Temp_One_Read.UP_Far.clear();
    Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size();
    Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
   short BP_End = Temp_One_Read.ReadLengthMinus; // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
   std::vector<UniquePoint> UP; // temporary container for unique far ends
    Temp_One_Read.MAX_SNP_ERROR = (short) (Temp_One_Read.UnmatchedSeq.size () * Seq_Error_Rate);
    Temp_One_Read.TOTAL_SNP_ERROR_CHECKED = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
    Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus = Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
   std::vector<unsigned int> PD_Plus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
   std::vector<unsigned int> PD_Minus[Temp_One_Read.TOTAL_SNP_ERROR_CHECKED];
    //std::cout << "In SearchFarEndAtPos " << Temp_One_Read.UnmatchedSeq << " " << Temp_One_Read.MAX_SNP_ERROR << " " << Temp_One_Read.TOTAL_SNP_ERROR_CHECKED << std::endl;
    //std::cout <<  chromosome.size() << std::endl;
   int Start = SearchCenter - Range - Temp_One_Read.ReadLength;
   int End = SearchCenter + Range + Temp_One_Read.ReadLength;

   for (int CheckIndex = 0; CheckIndex < Temp_One_Read.TOTAL_SNP_ERROR_CHECKED; CheckIndex++) {
      PD_Plus[CheckIndex].reserve(End - Start + 1);
      PD_Minus[CheckIndex].reserve(End - Start + 1);
   }

   char CurrentBase = Temp_One_Read.UnmatchedSeq[0];
   char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];

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
   short BP_Start = 10;
   // std::cout << " + " << PD_Plus[0].size() << " - " << PD_Minus[0].size() << std::endl;
   if (PD_Minus[0].size() + PD_Plus[0].size() > 0) { // skip all reads starting with 'N'
      CheckBoth(Temp_One_Read, chromosome, Temp_One_Read.UnmatchedSeq, PD_Plus, PD_Minus, BP_Start, BP_End, FirstBase, UP);
   }
    //std::cout << "after CheckBoth" << std::endl;
    //std::cout << "UP.size() " << UP.size() << " ReadLength: " << Temp_One_Read.UnmatchedSeq.size() 
    //          << " " << Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size() - 1].LengthStr << " " << UP[UP.size() - 1].LengthStr << " Sum: " << Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size() - 1].LengthStr + UP[UP.size() - 1].LengthStr << std::endl;

   if (UP.empty()) {}
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
         std::cout << "We shouldn't get here: farend_searcher.cpp at line 90" << std::endl;
      }
   }
   UP.clear();
   // std::cout << "end of SearchFarEndAtPos" << std::endl; 
}

FarEndSearcher::~FarEndSearcher()
{
}
