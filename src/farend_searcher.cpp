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

FarEndSearchPerRegion::FarEndSearchPerRegion( const Chromosome* Chromosome, unsigned short NumberOfErrors, unsigned size)
: CurrentChromosome ( Chromosome )
{
    PosVector emptyPosVector;
    emptyPosVector.reserve(size);
    PD_Plus.assign( NumberOfErrors, emptyPosVector);
    PD_Minus.assign( NumberOfErrors, emptyPosVector);
    
}

void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions )
{
    
    
   // step 1 find out which chromosomes in Regions: set? linear pass of regions
   // step 2 for each identified chromsome, for each regions on the chromosme, do the business.
   char CurrentBase = Temp_One_Read.getUnmatchedSeq()[0];
   char CurrentBaseRC = Convert2RC4N[(short) CurrentBase];
    
   if (CurrentBase == 'N' || Temp_One_Read.MaxLenCloseEnd() == 0) return;
   int CurrentReadLength = Temp_One_Read.getReadLength();
    
   std::vector <FarEndSearchPerRegion> WholeGenomeSearchResult;
   unsigned NumberOfHits = 0;
	/*if (Temp_One_Read.Name=="@read_6990/2" ) {
   	 std::cout << "6990/2Number of regions: " << Regions.size() << " " << Temp_One_Read.Name << " " << Temp_One_Read.FragName << " " << Temp_One_Read.MatchedD << " " << Temp_One_Read.MatchedRelPos << std::endl;
	}*/

   for (unsigned RegionIndex = 0; RegionIndex < Regions.size(); RegionIndex++) {
       
       FarEndSearchPerRegion CurrentRegion(Regions[RegionIndex].getChromosome(), Temp_One_Read.getTOTAL_SNP_ERROR_CHECKED(), Regions[RegionIndex].getSize());

       int Start = Regions[RegionIndex].getStart() - CurrentReadLength;
       int End = Regions[RegionIndex].getEnd() + CurrentReadLength;
       if (Start < 0) Start = End -1;
       const std::string & chromosome = Regions[RegionIndex].getChromosome()->getSeq();
       //std::cout << Regions[RegionIndex].getChromosome()->getName() << " " << chromosome.size() - g_SpacerBeforeAfter * 2 << " " << Start - g_SpacerBeforeAfter << " " << End - g_SpacerBeforeAfter << std::endl;
       for (int pos = Start; pos < End; pos++) {
           if (chromosome.at(pos) == CurrentBase) {
               CurrentRegion.PD_Plus[0].push_back(pos); // else
           }
           else if (chromosome.at(pos) == CurrentBaseRC) {
               CurrentRegion.PD_Minus[0].push_back(pos);
           }
       }
       NumberOfHits += CurrentRegion.PD_Plus[0].size() + CurrentRegion.PD_Minus[0].size();
       WholeGenomeSearchResult.push_back(CurrentRegion);
   }

	if (NumberOfHits>0) {
      short BP_Start = 10; // perhaps use global constant like "g_MinimumLengthToReportMatch"
      short BP_End = Temp_One_Read.getReadLengthMinus(); // matched far end should be between BP_Start and BP_End bases long (including BP_Start and End)
      SortedUniquePoints UP; // temporary container for unique far ends
      CheckBoth(Temp_One_Read, Temp_One_Read.getUnmatchedSeq(), WholeGenomeSearchResult, BP_Start, BP_End, 1, UP);
        
	  	if ( UP.MaxLen() > Temp_One_Read.MaxLenFarEnd() ) {
	 		Temp_One_Read.UP_Far.swap(UP);
      }
      UP.clear(); // may not be necessary as this is deleted from the stack anyway
	}
}



/*void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions )
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
        
		if ( UP.MaxLen() > Temp_One_Read.MaxLenFarEnd() ) {
			Temp_One_Read.UP_Far.swap(UP);
		}
      UP.clear(); // may not be necessary as this is deleted from the stack anyway
	}
}*/

