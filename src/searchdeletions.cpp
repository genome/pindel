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

// Pindel include files
#include "control_state.h"
#include "search_variant.h"
#include "searchdeletions.h"
#include "pindel.h"
#include "reporter.h"

SearchDeletions::SearchDeletions()
{
   typeOfVariant = "deletions";
}

SearchDeletions::~SearchDeletions()
{

}

bool SearchDeletions::decisionBranch1(const SPLIT_READ& read, unsigned int CloseIndex, int FarIndex)
{
   return ( read.UP_Far[FarIndex].LengthStr + read.UP_Close[CloseIndex].LengthStr == read.getReadLength() )
          && ( read.UP_Far[FarIndex].AbsLoc > (read.UP_Close[CloseIndex]. AbsLoc + 1 ));
}

bool SearchDeletions::decisionBranch2(const SPLIT_READ& read, unsigned int CloseIndex, int FarIndex)
{
   return read.UP_Close[CloseIndex].LengthStr + read.UP_Far[FarIndex].LengthStr == read.getReadLength()
          && read. UP_Close[CloseIndex].AbsLoc > read.UP_Far[FarIndex]. AbsLoc + 1;
}

unsigned int SearchDeletions::calculateIndelSize(const SPLIT_READ& read )
{
   return (read.Right - read.Left) - read.getReadLengthMinus();
}

std::string SearchDeletions::getInsertedStr1(const SPLIT_READ& read )
{
   return "";
}

std::string SearchDeletions::getInsertedStr2(const SPLIT_READ& read )
{
   return "";
}

void SearchDeletions::outputResults(BDData & g_bdData, ControlState& currentState, std::vector<unsigned> Vars[], const unsigned NumBoxes, const SearchWindow& currentWindow)
{
   std::ofstream DeletionOutf(userSettings->getDOutputFilename().c_str(), std::ios::app);
   SortOutputD(currentState, NumBoxes, currentWindow.getChromosome()->getSeq(), currentState.Reads_SR, Vars, DeletionOutf);
   DeletionOutf.close();
}
