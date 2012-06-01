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
#include "control_state.h"
#include "search_variant.h"
#include "searchshortinsertions.h"
#include "pindel.h"
#include "reporter.h"

SearchShortInsertions::SearchShortInsertions()
{
   typeOfVariant = "short insertions";
}

SearchShortInsertions::~SearchShortInsertions()
{

}

bool SearchShortInsertions::decisionBranch1(ControlState& currentState,
      unsigned ReadIndex, unsigned int CloseIndex, int FarIndex)
{
   return currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc
          == currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].AbsLoc + 1
          && currentState.Reads_SR[ReadIndex].UP_Close[CloseIndex].LengthStr
          + currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. LengthStr
          < currentState.Reads_SR[ReadIndex].ReadLength;
}

bool SearchShortInsertions::decisionBranch2(ControlState& currentState,
      unsigned ReadIndex, unsigned int CloseIndex, int FarIndex)
{
   return currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].AbsLoc
          == currentState.Reads_SR[ReadIndex].UP_Far[FarIndex]. AbsLoc + 1
          && currentState.Reads_SR[ReadIndex]. UP_Far[FarIndex].LengthStr
          + currentState.Reads_SR[ReadIndex]. UP_Close[CloseIndex].LengthStr
          < currentState.Reads_SR[ReadIndex].ReadLength;
}

unsigned int SearchShortInsertions::calculateIndelSize(
   ControlState& currentState, unsigned ReadIndex)
{
   return currentState.Reads_SR[ReadIndex].ReadLengthMinus
          - (currentState.Reads_SR[ReadIndex].Right
             - currentState.Reads_SR[ReadIndex].Left);
}

std::string SearchShortInsertions::getInsertedStr1(ControlState& currentState,
      unsigned ReadIndex)
{
   return ReverseComplement(currentState.Reads_SR[ReadIndex]. UnmatchedSeq). substr(
             currentState.Reads_SR[ReadIndex].BP + 1,
             currentState.Reads_SR[ReadIndex]. IndelSize);
}

std::string SearchShortInsertions::getInsertedStr2(ControlState& currentState,
      unsigned ReadIndex)
{
   return currentState.Reads_SR[ReadIndex].UnmatchedSeq. substr(
             currentState.Reads_SR[ReadIndex].BP + 1,
             currentState.Reads_SR[ReadIndex]. IndelSize);
}

void SearchShortInsertions::outputResults(ControlState& currentState,
      std::vector<unsigned> Vars[], const unsigned NumBoxes)
{
   std::ofstream SIoutputfile(currentState.SIOutputFilename.c_str(),
                              std::ios::app);
   SortOutputSI(NumBoxes, currentState.CurrentChrSeq, currentState.Reads_SR, Vars,
                SIoutputfile);
   SIoutputfile.close();
}
