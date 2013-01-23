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

// Pindel header files
#include "pindel.h"
#include "logdef.h"
#include "refreader.h"

RefReader::RefReader(std::ifstream & inf_Seq_in, std::string & TheInput_in)
{
   inf_Seq = &inf_Seq_in;
   TheInput = &TheInput_in;
}

RefReader::~RefReader()
{
}

void RefReader::ReadSeq(bool WhetherBuildUp)
{
   TheInput->clear ();
   std::string TempLine, TempChrName;

   if (WhetherBuildUp) {
      CopyThisSeq();
   }
   else {
      SkipThisSeq();
   }

   if (!TheInput->empty ()) {
      std::string Spacer = "";
      for (unsigned i = 0; i < g_SpacerBeforeAfter; i++) {
         Spacer += "N";
      }
      *TheInput = Spacer + *TheInput + Spacer;
   }
}

void RefReader::SkipThisSeq()
{
   char TempChar;
   // Skip until the next sequence start.
   while (*inf_Seq >> TempChar) {
      if (TempChar == '>') {
         break;
      }
   }
}

void RefReader::CopyThisSeq()
{
   char TempChar;
   short Diff2UpperCase = 'A' - 'a';
   while (*inf_Seq >> TempChar) {
      if (TempChar != '\n' && TempChar != '\r') {
         if (TempChar == '>') {
            break;
         }
         else {
            if ('a' <= TempChar && TempChar <= 'z') {
               TempChar = TempChar + Diff2UpperCase;
            }
            switch (TempChar) {
            case 'A':
               *TheInput += 'A';
               break;	// 00000000
            case 'C':
               *TheInput += 'C';
               break;	// 00010000
            case 'G':
               *TheInput += 'G';
               break;	// 00100000
            case 'T':
               *TheInput += 'T';
               break;	// 00110000
            default:
               *TheInput += 'N';	// 01000000
            }
         }						// else TempChar
      }
   }
}
