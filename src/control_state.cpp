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

ControlState::ControlState()
{
   lowerBinBorder = 0;
   upperBinBorder = 0;
   endRegionPlusBuffer = 0;

   CountFarEnd = 0;
   CountFarEndPlus = 0;
   CountFarEndMinus = 0;

	lineReader = 0;
	inf_Pindel_Reads = 0;
}

ControlState::~ControlState()
{
   // TODO Auto-generated destructor stub
}
