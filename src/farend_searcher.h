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

#ifndef FARENDSEARCHER_H_
#define FARENDSEARCHER_H_

typedef std::vector <unsigned> PosVector;

class FarEndSearchPerRegion {
public:
FarEndSearchPerRegion( const Chromosome* Chromosome, unsigned short NumberOfErrors, unsigned size)
: CurrentChromosome ( Chromosome )
{
    PosVector emptyPosVector;
    emptyPosVector.reserve(size);
    PD_Plus.assign( NumberOfErrors, emptyPosVector);
    PD_Minus.assign( NumberOfErrors, emptyPosVector);
    emptyPosVector.clear();
}

~FarEndSearchPerRegion() {
	//emptyPosVector.clear();
	PD_Plus.clear();
	PD_Minus.clear();
}
    

    const Chromosome * CurrentChromosome;
    std::vector< PosVector > PD_Plus;
    std::vector< PosVector > PD_Minus;
    //FarEndSearchPerRegion( const Chromosome* Chromosome, unsigned short NumberOfErrors, unsigned size);
};




void SearchFarEndAtPos(SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions );

void SearchFarEndAtPosPerfect(SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions );

#endif /* FARENDSEARCHER_H_ */
