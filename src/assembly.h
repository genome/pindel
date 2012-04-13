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

#include "control_state.h"
#include "pindel.h"

#ifndef ASSEMBLY_H
#define	ASSEMBLY_H

void doAssembly (ControlState & CurrentState, ParCollection & par);

struct Chromosome {
    Chromosome() {
        ChrName = "";
        Seq = ""; 
    }
    std::string ChrName;
    std::string Seq;
};

struct Assembly {
    Assembly() {
        Type = "";
        ChrA = "";
        PosA = 0;
        CI_A = 0;
        ChrB = "";
        PosB = 0;
        CI_B = 0;        
    }
    std::string Type;
    std::string ChrA;
    unsigned PosA;
    unsigned CI_A;
    std::string ChrB;
    unsigned PosB;
    unsigned CI_B;
};

#endif /* ASSEMBLY_H */
