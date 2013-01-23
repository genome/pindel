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
#include <vector>
#include <map>
#include <string>

#ifndef rdrp_H
#define	rdrp_H

struct ChrNameAndSize {
    std::string ChrName;
    unsigned ChrSize;
};

struct DiscordantRP {
    std::string Read1_ChrName;
    unsigned Read1_Pos;
    char Read1_D;
    std::string Read2_ChrName;
    unsigned Read2_Pos;
    char Read2_D;
    std::string SampleName;
};

struct BAM_Path_IS {
    BAM_Path_IS() {
        BamFile = "";
        InsertSize = 0;
    }
    std::string BamFile;
    int InsertSize;
};

struct SampleAndBAMFiles {
    std::string SampleName;
    std::vector <BAM_Path_IS> BAMs;
};

void do_rd_rp (ControlState & CurrentState, ParCollection & par);
short getWholeGenomeSize(ControlState & CurrentState, std::vector <ChrNameAndSize> & ChrSizes);

#endif /* rdrp_H */
