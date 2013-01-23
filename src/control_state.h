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

#ifndef CONTROLSTATE_H_
#define CONTROLSTATE_H_

// System header files
#include <iostream>
#include <fstream>

// Pindel header files
#include "pindel.h"
#include "pindel_read_reader.h"



struct BreakDancerCoordinate {
	static const unsigned int BREAKDANCER_WINDOWSPAN = 200;

	std::string chromosomeName;
	unsigned int position;

	BreakDancerCoordinate( const std::string& chrName, const unsigned int pos ) : chromosomeName( chrName ), position (pos) {};
	unsigned int startOfWindow() const;
	unsigned int endOfWindow() const; // NOTE: this is dangerous unless we also save the chromosome size somewhere...
	bool operator<(const BreakDancerCoordinate& other ) const;
	bool operator!=(const BreakDancerCoordinate& other ) const { return ( chromosomeName!=other.chromosomeName || position!=other.position ); }
};

struct BreakDancerEvent {
	BreakDancerCoordinate first, second;

	BreakDancerEvent( const BreakDancerCoordinate& bd1, const BreakDancerCoordinate& bd2 ) : first( bd1 ), second (bd2 ) {};
	
	friend bool sortOnFirstBDCoordinate( const BreakDancerEvent& event1, const BreakDancerEvent& event2 );
	friend bool sortOnSecondBDCoordinate( const BreakDancerEvent& event1, const BreakDancerEvent& event2 );
};

bool sortOnFirstBDCoordinate( const BreakDancerEvent& event1, const BreakDancerEvent& event2 );
bool sortOnSecondBDCoordinate( const BreakDancerEvent& event1, const BreakDancerEvent& event2 );

typedef std::vector< RP_READ> RPVector; 

class ControlState {
public:

	ControlState();
	virtual ~ControlState();

    std::ifstream inf_AssemblyInput;
    std::ifstream inf_GenotypingInput;
	LineReader *lineReader;
	PindelReadReader *inf_Pindel_Reads;
	std::vector<bam_info> bams_to_parse;
	std::vector<std::string> pindelfilesToParse;
    
	std::string CurrentChrSeq;
	std::string CurrentChrName;

	std::vector<SPLIT_READ> InputReads_SR, Reads_SR, FutureReads_SR;
    std::vector <RPVector> Reads_RP;
 
	int CountFarEnd, CountFarEndPlus, CountFarEndMinus;

	//TODO: explain what are stored in these two vectors.
	std::vector<BreakDancerEvent> All_BD_events_WG, All_BD_events;
};

#endif /* CONTROLSTATE_H_ */
