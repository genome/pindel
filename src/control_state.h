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
//#include "bddata.h"


struct BED {
	std::string ChrName;
	unsigned Start;
	unsigned End;
};

struct BreakDancerCoordinate {

private:
	const std::string m_chromosomeName;

public:
	static const unsigned int BREAKDANCER_WINDOWSPAN = 200;
	BreakDancerCoordinate& operator=( const BreakDancerCoordinate& other );

	unsigned int position;
    unsigned int position2;	//Han(2013.06.17)
	BreakDancerCoordinate( const BreakDancerCoordinate& other );
    BreakDancerCoordinate( const std::string& chromosomeName, const unsigned int pos, const unsigned int pos2 );	//Han(2013.06.17)
	BreakDancerCoordinate( const std::string& chromosomeName, const unsigned int pos );
	BreakDancerCoordinate( const Chromosome* const chromosome, const unsigned int pos );
	unsigned int startOfWindow() const;
	unsigned int endOfWindow() const; // NOTE: this is dangerous unless we also save the chromosome size somewhere...
	bool operator<(const BreakDancerCoordinate& other ) const;
	bool operator!=(const BreakDancerCoordinate& other ) const { return ( getChromosomeName()!=other.getChromosomeName() || position!=other.position ); }
	const std::string& getChromosomeName() const { return m_chromosomeName; }
	const Chromosome* getChromosome() const;
};

struct BreakDancerEvent {
	BreakDancerCoordinate first, second;
	BreakDancerEvent( const BreakDancerEvent& other );
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

	std::vector <SPLIT_READ> InputReads_SR, Reads_SR, FutureReads_SR, InterChromosome_SR, OneEndMappedReads;
	std::vector <RPVector> Reads_RP;
	std::vector <RP_READ> Reads_RP_Discovery, Reads_RP_Discovery_InterChr;
	std::vector <REF_READ> RefSupportingReads;
	//SearchWindow::SearchWindow CURRENT_WINDOW;
	unsigned RegionStart, RegionEnd;
	int CountFarEnd, CountFarEndPlus, CountFarEndMinus;

	//TODO: explain what are stored in these two vectors.
	std::vector<BreakDancerEvent> External_BD, All_BD_events_WG, All_BD_events;

	std::vector<BED> IncludeBed, ExcludeBed;
};

#endif /* CONTROLSTATE_H_ */
