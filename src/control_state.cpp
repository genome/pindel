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

//Han(2013.06.17)
BreakDancerCoordinate::BreakDancerCoordinate( const std::string& chromosomeName, const unsigned int pos, const unsigned int pos2 ) : m_chromosomeName ( chromosomeName )
{
    if (pos <= pos2) {
        position = pos;
        position2 = pos2; //Han(2013.06.17)
    }
    else {
        position = pos2;
        position2 = pos;
    }
}

BreakDancerCoordinate::BreakDancerCoordinate( const std::string& chromosomeName, const unsigned int pos ) : m_chromosomeName ( chromosomeName )
{
	position = pos;
    position2 = pos; //Han(2013.06.17)
}

BreakDancerCoordinate::BreakDancerCoordinate( const Chromosome* const chromosome, const unsigned int pos ): m_chromosomeName( chromosome->getName() )
{
	position = pos;
    position2 = pos; //Han(2013.06.17)
}

BreakDancerCoordinate::BreakDancerCoordinate( const BreakDancerCoordinate& other ) : m_chromosomeName ( other.m_chromosomeName ), position (other.position ), position2 (other.position2 )	//Han(2013.06.17)
{}


BreakDancerCoordinate& BreakDancerCoordinate::operator=(const BreakDancerCoordinate& other)
{
   if (this != &other) {
      this->BreakDancerCoordinate::~BreakDancerCoordinate(); // explicit non-virtual destructor
      new (this) BreakDancerCoordinate( other ); // placement new
   }
   return *this;
}

const Chromosome* BreakDancerCoordinate::getChromosome() const 
{ 
	const Chromosome* foundChromosome = g_genome.getChr(m_chromosomeName); 
	if (foundChromosome==NULL) {
		std::cout << "Error: chromosome with name : " << m_chromosomeName << " not yet loaded into memory. Aborting.\n";
		exit( EXIT_FAILURE );
	}
	return foundChromosome;
}

unsigned int BreakDancerCoordinate::startOfWindow() const
{	//Han(2013.06.17)
	unsigned int tmp_position = position;
	if ( position2 < position && position2 > 0 )  tmp_position = position2;
	if ( tmp_position>=BREAKDANCER_WINDOWSPAN) {return tmp_position - BREAKDANCER_WINDOWSPAN; }
	else { return 0; }
    
}

BreakDancerEvent::BreakDancerEvent( const BreakDancerEvent& other ) : first (other.first), second( other.second)
{
}

unsigned int BreakDancerCoordinate::endOfWindow() const // NOTE: this is dangerous unless we also save the chromosome size somewhere...
{	//Han(2013.06.17)
	unsigned int tmp_position = position;
	if ( position2 > position && position2 > 0 )  tmp_position = position2;
	return tmp_position + BREAKDANCER_WINDOWSPAN;
}

bool BreakDancerCoordinate::operator<(const BreakDancerCoordinate& other ) const
{
	if (getChromosomeName()!=other.getChromosomeName()) {
		return (getChromosomeName() < other.getChromosomeName() );
	}
	else if (position!=other.position ) {
		return (position < other.position );
	} 	
	else {
		return false; // are equal, so no <
	}
}

bool sortOnFirstBDCoordinate( const BreakDancerEvent& event1, const BreakDancerEvent& event2 )
{
	if (event1.first != event2.first ) {
		return (event1.first < event2.first );
	}
	else if (event1.second != event2.second ) {
		return (event1.second < event2.second );
	}
	else return false; // they're exactly equal so event1 is not 'smaller than' event2
}

bool sortOnSecondBDCoordinate( const BreakDancerEvent& event1, const BreakDancerEvent& event2 )
{
	if (event1.second != event2.second ) {
		return (event1.second < event2.second );
	}
	else if (event1.first != event2.first ) {
		return (event1.first < event2.first );
	}
	else return false; // they're exactly equal so event1 is not 'smaller than' event2
}



ControlState::ControlState()
{
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
