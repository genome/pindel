//
//  pindel_read_reader.cpp
//  pindel-gz
//
//  Created by Matthijs Moed on 4/19/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

#include "pindel_read_reader.h"
#include "line_reader.h"

#include <string>
#include <sstream>
using std::string;
using std::istringstream;


PindelReadReader::PindelReadReader(LineReader &_reader):reader(_reader)
{
	advance();
}


PindelReadReader::~PindelReadReader()
{
}


void PindelReadReader::Reset()
{
	reader.Reset();
	advance();
}


bool PindelReadReader::HasNext()
{
	return (reader.HasNext() || buffer.Name!="" );
}


SPLIT_READ PindelReadReader::NextRead()
{
	const SPLIT_READ tmp = buffer;
	advance();
	return tmp;
}


void PindelReadReader::advance()
{
	buffer = SPLIT_READ();
	
	// Read line by line
	 
	buffer.Name = reader.NextLine();
	buffer.setUnmatchedSeq( reader.NextLine());
	
	istringstream iss(reader.NextLine().c_str());
	iss >> buffer.MatchedD;
	iss >> buffer.FragName;
	iss >> buffer.MatchedRelPos;
	iss >> buffer.MS;
	iss >> buffer.InsertSize;
	iss >> buffer.Tag;
}
