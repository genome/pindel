//
//  gz_line_reader.cpp
//  pindel-gz
//
//  Created by Matthijs Moed on 4/19/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

#include <zlib.h>
#include <cstring>
#include "gz_line_reader.h"

using std::string;


const int GZLineReader::BUFFER_SIZE = 128;


void GZLineReader::Open()
{
	in = gzopen(filename.c_str(), "rb");
	cBuffer = new char[BUFFER_SIZE];
	eof = false;
	Advance();
}


void GZLineReader::Close()
{
	gzclose(in);
}


void GZLineReader::Advance()
{
	if (eof)
		return;
	oss.str("");
	
	do
	{
		memset(cBuffer, ' ', BUFFER_SIZE);
		const char *returnCode = gzgets(in, cBuffer, BUFFER_SIZE);
		if (returnCode == NULL)
		{
			eof = true;
			break;
		}
		
		// Strip newline if present
		
		newLinePtr = strchr(cBuffer, int('\n'));
		if (newLinePtr != NULL)
			oss << std::string(cBuffer, newLinePtr - cBuffer);
		else
			oss << cBuffer;
	} while (newLinePtr == NULL);
	buffer = oss.str();
}

GZLineReader::GZLineReader(const char *_filename):LineReader(), filename(_filename)
{
	Open();
	cBuffer = new char[BUFFER_SIZE];
	Advance();
}


GZLineReader::~GZLineReader()
{
	Close();
}


void GZLineReader::Reset()
{
	// Sinze gzseek can be extremely slow when a file is opened in read mode, we simply reopen the file
	
	Close();
	Open();
}


bool GZLineReader::HasNext()
{
	return !eof && buffer.length() > 0;
}


string GZLineReader::NextLine()
{
	std::string tmp = buffer;
	Advance();
	return tmp;
}