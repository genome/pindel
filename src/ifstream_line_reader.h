//
//  ifstream_line_reader.h
//  pindel-gz
//
//  Created by Matthijs Moed on 4/19/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

#ifndef __IFSTREAMLINEREADER_H__
#define __IFSTREAMLINEREADER_H__

#include <fstream>
#include <string>
#include <stdexcept>
#include "line_reader.h"



class IfstreamLineReader:public LineReader
{
	private:
	
		std::ifstream *in;
		const std::string filename;
		std::string buffer;
		
		
		void Advance();
		
	
	public:
	
		IfstreamLineReader(const char *_filename);
			
			
		~IfstreamLineReader();
			
			
		virtual void Reset();
			
			
		virtual bool HasNext();
			
			
		virtual std::string NextLine();
};

#endif
