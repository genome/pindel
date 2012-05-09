//
//  gz_line_reader.h
//  pindel-gz
//
//  Created by Matthijs Moed on 4/19/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

#ifndef __GZLINEREADER_H__
#define __GZLINEREADER_H__

#include <string>
#include <sstream>
#include "line_reader.h"


class GZLineReader:public LineReader
{
	private:
		
		std::string buffer;
		gzFile in;
		std::string filename;
		const static int BUFFER_SIZE;
		bool eof;
		
		char *cBuffer;
		std::ostringstream oss;
		char *newLinePtr;
		
		
		void Open();
			
			
		void Close();		
			
		void Advance();
		
	
	public:
		
		GZLineReader(const char *_filename);
			
			
		~GZLineReader();
			
		virtual void Reset();
			
			
		virtual bool HasNext();
			
			
		virtual std::string NextLine();
};



#endif
