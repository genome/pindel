//
//  pindel_read_reader.h
//  pindel-gz
//
//  Created by Matthijs Moed on 4/19/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

#ifndef __PINDELREADREADER_H__
#define __PINDELREADREADER_H__


#include "pindel.h" // For SPLIT_READ struct

class LineReader;


class PindelReadReader
{
	
	private:
		
		virtual void advance();
		
		LineReader &reader;
		SPLIT_READ buffer;
	
	public:
		
		PindelReadReader(LineReader &_reader);
		
		
		virtual ~PindelReadReader();
		
		
		virtual void Reset();
			
			
		virtual bool HasNext();
			
			
		SPLIT_READ NextRead();
};


#endif
