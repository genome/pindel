//
//  LineReader.h
//  pindel-gz
//
//  Created by Matthijs Moed on 4/19/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

#ifndef __LINEREADER_H__
#define __LINEREADER_H__

#include <string>


class LineReader
{
	public:
		
		LineReader();
		virtual ~LineReader();
	
		virtual void		Reset() = 0;
		virtual bool		HasNext() = 0;
		virtual std::string NextLine() = 0;
};


#endif
