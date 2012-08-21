
#include "ifstream_line_reader.h"
#include "pindel.h" // for safeGetline

using std::ifstream;
//using std::getline;
using std::string;


IfstreamLineReader::IfstreamLineReader(const char *_filename):LineReader(), filename(_filename)
{
	in = new ifstream(_filename);
	Advance();
}


IfstreamLineReader::~IfstreamLineReader()
{
	if (in -> is_open())
		in -> close();
	delete in;
}


void IfstreamLineReader::Reset()
{
	if (in -> is_open()) { //EW180612: add in->clear() needed for some compilers
		in -> clear();
		in -> seekg(0);
	}
	Advance();
}


bool IfstreamLineReader::HasNext()
{
	return !in -> eof(); 
}


string IfstreamLineReader::NextLine()
{
	string tmp = buffer;
	Advance();
	return tmp;
}

void IfstreamLineReader::Advance()
{
	safeGetline(*in, buffer);
}
