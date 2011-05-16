/* 'Parameter' stores an individual parameter; it is set by the command line parameters, and used by the program. */
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <omp.h>

#include "parameter.h"

Parameter::Parameter (const std::string & shortName, const std::string & longName,
											const std::string & description, const bool required)
{
	d_required = required;
	d_shortName = shortName;
	d_longName = longName;
	d_description = description;
	d_isSet = false;
}

/* 'getWord' removes the first word from the line, a bit like perl's unshift. */
const std::string
Parameter::getWord (std::string & line) const
{
	std::string head, tail;
	size_t endPos = line.find (" ");
	if (endPos == line.npos)
		{
			head = line;
			tail = "";
		}
	else
		{
			head = line.substr (0, endPos);
			tail = line.substr (endPos + 1);
		}

	line = tail;
	return head;
}

/* 'makeNiceLine' adds \n's on the right position of the description to get it fitting nicely in a window. */
const std::string
Parameter::makeNiceLine (const std::string & rawDescription) const
{
	std::string neatLine = describeTab ();
	std::string words = rawDescription;
	const size_t LIMIT = d_MAX_LINE_WIDTH - d_DESCRIBE_WIDTH;
	size_t lineSize = 0;
	while (words.size () > 0)
		{
			std::string newWord = getWord (words);
			if (newWord.size () + lineSize > LIMIT)
				{
					neatLine += '\n';
					neatLine += describeTab ();
					lineSize = 0;
				}
			neatLine += newWord + " ";
			lineSize += newWord.size () + 1;	// plus space!
		}

	return neatLine;
}

void
Parameter::describe () const
{
	for (int i = 0; i < d_DESCRIBE_WIDTH; i++)
		{
			std::cout << " ";
		}
	std::cout << d_shortName << "/";
	std::cout << d_longName << std::endl;

	//for (int i=0; i<d_DESCRIBE_WIDTH; i++)  { cout << " ";} // TODO: Ask Kai whether this can be removed
	std::cout << makeNiceLine (d_description);
	//if ( d_required ) { cout << " required parameter" ; } // TODO: Ask Kai whether this can be removed
	std::cout << std::endl << std::endl;
}

/* IntParameter */
IntParameter::IntParameter (int *par_ptr, const std::string & shortName,
														const std::string & longName,
														const std::string & description, const bool required,
														const int value):
Parameter (shortName, longName, description, required),
d_data_ptr (par_ptr)
{
	*d_data_ptr = value;
}

void
IntParameter::setValue (const std::string & value)
{
	setValue (atoi (value.c_str ()));
}

void
IntParameter::setValue (const int value)
{
	*d_data_ptr = value;
	set ();
}

/* BoolParameter */
BoolParameter::BoolParameter (bool * par_ptr, const std::string & shortName, const std::string & longName, const std::string & description, const bool required, const bool value):
Parameter (shortName, longName, description,
					 required)
{
	d_data_ptr = par_ptr;
	*d_data_ptr = value;
}

void
BoolParameter::setValue (const std::string & value)
{
	char firstChar = tolower (value[0]);
	setValue ((firstChar == 'f' || firstChar == '0') ? false : true);
}

void
BoolParameter::setValue (const bool value)
{
	*d_data_ptr = value;
	set ();
}

/* FloatParameter */
FloatParameter::FloatParameter (double *par_ptr, const std::string & shortName,
																const std::string & longName,
																const std::string & description,
																const bool required, const double value):
Parameter (shortName, longName, description, required)
{
	d_data_ptr = par_ptr;
	*d_data_ptr = value;
}

void
FloatParameter::setValue (const std::string & value)
{
	setValue (atof (value.c_str ()));
}

void
FloatParameter::setValue (const double value)
{
	*d_data_ptr = value;
	set ();
}

/* StringParameter */
StringParameter::StringParameter (std::string * par_ptr, const std::string & shortName,
																	const std::string & longName,
																	const std::string & description,
																	const bool required, const std::string & value):
Parameter (shortName, longName, description, required)
{
	d_data_ptr = par_ptr;
	*d_data_ptr = value;
}

void
StringParameter::setValue (const std::string & value)
{
	*d_data_ptr = value;
	set ();
}