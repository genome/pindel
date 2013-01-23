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

#ifndef PARAMETER_H
#define	PARAMETER_H

#include <iostream>
#include <string>

class Parameter
{
public:
	bool isRequired () const
	{
		return d_required;
	}
	void describe () const;
	std::string getDescription () const
	{
		return d_description;
	}
	std::string getShortName () const
	{
		return d_shortName;
	}
	std::string getLongName () const
	{
		return d_longName;
	}
	bool hasName (const std::string & name) const
	{
		return (d_shortName.compare (name) == 0
						|| d_longName.compare (name) == 0);
	}
	bool isSet () const
	{
		return d_isSet;
	}
	virtual void setValue (const std::string & value)
	{
		std::cout << "WHAT!" << std::endl;
	};
	virtual void setValue (const int value)
	{
	};
	virtual void setValue (const double value)
	{
	};
	virtual void setValue (const bool value)
	{
	};
	virtual int getIValue () const
	{
		return 0;
	};
	virtual bool getBValue () const
	{
		return false;
	};
	virtual std::string getSValue () const
	{
		return "";
	};
	virtual double getFValue () const
	{
		return 0.0f;
	};

	virtual bool isUnary () const
	{
		return false;
	}

	Parameter (const std::string & shortName, const std::string & longName,
						 const std::string & description, const bool required);
	virtual ~Parameter() {}

protected:
	void set ()
	{
		d_isSet = true;
	}															//cout << "setting " << d_shortName << endl; } // TODO: Ask Kai whether this can be removed


private:
	bool d_required;
	std::string d_shortName;
	std::string d_longName;
	std::string d_description;
	bool d_isSet;
	static const int d_DESCRIBE_WIDTH = 11;
	static const int d_MAX_LINE_WIDTH = 80;
  const std::string getWord (std::string & line, bool& getLineEnd) const;
	const std::string makeNiceLine (const std::string & rawDescription) const;
	const std::string describeTab () const
	{															//return string(' ',d_DESCRIBE_WIDTH);  // TODO: Ask Kai whether this can be removed
		std::string descriptionTab = "";
		for (int i = 0; i < d_DESCRIBE_WIDTH; i++)
			{
				descriptionTab += " ";
			} return descriptionTab;
	}
};

class IntParameter:public Parameter
{
public:
	IntParameter (int *par_ptr, const std::string & shortName,
								const std::string & longName, const std::string & description,
								const bool required, const int value);

	virtual int getIValue () const
	{
		return *d_data_ptr;
	}
	virtual void setValue (const std::string & value);
	virtual void setValue (const int value);
	virtual ~IntParameter() {}
private:
	int *d_data_ptr;
};

class UIntParameter:public Parameter
{
public:
	UIntParameter (unsigned int *par_ptr, const std::string & shortName,
								const std::string & longName, const std::string & description,
								const bool required, const unsigned int value);

	virtual unsigned int getUIValue () const
	{
		return *d_data_ptr;
	}
	virtual void setValue (const std::string & value);
	virtual void setValue (const unsigned int value);
	virtual ~UIntParameter() {}
private:
	unsigned int *d_data_ptr;
};

class BoolParameter:public Parameter
{
public:
	BoolParameter (bool * par_ptr, const std::string & shortName,
								 const std::string & longName, const std::string & description,
								 const bool required, const bool value);

	virtual bool getBValue () const
	{
		return *d_data_ptr;
	}
	virtual void setValue (const std::string & value);
	virtual void setValue (const bool value);
	virtual bool isUnary () const
	{
		return true;
	}
	virtual ~BoolParameter() {}
private:
	  bool * d_data_ptr;
};

/* FloatParameter */
class FloatParameter:public Parameter
{
public:
	FloatParameter (double *par_ptr, const std::string & shortName,
									const std::string & longName, const std::string & description,
									const bool required, const double value);

	double getFValue () const
	{
		return *d_data_ptr;
	}
	void setValue (const std::string & value);
	void setValue (const double value);
	virtual ~FloatParameter() {}
private:
	double *d_data_ptr;
};

class StringParameter:public Parameter
{
public:
	StringParameter (std::string * par_ptr, const std::string & shortName,
									 const std::string & longName, const std::string & description,
									 const bool required, const std::string & value);

	std::string getSValue () const
	{
		return *d_data_ptr;
	}
	void setValue (const std::string & value);
	virtual ~StringParameter() {}
private:
	  std::string * d_data_ptr;
};



#endif	/* PARAMETER_H */

