/* 
 * File:   parameter.h
 * Author: enckevort
 *
 * Created on 16 mei 2011, 15:57
 */

#ifndef PARAMETER_H
#define	PARAMETER_H

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
  const std::string getWord (std::string & line) const;
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
private:
	int *d_data_ptr;
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
private:
	  std::string * d_data_ptr;
};



#endif	/* PARAMETER_H */

