/* UserDefinedSettings' manages the parameters that have been set by the users or have predefined default values. */


// singleton pattern

#include <string>

#ifndef USER_DEFINED_SETTINGS_H
#define USER_DEFINED_SETTINGS_H

class UserDefinedSettings {

public:
   static UserDefinedSettings* Instance();

	std::string bamConfigFilename;
	std::string breakdancerFilename;
	std::string breakdancerOutputFilename;
	double FLOAT_WINDOW_SIZE;
   std::string inf_AssemblyInputFilename; 
	std::string inf_GenotypingInputFilename;
	std::string logFilename;	
	int MaxRangeIndex;
	int numThreads;
	std::string outputFilename;
	std::string pindelConfigFilename;
	std::string pindelFilename;
	std::string referenceFilename;
	bool reportOnlyCloseMappedReads;
	std::string SearchRegion;
	double Seq_Error_Rate;
	bool showHelp;

private:
   UserDefinedSettings();
	static UserDefinedSettings* m_instance;
};

#endif
