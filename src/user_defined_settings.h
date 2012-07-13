/* UserDefinedSettings' manages the parameters that have been set by the users or have predefined default values. */


// singleton pattern

#include <string>

#ifndef USER_DEFINED_SETTINGS_H
#define USER_DEFINED_SETTINGS_H

class UserDefinedSettings {

public:
   static UserDefinedSettings* Instance();

	int ADDITIONAL_MISMATCH; // user
	bool Analyze_BP; 
	bool Analyze_INV; 
	bool Analyze_LI; 
	bool Analyze_TD; 
	unsigned int BalanceCutoff;
	std::string bamConfigFilename;
	std::string breakdancerFilename;
	std::string breakdancerOutputFilename;
	double FLOAT_WINDOW_SIZE;
   std::string inf_AssemblyInputFilename; 
	std::string inf_GenotypingInputFilename;
	std::string logFilename;	
	double MaximumAllowedMismatchRate;

	int MaxRangeIndex;
	unsigned int minimalAnchorQuality;
	int MIN_IndelSize_Inversion;
	int Min_Num_Matched_Bases;
	int Min_Perfect_Match_Around_BP;
	unsigned int NumRead2ReportCutOff;
	int numThreads;
	std::string outputFilename;
	std::string pindelConfigFilename;
	std::string pindelFilename;
	std::string referenceFilename;
	bool ReportCloseMappedRead;
	bool reportOnlyCloseMappedReads;
	std::string SearchRegion;
	double Seq_Error_Rate;
	bool showHelp;

	bool reportCloseMappedReads() const { return ( ReportCloseMappedRead || reportOnlyCloseMappedReads ); };

	std::string getSIOutputFilename() const { return outputFilename + "_SI"; };
	std::string getDOutputFilename() const { return outputFilename + "_D"; };
	std::string getTDOutputFilename() const { return outputFilename + "_TD"; };
	std::string getINVOutputFilename() const { return outputFilename + "_INV"; };
	std::string getLIOutputFilename() const { return outputFilename + "_LI"; };
	std::string getBPOutputFilename() const { return outputFilename + "_BP"; };
	std::string getCloseEndOutputFilename() const { return outputFilename + "_CloseEndMapped"; };
	std::string getASMOutputFilename() const { return outputFilename + "_ASM"; };
	std::string getGTOutputFilename() const { return outputFilename + "_GT"; };

private:
   UserDefinedSettings();
	static UserDefinedSettings* m_instance;
};

#endif
