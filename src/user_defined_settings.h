/* UserDefinedSettings' manages the parameters that have been set by the users or have predefined default values. */


// singleton pattern

#include <string>

#ifndef USER_DEFINED_SETTINGS_H
#define USER_DEFINED_SETTINGS_H

class SearchRegion {

public:
	SearchRegion(const std::string& regionString);
	bool isStartDefined() const; 
	bool isEndDefined() const;
	bool isTargetChromosomeDefined() const;
	std::string getTargetChromosomeName() const;
	int getStart() const;
	int getEnd() const;

private:
	SearchRegion();
	bool m_startDefined;
	bool m_endDefined;
	std::string m_targetChromosomeName;
	int m_start;
	int m_end;
};

class UserDefinedSettings {

public:
   static UserDefinedSettings* Instance();

	int ADDITIONAL_MISMATCH; // user
	bool Analyze_BP; 
	bool Analyze_INV; 
	bool Analyze_LI; 
	bool Analyze_TD; 
    bool Analyze_MEI;
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
	std::string userDefinedRegion;
	double Seq_Error_Rate;
	bool showHelp;

	bool reportCloseMappedReads() const { return ( ReportCloseMappedRead || reportOnlyCloseMappedReads ); };
	bool singlePindelFileAsInput() const { return pindelFilename!=""; };
	bool pindelConfigFileAsInput() const { return pindelConfigFilename != ""; };
	bool bamFilesAsInput() const { return bamConfigFilename != ""; };
	bool pindelFilesAsInput() const { return ( singlePindelFileAsInput() || pindelConfigFileAsInput() ); };

	std::string getSIOutputFilename() const { return outputFilename + "_SI"; };
	std::string getDOutputFilename() const { return outputFilename + "_D"; };
	std::string getTDOutputFilename() const { return outputFilename + "_TD"; };
	std::string getINVOutputFilename() const { return outputFilename + "_INV"; };
	std::string getLIOutputFilename() const { return outputFilename + "_LI"; };
	std::string getBPOutputFilename() const { return outputFilename + "_BP"; };
	std::string getCloseEndOutputFilename() const { return outputFilename + "_CloseEndMapped"; };
	std::string getASMOutputFilename() const { return outputFilename + "_ASM"; };
	std::string getGTOutputFilename() const { return outputFilename + "_GT"; };

	bool loopOverAllChromosomes() { return ! getRegion()->isTargetChromosomeDefined(); };

	SearchRegion* getRegion();

private:
   UserDefinedSettings();
	static UserDefinedSettings* m_instance;
	SearchRegion* m_region; 
};

#endif
