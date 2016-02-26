/* UserDefinedSettings' manages the parameters that have been set by the users or have predefined default values. */


// singleton pattern

#include <string>

#ifndef USER_DEFINED_SETTINGS_H
#define USER_DEFINED_SETTINGS_H

class SearchRegion
{

public:
   SearchRegion(const std::string& regionString);
   void SetRegion(const std::string& ChrName, const unsigned int start, const unsigned int end);
   bool isStartDefined() const;
   bool isEndDefined() const;
   bool isTargetChromosomeDefined() const;
   std::string getTargetChromosomeName() const;
   unsigned int getStart() const;
   unsigned int getEnd() const;

private:
   SearchRegion();
   bool m_startDefined;
   bool m_endDefined;
   std::string m_targetChromosomeName;
   unsigned int m_start;
   unsigned int m_end;
};

class UserDefinedSettings
{

public:

   UserDefinedSettings()
   {
      m_instance = NULL;
      m_region = NULL;

      ADDITIONAL_MISMATCH = 2; // user
      Analyze_BP = true;
      Analyze_INV = true;
      Analyze_LI = true;
      Analyze_TD = true;
      Analyze_TD = false;
      BalanceCutoff = 0;
      bamConfigFilename = "";
      breakdancerFilename = "";
      breakdancerOutputFilename = "";

      inf_InclusiveBedFileName = "";
      inf_ExclusiveBedFileName = "";

      PloidyFileName = "";
      FLOAT_WINDOW_SIZE = 0.0;
      inf_AssemblyInputFilename = "";
      inf_GenotypingInputFilename = "";
      logFilename = "";
      MaximumAllowedMismatchRate = 0.0;
      sensitivity = 0.0; // -E
      minClose = 8;

      MaxRangeIndex = 5;
      minimalAnchorQuality = 30;
      MIN_IndelSize_Inversion = 30;
      Min_Num_Matched_Bases = 30;
      Min_Perfect_Match_Around_BP = 3;
      NumRead2ReportCutOff = 3;
      numThreads = 1;
      outputFilename = "";
      pindelConfigFilename = "";
      pindelFilename = "";
      referenceFilename = "";
      SearchDiscordantReadPair = true;
      ReportCloseMappedRead = false;
      reportOnlyCloseMappedReads = false;
      reportInterchromosomalEvents = true;
      IndelCorrection = false;
      NormalSamples = false;
      userDefinedRegion = "";
      Seq_Error_Rate = 0.0;
      showHelp = false;
      WhetherEnforeStrandBias = true;
      NM = 2;
   }

   static UserDefinedSettings* Instance();



   int ADDITIONAL_MISMATCH; // user
   bool Analyze_BP;
   bool Analyze_INV;
   bool Analyze_LI;
   bool Analyze_TD;
   bool Analyze_DD;
   unsigned int BalanceCutoff;
   std::string bamConfigFilename;
   std::string breakdancerFilename;
   std::string breakdancerOutputFilename;

   std::string inf_InclusiveBedFileName;
   std::string inf_ExclusiveBedFileName;

   std::string PloidyFileName;
   double FLOAT_WINDOW_SIZE;
   std::string inf_AssemblyInputFilename;
   std::string inf_GenotypingInputFilename;
   std::string logFilename;
   double MaximumAllowedMismatchRate;
   double sensitivity; // -E
   int minClose;

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
   bool SearchDiscordantReadPair;
   bool ReportCloseMappedRead;
   bool reportOnlyCloseMappedReads;
   bool reportInterchromosomalEvents;
   bool IndelCorrection;
   bool NormalSamples;
   std::string userDefinedRegion;
   double Seq_Error_Rate;
   bool showHelp;
   bool WhetherEnforeStrandBias;
   // DD detection: max distance of breakpoints for same event.
   int MAX_DD_BREAKPOINT_DISTANCE;

   // DD detection: max distance reads for assigning to same cluster.
   int MAX_DISTANCE_CLUSTER_READS;

   // DD detection: minimal cluster size needed to estimate breakpoint.
   int MIN_DD_CLUSTER_SIZE;

   // DD detection: minimal n/o split reads to call exact breakpoint.
   int MIN_DD_BREAKPOINT_SUPPORT;

   // DD detection: minimal mapping distance for reads to be discordant.
   int MIN_DD_MAP_DISTANCE;

   // DD detection: report discordant mates.
   bool DD_REPORT_DUPLICATION_READS;

   int NM;

   bool reportCloseMappedReads() const
   {
      return ( ReportCloseMappedRead || reportOnlyCloseMappedReads );
   };
   bool singlePindelFileAsInput() const
   {
      return pindelFilename!="";
   };
   bool pindelConfigFileAsInput() const
   {
      return pindelConfigFilename != "";
   };
   bool bamFilesAsInput() const
   {
      return bamConfigFilename != "";
   };
   bool pindelFilesAsInput() const
   {
      return ( singlePindelFileAsInput() || pindelConfigFileAsInput() );
   };

   std::string getRefFilename() const
   {
      return referenceFilename;
   };
   std::string getInclusiveBedFileName() const
   {
      return inf_InclusiveBedFileName;
   };
   std::string getExclusiveBedFileName() const
   {
      return inf_ExclusiveBedFileName;
   };

   std::string getSIOutputFilename() const
   {
      return outputFilename + "_SI";
   };
   std::string getDOutputFilename() const
   {
      return outputFilename + "_D";
   };
   std::string getTDOutputFilename() const
   {
      return outputFilename + "_TD";
   };
   std::string getINVOutputFilename() const
   {
      return outputFilename + "_INV";
   };
   std::string getLIOutputFilename() const
   {
      return outputFilename + "_LI";
   };
   std::string getBPOutputFilename() const
   {
      return outputFilename + "_BP";
   };
   std::string getCloseEndOutputFilename() const
   {
      return outputFilename + "_CloseEndMapped";
   };
   std::string getASMOutputFilename() const
   {
      return outputFilename + "_ASM";
   };
   std::string getGTOutputFilename() const
   {
      return outputFilename + "_GT";
   };
   std::string getMEIOutputFilename() const
   {
      return outputFilename + "_DD";
   };
   std::string getINTOutputFilename() const
   {
      return outputFilename + "_INT";
   };
   std::string getContigOutputFilename() const
   {
      return outputFilename + "_contig";
   };
   std::string getRPOutputFilename() const
   {
      return outputFilename + "_RP";
   };
   std::string getIndelConsensusOutputFilename() const
   {
      return outputFilename + "_CINDEL";
   };
   bool loopOverAllChromosomes()
   {
      return ! getRegion()->isTargetChromosomeDefined();
   };
   SearchRegion* getRegion();

private:
   //UserDefinedSettings();
   static UserDefinedSettings* m_instance;
   SearchRegion* m_region;
};

#endif
