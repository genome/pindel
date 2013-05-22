/* UserDefinedSettings' manages the parameters that have been set by the users or have predefined default values. */

#include <cstdlib> // for exit function
#include <iostream>

#include "logstream.h"
#include "user_defined_settings.h"

#include "pindel.h"

UserDefinedSettings* UserDefinedSettings::m_instance = NULL;

/*
UserDefinedSettings::UserDefinedSettings()
{
	m_region = NULL;
}
*/

UserDefinedSettings* UserDefinedSettings::Instance() {
	if (m_instance == NULL ) {
		m_instance = new UserDefinedSettings;
	}
	return m_instance;
}

SearchRegion* UserDefinedSettings::getRegion() 
{
	if (m_region == NULL) {
		m_region = new SearchRegion( userDefinedRegion );
	}
	return m_region;
}

bool SearchRegion::isStartDefined() const
{
	return m_startDefined;
} 

bool SearchRegion::isEndDefined() const
{
	return m_endDefined;
} 

std::string uppercase( const std::string& input )
{
    std::string output = input;
    for(unsigned int pos=0; pos<input.length(); pos++ ) {
        output[ pos] = toupper( input[pos] );
    }
    return output;
}

bool SearchRegion::isTargetChromosomeDefined() const
{
	return ( uppercase(m_targetChromosomeName) != "ALL" );
} 

std::string SearchRegion::getTargetChromosomeName() const
{
	return m_targetChromosomeName;
} 

unsigned int SearchRegion::getStart() const
{
	if (!m_startDefined) {
		std::cout << "Error, the region start is requested but has not been defined!\n";
		exit( EXIT_FAILURE );
	}
	return m_start;
} 

unsigned int SearchRegion::getEnd() const
{
	if (!m_endDefined) {
		std::cout << "Error, the region end is requested but has not been defined!\n";
		exit( EXIT_FAILURE );
	}
	return m_end;
} 

/** 'eliminate' eliminates a character from the input string. */
void eliminate(const char ch, // in: character to be eliminated from the string
               std::string & str // modif: string that needs to be modified
              )
{
    size_t eliminateCharPos = str.find(ch);
    while (eliminateCharPos != std::string::npos) {
        str.erase(eliminateCharPos, 1);
        eliminateCharPos = str.find(ch);
    }
}

void SearchRegion::SetRegion(const std::string& ChrName, const unsigned int start, const unsigned int end) {
	m_targetChromosomeName = ChrName;
	m_startDefined = true;
	m_endDefined = true;
	m_start = start;
	m_end = end;
}

SearchRegion::SearchRegion(const std::string& regionString )
{

	std::cout << "SearchRegion::SearchRegion" << std::endl;
	size_t separatorPos = regionString.find(":");
	bool correctParse = false;
	m_startDefined = false;
	m_endDefined = false;
	m_start = -1;
	m_end = -1;

    // found a separator
   if (separatorPos != std::string::npos) {
      m_targetChromosomeName = regionString.substr(0, separatorPos);
      std::string coordinates = regionString.substr(separatorPos + 1);
      eliminate(',', coordinates); // removes the ',' in 1,000 or 1,000,000 that users may add for readability but wreak havoc with atoi
      size_t startEndSeparatorPos = coordinates.find("-");

      // there are two coordinates
      if (startEndSeparatorPos != std::string::npos) {
         std::string secondPositionStr = coordinates.substr( startEndSeparatorPos + 1);
         m_end = atoi(secondPositionStr.c_str());
         m_endDefined = true;
      }
      m_start = atoi(coordinates.c_str());
      m_startDefined = true;

      //LOG_DEBUG(*logStream << "sor: " << m_start << "eor: " << m_end << std::endl);
      if (m_start < 0 || (m_endDefined && ( m_end < m_start ) )) {
         correctParse = false;
      }
      else {
         correctParse = true;
      }
   }
   // no separator found
   else {
      m_targetChromosomeName = regionString;
	m_start = 1;
	m_end = 1;
	for (unsigned index = 0; index < g_ChrNameAndSizeAndIndex.size(); index++) {
		if (g_ChrNameAndSizeAndIndex[index].ChrName == m_targetChromosomeName)
			m_end = g_ChrNameAndSizeAndIndex[index].ChrSize;
	}
	//getChrSize(m_targetChromosomeName);
	m_startDefined = true;
	m_endDefined = true;
      correctParse = true;
   }

   if (!correctParse) {
     *logStream << "I cannot parse the region '" << regionString
                  << "'. Please give region in the format -c ALL, -c <chromosome_name> "
                  "(for example -c 20) or -c <chromosome_name>:<start_position>[-<end_position>], for example -c II:1,000 or "
                  "-c II:1,000-50,000. If an end position is specified, it must be larger than the start position."
                  << std::endl;
      exit ( EXIT_FAILURE);
   }
}
