/* tag_collection.cpp

   Collects the sample tags, and returns them efficiently for calling functions.

   Created by Eric-Wubbo Lameijer, Leiden University Medical Center, e.m.w.lameijer@lumc.nl, 
   September 6th, 2011
*/

#include <map>
#include <string>
#include <vector>


#include "tag_collection.h"

void TagCollection::addTag( const std::string& tag ) 
{
	std::map<std::string,int>::iterator it = m_tagMap.find( tag );
	if (it==m_tagMap.end()) { // element not in map yet, so insert it
   	m_tagMap.insert( std::pair<std::string, int>(tag, -1 ));
   	m_tagVector.clear();
		int count = 0;
   	for (std::map< std::string, int>::iterator it=m_tagMap.begin(); it!=m_tagMap.end(); it++  ) {
      	m_tagVector.push_back( it->first );
     		it->second = count;
      	count++;
   	}
	}
}

int TagCollection::getIndexFromTag( const std::string& tag) 
{
   return m_tagMap[ tag ];
}

std::string TagCollection::getTagFromIndex( const int tagIndex ) const
{
	return m_tagVector[ tagIndex ];
}

int TagCollection::getSize() const
{
	return m_tagVector.size();
}
