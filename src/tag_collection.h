/* tag_collection.h

   Collects the sample tags, and returns them efficiently for calling functions.

   Created by Eric-Wubbo Lameijer, Leiden University Medical Center, e.m.w.lameijer@lumc.nl, 
   September 6th, 2011
*/

#ifndef TAG_COLLECTION_H
#define TAG_COLLECTION_H


#include <map>
#include <string>
#include <vector>

class TagCollection {

public:
   void addTag( const std::string& tag );
   int getIndexFromTag( const std::string& tag );
   std::string getTagFromIndex( const int tagIndex ) const;
	int getSize() const;

private:
   std::map< std::string, int >  m_tagMap; // for tag to index
   std::vector< std::string > m_tagVector; // double representation of same values; inelegant, but speeds up searching
};




#endif // TAG_COLLECTION_H
