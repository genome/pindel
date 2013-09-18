/* 'ReadBuffer' stores raw reads, and when it is full, invokes close-end detection to save 
	those whose close ends can be mapped into the Reads() vector.

	By Eric-Wubbo Lameijer, Leiden University Medical Center, e.m.w.lameijer@lumc.nl
	30-08-11	
*/

#ifndef READ_BUFFER_H
#define	READ_BUFFER_H

#include <vector>

#include "pindel.h"

class ReadBuffer {
   public:
      ReadBuffer( int capacity, std::vector<SPLIT_READ> &referenceToFilteredReads, std::vector<SPLIT_READ> &referenceToOneEndMappedReads, const std::string& chromosome, const Chromosome& chr);
      void flush();
      void addRead( const SPLIT_READ& newRead );
      ~ReadBuffer();
   private:
      void close();

		std::vector<SPLIT_READ> &m_filteredReads;
    std::vector<SPLIT_READ> &m_OneEndMappedReads;

      const int m_CAPACITY;
		const std::string& m_CHROMOSOME;
      const Chromosome& chr;
      int m_currentsize;
      std::vector<SPLIT_READ> m_rawreads;
};

#endif // READ_BUFFER_H


