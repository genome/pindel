/* 'ReadBuffer' stores raw reads, and when it is full, invokes close-end detection to save
	those whose close ends can be mapped into the Reads() vector.

	By Eric-Wubbo Lameijer, Leiden University Medical Center, e.m.w.lameijer@lumc.nl
	31-08-11
*/

#include "read_buffer.h"
#include "pindel.h"

/* 'ReadBuffer' constructor: sets the capacity of the buffer */
ReadBuffer::ReadBuffer( int capacity, std::vector<SPLIT_READ> &referenceToFilteredReads, const std::string& chromosome)
   : m_filteredReads(referenceToFilteredReads), m_CAPACITY( capacity), m_CHROMOSOME( chromosome )
{
   m_currentsize = 0;
   m_rawreads.reserve( m_CAPACITY );
}

/* 'ReadBuffer::close' clears the buffer, and closes the output stream. */
void ReadBuffer::close()
{
   flush();
   m_rawreads.clear();
}

/* WriteBuffer' destructor clears the buffer, and closes the output stream via
   'WriteBuffer::close'. */
ReadBuffer::~ReadBuffer()
{
   close();
}

/* 'ReadBuffer::flush' writes the current contents of the buffer to the
   designated output file, and then clears the contents to be ready to receive
   the next reads. */
void ReadBuffer::flush()
{
   #pragma omp parallel for
   for (int i=0; i<m_currentsize ; i++ ) {
      GetCloseEnd(m_CHROMOSOME, m_rawreads[i]);
      if (m_rawreads[i].hasCloseEnd()) {
         updateReadAfterCloseEndMapping(m_rawreads[i]);
         #pragma omp critical
         m_filteredReads.push_back(m_rawreads[i]);
      }
   }
   m_rawreads.clear();
   m_currentsize = 0;
}

/* 'ReadBuffer::addRawRead' adds the read "newRead" to the buffer, and invokes
   the flush function when the buffer reaches its capacity. */
void ReadBuffer::addRead( const SPLIT_READ& newRead )
{
   if (m_currentsize == m_CAPACITY ) {
      flush();
   }
   m_rawreads.push_back(newRead);
   m_currentsize++;
}

