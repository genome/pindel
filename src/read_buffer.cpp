/* 'ReadBuffer' stores raw reads, and when it is full, invokes close-end detection to save
	those whose close ends can be mapped into the Reads() vector.

	By Eric-Wubbo Lameijer, Leiden University Medical Center, e.m.w.lameijer@lumc.nl
	31-08-11
*/

#include "read_buffer.h"
#include "pindel.h"

/* 'ReadBuffer' constructor: sets the capacity of the buffer */
ReadBuffer::ReadBuffer( int capacity, std::vector<SPLIT_READ> &referenceToFilteredReads, std::vector<SPLIT_READ> &referenceToOneEndMappedReads, const std::string& chromosome, const Chromosome& chr)
   : m_filteredReads(referenceToFilteredReads), m_OneEndMappedReads(referenceToOneEndMappedReads), m_CAPACITY( capacity), m_CHROMOSOME( chromosome ), chr(chr)
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
	// std::cout << "in flush " << std::endl;
	#pragma omp parallel for
	for (int i=0; i<m_currentsize ; i++ ) {
		// std::cout << "before GetCloseEnd " << std::endl;
		//std::map<std::string, unsigned>::iterator it = g_ReadSeq2Index.find(m_rawreads[i].UnmatchedSeq);
       
		//if (it == g_ReadSeq2Index.end()) 
		{
		//if (m_rawreads[i].MapperSplit ) {
		//	std::cout << "skip close end search" << std::endl;
		//}
            if (m_rawreads[i].MapperSplit == false)
                GetCloseEnd(/*m_CHROMOSOME*/ chr, m_rawreads[i]);
		
			if (m_rawreads[i].hasCloseEnd()) {
 				updateReadAfterCloseEndMapping(m_rawreads[i]);
               
 				#pragma omp critical 
 				{
		//			g_ReadSeq2Index.insert(std::pair<std::string, unsigned> (m_rawreads[i].UnmatchedSeq, m_filteredReads.size()));
					m_rawreads[i].SampleName2Number.insert(std::pair <std::string, unsigned> (m_rawreads[i].Tag, 1));
					m_filteredReads.push_back(m_rawreads[i]);
				}
               
			}
           	//	else {
               			//if (m_rawreads[i].Name == "@DD7DT8Q1:4:1106:17724:13906#GTACCT/1") {
              		 	//    std::cout << "m_rawreads[i] no close end" << std::endl;
               			//}
               			//#pragma omp critical
               			//m_OneEndMappedReads.push_back(m_rawreads[i]);
           	//	}
       		}
       		//else { // SampleName2Number std::map <std::string, unsigned> SampleName2Number;
           	//	#pragma omp critical 
           	//	{
		//		unsigned ReadIndex = it -> second; // m_filteredReads[ReadIndex]
               	//		std::map <std::string, unsigned>::iterator it_SampleName = m_filteredReads[ReadIndex].SampleName2Number.find(m_rawreads[i].Tag);
//
  //             			if (it_SampleName == m_filteredReads[ReadIndex].SampleName2Number.end()) {
//					//std::cout << "adding " << m_rawreads[i].Tag << "\t1" << std::endl;
  //                 			m_filteredReads[ReadIndex].SampleName2Number.insert(std::pair <std::string, unsigned> (m_rawreads[i].Tag, 1));
    //           			}
      //         			else {
	//				
	//				it_SampleName -> second++;
	//				//std::cout << "increasing " << m_rawreads[i].Tag << "\t" << it_SampleName -> second << std::endl;
	//			}
	//		}
       		//}

      		// std::cout << "after GetCloseEnd " << std::endl;

   	}
    	//std::cout << "end of flush " << std::endl;
   	m_rawreads.clear();
   	m_currentsize = 0;
    	//std::cout << "existing flush " << std::endl;
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

