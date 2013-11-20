/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// System header files
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>

// Samtools header files
#include "bam.h"
#include "sam.h"
#include "kstring.h"
#include "kseq.h"
#include "khash.h"
#include "ksort.h"

// Pindel header files
#include "logstream.h"
#include "pindel.h"
#include "reader.h"
#include "logdef.h"
#include "refreader.h"
#include "control_state.h"
#include "read_buffer.h"
#include "line_reader.h"
#include "pindel_read_reader.h"
#include "user_defined_settings.h"
// Static function declaration

static int fetch_func_SR (const bam1_t * b1, void *data);
static int fetch_func_RP (const bam1_t * b1, void *data);
static int fetch_func_RP_Discovery (const bam1_t * b1, void *data);
int32_t bam_cigar2mismatch( const bam1_core_t *readCore, const uint32_t *cigar);

const int BUFFER_SIZE = 50000;

unsigned int g_NumReadInWindow = 0; // #################
unsigned int g_NumReadScanned = 0;
unsigned int g_InWinPlus = 0; // #################
unsigned int g_InWinMinus = 0; // #################
unsigned int g_CloseMappedPlus = 0; // #################
unsigned int g_CloseMappedMinus = 0; // #################

// Reader (BamReader and PindelReader)

//init hash/maps for read pairing on the fly
KSORT_INIT_GENERIC (uint32_t) KHASH_MAP_INIT_STR (read_name, bam1_t *)

void GetReadSeq (const bam1_t* bamOfRead, std::string & c_sequence) {
	//std::string c_sequence;
	const bam1_core_t *bamCore = &bamOfRead->core;
	uint8_t *s = bam1_seq (bamOfRead);
	for (int i = 0; i < bamCore -> l_qseq; ++i) {
		c_sequence.append (1, bam_nt16_rev_table[bam1_seqi (s, i)]);
    	}	
}

struct fetch_func_data_SR {
   fetch_func_data_SR () : CurrentChrSeq( NULL ){
      LeftReads = NULL;
       OneEndMappedReads = NULL;
       RefSupportingReads = NULL;
      read_to_map_qual = NULL;
      header = NULL;
      b1_flags = NULL;
      b2_flags = NULL;   
      Tag = "";
      InsertSize = 0;
   }
   fetch_func_data_SR (const std::string* chromosomeSeq) : CurrentChrSeq( chromosomeSeq ){
      LeftReads = NULL;
       OneEndMappedReads = NULL;
       RefSupportingReads = NULL;
      read_to_map_qual = NULL;
      header = NULL;
      b1_flags = NULL;
      b2_flags = NULL;   
      Tag = "";
      InsertSize = 0;
   }
   ReadBuffer *readBuffer;
   std::vector < SPLIT_READ > *LeftReads;
   std::vector < SPLIT_READ > *OneEndMappedReads;
   std::vector <REF_READ> *RefSupportingReads;
   khash_t (read_name) * read_to_map_qual;
   bam_header_t *header;
   flags_hit *b1_flags;
   flags_hit *b2_flags;
   const std::string * CurrentChrSeq;
   std::string Tag;
   int InsertSize;
};

struct fetch_func_data_RP {
    fetch_func_data_RP () : CurrentChrSeq( NULL ){
        LeftReads = NULL;
        LeftReads_InterChr = NULL;
        read_to_map_qual = NULL;
        header = NULL;
        b1_flags = NULL;
        b2_flags = NULL;
        CurrentChrSeq = NULL;
        Tag = "";
        InsertSize = 0;
    }
   fetch_func_data_RP (const std::string* chromosomeSeq) : CurrentChrSeq( chromosomeSeq ){
      LeftReads = NULL;
       LeftReads_InterChr = NULL;

      read_to_map_qual = NULL;
      header = NULL;
      b1_flags = NULL;
      b2_flags = NULL;   
      Tag = "";
      InsertSize = 0;
   }
    //ReadBuffer *readBuffer;
    std::vector < RP_READ > *LeftReads;
    std::vector < RP_READ > *LeftReads_InterChr;
    khash_t (read_name) * read_to_map_qual;
    bam_header_t *header;
    flags_hit *b1_flags;
    flags_hit *b2_flags;
    const std::string * CurrentChrSeq;
    std::string Tag;
    int InsertSize;
};



// Return sample name for given alignment, according to known read group to sample mapping.
void get_read_group(const bam1_t* alignment, std::string& read_group) {
    uint8_t* s = bam_aux_get(alignment, "RG");
    if (s != NULL) {
        read_group = bam_aux2Z(s);
    }
}


void GetOneChrSeq (std::ifstream & fastaFile, std::string & chromosomeSequence, bool WhetherBuildUp)
{
   RefReader* rr = new RefReader(fastaFile, chromosomeSequence );
   rr->ReadSeq(WhetherBuildUp);
   delete rr;
}


double safeDivide( int dividend, int divisor )
{
   if (divisor == 0) {
      return 0;
   }
   else {
      return ((double)dividend / divisor );
   }
}

void showReadStats(const std::vector<SPLIT_READ>& Reads, const std::vector<SPLIT_READ>& OneEndMappedReads)
{
   LOG_INFO(*logStream << "Number of problematic reads in current window:            \t" << g_NumReadInWindow <<
            ", + " << g_InWinPlus << " - " << g_InWinMinus << std::endl);
   LOG_INFO(*logStream << "Number of split-reads where the close end could be mapped:\t" << Reads.size () <<
            ", + " << g_CloseMappedPlus << " - " << g_CloseMappedMinus << std::endl);
   //LOG_INFO(*logStream << "Number of hanging reads (no close end mapped):            \t" << OneEndMappedReads.size () << ", + " << g_InWinPlus - g_CloseMappedPlus << " - " << g_InWinMinus - g_CloseMappedMinus << std::endl);
   LOG_INFO(*logStream << "Percentage of problematic reads with close end mapped:    \t+ " << std::setprecision(2) << std::fixed << safeDivide( (int)(g_CloseMappedPlus * 100.0) , g_InWinPlus ) <<
            "% - " << safeDivide( (int)(g_CloseMappedMinus * 100.0) , g_InWinMinus ) << "%");
   *logStream << std::endl;
}

/** 'ReadInRead' reads in reads from Pindel input file. */
short ReadInRead (PindelReadReader & inf_ReadSeq, const std::string & FragName,
            const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads,
            const SearchWindow& currentWindow)
{
	LOG_INFO(*logStream << "Scanning and processing reads anchored in " << FragName << std::endl);
	SPLIT_READ Temp_One_Read;
	std::vector < SPLIT_READ > BufferReads;
	std::string TempQC, TempLine, TempStr, TempFragName;

	inf_ReadSeq.Reset();
	
	int UPCLOSE_COUNTER = 0;
	//loop over all reads in the file
	while (inf_ReadSeq.HasNext()) {
		Temp_One_Read = inf_ReadSeq.NextRead();
		// first character of readname should be '@'
		if (Temp_One_Read.Name[0] != FirstCharReadName) { 
			LOG_WARN(*logStream << "Something wrong with the read name: " << Temp_One_Read.Name << std::endl);
			Reads.clear ();
			return EXIT_FAILURE;
		}
		g_NumReadScanned++;

		if (Temp_One_Read.MatchedD != Minus && Temp_One_Read.MatchedD != Plus) {
			LOG_INFO(*logStream << Temp_One_Read.Name << std::endl
						<< Temp_One_Read.getUnmatchedSeq() << std::endl
						<< Temp_One_Read.MatchedD << " ..." << std::endl);
			LOG_INFO(*logStream << "+/-" << std::endl);
			return EXIT_FAILURE;
		}

		if (Temp_One_Read.MatchedRelPos > g_maxPos) {
			g_maxPos = Temp_One_Read.MatchedRelPos;
		}

		if (Temp_One_Read.FragName == FragName && Temp_One_Read.MatchedRelPos >= currentWindow.getStart()
			&& Temp_One_Read.MatchedRelPos < currentWindow.getEnd()) {
			g_NumReadInWindow++;

			if (Temp_One_Read.MatchedD == Plus) {
				g_InWinPlus++;
			}
			else {
				g_InWinMinus++; 
			}
			if (Temp_One_Read.MatchedRelPos > currentWindow.getChromosome()->getBiolSize()) { // perhaps make setter for MatchedRelpos, so these details can be factored out?
				Temp_One_Read.MatchedRelPos = currentWindow.getChromosome()->getBiolSize();
			}
			if (Temp_One_Read.MatchedRelPos < 0) {
				Temp_One_Read.MatchedRelPos = 0;
			}
			BufferReads.push_back(Temp_One_Read);
			if (BufferReads.size () >= NumberOfReadsPerBuffer) {
				#pragma omp parallel default(shared)
				{
					#pragma omp for
					// openMP 2.5 requires signed loop index
					for (int BufferReadsIndex = 0;  BufferReadsIndex < (int)NumberOfReadsPerBuffer; BufferReadsIndex++) {
						GetCloseEnd( *currentWindow.getChromosome() /*CurrentChr*/, BufferReads[BufferReadsIndex] );
					}
				}
				// EW: would it be useful to fuse this loop with the previous one?
				for (unsigned int BufferReadsIndex = 0; BufferReadsIndex < NumberOfReadsPerBuffer; BufferReadsIndex++) {
					SPLIT_READ& currentRead = BufferReads[BufferReadsIndex];
					if (currentRead.UP_Close.size ()) {
					// remove UPCLOSE_COUNTER? It seems to duplicate the effort of g_CloseMappedPlus/Minus
						UPCLOSE_COUNTER++;
						if (g_reportLength < currentRead.getReadLength()) {
							g_reportLength = currentRead.getReadLength();
						}
						// are the next two default settings? If so, we can put them in the SPLIT_READ constructor.
						currentRead.Used = false;
						currentRead.UniqueRead = true; 
	
						//                  LOG_DEBUG(*logStream << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t");
						CleanUniquePoints (currentRead.UP_Close);
						//                  LOG_DEBUG(*logStream << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << std::endl);
								
						// CloseEndLength duplicates MaxLenCloseEnd() (and is only used in reader, so does not seem time-saving) remove?
						currentRead.CloseEndLength = currentRead.UP_Close[currentRead.UP_Close.size() - 1].LengthStr;
						if (currentRead.MatchedD == Plus) {
							currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size() - 1].AbsLoc + 1 - currentRead.CloseEndLength;
						}
				                else {
							currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size() - 1].AbsLoc + currentRead.CloseEndLength -                      											currentRead.getReadLength();
						}
						// we may want to remove some of the LOG_DEBUG statements; these mainly increase code-reading time
				                LOG_DEBUG(*logStream << "Pushing back!" << std::endl);
						currentRead.SampleName2Number.insert(std::pair <std::string, unsigned> (currentRead.Tag, 1));
				                Reads.push_back (currentRead);
						if (currentRead.MatchedD == Plus) {
							g_CloseMappedPlus++;
						}
						else {
							g_CloseMappedMinus++;
						}
						g_sampleNames.insert(currentRead.Tag);
					}
					//else OneEndMappedReads.push_back(currentRead); // OneEndMappedReads
			        }
				BufferReads.clear ();
			}										// if buffer-reads threatens to overflow
		}												// if the read is in the correct bin
	}														// while loop over each read
	LOG_INFO(*logStream << "last one: " << BufferReads.size () << " and UPCLOSE= " << UPCLOSE_COUNTER << std::endl);
	// the below can be fused with the above
	#pragma omp parallel default(shared)
	{
	#pragma omp for
		for (int BufferReadsIndex = 0; BufferReadsIndex < (int)BufferReads.size (); BufferReadsIndex++) { // signed type required by OpenMP 2.5
			GetCloseEnd ( *currentWindow.getChromosome() /*CurrentChr*/, BufferReads[BufferReadsIndex]);
		}
	}
	for (unsigned int BufferReadsIndex = 0; BufferReadsIndex < BufferReads.size (); BufferReadsIndex++) {
		SPLIT_READ& currentRead = BufferReads[BufferReadsIndex];
		if (currentRead.UP_Close.size ()) {
			UPCLOSE_COUNTER++;
			if (g_reportLength < currentRead.getReadLength()) {
				g_reportLength = currentRead.getReadLength();
			}
			currentRead.Used = false;
			currentRead.UniqueRead = true; 
			//         LOG_DEBUG(*logStream << Temp_One_Read.MatchedD  << "\t" << Temp_One_Read.UP_Close.size() << "\t");

			CleanUniquePoints (currentRead.UP_Close);

			//         LOG_DEBUG(*logStream << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << std::endl);

			currentRead.CloseEndLength = currentRead.UP_Close[currentRead.UP_Close.size () - 1].LengthStr;
			if (currentRead.MatchedD == Plus) {
				currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size () - 1].AbsLoc + 1 - currentRead.CloseEndLength;
			}
			else {
				currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size () - 1].AbsLoc + currentRead.CloseEndLength - currentRead.getReadLength();
			}
			currentRead.SampleName2Number.insert(std::pair <std::string, unsigned> (currentRead.Tag, 1));
			Reads.push_back (currentRead);
			if (currentRead.MatchedD == Plus) {
				g_CloseMappedPlus++;
			}
			else {
				g_CloseMappedMinus++;
			}
			g_sampleNames.insert(currentRead.Tag);
		}
		//else OneEndMappedReads.push_back(currentRead); // OneEndMappedReads
	}
	BufferReads.clear ();

	// what do we want with this? If we want to keep it like this, we would not want to keep it a global boolean, but a static function boolean. Or for none or all chromosomes?
	if (FirstChr) {
		LOG_INFO(*logStream << std::endl << "The last read Pindel scanned: " << std::endl);
		// we may want to factor this out into a function "ReportBasicReadData" or such...
		LOG_INFO(*logStream << Temp_One_Read.Name << std::endl << Temp_One_Read.getUnmatchedSeq() << std::endl << Temp_One_Read.MatchedD << "\t"
				<< Temp_One_Read.FragName << "\t"
				<< Temp_One_Read.MatchedRelPos << "\t"
				<< Temp_One_Read.MS << "\t"
				<< Temp_One_Read.InsertSize << "\t" << Temp_One_Read.Tag << std::endl << std::endl);
		FirstChr = false;
	}
	std::vector < SPLIT_READ > HangingReads;
	showReadStats(Reads, HangingReads);

	// 0 means "success?" basically just let caller find out that there are no reads.
	if (Reads.size() == 0) {
		return 0;
	}
	LOG_INFO(*logStream << " finished!" << std::endl);
	return 0;
}

bool ReadInBamReads_RP (const char *bam_path, const std::string & FragName,
                        const std::string * CurrentChrSeq, std::vector <RP_READ> &LeftReads, 
                        int InsertSize, std::string Tag, const SearchWindow& currentWindow ) 
{ 
    bamFile fp;
    fp = bam_open (bam_path, "r");
    assert (fp);
    bam_index_t *idx;
    idx = bam_index_load (bam_path);	// load BAM index
    assert (idx);
    bam_header_t *header = bam_header_read (fp);
    bam_init_header_hash (header);
    assert (header);
    //need thing that converts "tid" to "chromosome name"
    int tid;
    tid = bam_get_tid (header, FragName.c_str ());

    //kai does the below line in readinreads. dunno why yet  
    fetch_func_data_RP data( CurrentChrSeq );
    data.header = header;
    //data.CurrentChrSeq = CurrentChrSeq;
    data.LeftReads = &LeftReads;
    //data.read_to_map_qual = NULL;
    //data.read_to_map_qual = kh_init (read_name);
    flags_hit b1_flags;//, b2_flags;
    data.b1_flags = &b1_flags;
    //data.b2_flags = &b2_flags;
    data.InsertSize = InsertSize;
    data.Tag = Tag;
    bam_fetch (fp, idx, tid, currentWindow.getStart(), currentWindow.getEnd(), &data, fetch_func_RP);
    /*
    khint_t key;
    if (kh_size (data.read_to_map_qual) > 0) {
        for (key = kh_begin (data.read_to_map_qual);
             key != kh_end (data.read_to_map_qual); ++key) {
            if (kh_exist (data.read_to_map_qual, key)) {
                bam_destroy1 (kh_value (data.read_to_map_qual, key));
                free ((char *) kh_key (data.read_to_map_qual, key));
            }
        }
    }
	*/
    kh_clear (read_name, data.read_to_map_qual);
    kh_destroy (read_name, data.read_to_map_qual);
    
    bam_header_destroy (header);
    bam_index_destroy (idx);
    bam_close (fp);
    return true;
}

bool ReadInBamReads_RP_Discovery (const char *bam_path, 
				const std::string & FragName, 
				const std::string * CurrentChrSeq, 
				std::vector <RP_READ> &LeftReads, 
				std::vector <RP_READ> &LeftReads_InterChr, 
				int InsertSize, std::string Tag, 
				const SearchWindow& currentWindow ){

	//std::cout << "Entering ReadInBamReads_RP_Discovery" << std::endl;
	bamFile fp;
	fp = bam_open (bam_path, "r");
	assert (fp);
	bam_index_t *idx;
	idx = bam_index_load (bam_path);	// load BAM index
	assert (idx);
	bam_header_t *header = bam_header_read (fp);
	bam_init_header_hash (header);
	assert (header);
	//need thing that converts "tid" to "chromosome name"
	int tid;
	tid = bam_get_tid (header, FragName.c_str ());
    
	fetch_func_data_RP data( CurrentChrSeq );
	data.header = header;
	//data.CurrentChrSeq = CurrentChrSeq;
	data.LeftReads = &LeftReads;
	data.LeftReads_InterChr = &LeftReads_InterChr;
	//data.read_to_map_qual = NULL;
	//data.read_to_map_qual = kh_init (read_name);
	flags_hit b1_flags;//, b2_flags;
	data.b1_flags = &b1_flags;
	//data.b2_flags = &b2_flags;
	data.InsertSize = InsertSize;
	data.Tag = Tag;
	//std::cout << "Before bam_fetch" << std::endl;
	bam_fetch (fp, idx, tid, currentWindow.getStart(), currentWindow.getEnd(), &data, fetch_func_RP_Discovery);
	//std::cout << "After bam_fetch" << std::endl;
	/*
	khint_t key;
	if (kh_size (data.read_to_map_qual) > 0) {
		for (key = kh_begin (data.read_to_map_qual); key != kh_end (data.read_to_map_qual); ++key) {
			if (kh_exist (data.read_to_map_qual, key)) {
				bam_destroy1 (kh_value (data.read_to_map_qual, key));
				free ((char *) kh_key (data.read_to_map_qual, key));
			}
		}
	}
	*/
	kh_clear (read_name, data.read_to_map_qual);
	kh_destroy (read_name, data.read_to_map_qual);
	
	bam_header_destroy (header);
	bam_index_destroy (idx);
	bam_close (fp);
	//std::cout << "Leaving ReadInBamReads_RP_Discovery" << std::endl;
	return true;
}


bool ReadInBamReads_SR (const char *bam_path, const std::string & FragName,
                const std::string * CurrentChrSeq,
                std::vector < SPLIT_READ > &LeftReads,
                std::vector < SPLIT_READ > &OneEndMappedReads,
                std::vector <REF_READ> &RefSupportingReads,
                int InsertSize,
                std::string Tag,
                const SearchWindow& window,
                ReadBuffer& readBuffer,
                bool verbose)
{
   //std:: cout << " in ReadInBamReads_SR " << std::endl;
   bamFile fp;
   fp = bam_open (bam_path, "r");
   assert (fp);
   bam_index_t *idx;
   idx = bam_index_load (bam_path);	// load BAM index
   assert (idx);
   bam_header_t *header = bam_header_read (fp);
   bam_init_header_hash (header);
   assert (header);
   //need thing that converts "tid" to "chromosome name"
   int tid;
   tid = bam_get_tid (header, FragName.c_str ());

   //kai does the below line in readinreads. dunno why yet
   fetch_func_data_SR data( CurrentChrSeq );
   data.header = header;
   //data.CurrentChrSeq = CurrentChrSeq;
   data.LeftReads = &LeftReads;
    data.OneEndMappedReads = &OneEndMappedReads;
    data.RefSupportingReads = &RefSupportingReads;
   data.read_to_map_qual = NULL;
   data.read_to_map_qual = kh_init (read_name);
   flags_hit b1_flags, b2_flags;
   data.b1_flags = &b1_flags;
   data.b2_flags = &b2_flags;
   data.InsertSize = InsertSize;
   data.Tag = Tag;
   data.readBuffer=&readBuffer;
   // std:: cout << " before bam_fetch " << std::endl;
	//g_ReadSeq2Index.clear();
   bam_fetch (fp, idx, tid, window.getStart(), window.getEnd(), &data, fetch_func_SR);
   // std:: cout << " after bam_fetch " << std::endl;
   readBuffer.flush();
	//g_ReadSeq2Index.clear();
	
   // std:: cout << " after flush " << std::endl;
   if (verbose) {
       showReadStats(LeftReads, OneEndMappedReads);
   }
	//std:: cout << "1 " << std::endl;
   khint_t key;
   if (kh_size (data.read_to_map_qual) > 0) {
      for (key = kh_begin (data.read_to_map_qual);
            key != kh_end (data.read_to_map_qual); ++key) {
         if (kh_exist (data.read_to_map_qual, key)) {
            bam_destroy1 (kh_value (data.read_to_map_qual, key));
            free ((char *) kh_key (data.read_to_map_qual, key));
         }
      }
   }
//std:: cout << "2 " << std::endl;
   kh_clear (read_name, data.read_to_map_qual);
   kh_destroy (read_name, data.read_to_map_qual);

   bam_header_destroy (header);
   bam_index_destroy (idx);
   bam_close (fp);
   //std:: cout << " existing ReadInBamReads_SR " << std::endl;
   return true;
}

bool isGoodAnchor( const flags_hit *read, const bam1_t * bamOfRead ) //bam1_qname
{
//return true;
	const bam1_core_t *bamCore = &bamOfRead->core;
	//std::cout << "1";
	if (bamCore->flag & BAM_FSECONDARY || bamCore->flag & BAM_FQCFAIL || bamCore->flag & BAM_FDUP) return false;
	//std::cout << "2";
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
	if (userSettings->minimalAnchorQuality == 0) return true;
	//std::cout << "3";
		//std::string NR = bam1_qname(bamOfRead);
		//std::string seq;
		//GetReadSeq(bamOfRead, seq);
		//if (NR == "HWI-ST568:267:C1BD9ACXX:4:2101:18037:80445")
		//	std::cout << "isGoodAnchor HWI-ST568:267:C1BD9ACXX:4:2101:18037:80445 " << seq << std::endl;

//return true;
	
	
	//int maxEdits = int (bamCore->l_qseq * userSettings->MaximumAllowedMismatchRate) + 1;
/*      
	const uint8_t *nm = bam_aux_get(bamOfRead, "NM");
	if (nm) {
		int32_t nm_value = bam_aux2i(nm);
		//std::cout << bam1_qname(bamOfRead) << std::endl;
		std::string NR = bam1_qname(bamOfRead);
		std::string seq;
		GetReadSeq(bamOfRead, seq);
		if (NR == "HWI-ST568:267:C1BD9ACXX:4:2101:18037:80445")
			std::cout << "isGoodAnchor " << nm_value << " " << userSettings->NM << " " << seq << std::endl;
		//std::cout << "isGoodAnchor " << (*nm) << " " << userSettings->NM << std::endl;
		if ((nm_value > maxEdits) || (nm_value > userSettings->NM)) {
			return false;
		}
	}

*/
	//unsigned int mappingQuality = bamCore->qual;

	if (!read->mapped) return false;
	//std::cout << "4";
	if (bamCore->qual < userSettings->minimalAnchorQuality) return false;
	//std::cout << "5";
	return true;
/*
	if (read->unique || read->sw) return false;

	

	if (read->edits > maxEdits) return false;
	if (read->suboptimal) return false;
*/
}

bool isRefRead ( const flags_hit *read, const bam1_t * bamOfRead )
{
const bam1_core_t *bamCore = &bamOfRead->core;

if (bamCore->flag & BAM_FSECONDARY || bamCore->flag & BAM_FQCFAIL || bamCore->flag & BAM_FDUP) return false;

	//std::cout << "isRefRead 1" << std::endl;
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
	//std::cout << "isRefRead 2" << std::endl;
	const uint8_t *nm = bam_aux_get(bamOfRead, "NM");
	//std::cout << "isRefRead 3" << std::endl;
	//const bam1_core_t *bamCore = &bamOfRead->core;
	int maxEdits = int (bamCore->l_qseq * userSettings->MaximumAllowedMismatchRate) + 1;
	//std::cout << "isRefRead 4" << std::endl;
	uint32_t *cigar_pointer = bam1_cigar (bamOfRead);
	int cigarMismatchedBases = bam_cigar2mismatch (bamCore, cigar_pointer);
    	//std::cout << "isRefRead 5" << std::endl;

	if (nm) {
		int32_t nm_value = bam_aux2i(nm);
		//std::string NR = bam1_qname(bamOfRead);
		//std::string seq;
		//GetReadSeq(bamOfRead, seq);
		//if (NR == "HWI-ST568:267:C1BD9ACXX:4:2101:18037:80445")
		//	std::cout << "isRefRead " << nm_value << " " << userSettings->NM << " " << seq << std::endl;
		if ((nm_value > userSettings->NM) || (nm_value > maxEdits)) {
			return false;
		}
	}

	if (read->mapped && read->edits <= 2 && cigarMismatchedBases <= 2) {
		return true;
	}
	else {
		return false;
	}
}

bool isWeirdRead( const flags_hit *read, const bam1_t * bamOfRead )
{

const bam1_core_t *bamCore = &bamOfRead->core;

if (bamCore->flag & BAM_FSECONDARY || bamCore->flag & BAM_FQCFAIL || bamCore->flag & BAM_FDUP) return false;

	if (!(read->mapped)) return true;

	const uint8_t *nm = bam_aux_get(bamOfRead, "NM");
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();

	if (nm) {
		int32_t nm_value = bam_aux2i(nm);
		//std::string NR = bam1_qname(bamOfRead);
		//std::string seq;
		//GetReadSeq(bamOfRead, seq);
		//if (NR == "HWI-ST568:267:C1BD9ACXX:4:2101:18037:80445")
		//	std::cout << "isWeirdRead " << nm_value << " " << userSettings->NM << " " << seq << std::endl;
		if (nm_value) return true;
	}
	
	
	//const bam1_core_t *bamCore;
	bamCore = &bamOfRead->core; 

	//int maxEdits = int (bamCore->l_qseq * userSettings->MaximumAllowedMismatchRate) + 1;

	uint32_t *cigar_pointer = bam1_cigar (bamOfRead);
	int cigarMismatchedBases = bam_cigar2mismatch (&bamOfRead->core, cigar_pointer);

	if ( read->edits + cigarMismatchedBases > 0) {
		return true;
	}

	return false;
}

bool WhetherSimpleCigar(const bam1_t * unmapped_read, const bam1_core_t * c, const uint32_t * cigar, SPLIT_READ & SR_Read, int & indelsize) {
	//std::cout << "entering WhetherSimpleCigar" << std::endl;
	const uint8_t *nm = bam_aux_get(unmapped_read, "NM");
	if (nm) {
		int32_t nm_value = bam_aux2i(nm);
		if (nm_value) return false;
	}
        uint32_t k;
        unsigned CountNonM = 0;
        unsigned CountIndel = 0;
        unsigned CountM = 0;
//	if (SR_Read.Name == "@1_13911_14618_0_1_0_0_1:0:0_1:0:0_2e3d0/1" || SR_Read.Name == "@1_13911_14618_0_1_0_0_1:0:0_1:0:0_2e3d0/2")
//		std::cout << "Standard: \tMatch " << BAM_CMATCH << "\tINS " << BAM_CINS << "\tDEL " << BAM_CDEL << std::endl;
        
        for (k = 0; k < c->n_cigar; ++k) {
		
            int op = cigar[k] & BAM_CIGAR_MASK;
//		if (SR_Read.Name == "@1_13911_14618_0_1_0_0_1:0:0_1:0:0_2e3d0/1" || SR_Read.Name == "@1_13911_14618_0_1_0_0_1:0:0_1:0:0_2e3d0/2")
//			std::cout << k << "\t" << op << "\t" << (cigar[k] >> BAM_CIGAR_SHIFT) << std::endl;
            if (op == BAM_CMATCH) CountM++;
            else CountNonM++;
            if (op == BAM_CINS || op == BAM_CDEL) CountIndel++;
        }
//	if (SR_Read.Name == "@1_13911_14618_0_1_0_0_1:0:0_1:0:0_2e3d0/1" || SR_Read.Name == "@1_13911_14618_0_1_0_0_1:0:0_1:0:0_2e3d0/2")
//		std::cout << "##################\t" << CountM << " " << CountNonM << " " << CountIndel << std::endl;
        if (CountM == 2 && CountNonM == 1 && CountIndel == 1) {
//		std::cout << "good splitmapper" << std::endl;
		int op = cigar[1] & BAM_CIGAR_MASK;
		if (op == BAM_CDEL)
			indelsize = (cigar[1] >> BAM_CIGAR_SHIFT) * (-1);
		else if (op == BAM_CINS)
			indelsize = (cigar[1] >> BAM_CIGAR_SHIFT);
		return true; // this read just contains one indel mapped by the aligner
	}
        else return false;
}

void AddUniquePoint(const bam1_t * unmapped_read, const bam1_core_t * c, const uint32_t * cigar, SPLIT_READ & SR_Read, const int & IndelSize) {
    //unsigned LeftLength, RightLength, MiddleLength;
    std::string IndelType;
    if (c->n_cigar != 3) return;
    int32_t LeftLength = cigar[0] >> BAM_CIGAR_SHIFT;
    //int32_t VariantLength = cigar[1] >> BAM_CIGAR_SHIFT;
    int32_t RightLength = cigar[2] >> BAM_CIGAR_SHIFT;
    int op = cigar[1] & BAM_CIGAR_MASK;
    if (op == BAM_CINS) IndelType = "I";
    else if (op == BAM_CDEL) IndelType = "D";
    else return;
	if (c -> qual < (short)userSettings->minimalAnchorQuality) return;
    if (SR_Read.MatchedD == '+') {
        UniquePoint TempOneClose(g_genome.getChr(SR_Read.FragName), LeftLength, c->pos + LeftLength + g_SpacerBeforeAfter - 1, FORWARD, ANTISENSE, 0);
        UniquePoint TempOneFar(g_genome.getChr(SR_Read.FragName), RightLength, c->pos + SR_Read.getReadLength() - RightLength + g_SpacerBeforeAfter - IndelSize, BACKWARD, ANTISENSE, 0);
        SR_Read.UP_Close.push_back(TempOneClose);
        SR_Read.UP_Far.push_back(TempOneFar);
    }
    else { // -
        UniquePoint TempOneFar(g_genome.getChr(SR_Read.FragName), LeftLength, c->pos + LeftLength + g_SpacerBeforeAfter - 1, FORWARD, SENSE, 0);
        UniquePoint TempOneClose(g_genome.getChr(SR_Read.FragName), RightLength, c->pos + SR_Read.getReadLength() - RightLength + g_SpacerBeforeAfter - IndelSize, BACKWARD, SENSE, 0);
        SR_Read.UP_Close.push_back(TempOneClose);
        SR_Read.UP_Far.push_back(TempOneFar);
    }

	//if (SR_Read.Name == "@GAIIX-599_0004:5:118:12447:21053/2") {
	//	std::cout << "IndelSize " << IndelSize << std::endl;
	//	std::cout << SR_Read;
	//}


    SR_Read.MapperSplit = true;
	g_NumberOfGapAlignedReads++;
	//std::cout << "mapper split" << std::endl;
    //UniquePoint
    //SortedUniquePoints UP_Close
    //SortedUniquePoints UP_Far;
}

void build_record_SR (const bam1_t * mapped_read, const bam1_t * unmapped_read, void *data)
{   // userSettings->minimalAnchorQuality
    //UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
    SPLIT_READ Temp_One_Read;
    fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
    bam_header_t *header = (bam_header_t *) data_for_bam->header;
    //std::string & CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq; //
    std::string Tag = (std::string) data_for_bam->Tag;
    int InsertSize = (int) data_for_bam->InsertSize;
    
    const bam1_core_t *mapped_core;
    const bam1_core_t *unmapped_core;
    mapped_core = &mapped_read->core;
    unmapped_core = &unmapped_read->core;
    Temp_One_Read.MS = mapped_core->qual;
    if ((short)Temp_One_Read.MS < (short)userSettings->minimalAnchorQuality) return;
    Temp_One_Read.Name = "@";
    Temp_One_Read.Name.append ((const char *) bam1_qname (unmapped_read));
    if (unmapped_core->flag & BAM_FREAD1) {
        Temp_One_Read.Name.append ("/1");
    }
    else if (unmapped_core->flag & BAM_FREAD2) {
        Temp_One_Read.Name.append ("/2");
    }
    
    // Determine sample name for read.
    get_read_group(mapped_read, Temp_One_Read.read_group);
    
    //if (Temp_One_Read.Name == "@DD7DT8Q1:4:1106:17724:13906#GTACCT/1") {
    //    std::cout << "I am here." << std::endl;
    //}
    std::string c_sequence;
    int i;
    uint8_t *s = bam1_seq (unmapped_read);
    for (i = 0; i < unmapped_core->l_qseq; ++i) {
        c_sequence.append (1, bam_nt16_rev_table[bam1_seqi (s, i)]);
    }
    //rudimentary n filter
    int length = unmapped_core->l_qseq;
    while (c_sequence[0] == 'N') {
        c_sequence.erase (0, 1);
        length--;
    }
    if (c_sequence.size () > 0) {
        while (c_sequence[length - 1] == 'N') {
            c_sequence.erase (length - 1, 1);
            length--;
        }
    }
    int n_count = 0;
    size_t found = c_sequence.find ('N', 0);
    int max_ns = (int)(length * .10);
    while (found != std::string::npos) {
        n_count++;
        found = c_sequence.find ('N', found + 1);
    }
    if (n_count > max_ns || length < 22) {
        return;
    }

    if (unmapped_core->flag & BAM_FREVERSE) {
        Temp_One_Read.setUnmatchedSeq( ReverseComplement (c_sequence) );
    }
    else {
        Temp_One_Read.setUnmatchedSeq( c_sequence );
    }
    Temp_One_Read.MatchedRelPos = mapped_core->pos;
    uint32_t *cigar_pointer_mapped = bam1_cigar (mapped_read);
    uint32_t *cigar_pointer_unmapped = bam1_cigar (unmapped_read);
    if (mapped_core->flag & BAM_FREVERSE) {
        Temp_One_Read.MatchedD = '-';
        
        int Rlength = bam_cigar2len (mapped_core, cigar_pointer_mapped);
        Temp_One_Read.MatchedRelPos += Rlength + InsertSize;
    }
    else {
        Temp_One_Read.MatchedD = '+';
        Temp_One_Read.MatchedRelPos -= InsertSize;
    }

    //FIXME pass these through from the command line with a struct
    Temp_One_Read.InsertSize = InsertSize;
	if (InsertSize <= length ) {
		*logStream << "Error: the insert size is only " << InsertSize << " while the read length is " << length << " " << Temp_One_Read.getReadLength() << std::endl;
        *logStream << "in paired end sequencing, the insert size is the total size of the fragment to be sequenced, with a read length of 100 bp the entire fragment may for example look like\n\n";
		*logStream << "|----100 bp: first read of the read pair--|-------------------300 bp: unsequenced DNA---------------------|----100 bp: second read of the read pair--|\n";
        *logStream << "<-----------------------------------------------------------insert size=500 ------------------------------------------------------------------------->\n\n";
		*logStream << "In the configuration file (the -i option) please check/correct the insert size (second item on each line). If you continue to have problems, please contact us (Kai Ye, kye@genome.wustl.edu)\n";
		exit( EXIT_FAILURE );
	}
    Temp_One_Read.Tag = Tag;

    Temp_One_Read.FragName = header->target_name[mapped_core->tid];
    //Temp_One_Read.FragName = FragName;
//    LOG_DEBUG(*logStream << Temp_One_Read.Name << std::endl
//              << Temp_One_Read.getUnmatchedSeq() << std::endl);
//    LOG_DEBUG(*logStream << Temp_One_Read.MatchedD <<
//              "\t" << Temp_One_Read.FragName <<
//              "\t" << Temp_One_Read.MatchedRelPos <<
//              "\t" << Temp_One_Read.MS <<
//              "\t" << Temp_One_Read.InsertSize <<
//              "\t" << Temp_One_Read.Tag << std::endl);
    
    
    g_NumReadInWindow++;
    
    if (Temp_One_Read.MatchedD == Plus) {
        g_InWinPlus++;
    }
    else {
        g_InWinMinus++;
    }
    if (Temp_One_Read.MatchedRelPos > data_for_bam->CurrentChrSeq->size()-2 * g_SpacerBeforeAfter) {
        Temp_One_Read.MatchedRelPos = data_for_bam->CurrentChrSeq->size()-2 * g_SpacerBeforeAfter;
    }
    if (Temp_One_Read.MatchedRelPos < 1) {
        Temp_One_Read.MatchedRelPos = 0;
    }
    //if (Temp_One_Read.Name == "@DD7DT8Q1:4:1106:17724:13906#GTACCT/1") {
    //    std::cout << "I am there." << std::endl;
    //}
	int IndelSize;
    if (WhetherSimpleCigar(unmapped_read, unmapped_core, cigar_pointer_unmapped, Temp_One_Read, IndelSize)) {
        AddUniquePoint(unmapped_read, unmapped_core, cigar_pointer_unmapped, Temp_One_Read, IndelSize);
    }

    data_for_bam->readBuffer->addRead(Temp_One_Read);
    return;
}


void build_record_RefRead (const bam1_t * mapped_read, const bam1_t * ref_read, void *data)
{   // std::vector <REF_READ> *RefSupportingReads;
   // UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
    REF_READ One_RefRead;
    fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
    bam_header_t *header = (bam_header_t *) data_for_bam->header;
    //std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
    //std::string Tag = (std::string) data_for_bam->Tag;
    //int InsertSize = (int) data_for_bam->InsertSize;
    
    const bam1_core_t *mapped_core;
    const bam1_core_t *Ref_read_core;
    mapped_core = &mapped_read->core;
    Ref_read_core = &ref_read->core;
    
    
    if (Ref_read_core->qual < userSettings->minimalAnchorQuality) return;

    One_RefRead.Pos = Ref_read_core->pos; // mapped_core->pos;
    One_RefRead.Tag = (std::string) data_for_bam->Tag;
    One_RefRead.MQ = Ref_read_core->qual;
    One_RefRead.FragName = header->target_name[Ref_read_core->tid];
    One_RefRead.ReadLength = Ref_read_core->l_qseq;

    data_for_bam->RefSupportingReads->push_back(One_RefRead);// addRead(Temp_One_Read);
    return;
}

void build_record_RP (const bam1_t * r1, void *data)
{
    
    const bam1_core_t * r1_core;
    //const bam1_core_t * r2_core;
    r1_core = &r1->core;
    //r2_core = &r2->core;
    


    RP_READ Temp_One_Read;
    fetch_func_data_RP *data_for_bam = (fetch_func_data_RP *) data;
    bam_header_t *header = (bam_header_t *) data_for_bam->header;
    //std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
    std::string Tag = (std::string) data_for_bam->Tag;
    
    if (!(r1_core->flag & BAM_FUNMAP || r1_core->flag & BAM_FMUNMAP)) { // both reads are mapped.
        if ((r1_core->tid != r1_core->mtid) || abs(r1_core->isize) > r1_core->l_qseq + 2 * data_for_bam->InsertSize) {
            Temp_One_Read.ReadName = "";
            Temp_One_Read.ReadName.append ((const char *) bam1_qname (r1));
            if (r1_core->flag & BAM_FREVERSE) {
                Temp_One_Read.DA = '-';
            }
            else Temp_One_Read.DA = '+';
            if (r1_core->flag & BAM_FMREVERSE) {
                Temp_One_Read.DB = '-';
            }
            else Temp_One_Read.DB = '+';
            
            Temp_One_Read.PosA = r1_core->pos;
            if (Temp_One_Read.DA == '+') Temp_One_Read.PosA = Temp_One_Read.PosA;// + r1_core->l_qseq;
            //else Temp_One_Read.PosA = Temp_One_Read.PosA;
            Temp_One_Read.PosB = r1_core->mpos;
            if (Temp_One_Read.DB == '+') Temp_One_Read.PosB = Temp_One_Read.PosB;// + r1_core->l_qseq;
            //else Temp_One_Read.PosB = Temp_One_Read.PosB;
            //std::cout << Temp_One_Read.ReadName << " " << Temp_One_Read.DA << " " << Temp_One_Read.PosA << " " << Temp_One_Read.DB << " " << Temp_One_Read.PosB << std::endl;
            Temp_One_Read.MQA = r1_core->qual;
            Temp_One_Read.MQB = r1_core->qual;
            Temp_One_Read.ChrNameA = header->target_name[r1_core->tid];
            Temp_One_Read.ChrNameB = header->target_name[r1_core->mtid];
            Temp_One_Read.InsertSize = data_for_bam->InsertSize;
            //FIXME pass these through from the command line with a struct
            Temp_One_Read.Tag = Tag;
            
            data_for_bam->LeftReads->push_back(Temp_One_Read);
        }
    }
    return;
}

void build_record_RP_Discovery (const bam1_t * r1, void *data) {

	//UserDefinedSettings *userSettings = UserDefinedSettings::Instance();
	//std::cout << "entering build_record_RP_Discovery" << std::endl;
	const bam1_core_t * r1_core;
	//const bam1_core_t * r2_core;
	r1_core = &r1->core;
	//r2_core = &r2->core;
    
	RP_READ Temp_One_Read;
	fetch_func_data_RP *data_for_bam = (fetch_func_data_RP *) data;
	bam_header_t *header = (bam_header_t *) data_for_bam->header;
	//std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
	std::string Tag = (std::string) data_for_bam->Tag;
	//std::string read_name = bam1_qname (r1);
    
	if (!(r1_core->flag & BAM_FPAIRED)) return;

	if (r1_core->qual < userSettings->minimalAnchorQuality) return;
    
	if (!(r1_core->flag & BAM_FUNMAP || r1_core->flag & BAM_FMUNMAP)) { // both reads are mapped.
		/*std::cout << "Read Insertsize=" << abs(r1_core->isize) << ", qseq=" << r1_core->l_qseq << " BAM ISize=" << data_for_bam->InsertSize 
			<< "Total is " << r1_core->l_qseq *2 + 2 * data_for_bam->InsertSize << "\n";*/
		//if ("read_10038" == read_name || "read_100478" ==  read_name) {
		//    std::cout << "########print " << read_name << " is on line 854." << std::endl;
		//    std::cout << r1_core->tid << " " << r1_core->mtid << " " << abs(r1_core->isize) << " "
		//    << (r1_core->flag & BAM_FREVERSE) << " "
		//    << (r1_core->flag & BAM_FMREVERSE) << std::endl;
		//}
		if ((r1_core->tid != r1_core->mtid) || abs(r1_core->isize) > 3 * data_for_bam->InsertSize + 1000 || ((bool)(r1_core->flag & BAM_FREVERSE) == (bool)(r1_core->flag & BAM_FMREVERSE))) { // different chr or same strand or insert is too large
			//std::cout << "passed the test\n";
			//if ("read_10038" == read_name) {
			//    std::cout << "print " << read_name << " is on line 859." << std::endl;
			//}
			if (r1_core->flag & BAM_FREVERSE) {
				Temp_One_Read.DA = '-';
			}
			else Temp_One_Read.DA = '+';
			if (r1_core->flag & BAM_FMREVERSE) {
				Temp_One_Read.DB = '-';
			}
			else Temp_One_Read.DB = '+';
	            
			Temp_One_Read.PosA = r1_core->pos;
			//if (Temp_One_Read.DA == '+') Temp_One_Read.PosA = Temp_One_Read.PosA;// + r1_core->l_qseq;
			//else Temp_One_Read.PosA = Temp_One_Read.PosA;
			Temp_One_Read.PosB = r1_core->mpos;
			Temp_One_Read.OriginalPosA = Temp_One_Read.PosA;
			Temp_One_Read.OriginalPosB = Temp_One_Read.PosB;
			//if (Temp_One_Read.DB == '+') Temp_One_Read.PosB = Temp_One_Read.PosB;// + r1_core->l_qseq;
	            
			//if (r1_core->tid == r1_core->mtid && Temp_One_Read.PosA < Temp_One_Read.PosB) return;
	            
			Temp_One_Read.Experimental_InsertSize = data_for_bam->InsertSize;
			Temp_One_Read.ReadName = "";
			Temp_One_Read.ReadName.append ((const char *) bam1_qname (r1));
			//if ("read_10038" == Temp_One_Read.ReadName) {
			//    std::cout << "print " << Temp_One_Read.ReadName << " is on line 881." << std::endl;
			//}
			//else Temp_One_Read.PosB = Temp_One_Read.PosB;
			//std::cout << Temp_One_Read.ReadName << " " << Temp_One_Read.DA << " " << Temp_One_Read.PosA << " " << Temp_One_Read.DB << " " << Temp_One_Read.PosB << std::endl;
			Temp_One_Read.MQA = r1_core->qual;
			Temp_One_Read.MQB = r1_core->qual;
			Temp_One_Read.ChrNameA = header->target_name[r1_core->tid];
			Temp_One_Read.ChrNameB = header->target_name[r1_core->mtid];
			Temp_One_Read.InsertSize = data_for_bam->InsertSize;
			//FIXME pass these through from the command line with a struct
			Temp_One_Read.Tag = Tag;
			Temp_One_Read.ReadLength = r1_core->l_qseq;
			//std::cout << Temp_One_Read.ReadName << " " << Temp_One_Read.DA << " " << Temp_One_Read.PosA << " " << Temp_One_Read.DB << " " << Temp_One_Read.PosB << std::endl;
			if (r1_core->tid == r1_core->mtid) {
				if (Temp_One_Read.PosA < Temp_One_Read.PosB) 
					data_for_bam->LeftReads->push_back(Temp_One_Read);
				else {
					RP_READ Temp_One_Read_another = Temp_One_Read;
					Temp_One_Read_another.DA = Temp_One_Read.DB;
					Temp_One_Read_another.DB = Temp_One_Read.DA;
					Temp_One_Read_another.PosA = Temp_One_Read.PosB;
					Temp_One_Read_another.PosB = Temp_One_Read.PosA;
					Temp_One_Read_another.OriginalPosA = Temp_One_Read.OriginalPosB;
					Temp_One_Read_another.OriginalPosB = Temp_One_Read.OriginalPosA;
					Temp_One_Read_another.ChrNameA = Temp_One_Read.ChrNameB;
					Temp_One_Read_another.ChrNameB = Temp_One_Read.ChrNameA;
					data_for_bam->LeftReads->push_back(Temp_One_Read_another);
				}	
				/*
				RP_READ Temp_One_Read_another = Temp_One_Read;
				Temp_One_Read_another.DA = Temp_One_Read.DB;
				Temp_One_Read_another.DB = Temp_One_Read.DA;
				Temp_One_Read_another.PosA = Temp_One_Read.PosB;
				Temp_One_Read_another.PosB = Temp_One_Read.PosA;
				Temp_One_Read_another.OriginalPosA = Temp_One_Read.OriginalPosB;
				Temp_One_Read_another.OriginalPosB = Temp_One_Read.OriginalPosA;
				Temp_One_Read_another.ChrNameA = Temp_One_Read.ChrNameB;
				Temp_One_Read_another.ChrNameB = Temp_One_Read.ChrNameA;
				data_for_bam->LeftReads->push_back(Temp_One_Read_another);
				*/
			}
			else    data_for_bam->LeftReads_InterChr->push_back(Temp_One_Read);
			//	std::cout << "finished the procedure\n";
		}

	}
	//std::cout << "leaving build_record_RP_Discovery" << std::endl;
	return;
}

static int fetch_func_SR (const bam1_t * b1, void *data)
{
//	std::cout << "in fetch_func_SR " << std::endl;
   g_NumReadScanned++;

   fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
   khash_t (read_name) * read_to_map_qual =
      (khash_t (read_name) *) data_for_bam->read_to_map_qual;
   flags_hit *b1_flags = data_for_bam->b1_flags;
   flags_hit *b2_flags = data_for_bam->b2_flags;
   //const std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
   //SPLIT_READ Temp_One_Read;
   const bam1_core_t *b1_core;
   bam1_t *b2;
   bam1_core_t *b2_core;
   b1_core = &b1->core;
   //std::string read_name = bam1_qname (b1);
   // if (read_name == "DD7DT8Q1:4:1106:17724:13906#GTACCT") std::cout << " +++++++ here ++++++++ " << std::endl;
   khint_t key = kh_get (read_name, read_to_map_qual, bam1_qname (b1));

   if (key == kh_end (read_to_map_qual)) {
      int ret=0;
      key = kh_put (read_name, read_to_map_qual, strdup (bam1_qname (b1)), &ret);
	//	key = kh_put (read_name, read_to_map_qual, bam1_qname (b1), &ret);
      kh_value (read_to_map_qual, key) = bam_dup1 (b1);
      return 0;
   }
   else {
      b2 = bam_dup1 (kh_value (read_to_map_qual, key));
      bam_destroy1 (kh_value (read_to_map_qual, key));
      b2_core = &b2->core;
      //this seems stupid, but in order to manage the read names, necessary
      free ((char *) kh_key (read_to_map_qual, key));
      kh_del (read_name, read_to_map_qual, key);
      //std::string c_sequence;
   }
	/*
	std::string RN = bam1_qname(b1);
	//std::cout << "|" << RN << "|" << std::endl;
	std::string query = "1_112673_113350_0_1_0_0_0:0:0_0:0:0_c9d04"; 
	if (RN == query)
			std::cout << "################################################## found " << query << std::endl;
	*/

   parse_flags_and_tags (b1, b1_flags);
   parse_flags_and_tags (b2, b2_flags);

//std::cout << "isGoodAnchor( b1_flags, b1)" << std::endl;

	if (isGoodAnchor( b1_flags, b1)) {
	//	if (RN == query) std::cout << "isGoodAnchor" << std::endl;
		if (isWeirdRead( b2_flags, b2 )) {
	//		if (RN == query) std::cout << "################################################## found b1 b2 " << query << std::endl;
			build_record_SR (b1, b2, data);
		}
		if (isRefRead( b2_flags, b2 )) {
	//		if (RN == query) std::cout << "################################################## found ref b2 b1 " << query << std::endl;
	        	build_record_RefRead (b1, b2, data);
			//std::cout << "refread 1" << std::endl;
		}
	}



//std::cout << "isGoodAnchor( b2_flags, b2)" << std::endl;
	if (isGoodAnchor( b2_flags, b2) ) {
	//	if (RN == query) std::cout << "isGoodAnchor" << std::endl;
		if (isWeirdRead( b1_flags, b1 )) {
	//		if (RN == query) std::cout << "################################################## found b2 b1 " << query << std::endl;
			build_record_SR (b2, b1, data);
		}
		if (isRefRead(b1_flags, b1 )) {
	//		if (RN == query) std::cout << "################################################## found ref b1 b2 " << query << std::endl;
			build_record_RefRead(b2, b1, data);
		}
	}


//std::cout << "before bam_destroy" << std::endl;
	bam_destroy1 (b2);
	//std::cout << "existing fetch_func_SR " << std::endl;
	return 0;
}

static int fetch_func_RP (const bam1_t * b1, void *data)
{
    g_NumReadScanned++;
    fetch_func_data_RP *data_for_bam = (fetch_func_data_RP *) data;
    //khash_t (read_name) * read_to_map_qual =
    //(khash_t (read_name) *) data_for_bam->read_to_map_qual;
    flags_hit *b1_flags = data_for_bam->b1_flags;
    //flags_hit *b2_flags = data_for_bam->b2_flags;
    //const std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;

    RP_READ Temp_One_Read;
    const bam1_core_t *b1_core;
    //bam1_t *b2;
    //bam1_core_t *b2_core;
    b1_core = &b1->core;
    //std::string read_name = bam1_qname (b1);
    /*
    khint_t key = kh_get (read_name, read_to_map_qual, bam1_qname (b1));
    if (key == kh_end (read_to_map_qual)) {
        int ret=0;
        key = kh_put (read_name, read_to_map_qual, strdup (bam1_qname (b1)), &ret);
        kh_value (read_to_map_qual, key) = bam_dup1 (b1);
        return 0;
    }
    else {
        b2 = bam_dup1 (kh_value (read_to_map_qual, key));
        bam_destroy1 (kh_value (read_to_map_qual, key));
        b2_core = &b2->core;
        //this seems stupid, but in order to manage the read names, necessary
        free ((char *) kh_key (read_to_map_qual, key));
        kh_del (read_name, read_to_map_qual, key);
    }
    */ 

    parse_flags_and_tags (b1, b1_flags);
    //parse_flags_and_tags (b2, b2_flags);

        build_record_RP (b1, data);

    //bam_destroy1 (b2);
    
    return 0;
}

static int fetch_func_RP_Discovery (const bam1_t * b1, void *data)
{   //std::cout << "entering fetch_func_RP_Discovery" << std::endl;
    g_NumReadScanned++;
    fetch_func_data_RP *data_for_bam = (fetch_func_data_RP *) data;
    //khash_t (read_name) * read_to_map_qual =
    //(khash_t (read_name) *) data_for_bam->read_to_map_qual;
    flags_hit *b1_flags = data_for_bam->b1_flags;
    //flags_hit *b2_flags = data_for_bam->b2_flags;
    //const std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
    
    RP_READ Temp_One_Read;
    const bam1_core_t *b1_core;
    //bam1_t *b2;
    //bam1_core_t *b2_core;
    b1_core = &b1->core;
    //std::string read_name = bam1_qname (b1);
    //if ("read_10038" == read_name) {
    //    std::cout << "see read_10038 in fetch_func_RP_Discovery" << std::endl;
    //}
    /*
     khint_t key = kh_get (read_name, read_to_map_qual, bam1_qname (b1));
     if (key == kh_end (read_to_map_qual)) {
     int ret=0;
     key = kh_put (read_name, read_to_map_qual, strdup (bam1_qname (b1)), &ret);
     kh_value (read_to_map_qual, key) = bam_dup1 (b1);
     return 0;
     }
     else {
     b2 = bam_dup1 (kh_value (read_to_map_qual, key));
     bam_destroy1 (kh_value (read_to_map_qual, key));
     b2_core = &b2->core;
     //this seems stupid, but in order to manage the read names, necessary
     free ((char *) kh_key (read_to_map_qual, key));
     kh_del (read_name, read_to_map_qual, key);
     }
     */
    
    parse_flags_and_tags (b1, b1_flags);
    //parse_flags_and_tags (b2, b2_flags);
    //std::cout << "before build_record_RP_Discovery" << std::endl;
    
    build_record_RP_Discovery (b1, data);
    //std::cout << "after build_record_RP_Discovery" << std::endl;
    
    //bam_destroy1 (b2);
    //std::cout << "leaving fetch_func_RP_Discovery" << std::endl;
    return 0;
}

/* 'isInBin' returns whether the read "read" is in the designated bin. */
/*bool isInBin (const SPLIT_READ & read)
{
   if (read.MatchedRelPos > g_maxPos) {
      g_maxPos = read.MatchedRelPos;
   }
   return (((signed int)read.MatchedRelPos >= (signed int)(g_binIndex * WINDOW_SIZE)) &&
           ((signed int)read.MatchedRelPos < (signed int)((g_binIndex + 1) * WINDOW_SIZE)));
}*/


void parse_flags_and_tags (const bam1_t * b, flags_hit * flags)
{
   const bam1_core_t *c = &b->core;
   char xt_code = 0;
   int mf_code = 0, nm_code = 0, best_hits = 0;
   flags->unique = 0;
   flags->mapped = !(c->flag & BAM_FUNMAP);

   uint8_t *s = bam_aux_get (b, "XT");
   if (s != 0) {
      xt_code = bam_aux2A (s);
      if (xt_code == 'U') {
         flags->unique = 1;
      }
      else {
         flags->unique = 0;
      }

   }
   else {
      flags->unique = 1;
   }

   s = NULL;
   s = bam_aux_get (b, "X0");
   if (s != 0) {
      best_hits = bam_aux2i (s);
   }
   s = NULL;
   s = bam_aux_get (b, "X1");

   if (s != 0) {
      int sub_hits = bam_aux2i (s);

      if (best_hits + sub_hits == 1) {
         flags->suboptimal = 0;
      }
      else {
         flags->suboptimal = 1;
      }
   }
   else {
      flags->suboptimal = 0;
   }

   if (xt_code == 'M' || mf_code == 130) {
      flags->sw = 1;
      //short term fix to unset unique if the maq read was s-w mapped. bwa can't set U and M at once.
   }
   else {
      flags->sw = 0;
   }
   s = NULL;
   s = bam_aux_get (b, "NM");
   if (s != 0) {
      nm_code = bam_aux2i (s);
      flags->edits = nm_code;
   }
   else {
      nm_code = 0;
      flags->edits = nm_code;
   }

   return;
}


int32_t bam_cigar2len (const bam1_core_t * c, const uint32_t * cigar)
{
   uint32_t k;
   int32_t l = 0;
   for (k = 0; k < c->n_cigar; ++k) {
      int op = cigar[k] & BAM_CIGAR_MASK;
      if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP) {
         l += cigar[k] >> BAM_CIGAR_SHIFT;
      }
      if (op == BAM_CDEL) {
         l -= cigar[k] >> BAM_CIGAR_SHIFT;
      }
   }
   return l;
}

int32_t bam_cigar2mismatch( const bam1_core_t *readCore, const uint32_t *cigar) {
	uint32_t cigarIndex;
	int32_t numberOfMismatches = 0;
	for (cigarIndex = 0; cigarIndex < readCore->n_cigar; ++cigarIndex) {
		int elementType = cigar[cigarIndex] & BAM_CIGAR_MASK;
		if (elementType != BAM_CMATCH ) {
			numberOfMismatches += cigar[cigarIndex] >> BAM_CIGAR_SHIFT;
		}
	}
	return numberOfMismatches;
}

short get_RP_Reads(ControlState& currentState, const SearchWindow& currentWindow ) 
{
	short ReturnFromReadingReads;
   RPVector TempOneRPVector;

   if (userSettings->bamFilesAsInput()) {
       ReturnFromReadingReads = 0;
       for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
          currentState.Reads_RP.push_back(TempOneRPVector);
          ReturnFromReadingReads = ReadInBamReads_RP(
                                                       currentState.bams_to_parse[i].BamFile.c_str(),
                                                       currentWindow.getChromosome()->getName(), 
																		 &currentWindow.getChromosome()->getSeq(),
                                                       //&currentState.CurrentChrSeq,
                                                       currentState.Reads_RP[i],
                                                       currentState.bams_to_parse[i].InsertSize,
                                                       currentState.bams_to_parse[i].Tag,
                                                       currentWindow);
          if (ReturnFromReadingReads == 0) {
              LOG_ERROR(*logStream << "Bam read failed: " << currentState.bams_to_parse[i].BamFile << std::endl);
              return 1;
          }
          else if (currentState.Reads_RP.size() == 0) {
          }
          else std::cout << currentState.bams_to_parse[i].BamFile << " RP " << currentState.Reads_RP[i].size() << std::endl;
       }
    }
    return 0;
}

short get_RP_Reads_Discovery(ControlState& currentState, const SearchWindow& currentWindow )
{
	short ReturnFromReadingReads;
	RPVector TempOneRPVector;
    
	if (userSettings->bamFilesAsInput()) {
		ReturnFromReadingReads = 0;
		for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
			//currentState.Reads_RP.push_back(TempOneRPVector);
			//std::cout << "Current BAM file: " << currentState.bams_to_parse[i].BamFile.c_str() << std::endl;
			ReturnFromReadingReads = ReadInBamReads_RP_Discovery(
                                                       currentState.bams_to_parse[i].BamFile.c_str(),
                                                       currentWindow.getChromosome()->getName(),
                                                       &currentWindow.getChromosome()->getSeq(),
                                                       //&currentState.CurrentChrSeq,
                                                       currentState.Reads_RP_Discovery,
                                                       currentState.Reads_RP_Discovery_InterChr,
                                                       currentState.bams_to_parse[i].InsertSize,
                                                       currentState.bams_to_parse[i].Tag,
                                                       currentWindow);
			if (ReturnFromReadingReads == 0) {
				*logStream << "Bamfile " << currentState.bams_to_parse[i].BamFile << " has some kind of error. Aborting Pindel\n";
				return 1;
			}
			if (currentState.Reads_RP_Discovery.size() == 0) {
				*logStream << "No discordant RP reads in Bamfile " << currentState.bams_to_parse[i].BamFile.c_str() << std::endl;
			}
			else {
				std::cout << currentState.bams_to_parse[i].BamFile << " RP " << currentState.Reads_RP_Discovery.size() << std::endl;
			}
		}
	}
	return 0;
}

void readInPindelReads(PindelReadReader &reader, const std::string& pindelFilename, ControlState& currentState, const SearchWindow& currentWindow )
{
	int ReturnFromReadingReads = ReadInRead(  reader, currentWindow.getChromosome()->getName(),
                                             currentWindow.getChromosome()->getSeq(), 
//															currentState.CurrentChrSeq,
															currentState.Reads_SR,
                                             currentWindow);
	if (ReturnFromReadingReads == 1) {
		LOG_ERROR(*logStream << "malformed record detected in " << pindelFilename << std::endl);
		exit( EXIT_FAILURE );
	}
	else if (currentState.Reads_SR.size() == 0) {
		LOG_ERROR(*logStream << "No reads found in " << pindelFilename << std::endl);
	}
}


short get_SR_Reads(ControlState& currentState, const SearchWindow& currentWindow )
{
   // std::cout << "in get_SR_Reads " << std::endl;
    
	g_NumReadInWindow = 0; // #################
	g_NumReadScanned = 0;
	g_InWinPlus = 0; // #################
	g_InWinMinus = 0; // #################
	g_CloseMappedPlus = 0; // #################
	g_CloseMappedMinus = 0; // #################
  // std::cout << "getReads " << currentWindow.getChromosome()->getName() << " " << currentWindow.getChromosome()->getSeq().size() << std::endl;
   short ReturnFromReadingReads;
   ReadBuffer readBuffer(BUFFER_SIZE, currentState.Reads_SR, currentState.OneEndMappedReads, currentWindow.getChromosome()->getSeq(), *currentWindow.getChromosome());
	//UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
   if (userSettings->bamFilesAsInput()) {
      ReturnFromReadingReads = 0;
      for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
         *logStream << "\nInsertsize in config: " << currentState.bams_to_parse[i].InsertSize << std::endl;
         //std::cout << "before ReadInBamReads_SR " << std::endl;
         ReturnFromReadingReads = ReadInBamReads_SR(
                                                    currentState.bams_to_parse[i].BamFile.c_str(),
                                                    currentWindow.getChromosome()->getName(), 
                                                    &currentWindow.getChromosome()->getSeq(),
																	//	&currentState.CurrentChrSeq, 
                                                    currentState.Reads_SR,
                                                    currentState.OneEndMappedReads,
                                                    currentState.RefSupportingReads,
                                                    currentState.bams_to_parse[i].InsertSize,
                                                    currentState.bams_to_parse[i].Tag,
                                                    currentWindow, readBuffer );
         //std::cout << "after ReadInBamReads_SR " << ReturnFromReadingReads << std::endl;
         if (ReturnFromReadingReads == 0) { // perhaps 'false'? ReadInBamReads returns a boolean...
            LOG_ERROR(*logStream << "Bam read failed: " << currentState.bams_to_parse[i].BamFile << std::endl);
            return EXIT_FAILURE;
         }
         else if (currentState.Reads_SR.size() == 0) {
            LOG_ERROR(*logStream << "No currentState.Reads for " << currentWindow.getChromosome()->getName() << " found in " << currentState.bams_to_parse[i].BamFile << std::endl);
         }
         *logStream << "BAM file index\t" << i << "\nBam file name\t" << currentState.bams_to_parse[i].BamFile.c_str() << "\nNumber of split-reads so far\t" << currentState.Reads_SR.size() << "\n" << std::endl;
      }  
   }

   if (userSettings->pindelConfigFileAsInput()) {	
		for (unsigned int fileIndex=0; fileIndex<currentState.pindelfilesToParse.size(); fileIndex++ ) {
			const std::string &filename = currentState.pindelfilesToParse[fileIndex];
			currentState.lineReader = getLineReaderByFilename(filename.c_str());
			currentState.inf_Pindel_Reads = new PindelReadReader(*currentState.lineReader);
			readInPindelReads(*currentState.inf_Pindel_Reads, filename, currentState, currentWindow );

			delete currentState.inf_Pindel_Reads;
			delete currentState.lineReader;
		}
	}

	if (userSettings->singlePindelFileAsInput()) {
		readInPindelReads(*currentState.inf_Pindel_Reads, userSettings->pindelFilename, currentState, currentWindow );
	}
   //std::cout << "existing get_SR_Reads " << std::endl;
   return 0;
}

/*

short get_RP_Reads(ControlState& currentState, const SearchWindow& currentWindow )
{
	g_NumReadInWindow = 0; // #################
	g_NumReadScanned = 0;
	g_InWinPlus = 0; // #################
	g_InWinMinus = 0; // #################
	g_CloseMappedPlus = 0; // #################
	g_CloseMappedMinus = 0; // #################
    std::cout << "getReads " << currentWindow.getChromosome()->getName() << " " << currentWindow.getChromosome()->getSeq().size() << std::endl;
    short ReturnFromReadingReads;
    //ReadBuffer readBuffer(BUFFER_SIZE, currentState.Reads_RP, currentWindow.getChromosome()->getSeq());
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
    if (userSettings->bamFilesAsInput()) {
        ReturnFromReadingReads = 0;
        for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
            *logStream << "Insertsize in bamreads: " << currentState.bams_to_parse[i].InsertSize << std::endl;
            ReturnFromReadingReads = ReadInBamReads_RP(
                                                       currentState.bams_to_parse[i].BamFile.c_str(),
                                                       currentWindow.getChromosome()->getName(),
                                                       &currentWindow.getChromosome()->getSeq(),
                                                       //	&currentState.CurrentChrSeq,
                                                       currentState.Reads_RP,
                                                       currentState.bams_to_parse[i].InsertSize,
                                                       currentState.bams_to_parse[i].Tag,
                                                       currentWindow);
            if (ReturnFromReadingReads == 0) { // perhaps 'false'? ReadInBamReads returns a boolean...
                LOG_ERROR(*logStream << "Bam read failed: " << currentState.bams_to_parse[i].BamFile << std::endl);
                return EXIT_FAILURE;
            }
            else if (currentState.Reads_RP.size() == 0) {
                LOG_ERROR(*logStream << "No currentState.Reads for " << currentWindow.getChromosome()->getName() << " found in " << currentState.bams_to_parse[i].BamFile << std::endl);
            }
            *logStream << "BAM file index\t" << i << "\t" << currentState.Reads_RP.size() << std::endl;
        }
    }    
    return 0;
}
*/


/* 'updateReadAfterCloseEndMapping' (EWL, 31thAug2011) */
void updateReadAfterCloseEndMapping( SPLIT_READ& Temp_One_Read )
{
	// do we define g_reportLength if we read in BAM files? =>EW: well, this is only called by ReadBuffer::Flush
    if (g_reportLength < Temp_One_Read.getReadLength()) {
        g_reportLength = Temp_One_Read.getReadLength();
    }
    Temp_One_Read.Used = false;
    Temp_One_Read.UniqueRead = true;
//    LOG_DEBUG(*logStream << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t" << std::endl);

    CleanUniquePoints (Temp_One_Read.UP_Close);

//    LOG_DEBUG(*logStream << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << std::endl);

    Temp_One_Read.CloseEndLength = Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size () - 1].LengthStr; // possibly remove CloseEndLength

    if (Temp_One_Read.MatchedD == Plus) {
        Temp_One_Read.LeftMostPos = Temp_One_Read.UP_Close[0].AbsLoc + 1 - Temp_One_Read.UP_Close[0].LengthStr;
        g_CloseMappedPlus++;
    }
    else {
        Temp_One_Read.LeftMostPos = Temp_One_Read.UP_Close[0].AbsLoc +	Temp_One_Read.UP_Close[0].LengthStr - Temp_One_Read.getReadLength();
        g_CloseMappedMinus++;
    }
}
