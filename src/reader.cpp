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
struct fetch_func_data_SR {
   fetch_func_data_SR () {
      LeftReads = NULL;
      read_to_map_qual = NULL;
      header = NULL;
      b1_flags = NULL;
      b2_flags = NULL;
      CurrentChrSeq = NULL;
      Tag = "";
      InsertSize = 0;
   }
   ReadBuffer *readBuffer;
   std::vector < SPLIT_READ > *LeftReads;
   khash_t (read_name) * read_to_map_qual;
   bam_header_t *header;
   flags_hit *b1_flags;
   flags_hit *b2_flags;
   std::string * CurrentChrSeq;
   std::string Tag;
   int InsertSize;
};

struct fetch_func_data_RP {
    fetch_func_data_RP () {
        LeftReads = NULL;
        read_to_map_qual = NULL;
        header = NULL;
        b1_flags = NULL;
        b2_flags = NULL;
        CurrentChrSeq = NULL;
        Tag = "";
        InsertSize = 0;
    }
    //ReadBuffer *readBuffer;
    std::vector < RP_READ > *LeftReads;
    khash_t (read_name) * read_to_map_qual;
    bam_header_t *header;
    flags_hit *b1_flags;
    flags_hit *b2_flags;
    std::string * CurrentChrSeq;
    std::string Tag;
    int InsertSize;
};

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

void showReadStats(const std::vector<SPLIT_READ>& Reads)
{
   LOG_INFO(*logStream << "Number of reads in current window:                  \t" << g_NumReadInWindow <<
            ", + " << g_InWinPlus << " - " << g_InWinMinus << std::endl);
   LOG_INFO(*logStream << "Number of reads where the close end could be mapped:\t" << Reads.size () <<
            ", + " << g_CloseMappedPlus << " - " << g_CloseMappedMinus << std::endl);
   LOG_INFO(*logStream << "Percentage of reads which could be mapped: + " << std::setprecision(2) << std::fixed << safeDivide( (int)(g_CloseMappedPlus * 100.0) , g_InWinPlus ) <<
            "% - " << safeDivide( (int)(g_CloseMappedMinus * 100.0) , g_InWinMinus ) << "%\n");
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
	// LeftReads does not seem to be defined here... remove line?
   LOG_DEBUG(*logStream << LeftReads.size() << std::endl);
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
         if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size) { // perhaps make setter for MatchedRelpos, so these details can be factored out?
            Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
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
                  GetCloseEnd (CurrentChr, BufferReads[BufferReadsIndex]);
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

                  LOG_DEBUG(*logStream << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t");
                  CleanUniquePoints (currentRead.UP_Close);
                  LOG_DEBUG(*logStream << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << std::endl);
						
						// CloseEndLength duplicates MaxLenCloseEnd() (and is only used in reader, so does not seem time-saving) remove?
                  currentRead.CloseEndLength = currentRead.UP_Close[currentRead.UP_Close.size() - 1].LengthStr;
                  if (currentRead.MatchedD == Plus) {
							currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size() - 1].AbsLoc + 1 - currentRead.CloseEndLength;
						}
                  else {
                     currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size() - 1].AbsLoc + currentRead.CloseEndLength -                      										currentRead.getReadLength();
						}
						// we may want to remove some of the LOG_DEBUG statements; these mainly increase code-reading time
                  LOG_DEBUG(*logStream << "Pushing back!" << std::endl);
                  Reads.push_back (currentRead);
                  if (currentRead.MatchedD == Plus) {
                     g_CloseMappedPlus++;
                  }
                  else {
                     g_CloseMappedMinus++;
                  }
                  g_sampleNames.insert(currentRead.Tag);
               }
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
         GetCloseEnd (CurrentChr, BufferReads[BufferReadsIndex]);
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
         LOG_DEBUG(*logStream << Temp_One_Read.MatchedD  << "\t" << Temp_One_Read.UP_Close.size() << "\t");

         CleanUniquePoints (currentRead.UP_Close);

         LOG_DEBUG(*logStream << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << std::endl);

         currentRead.CloseEndLength = currentRead.UP_Close[currentRead.UP_Close.size () - 1].LengthStr;
         if (currentRead.MatchedD == Plus) {
            currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size () - 1].AbsLoc + 1 - currentRead.CloseEndLength;
			}
         else {
            currentRead.LeftMostPos = currentRead.UP_Close[currentRead.UP_Close.size () - 1].AbsLoc + currentRead.CloseEndLength - currentRead.getReadLength();
			}
         Reads.push_back (currentRead);
         if (currentRead.MatchedD == Plus) {
            g_CloseMappedPlus++;
         }
         else {
            g_CloseMappedMinus++;
         }
         g_sampleNames.insert(currentRead.Tag);
      }
   }
   BufferReads.clear ();

	// what do we want with this? If we want to keep it like this, we would not want to keep it a global boolean, but a static function boolean. Or for none or all chromosomes?
   if (FirstChr) {
      LOG_INFO(*logStream << std::endl << "The last read Pindel scanned: " << std::endl);
		// we may want to factor this out into a function "ReportBasicReadData" or such...
      LOG_INFO(*logStream << Temp_One_Read.Name << std::endl
               << Temp_One_Read.getUnmatchedSeq() << std::endl
               << Temp_One_Read.MatchedD << "\t"
               << Temp_One_Read.FragName << "\t"
               << Temp_One_Read.MatchedRelPos << "\t"
               << Temp_One_Read.MS << "\t"
               << Temp_One_Read.InsertSize << "\t" << Temp_One_Read.Tag << std::endl << std::endl);
      FirstChr = false;
   }
   showReadStats(Reads);

	// 0 means "success?" basically just let caller find out that there are no reads.
   if (Reads.size() == 0) {
      return 0;
   }
   LOG_DEBUG(*logStream << LeftReads.size() << std::endl);
   LOG_INFO(*logStream << " finished!" << std::endl);
   return 0;
}

bool ReadInBamReads_RP (const char *bam_path, const std::string & FragName,
                        std::string * CurrentChrSeq, std::vector <RP_READ> &LeftReads, 
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
    fetch_func_data_RP data;
    data.header = header;
    data.CurrentChrSeq = CurrentChrSeq;
    data.LeftReads = &LeftReads;
    data.read_to_map_qual = NULL;
    data.read_to_map_qual = kh_init (read_name);
    flags_hit b1_flags, b2_flags;
    data.b1_flags = &b1_flags;
    data.b2_flags = &b2_flags;
    data.InsertSize = InsertSize;
    data.Tag = Tag;
    bam_fetch (fp, idx, tid, currentWindow.getStart(), currentWindow.getEnd(), &data, fetch_func_RP);
    
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
    kh_clear (read_name, data.read_to_map_qual);
    kh_destroy (read_name, data.read_to_map_qual);
    
    bam_header_destroy (header);
    bam_index_destroy (idx);
    bam_close (fp);
    return true;
}


bool ReadInBamReads_SR (const char *bam_path, const std::string & FragName,
                std::string * CurrentChrSeq,
                std::vector < SPLIT_READ > &LeftReads,
                int InsertSize,
                std::string Tag,
                const SearchWindow& window,
                ReadBuffer& readBuffer)
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
   fetch_func_data_SR data;
   data.header = header;
   data.CurrentChrSeq = CurrentChrSeq;
   data.LeftReads = &LeftReads;
   data.read_to_map_qual = NULL;
   data.read_to_map_qual = kh_init (read_name);
   flags_hit b1_flags, b2_flags;
   data.b1_flags = &b1_flags;
   data.b2_flags = &b2_flags;
   data.InsertSize = InsertSize;
   data.Tag = Tag;
   data.readBuffer=&readBuffer;
   bam_fetch (fp, idx, tid, window.getStart(), window.getEnd(), &data, fetch_func_SR);
   readBuffer.flush(); 
   showReadStats(LeftReads);

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
   kh_clear (read_name, data.read_to_map_qual);
   kh_destroy (read_name, data.read_to_map_qual);

   bam_header_destroy (header);
   bam_index_destroy (idx);
   bam_close (fp);
   return true;
}

bool isGoodAnchor( const flags_hit *read, const bam1_core_t *bamCore )
{
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
   int maxEdits = int (bamCore->l_qseq * userSettings->MaximumAllowedMismatchRate) + 1;
   unsigned int mappingQuality = bamCore->qual;

   return ( read->mapped &&
            ( mappingQuality >= userSettings->minimalAnchorQuality ) &&
            ( read->unique || read->sw ) &&
            ( ! read->suboptimal ) &&
            ( read->edits <= maxEdits )
          );
}

bool isWeirdRead( const flags_hit *read, const bam1_t * bamOfRead )
{
   if ( ! read->mapped ) {
      return true;
   }

   uint32_t *cigar_pointer = bam1_cigar (bamOfRead);
   int cigarMismatchedBases = bam_cigar2mismatch (&bamOfRead->core, cigar_pointer);

   if ( read->edits + cigarMismatchedBases > 0 ) {
      return true;
   }
   else {
      return false;
   }
}


void build_record_SR (const bam1_t * mapped_read, const bam1_t * unmapped_read, void *data)
{
    SPLIT_READ Temp_One_Read;
    fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
    bam_header_t *header = (bam_header_t *) data_for_bam->header;
    std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
    std::string Tag = (std::string) data_for_bam->Tag;
    int InsertSize = (int) data_for_bam->InsertSize;
    
    const bam1_core_t *mapped_core;
    const bam1_core_t *unmapped_core;
    mapped_core = &mapped_read->core;
    unmapped_core = &unmapped_read->core;
    Temp_One_Read.Name = "@";
    Temp_One_Read.Name.append ((const char *) bam1_qname (unmapped_read));
    if (unmapped_core->flag & BAM_FREAD1) {
        Temp_One_Read.Name.append ("/1");
    }
    else if (unmapped_core->flag & BAM_FREAD2) {
        Temp_One_Read.Name.append ("/2");
    }
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
    if (mapped_core->flag & BAM_FREVERSE) {
        Temp_One_Read.MatchedD = '-';
        uint32_t *cigar_pointer = bam1_cigar (mapped_read);
        //reusing length to be something else now. eat it.
        length = bam_cigar2len (mapped_core, cigar_pointer);
        Temp_One_Read.MatchedRelPos += length;
    }
    else {
        Temp_One_Read.MatchedD = '+';
    }
    Temp_One_Read.MS = mapped_core->qual;
    //FIXME pass these through from the command line with a struct
    Temp_One_Read.InsertSize = InsertSize;
	if (InsertSize <= length ) {
		*logStream << "Error: the insert size is only " << InsertSize << " while the read length is " << Temp_One_Read.getReadLength() << std::endl;
        *logStream << "in paired end sequencing, the insert size is the total size of the fragment to be sequenced, with a read length of 100 bp the entire fragment may for example look like\n\n";
		*logStream << "|----100 bp: first read of the read pair--|-------------------300 bp: unsequenced DNA---------------------|----100 bp: second read of the read pair--|\n";
        *logStream << "<-----------------------------------------------------------insert size=500 ------------------------------------------------------------------------->\n\n";
		*logStream << "In the configuration file (the -i option) please check/correct the insert size (second item on each line). If you continue to have problems, please contact us (Kai Ye, k.ye@lumc.nl)\n";
		exit( EXIT_FAILURE );
	}
    Temp_One_Read.Tag = Tag;
    if (((mapped_core->tid == unmapped_core->tid) &&
         (mapped_core->pos != unmapped_core->pos)) &&
        (abs (mapped_core->isize) < Temp_One_Read.InsertSize)) {
        if (Temp_One_Read.MatchedD == '+') {
            Temp_One_Read.MatchedRelPos -= Temp_One_Read.InsertSize;
        }
        else {
            Temp_One_Read.MatchedRelPos += Temp_One_Read.InsertSize;
        }
    }
    
    Temp_One_Read.FragName = header->target_name[mapped_core->tid];
    //Temp_One_Read.FragName = FragName;
    LOG_DEBUG(cout << Temp_One_Read.Name << std::endl
              << Temp_One_Read.getUnmatchedSeq() << std::endl);
    LOG_DEBUG(cout << Temp_One_Read.MatchedD <<
              "\t" << Temp_One_Read.FragName <<
              "\t" << Temp_One_Read.MatchedRelPos <<
              "\t" << Temp_One_Read.MS <<
              "\t" << Temp_One_Read.InsertSize <<
              "\t" << Temp_One_Read.Tag << std.endl);
    
    
    g_NumReadInWindow++;
    
    if (Temp_One_Read.MatchedD == Plus) {
        g_InWinPlus++;
    }
    else {
        g_InWinMinus++;
    }
    if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size) {
        Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
    }
    if (Temp_One_Read.MatchedRelPos < 1) {
        Temp_One_Read.MatchedRelPos = 0;
    }
    data_for_bam->readBuffer->addRead(Temp_One_Read);
    return;
}

void build_record_RP (const bam1_t * r1, const bam1_t * r2, void *data)
{
    
    const bam1_core_t * r1_core;
    const bam1_core_t * r2_core;
    r1_core = &r1->core;
    r2_core = &r2->core;

    RP_READ Temp_One_Read;
    fetch_func_data_RP *data_for_bam = (fetch_func_data_RP *) data;
    bam_header_t *header = (bam_header_t *) data_for_bam->header;
    std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
    std::string Tag = (std::string) data_for_bam->Tag;

    Temp_One_Read.ReadName = "";
    Temp_One_Read.ReadName.append ((const char *) bam1_qname (r2));
    if (r1_core->flag & BAM_FREVERSE) {
        Temp_One_Read.DA = '+';
    }
    else Temp_One_Read.DA = '-';
    
    if (r2_core->flag & BAM_FREVERSE) {
        Temp_One_Read.DB = '+';
    }
    else Temp_One_Read.DB = '-';
    
    Temp_One_Read.PosA = r1_core->pos;
    Temp_One_Read.PosB = r2_core->pos;

    Temp_One_Read.MQA = r1_core->qual;
    Temp_One_Read.MQB = r2_core->qual;
    Temp_One_Read.ChrNameA = header->target_name[r1_core->tid];
    Temp_One_Read.ChrNameB = header->target_name[r2_core->tid];
    //FIXME pass these through from the command line with a struct
    Temp_One_Read.Tag = Tag;
    
    data_for_bam->LeftReads->push_back(Temp_One_Read);
    return;
}


static int fetch_func_SR (const bam1_t * b1, void *data)
{

   g_NumReadScanned++;

   fetch_func_data_SR *data_for_bam = (fetch_func_data_SR *) data;
   khash_t (read_name) * read_to_map_qual =
      (khash_t (read_name) *) data_for_bam->read_to_map_qual;
   flags_hit *b1_flags = data_for_bam->b1_flags;
   flags_hit *b2_flags = data_for_bam->b2_flags;
   const std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;
   SPLIT_READ Temp_One_Read;
   const bam1_core_t *b1_core;
   bam1_t *b2;
   bam1_core_t *b2_core;
   b1_core = &b1->core;
   std::string read_name = bam1_qname (b1);

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
      std::string c_sequence;
   }

   parse_flags_and_tags (b1, b1_flags);
   parse_flags_and_tags (b2, b2_flags);

   if (isGoodAnchor( b1_flags, b1_core ) && isWeirdRead( b2_flags, b2 ) ) {
      build_record_SR (b1, b2, data);
   }
   if (isGoodAnchor( b2_flags, b2_core ) && isWeirdRead( b1_flags, b1 ) ) {
      build_record_SR (b2, b1, data);
   }
   bam_destroy1 (b2);
   return 0;
}

static int fetch_func_RP (const bam1_t * b1, void *data)
{
    g_NumReadScanned++;
    fetch_func_data_RP *data_for_bam = (fetch_func_data_RP *) data;
    khash_t (read_name) * read_to_map_qual =
    (khash_t (read_name) *) data_for_bam->read_to_map_qual;
    flags_hit *b1_flags = data_for_bam->b1_flags;
    flags_hit *b2_flags = data_for_bam->b2_flags;
    const std::string CurrentChrSeq = *(std::string *) data_for_bam->CurrentChrSeq;

    RP_READ Temp_One_Read;
    const bam1_core_t *b1_core;
    bam1_t *b2;
    bam1_core_t *b2_core;
    b1_core = &b1->core;
    std::string read_name = bam1_qname (b1);

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

    parse_flags_and_tags (b1, b1_flags);
    parse_flags_and_tags (b2, b2_flags);

        build_record_RP (b1, b2, data);

    bam_destroy1 (b2);
    
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

int32_t bam_cigar2mismatch( const bam1_core_t *readCore, const uint32_t *cigar)
{
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

   if (UserDefinedSettings::Instance()->bamFilesAsInput()) {
      ReturnFromReadingReads = 0;
      for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
         currentState.Reads_RP.push_back(TempOneRPVector);
         ReturnFromReadingReads = ReadInBamReads_RP(
                                                       currentState.bams_to_parse[i].BamFile.c_str(),
                                                       currentState.CurrentChrName, 
                                                       &currentState.CurrentChrSeq,
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
    	}  
    }
    return 0;
}


short get_SR_Reads(ControlState& currentState, const SearchWindow& currentWindow )
{
	g_NumReadInWindow = 0; // #################
	g_NumReadScanned = 0;
	g_InWinPlus = 0; // #################
	g_InWinMinus = 0; // #################
	g_CloseMappedPlus = 0; // #################
	g_CloseMappedMinus = 0; // #################
   std::cout << "getReads " << currentState.CurrentChrName << " " << currentState.CurrentChrSeq.size() << std::endl;
   short ReturnFromReadingReads;
   ReadBuffer readBuffer(BUFFER_SIZE, currentState.Reads_SR, currentState.CurrentChrSeq);
	UserDefinedSettings* userSettings = UserDefinedSettings::Instance();
   if (userSettings->bamFilesAsInput()) {
      ReturnFromReadingReads = 0;
      for (unsigned int i = 0; i < currentState.bams_to_parse.size(); i++) {
         *logStream << "Insertsize in bamreads: " << currentState.bams_to_parse[i].InsertSize << std::endl;
         ReturnFromReadingReads = ReadInBamReads_SR(
                                                    currentState.bams_to_parse[i].BamFile.c_str(),
                                                    currentState.CurrentChrName, 
                                                    &currentState.CurrentChrSeq,
                                                    currentState.Reads_SR,
                                                    currentState.bams_to_parse[i].InsertSize,
                                                    currentState.bams_to_parse[i].Tag,
                                                    currentWindow, readBuffer );
         if (ReturnFromReadingReads == 0) { // perhaps 'false'? ReadInBamReads returns a boolean...
            LOG_ERROR(*logStream << "Bam read failed: " << currentState.bams_to_parse[i].BamFile << std::endl);
            return EXIT_FAILURE;
         }
         else if (currentState.Reads_SR.size() == 0) {
            LOG_ERROR(*logStream << "No currentState.Reads for " << currentState.CurrentChrName << " found in " << currentState.bams_to_parse[i].BamFile << std::endl);
         }
         *logStream << "BAM file index\t" << i << "\t" << currentState.Reads_SR.size() << std::endl;
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
   return 0;
}

void readInPindelReads(PindelReadReader &reader, const std::string& pindelFilename, ControlState& currentState, const SearchWindow& currentWindow )
{
	int ReturnFromReadingReads = ReadInRead(  reader, currentState.CurrentChrName,
                                             currentState.CurrentChrSeq, currentState.Reads_SR,
                                             currentWindow);
	if (ReturnFromReadingReads == 1) {
		LOG_ERROR(*logStream << "malformed record detected in " << pindelFilename << std::endl);
		exit( EXIT_FAILURE );
	}
	else if (currentState.Reads_SR.size() == 0) {
		LOG_ERROR(*logStream << "No reads found in " << pindelFilename << std::endl);
	}
}

/* 'updateReadAfterCloseEndMapping' (EWL, 31thAug2011) */
void updateReadAfterCloseEndMapping( SPLIT_READ& Temp_One_Read )
{
	// do we define g_reportLength if we read in BAM files? =>EW: well, this is only called by ReadBuffer::Flush
    if (g_reportLength < Temp_One_Read.getReadLength()) {
        g_reportLength = Temp_One_Read.getReadLength();
    }
    Temp_One_Read.Used = false;
    Temp_One_Read.UniqueRead = true;
    LOG_DEBUG(cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t");

    CleanUniquePoints (Temp_One_Read.UP_Close);

    LOG_DEBUG(cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl);

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
