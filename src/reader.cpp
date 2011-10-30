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
#include "pindel.h"
#include "reader.h"
#include "logdef.h"
#include "refreader.h"

// Static function declaration

static int fetch_func (const bam1_t * b1, void *data);
int32_t bam_cigar2mismatch( const bam1_core_t *readCore, const uint32_t *cigar);

// Reader (BamReader and PindelReader)

//init hash/maps for read pairing on the fly
KSORT_INIT_GENERIC (uint32_t) KHASH_MAP_INIT_STR (read_name, bam1_t *)
		 struct fetch_func_data
		 {
			 fetch_func_data ()
			 {
				 LeftReads = NULL;
				 read_to_map_qual = NULL;
				 header = NULL;
				 b1_flags = NULL;
				 b2_flags = NULL;
				 CurrentChr = NULL;
				 Tag = "";
				 InsertSize = 0;
			 }
			 ReadBuffer *readBuffer;
			 std::vector < SPLIT_READ > *LeftReads;
			 khash_t (read_name) * read_to_map_qual;
			 bam_header_t *header;
			 flagshit *b1_flags;
			 flagshit *b2_flags;
			 std::string * CurrentChr;
			 std::string Tag;
			 int InsertSize;
		 };

void
GetOneChrSeq (std::ifstream & fastaFile, std::string & chromosomeSequence, bool WhetherBuildUp)
{
	RefReader* rr = new RefReader(fastaFile, chromosomeSequence );
	rr->ReadSeq(WhetherBuildUp);
	delete rr;
}


double safeDivide( int dividend, int divisor )
{
	if (divisor == 0) return 0;
	else return ((double)dividend / divisor );
}

void showReadStats(const std::vector<SPLIT_READ>& Reads)
{
	//LOG_INFO(std::cout << "NumReadScanned:\t" << g_NumReadScanned << std::endl);
	LOG_INFO(std::cout << "Number of reads in current window:                  \t" << g_NumReadInWindow << 
								 ", + " << g_InWinPlus << " - " << g_InWinMinus << std::endl);
	LOG_INFO(std::cout << "Number of reads where the close end could be mapped:\t" << Reads.size () << 
								 ", + " << g_CloseMappedPlus << " - " << g_CloseMappedMinus << std::endl);
	LOG_INFO(std::cout << "Percentage of reads which could be mapped: + " << std::setprecision(2) << std::fixed << safeDivide( g_CloseMappedPlus * 100.0 , g_InWinPlus ) <<
								 "% - " << safeDivide( g_CloseMappedMinus * 100.0 , g_InWinMinus ) << "%\n");
	std::cout << std::endl;
	/*LOG_INFO(std::cout << "NumReadStored / NumReadInChr = " << 
	   	       safeDivide( Reads.size () * 100.0 , g_NumReadInChr ) << 
                " %" << std::endl << "InChrPlus \t" << g_InChrPlus << "\tGetPlus \t" <<
					 g_GetPlus << "\t" << 
                safeDivide( GetPlus * 100.0 , g_InChrPlus ) << " %" << std::endl << "InChrMinus\t" << g_InChrMinus << "\tGetMinus\t" <<
					 g_GetMinus << "\t" << safeDivide (g_GetMinus * 100.0 , g_InChrMinus )<< " %" << std::endl << std::endl);*/
}

short
ReadInRead (std::ifstream & inf_ReadSeq, const std::string & FragName,
						const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads,
						const unsigned int lowerBinBorder, const unsigned int upperBinBorder)
{
	LOG_INFO(std::cout << "Scanning and processing reads anchored in " << FragName << std::endl);
	SPLIT_READ Temp_One_Read;
	std::vector < SPLIT_READ > BufferReads;
	g_reportLength = 0;
	LOG_DEBUG(std::cout << LeftReads.size() << std::endl);
	std::string TempQC, TempLine, TempStr, TempFragName;

	inf_ReadSeq.clear ();
	inf_ReadSeq.seekg (0);
	LOG_DEBUG(std::cout << "MINUSEXTRA!" << std::endl);
	int UPCLOSE_COUNTER = 0;
	//loop over all reads in the file
	while (inf_ReadSeq >> Temp_One_Read.Name)
		{
			if (Temp_One_Read.Name[0] != FirstCharReadName)
				{												// !='@'
					LOG_WARN(std::cout << "Something wrong with the read name: " << Temp_One_Read.
									 Name << std::endl);
					Reads.clear ();
					return 1;
				}
			g_NumReadScanned++;
			// get (useless) rest of first line
			std::getline (inf_ReadSeq, TempLine);
			inf_ReadSeq >> Temp_One_Read.UnmatchedSeq;
			std::getline (inf_ReadSeq, TempLine);

			inf_ReadSeq >> Temp_One_Read.MatchedD;
			if (Temp_One_Read.MatchedD != Minus && Temp_One_Read.MatchedD != Plus)
				{
					LOG_INFO(std::cout << Temp_One_Read.Name << std::endl
									 << Temp_One_Read.UnmatchedSeq << std::endl
									 << Temp_One_Read.MatchedD << " ..." << std::endl);
					LOG_INFO(std::cout << "+/-" << std::endl);
					return 1;
				}
			//   >> TempInt 
			inf_ReadSeq >> Temp_One_Read.FragName
				>> Temp_One_Read.MatchedRelPos
				>> Temp_One_Read.MS >> Temp_One_Read.InsertSize >> Temp_One_Read.Tag;
			std::getline (inf_ReadSeq, TempLine);


			if ((signed int)Temp_One_Read.MatchedRelPos > g_maxPos)
				{
					g_maxPos = Temp_One_Read.MatchedRelPos;
				}

			if (Temp_One_Read.FragName == FragName
					&& Temp_One_Read.MatchedRelPos >= lowerBinBorder
					&& Temp_One_Read.MatchedRelPos < upperBinBorder)
				{
					Temp_One_Read.ReadLength = Temp_One_Read.UnmatchedSeq.size ();
					Temp_One_Read.ReadLengthMinus = Temp_One_Read.ReadLength - 1;
					g_NumReadInWindow++;

					Temp_One_Read.MAX_SNP_ERROR =
						(short) (Temp_One_Read.UnmatchedSeq.size () * Seq_Error_Rate);
					Temp_One_Read.TOTAL_SNP_ERROR_CHECKED =
						Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
					Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus =
						Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
					Temp_One_Read.MinClose = 8;
					Temp_One_Read.Found = false;
					if (Temp_One_Read.MatchedD == Plus) {
						g_InWinPlus++;
					}
					else
						g_InWinMinus++;				// this seems to be going correctly
					if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size)
						Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
					if (Temp_One_Read.MatchedRelPos < 1)
						Temp_One_Read.MatchedRelPos = 0;
					BufferReads.push_back (Temp_One_Read);
					if (BufferReads.size () >= NumberOfReadsPerBuffer)
						{
							LOG_DEBUG(std::cout << "here" << std::endl);
#pragma omp parallel default(shared)
							{
#pragma omp for
								for (int BufferReadsIndex = 0; // openMP 2.5 requires signed
									 BufferReadsIndex < (int)NumberOfReadsPerBuffer;
										 BufferReadsIndex++)
									{
										GetCloseEnd (CurrentChr, BufferReads[BufferReadsIndex]);
									}
							}
							for (unsigned int BufferReadsIndex = 0;
								 BufferReadsIndex < NumberOfReadsPerBuffer;
									 BufferReadsIndex++)
								{
									if (BufferReads[BufferReadsIndex].UP_Close.size ())
										{
											UPCLOSE_COUNTER++;
											if (g_reportLength <
													BufferReads[BufferReadsIndex].ReadLength)
												g_reportLength =
													BufferReads[BufferReadsIndex].ReadLength;
											BufferReads[BufferReadsIndex].Used = false;
											BufferReads[BufferReadsIndex].Unique = true;
											LOG_DEBUG(std::cout << Temp_One_Read.MatchedD << "\t" 
																<< Temp_One_Read.UP_Close.size() << "\t");
											CleanUniquePoints (BufferReads[BufferReadsIndex].
																				 UP_Close);
											LOG_DEBUG(std::cout << Temp_One_Read.UP_Close.size() << "\t" 
																<< Temp_One_Read.UP_Close[0].Direction << std::endl);
											BufferReads[BufferReadsIndex].CloseEndLength =
												BufferReads[BufferReadsIndex].
												UP_Close[BufferReads[BufferReadsIndex].UP_Close.
																 size () - 1].LengthStr;
											if (BufferReads[BufferReadsIndex].MatchedD == Plus)
												BufferReads[BufferReadsIndex].LeftMostPos =
													BufferReads[BufferReadsIndex].
													UP_Close[BufferReads[BufferReadsIndex].UP_Close.
																	 size () - 1].AbsLoc + 1 -
													BufferReads[BufferReadsIndex].CloseEndLength;
											else
												BufferReads[BufferReadsIndex].LeftMostPos =
													BufferReads[BufferReadsIndex].
													UP_Close[BufferReads[BufferReadsIndex].UP_Close.
																	 size () - 1].AbsLoc +
													BufferReads[BufferReadsIndex].CloseEndLength -
													BufferReads[BufferReadsIndex].ReadLength;
											LOG_DEBUG(std::cout << "Pushing back!" << std::endl);
											Reads.push_back (BufferReads[BufferReadsIndex]);
											if (BufferReads[BufferReadsIndex].MatchedD == Plus)
												g_CloseMappedPlus++;
											else
												g_CloseMappedMinus++;
											g_sampleNames.insert(BufferReads[BufferReadsIndex].Tag);
										}
								}
							BufferReads.clear ();
						}										// if buffer-reads threatens to overflow
				}												// if the read is in the correct bin
		}														// while loop over each read
	LOG_INFO(std::cout << "last one: " << BufferReads.
					 size () << " and UPCLOSE= " << UPCLOSE_COUNTER << std::endl);
#pragma omp parallel default(shared)
	{
#pragma omp for
		for (int BufferReadsIndex = 0; BufferReadsIndex < (int)BufferReads.size ();
				 BufferReadsIndex++) // signed type required by OpenMP 2.5
			{
				GetCloseEnd (CurrentChr, BufferReads[BufferReadsIndex]);
			}
	}
	for (unsigned int BufferReadsIndex = 0; BufferReadsIndex < BufferReads.size ();
			 BufferReadsIndex++)
		{
			if (BufferReads[BufferReadsIndex].UP_Close.size ())
				{
					UPCLOSE_COUNTER++;
					if (g_reportLength < BufferReads[BufferReadsIndex].ReadLength)
						g_reportLength = BufferReads[BufferReadsIndex].ReadLength;
					BufferReads[BufferReadsIndex].Used = false;
					BufferReads[BufferReadsIndex].Unique = true;

					LOG_DEBUG(std::cout << Temp_One_Read.MatchedD 
										<< "\t" << Temp_One_Read.UP_Close.size() << "\t");

					CleanUniquePoints (BufferReads[BufferReadsIndex].UP_Close);

					LOG_DEBUG(std::cout << Temp_One_Read.UP_Close.size()
										<< "\t" << Temp_One_Read.UP_Close[0].Direction << std::endl);

					BufferReads[BufferReadsIndex].CloseEndLength =
						BufferReads[BufferReadsIndex].
						UP_Close[BufferReads[BufferReadsIndex].UP_Close.size () -
										 1].LengthStr;
					if (BufferReads[BufferReadsIndex].MatchedD == Plus)
						BufferReads[BufferReadsIndex].LeftMostPos =
							BufferReads[BufferReadsIndex].
							UP_Close[BufferReads[BufferReadsIndex].UP_Close.size () -
											 1].AbsLoc + 1 -
							BufferReads[BufferReadsIndex].CloseEndLength;
					else
						BufferReads[BufferReadsIndex].LeftMostPos =
							BufferReads[BufferReadsIndex].
							UP_Close[BufferReads[BufferReadsIndex].UP_Close.size () -
											 1].AbsLoc +
							BufferReads[BufferReadsIndex].CloseEndLength -
							BufferReads[BufferReadsIndex].ReadLength;
					Reads.push_back (BufferReads[BufferReadsIndex]);
					if (BufferReads[BufferReadsIndex].MatchedD == Plus)
						g_CloseMappedPlus++;
					else {
						g_CloseMappedMinus++;
					}
					g_sampleNames.insert(BufferReads[BufferReadsIndex].Tag);
				}
		}
	BufferReads.clear ();

	if (FirstChr)
		{
      LOG_INFO(std::cout << std::endl << "The last read Pindel scanned: " << std::endl);
      LOG_INFO(std::cout << Temp_One_Read.Name << std::endl
                         << Temp_One_Read.UnmatchedSeq << std::endl
                         << Temp_One_Read.MatchedD << "\t"
                         << Temp_One_Read.FragName << "\t"
                         << Temp_One_Read.MatchedRelPos << "\t"
                         << Temp_One_Read.MS << "\t"
                         << Temp_One_Read.InsertSize << "\t" << Temp_One_Read.Tag << std::endl << std::endl);
		FirstChr = false;
		}
	showReadStats(Reads);	


	if (Reads.size () == 0)
		return 0;
	LOG_DEBUG(std::cout << LeftReads.size() << std::endl);
	LOG_INFO(std::cout << " finished!" << std::endl);
	return 0;
}






bool
ReadInBamReads (const char *bam_path, const std::string & FragName,
                std::string * CurrentChr,
                std::vector < SPLIT_READ > &LeftReads, 
                int InsertSize, 
                std::string Tag, 
                int binStart, 
                int binEnd, 
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


	fetch_func_data data;
	data.header = header;
	data.CurrentChr = CurrentChr;
	data.LeftReads = &LeftReads;
	data.read_to_map_qual = NULL;
	data.read_to_map_qual = kh_init (read_name);
	flagshit b1_flags, b2_flags;
	data.b1_flags = &b1_flags;
	data.b2_flags = &b2_flags;
	data.InsertSize = InsertSize;
	data.Tag = Tag;
	data.readBuffer=&readBuffer;
    //std::cout << "before bam fetch " << LeftReads.size() << std::endl;
	bam_fetch (fp, idx, tid, binStart, binEnd, &data, fetch_func);
    //std::cout << "after bam fetch " << LeftReads.size() << std::endl;
	showReadStats(LeftReads);
    //std::cout << "after showReadStats " << LeftReads.size() << std::endl;
    
	khint_t key;
	if (kh_size (data.read_to_map_qual) > 0)
		{
			for (key = kh_begin (data.read_to_map_qual);
					 key != kh_end (data.read_to_map_qual); ++key)
				{
					if (kh_exist (data.read_to_map_qual, key))
						{
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

bool isGoodAnchor( const flagshit *read, const bam1_core_t *bamCore )
{
    //std::cout << "isGoodAnchor" << std::endl;
	int maxEdits = int (bamCore->l_qseq * MaximumAllowedMismatchRate) + 1;
	unsigned int mappingQuality = bamCore->qual;
    /*
    std::cout << "mapped " << read->mapped << "\t"
              << "mappingQuality " << mappingQuality << "\t"
              << "g_minimalAnchorQuality " << g_minimalAnchorQuality << "\t"
              << "unique " << read->unique << "\t"
              << "sw " << read->sw << "\t"
              << "read->unique || read->sw " << (read->unique || read->sw) << "\t"
              << "suboptimal " << read->suboptimal << "\t"
              << "edits " << read->edits << "\t"
              << "maxEdits " << maxEdits << std::endl;
     */
	return ( read->mapped &&
				( mappingQuality >= g_minimalAnchorQuality ) &&
            //( read->unique || read->sw ) &&
				( ! read->suboptimal ) &&
				( read->edits <= maxEdits ) 
          );
}

bool isWeirdRead( const flagshit *read, const bam1_t * bamOfRead )
{
    //std::cout << "isWeirdRead" << std::endl;
	if ( ! read->mapped ) {
		return true;
	}
    
	uint32_t *cigar_pointer = bam1_cigar (bamOfRead);
	int cigarMismatchedBases = bam_cigar2mismatch (&bamOfRead->core, cigar_pointer);

	if ( read->edits + cigarMismatchedBases > 0 ) {
		return true;
	}
	else return false;
		
}

static int
fetch_func (const bam1_t * b1, void *data)
{

	g_NumReadScanned++;
    //std::cout << "g_NumReadScanned " << g_NumReadScanned << std::endl;
	fetch_func_data *data_for_bam = (fetch_func_data *) data;
	khash_t (read_name) * read_to_map_qual =
		(khash_t (read_name) *) data_for_bam->read_to_map_qual;
	flagshit *b1_flags = data_for_bam->b1_flags;
	flagshit *b2_flags = data_for_bam->b2_flags;
	const std::string CurrentChr = *(std::string *) data_for_bam->CurrentChr;
    //std::cout << "1" << std::endl;
	SPLIT_READ Temp_One_Read;
	const bam1_core_t *b1_core;
	bam1_t *b2;
	bam1_core_t *b2_core;
	b1_core = &b1->core;
	std::string read_name = bam1_qname (b1);
	//    if(!(b1_core->flag & BAM_FPROPER_PAIR)) { 
	//            return 0;
	//NO BUENO!        }
    //std::cout << "2" << std::endl;
	khint_t key = kh_get (read_name, read_to_map_qual, bam1_qname (b1));
    //std::cout << "2a" << std::endl;
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
    //std::cout << "3" << std::endl;
	parse_flags_and_tags (b1, b1_flags);
	parse_flags_and_tags (b2, b2_flags);
	//read_name = bam1_qname(b1); 
    //std::cout << "4" << std::endl;
	if (isGoodAnchor( b1_flags, b1_core ) && isWeirdRead( b2_flags, b2 ) ) {
        //std::cout << "condition 1" << std::endl;
		build_record (b1, b2, data);
	}
	if (isGoodAnchor( b2_flags, b2_core ) && isWeirdRead( b1_flags, b1 ) ) {
        //std::cout << "condition 2" << std::endl;
		build_record (b2, b1, data);
	}
	bam_destroy1 (b2);
    //std::cout << "5" << std::endl;
	return 0;
}

/* 'isInBin' returns whether the read "read" is in the designated bin. */
bool
isInBin (const SPLIT_READ & read)
{
	if ((int)read.MatchedRelPos > g_maxPos)
		{
			g_maxPos = (int)read.MatchedRelPos;
		}
	return (((int)read.MatchedRelPos >= (g_binIndex * WINDOW_SIZE)) &&
					((int)read.MatchedRelPos < ((g_binIndex + 1) * WINDOW_SIZE)));
}

void
build_record (const bam1_t * mapped_read, const bam1_t * unmapped_read,
							void *data)
{
    //std::cout << "build_record g_NumReadScanned " << g_NumReadScanned << std::endl;
	SPLIT_READ Temp_One_Read;
	fetch_func_data *data_for_bam = (fetch_func_data *) data;
	bam_header_t *header = (bam_header_t *) data_for_bam->header;
	std::string CurrentChr = *(std::string *) data_for_bam->CurrentChr;
	std::string Tag = (std::string) data_for_bam->Tag;
	int InsertSize = (int) data_for_bam->InsertSize;

	const bam1_core_t *mapped_core;
	const bam1_core_t *unmapped_core;
	mapped_core = &mapped_read->core;
	unmapped_core = &unmapped_read->core;
	Temp_One_Read.Name = "@";
	Temp_One_Read.Name.append ((const char *) bam1_qname (unmapped_read));
	if (unmapped_core->flag & BAM_FREAD1)
		{
			Temp_One_Read.Name.append ("/1");
		}
	else if (unmapped_core->flag & BAM_FREAD2)
		{
			Temp_One_Read.Name.append ("/2");
		}
	std::string c_sequence;
	int i;
	uint8_t *s = bam1_seq (unmapped_read);
	for (i = 0; i < unmapped_core->l_qseq; ++i)
		c_sequence.append (1, bam_nt16_rev_table[bam1_seqi (s, i)]);
	//rudimentary n filter
	int length = unmapped_core->l_qseq;
	while (c_sequence[0] == 'N')
		{
			c_sequence.erase (0, 1);
			length--;
		}
	if (c_sequence.size () > 0)
		{
			while (c_sequence[length - 1] == 'N')
				{
					c_sequence.erase (length - 1, 1);
					length--;
				}
		}
	int n_count = 0;
	size_t found = c_sequence.find ('N', 0);
	int max_ns = length * .10;
	while (found != std::string::npos)
		{
			n_count++;
			found = c_sequence.find ('N', found + 1);
		}
	if (n_count > max_ns || length < 22)
		{
			return;
		}
	//rudimentary n filter end
	Temp_One_Read.ReadLength = length;
	Temp_One_Read.ReadLengthMinus = length - 1;
	if (unmapped_core->flag & BAM_FREVERSE)
		{
			Temp_One_Read.UnmatchedSeq = ReverseComplement (c_sequence);
		}
	else
		{
			Temp_One_Read.UnmatchedSeq = c_sequence;
		}
	Temp_One_Read.MatchedRelPos = mapped_core->pos;
	if (mapped_core->flag & BAM_FREVERSE)
		{
			Temp_One_Read.MatchedD = '-';
			uint32_t *cigar_pointer = bam1_cigar (mapped_read);
			//reusing length to be something else now. eat it.
			length = bam_cigar2len (mapped_core, cigar_pointer);
			Temp_One_Read.MatchedRelPos += length;
		}
	else
		{
			Temp_One_Read.MatchedD = '+';
		}
	Temp_One_Read.MS = mapped_core->qual;
	//FIXME pass these through from the command line with a struct
	Temp_One_Read.InsertSize = InsertSize;
	Temp_One_Read.Tag = Tag;
	if (((mapped_core->tid == unmapped_core->tid) &&
			 (mapped_core->pos != unmapped_core->pos)) &&
			(abs (mapped_core->isize) < Temp_One_Read.InsertSize))
		{
			if (Temp_One_Read.MatchedD == '+')
				{
					Temp_One_Read.MatchedRelPos -= Temp_One_Read.InsertSize;
				}
			else
				{
					Temp_One_Read.MatchedRelPos += Temp_One_Read.InsertSize;
				}
		}

	std::string FragName = header->target_name[mapped_core->tid];
	Temp_One_Read.FragName = FragName;
	LOG_DEBUG(cout << Temp_One_Read.Name << std::endl 
						<< Temp_One_Read.UnmatchedSeq << std::endl);
	LOG_DEBUG(cout << Temp_One_Read.MatchedD << 
						"\t" << Temp_One_Read.FragName << 
						"\t" << Temp_One_Read.MatchedRelPos << 
						"\t" << Temp_One_Read.MS << 
						"\t" << Temp_One_Read.InsertSize << 
						"\t" << Temp_One_Read.Tag << std.endl);
	g_NumReadInWindow++;
	Temp_One_Read.MAX_SNP_ERROR =
		(short) trunc((double)0.5+Temp_One_Read.UnmatchedSeq.size () * Seq_Error_Rate);

	Temp_One_Read.TOTAL_SNP_ERROR_CHECKED =
		Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
	Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus =
		Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
	Temp_One_Read.MinClose = 8;

	if (Temp_One_Read.MatchedD == Plus) {
		g_InWinPlus++;
	}
	else
		g_InWinMinus++;
	if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size)
		Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
	if (Temp_One_Read.MatchedRelPos < 1)
		Temp_One_Read.MatchedRelPos = 0;
    //std::cout << "before data_for_bam" << std::endl;
	data_for_bam->readBuffer->addRead(Temp_One_Read);
    //std::cout << "after data_for_bam" << std::endl;
	return;
}

void
parse_flags_and_tags (const bam1_t * b, flagshit * flags)
{
	const bam1_core_t *c = &b->core;
	char xt_code = 0;
	int mf_code = 0, nm_code = 0, best_hits = 0;
	flags->unique = 0;
	flags->mapped = !(c->flag & BAM_FUNMAP);

	uint8_t *s = bam_aux_get (b, "XT");
	if (s != 0)
		{
			xt_code = bam_aux2A (s);
			if (xt_code == 'U')
				{
					flags->unique = 1;
				}
			else
				{
					flags->unique = 0;
				}

		}
	s = NULL;
	s = bam_aux_get (b, "MF");
	if (s != 0)
		{
			mf_code = bam_aux2i (s);
			//the below paradigm doesn't exist in maq i think we should assume the read is unique so this hack is here
			if (mf_code != 130)
				{
					flags->unique = 1;
				}
			else
				{
					flags->unique = 0;
				}
			flags->suboptimal = 0;
		}
	s = NULL;
	s = bam_aux_get (b, "X0");
	if (s != 0)
		{
			best_hits = bam_aux2i (s);
		}
	s = NULL;
	s = bam_aux_get (b, "X1");
	if (s != 0)
		{
			int sub_hits = bam_aux2i (s);

			if (best_hits + sub_hits == 1)
				{
					flags->suboptimal = 0;
				}
			else
				{
					flags->suboptimal = 1;
				}
		}
	if (xt_code == 'M' || mf_code == 130)
		{
			flags->sw = 1;
			//short term fix to unset unique if the maq read was s-w mapped. bwa can't set U and M at once.
		}
	else
		{
			flags->sw = 0;
		}
	s = NULL;
	s = bam_aux_get (b, "NM");
	if (s != 0)
		{
			nm_code = bam_aux2i (s);
			flags->edits = nm_code;
		}
	else
		{
			nm_code = 0;
			flags->edits = nm_code;
		}

	return;
}


int32_t
bam_cigar2len (const bam1_core_t * c, const uint32_t * cigar)
{
	uint32_t k;
	int32_t l = 0;
	for (k = 0; k < c->n_cigar; ++k)
		{
			int op = cigar[k] & BAM_CIGAR_MASK;
			if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP)
				l += cigar[k] >> BAM_CIGAR_SHIFT;
			if (op == BAM_CDEL)
				l -= cigar[k] >> BAM_CIGAR_SHIFT;
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
