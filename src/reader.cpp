#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
// Samtools
#include "bam.h"
#include "sam.h"
#include "kstring.h"
#include "kseq.h"
#include "khash.h"
#include "ksort.h"
// Pindel
#include "pindel.h"
#include "reader.h"

// Static function declaration

static int fetch_func (const bam1_t * b1, void *data);

// Reader (BamReader and PindelReader)

//init hash/maps for read pairing on the fly
KSORT_INIT_GENERIC (uint32_t) KHASH_MAP_INIT_STR (read_name, bam1_t *)
		 struct fetch_func_data
		 {
			 fetch_func_data ()
			 {
				 read_to_map_qual = NULL;
				 header = NULL;
				 b1_flags = NULL;
				 b2_flags = NULL;
				 CurrentChr = NULL;
				 Tag = "";
				 InsertSize = 0;
			 }
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
ReadInOneChr (std::ifstream & inf_Seq, std::string & TheInput, const std::string & ChrName)
{
	TheInput.clear ();
	char TempChar;
	std::string TempLine, TempChrName;
	inf_Seq >> TempChar;
	if (TempChar != '>')
		{
			std::cout << "Please use fasta format for the reference file." << std::endl;
			TheInput.clear ();
			return;
		}
	while (inf_Seq >> TempChrName)
		{
			std::cout << "Processing chromosome " << TempChrName << std::endl;
			if (!TheInput.empty ())
				{
					std::cout << "Skip the rest of chromosomes.\n";
					break;
				}
			if (TempChrName == ChrName)
				{
					std::getline (inf_Seq, TempLine);
					while (inf_Seq >> TempChar)
						{
							if (TempChar != '\n' && TempChar != '\r')
								{
									if (TempChar == '>')
										{
											break;
										}
									else
										{
											if ('a' <= TempChar && TempChar <= 'z')
												TempChar = TempChar + ('A' - 'a');
											switch (TempChar)
												{
												case 'A':
													TheInput += 'A';
													break;	// 00000000
												case 'C':
													TheInput += 'C';
													break;	// 00010000
												case 'G':
													TheInput += 'G';
													break;	// 00100000
												case 'T':
													TheInput += 'T';
													break;	// 00110000
												default:
													TheInput += 'N';	// 01000000
												}
										}						// else TempChar
								}
						}
				}
			else
				{
					std::getline (inf_Seq, TempLine);
					while (inf_Seq >> TempChar)
						{
							if (TempChar != '\n' && TempChar != '\r')
								{
									if (TempChar == '>')
										{
											break;
										}
									else
										{

										}						// else TempChar
								}
						}
				}
		}
	std::cout << ChrName << "\t" << TheInput.size () << "\t";
	if (!TheInput.empty ())
		{
			std::string Spacer = "";
			for (unsigned i = 0; i < SpacerBeforeAfter; i++)
				Spacer += "N";
			TheInput = Spacer + TheInput + Spacer;
		}
	std::cout << TheInput.size () << std::endl;

	// EWL280311: reset file for to allow re-reading for the next chromosome
	inf_Seq.clear ();
	inf_Seq.seekg (0);

	return;
}

void
GetOneChrSeq (std::ifstream & inf_Seq, std::string & TheInput, bool WhetherBuildUp)
{
	TheInput.clear ();
	char TempChar;
	std::string TempLine, TempChrName;

	if (WhetherBuildUp)
		{
			//getline(inf_Seq, TempLine);
			while (inf_Seq >> TempChar)
				{
					if (TempChar != '\n' && TempChar != '\r')
						{
							if (TempChar == '>')
								{
									break;
								}
							else
								{
									if ('a' <= TempChar && TempChar <= 'z')
										TempChar = TempChar + ('A' - 'a');
									switch (TempChar)
										{
										case 'A':
											TheInput += 'A';
											break;		// 00000000
										case 'C':
											TheInput += 'C';
											break;		// 00010000
										case 'G':
											TheInput += 'G';
											break;		// 00100000
										case 'T':
											TheInput += 'T';
											break;		// 00110000
										default:
											TheInput += 'N';	// 01000000
										}
								}								// else TempChar
						}
				}
		}
	else
		{
			//getline(inf_Seq, TempLine);
			while (inf_Seq >> TempChar)
				{
					if (TempChar != '\n' && TempChar != '\r')
						{
							if (TempChar == '>')
								{
									break;
								}
							else
								{

								}								// else TempChar
						}
				}
		}

	if (!TheInput.empty ())
		{
			std::string Spacer = "";
			for (unsigned i = 0; i < SpacerBeforeAfter; i++)
				Spacer += "N";
			TheInput = Spacer + TheInput + Spacer;
		}
	return;
}


short
ReadInRead (std::ifstream & inf_ReadSeq, const std::string & FragName,
						const std::string & CurrentChr, std::vector < SPLIT_READ > &Reads,
						const unsigned int lowerBinBorder, const unsigned int upperBinBorder)
{
	std::cout << "Scanning and processing reads anchored in " << FragName << std::endl;
	//short ADDITIONAL_MISMATCH = 1;
	SPLIT_READ Temp_One_Read;
	unsigned int NumReadScanned = 0;
	unsigned int NumReadInChr = 0;
	unsigned int InChrPlus = 0;
	unsigned int InChrMinus = 0;
	unsigned int GetPlus = 0;
	unsigned int GetMinus = 0;
	//NumberOfReadsPerBuffer;
	std::vector < SPLIT_READ > BufferReads;
	ReportLength = 0;
	//cout << LeftReads.size() << endl;
	std::string TempQC, TempLine, TempStr, TempFragName;
	//int TempInt;

	inf_ReadSeq.clear ();
	inf_ReadSeq.seekg (0);
	VectorTag.clear ();
//cout << "MINUSEXTRA!\n";
	int UPCLOSE_COUNTER = 0;
	//loop over all reads in the file
	while (inf_ReadSeq >> Temp_One_Read.Name)
		{
			if (Temp_One_Read.Name[0] != FirstCharReadName)
				{												// !='@'
					std::cout << "Something wrong with the read name: " << Temp_One_Read.
						Name << std::endl;
					Reads.clear ();
					return 1;
				}
			NumReadScanned++;
			// get (useless) rest of first line
			std::getline (inf_ReadSeq, TempLine);
			inf_ReadSeq >> Temp_One_Read.UnmatchedSeq;
			std::getline (inf_ReadSeq, TempLine);

			inf_ReadSeq >> Temp_One_Read.MatchedD;
			if (Temp_One_Read.MatchedD != Minus && Temp_One_Read.MatchedD != Plus)
				{
					std::cout << Temp_One_Read.Name << "\n"
						<< Temp_One_Read.UnmatchedSeq << "\n"
						<< Temp_One_Read.MatchedD << " ...\n";
					std::cout << "+/-" << std::endl;
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
					NumReadInChr++;

					Temp_One_Read.MAX_SNP_ERROR =
						(short) (Temp_One_Read.UnmatchedSeq.size () * Seq_Error_Rate);


					Temp_One_Read.TOTAL_SNP_ERROR_CHECKED =
						Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
					Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus =
						Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
					Temp_One_Read.MinClose =
						short (log ((double) (Temp_One_Read.InsertSize * 3)) / log (4.0) +
									 0.8) + 3;
					//Temp_One_Read.IndelSize = 0;
					Temp_One_Read.Found = false;
					if (Temp_One_Read.MatchedD == Plus)
						{
							InChrPlus++;
							//Temp_One_Read.UnmatchedSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
						}
					else
						InChrMinus++;				// this seems to be going correctly
					if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size)
						Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
					if (Temp_One_Read.MatchedRelPos < 1)
						Temp_One_Read.MatchedRelPos = 0;
					BufferReads.push_back (Temp_One_Read);
					if (BufferReads.size () >= NumberOfReadsPerBuffer)
						{
							//cout << "here" << endl;
#pragma omp parallel default(shared)
							{
#pragma omp for
								for (unsigned int BufferReadsIndex = 0;
										 BufferReadsIndex < NumberOfReadsPerBuffer;
										 BufferReadsIndex++)
									{
										GetCloseEnd (CurrentChr, BufferReads[BufferReadsIndex]);
									}
							}
							for (int BufferReadsIndex = 0;
									 BufferReadsIndex < NumberOfReadsPerBuffer;
									 BufferReadsIndex++)
								{
									if (BufferReads[BufferReadsIndex].UP_Close.size ())
										{
											UPCLOSE_COUNTER++;
											if (ReportLength <
													BufferReads[BufferReadsIndex].ReadLength)
												ReportLength =
													BufferReads[BufferReadsIndex].ReadLength;
											BufferReads[BufferReadsIndex].Used = false;
											BufferReads[BufferReadsIndex].Unique = true;
											//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t";
											CleanUniquePoints (BufferReads[BufferReadsIndex].
																				 UP_Close);
											//cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl;
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
											//cout << "Pushing back!\n";
											Reads.push_back (BufferReads[BufferReadsIndex]);
											if (BufferReads[BufferReadsIndex].MatchedD == Plus)
												GetPlus++;
											else
												GetMinus++;
											if (NotInVector
													(BufferReads[BufferReadsIndex].Tag, VectorTag))
												{
													VectorTag.push_back (BufferReads[BufferReadsIndex].
																							 Tag);
												}
										}
								}
							BufferReads.clear ();
						}										// if buffer-reads threatens to overflow
				}												// if the read is in the correct bin
		}														// while loop over each read
	std::cout << "last one: " << BufferReads.
		size () << " and UPCLOSE= " << UPCLOSE_COUNTER << std::endl;
#pragma omp parallel default(shared)
	{
#pragma omp for
		for (unsigned int BufferReadsIndex = 0; BufferReadsIndex < BufferReads.size ();
				 BufferReadsIndex++)
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
					if (ReportLength < BufferReads[BufferReadsIndex].ReadLength)
						ReportLength = BufferReads[BufferReadsIndex].ReadLength;
					BufferReads[BufferReadsIndex].Used = false;
					BufferReads[BufferReadsIndex].Unique = true;
					//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t";
					CleanUniquePoints (BufferReads[BufferReadsIndex].UP_Close);
					//cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl;
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
						GetPlus++;
					else
						GetMinus++;
					if (NotInVector (BufferReads[BufferReadsIndex].Tag, VectorTag))
						{
							VectorTag.push_back (BufferReads[BufferReadsIndex].Tag);
						}
				}
		}
	BufferReads.clear ();

	if (FirstChr)
		{
			std::cout << "\nThe last read Pindel scanned: \n";
			std::cout << Temp_One_Read.Name << "\n"
				<< Temp_One_Read.UnmatchedSeq << "\n"
				<< Temp_One_Read.MatchedD << "\t"
				<< Temp_One_Read.FragName << "\t"
				<< Temp_One_Read.MatchedRelPos << "\t"
				<< Temp_One_Read.MS << "\t"
				<< Temp_One_Read.InsertSize << "\t" << Temp_One_Read.Tag << "\n\n";
			FirstChr = false;
			//ReportLength = Temp_One_Read.UnmatchedSeq.size();
		}
	std::cout << "NumReadScanned:\t" << NumReadScanned << std::endl;
	std::cout << "NumReadInChr:\t" << NumReadInChr << std::endl;
	std::cout << "NumReadStored:\t" << Reads.size () << std::endl;
	std::cout << "NumReadStored / NumReadInChr = " << Reads.size () * 100.0 /
		NumReadInChr << " %\n" << "InChrPlus \t" << InChrPlus << "\tGetPlus \t" <<
		GetPlus << "\t" << GetPlus * 100.0 /
		InChrPlus << " %\n" << "InChrMinus\t" << InChrMinus << "\tGetMinus\t" <<
		GetMinus << "\t" << GetMinus * 100.0 / InChrMinus << " %\n" << std::endl;
	//inf_ReadSeq.close();
	if (Reads.size () == 0)
		return 0;
	//cout << LeftReads.size() << endl;
	std::cout << "sorting tags ... ";
	std::string Str4Exchange;
	for (unsigned short i = 0; i < VectorTag.size () - 1; i++)
		{
			for (unsigned short j = 1; j < VectorTag.size (); j++)
				{
					if (VectorTag[i].size () > VectorTag[j].size ())
						{
							Str4Exchange = VectorTag[i];
							VectorTag[i] = VectorTag[j];
							VectorTag[j] = Str4Exchange;
						}
					else if (VectorTag[i].size () == VectorTag[j].size ())
						{
							for (unsigned short k = 0; k < VectorTag[i].size (); k++)
								{
									if ((short) VectorTag[i][k] > (short) VectorTag[j][k])
										{
											Str4Exchange = VectorTag[i];
											VectorTag[i] = VectorTag[j];
											VectorTag[j] = Str4Exchange;
											break;
										}
								}
						}
				}
		}
	std::cout << " finished!" << std::endl;
	return 0;
}






bool
ReadInBamReads (const char *bam_path, const std::string & FragName,
								std::string * CurrentChr, std::vector < SPLIT_READ > &LeftReads,
								int InsertSize, std::string Tag, int binStart, int binEnd)
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
	//VectorTag.clear();


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
	//ADDITIONAL_MISMATCH = 1;
	//Seq_Error_Rate = 0.05;
	bam_fetch (fp, idx, tid, binStart, binEnd, &data, fetch_func);
	std::cout << "Bam:\t" << bam_path << "\tTag:\t" << Tag << std::endl;
	std::cout << "NumReadScanned:\t" << NumReadScanned << std::endl;
	std::cout << "NumReadInChr:\t" << NumReadInChr << std::endl;
	std::cout << "NumReadStored:\t" << LeftReads.size () << std::endl;
	std::cout << "NumReadStored / NumReadInChr = " << LeftReads.size () * 100.0 /
		NumReadInChr << " %\n" << "InChrPlus \t" << InChrPlus << "\tGetPlus \t" <<
		GetPlus << "\t" << GetPlus * 100.0 /
		InChrPlus << " %\n" << "InChrMinus\t" << InChrMinus << "\tGetMinus\t" <<
		GetMinus << "\t" << GetMinus * 100.0 / InChrMinus << " %\n" << std::endl;
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

	NumReadScanned = 0;
	NumReadInChr = 0;
	InChrPlus = 0;
	InChrMinus = 0;
	GetPlus = 0;
	GetMinus = 0;


	//kai does the below in read inreads dunno why yet
	std::string Str4Exchange;
	if (VectorTag.size () > 0)
		{
			for (unsigned short i = 0; i < VectorTag.size () - 1; i++)
				{
					for (unsigned short j = 1; j < VectorTag.size (); j++)
						{
							if (VectorTag[i].size () > VectorTag[j].size ())
								{
									Str4Exchange = VectorTag[i];
									VectorTag[i] = VectorTag[j];
									VectorTag[j] = Str4Exchange;
								}
							else if (VectorTag[i].size () == VectorTag[j].size ())
								{
									for (unsigned short k = 0; k < VectorTag[i].size (); k++)
										{
											if ((short) VectorTag[i][k] > (short) VectorTag[j][k])
												{
													Str4Exchange = VectorTag[i];
													VectorTag[i] = VectorTag[j];
													VectorTag[j] = Str4Exchange;
													break;
												}
										}
								}
						}
				}
		}
	//end kai
	bam_header_destroy (header);
	bam_index_destroy (idx);
	bam_close (fp);
	return true;
}

static int
fetch_func (const bam1_t * b1, void *data)
{

	NumReadScanned++;
	fetch_func_data *data_for_bam = (fetch_func_data *) data;
	khash_t (read_name) * read_to_map_qual =
		(khash_t (read_name) *) data_for_bam->read_to_map_qual;
	flagshit *b1_flags = data_for_bam->b1_flags;
	flagshit *b2_flags = data_for_bam->b2_flags;
	const std::string CurrentChr = *(std::string *) data_for_bam->CurrentChr;

	SPLIT_READ Temp_One_Read;
	const bam1_core_t *b1_core;
	bam1_t *b2;
	bam1_core_t *b2_core;
	b1_core = &b1->core;
	std::string read_name = bam1_qname (b1);
	//    if(!(b1_core->flag & BAM_FPROPER_PAIR)) { 
	//            return 0;
	//NO BUENO!        }
	khint_t key = kh_get (read_name, read_to_map_qual, bam1_qname (b1));
	if (key == kh_end (read_to_map_qual))
		{
			int ret;
			key =
				kh_put (read_name, read_to_map_qual, strdup (bam1_qname (b1)), &ret);
			kh_value (read_to_map_qual, key) = bam_dup1 (b1);
			return 0;
		}
	else
		{
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
	//read_name = bam1_qname(b1); 
	if (!(b1_flags->mapped) && !(b2_flags->mapped))
		{
			bam_destroy1 (b2);
			return 0;
		}

	if (!b1_flags->unique && !b1_flags->sw && !b2_flags->unique
			&& !b2_flags->sw)
		{
			bam_destroy1 (b2);
			return 0;
		}

	if (b1_flags->unique && b2_flags->unique && b1_flags->edits == 0
			&& b2_flags->edits == 0)
		{

			bam_destroy1 (b2);

			return 0;
		}
	if (b1_flags->suboptimal && b2_flags->suboptimal)
		{
			bam_destroy1 (b2);
			return 0;
		}
	int max_edits1 = int (b1_core->l_qseq * .05) + 1;
	int max_edits2 = int (b2_core->l_qseq * .05) + 1;

	if (b1_flags->edits > max_edits1 && b2_flags->edits > max_edits2)
		{
			bam_destroy1 (b2);
			return 0;
		}

	int which_read_is_mapped = 0;
	if (b1_flags->mapped && b2_flags->mapped)
		{
			if ((b1_flags->unique || b1_flags->sw) && b1_flags->edits <= max_edits1
					&& b2_flags->edits > 0)
				{
					which_read_is_mapped = 1;
				}
			else if ((b2_flags->unique || b2_flags->sw)
							 && b2_flags->edits <= max_edits1 && b1_flags->edits > 0)
				{
					which_read_is_mapped = 2;
				}
			else
				{
					bam_destroy1 (b2);
					return 0;
				}
		}
	else if (!b1_flags->mapped && b2_flags->edits <= max_edits2)
		{
			which_read_is_mapped = 2;
		}
	else if (!b2_flags->mapped && b1_flags->edits <= max_edits1)
		{
			which_read_is_mapped = 1;
		}
	else
		{
			bam_destroy1 (b2);
			return 0;
		}
	if (which_read_is_mapped == 1)
		{
			build_record (b1, b2, data);
			if ((b2_flags->unique || b2_flags->sw) && b2_flags->edits <= max_edits2
					&& b1_flags->edits > 0)
				{
					build_record (b2, b1, data);
				}
		}
	else if (which_read_is_mapped == 2)
		{
			build_record (b2, b1, data);
		}
	else
		{
			printf ("Should not get here\n");
			exit (1);
		}
	bam_destroy1 (b2);
	return 0;
}

/* 'isInBin' returns whether the read "read" is in the designated bin. */
bool
isInBin (const SPLIT_READ & read)
{
	if (read.MatchedRelPos > g_maxPos)
		{
			g_maxPos = read.MatchedRelPos;
		}
	return ((read.MatchedRelPos >= (g_binIndex * WINDOW_SIZE)) &&
					(read.MatchedRelPos < ((g_binIndex + 1) * WINDOW_SIZE)));
}

void
build_record (const bam1_t * mapped_read, const bam1_t * unmapped_read,
							void *data)
{

	SPLIT_READ Temp_One_Read;
	fetch_func_data *data_for_bam = (fetch_func_data *) data;
	std::vector < SPLIT_READ > *LeftReads =
		(std::vector < SPLIT_READ > *)data_for_bam->LeftReads;
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
	int n_count = 0, found = c_sequence.find ('N', 0);
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
	//cout << Temp_One_Read.Name << "\n" << Temp_One_Read.UnmatchedSeq << "\n";
	//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.FragName << "\t" << Temp_One_Read.MatchedRelPos << "\t" << Temp_One_Read.MS << "\t" << Temp_One_Read.InsertSize << "\t" << Temp_One_Read.Tag << "\n";
	NumReadInChr++;
	Temp_One_Read.MAX_SNP_ERROR =
		(short) (Temp_One_Read.UnmatchedSeq.size () * Seq_Error_Rate);
	Temp_One_Read.TOTAL_SNP_ERROR_CHECKED =
		Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH + 1;
	Temp_One_Read.TOTAL_SNP_ERROR_CHECKED_Minus =
		Temp_One_Read.MAX_SNP_ERROR + ADDITIONAL_MISMATCH;
	Temp_One_Read.MinClose =
		short (log ((double) (Temp_One_Read.InsertSize * 3)) / log (4.0) + 0.8) +
		3;
	//MinClose = short(log((double)(Temp_One_Read.InsertSize * 3))/log(4.0) + 0.8) + 3;// + MAX_SNP_ERROR;//atoi(argv[1]);
	//MinFar_I = MinClose + 1;//atoi(argv[2]);
	if (Temp_One_Read.MatchedD == Plus)
		{
			InChrPlus++;
			//Temp_One_Read.UnmatchedSeq = ReverseComplement(Temp_One_Read.UnmatchedSeq);
		}
	else
		InChrMinus++;
	if (Temp_One_Read.MatchedRelPos > CONS_Chr_Size)
		Temp_One_Read.MatchedRelPos = CONS_Chr_Size;
	if (Temp_One_Read.MatchedRelPos < 1)
		Temp_One_Read.MatchedRelPos = 0;
	GetCloseEnd (CurrentChr, Temp_One_Read);
	if (Temp_One_Read.UP_Close.size ())
		{
			if (ReportLength < Temp_One_Read.ReadLength)
				ReportLength = Temp_One_Read.ReadLength;
			Temp_One_Read.Used = false;
			Temp_One_Read.Unique = true;
			//cout << Temp_One_Read.MatchedD << "\t" << Temp_One_Read.UP_Close.size() << "\t";
			CleanUniquePoints (Temp_One_Read.UP_Close);
			//cout << Temp_One_Read.UP_Close.size() << "\t" << Temp_One_Read.UP_Close[0].Direction << endl;
			Temp_One_Read.CloseEndLength =
				Temp_One_Read.UP_Close[Temp_One_Read.UP_Close.size () - 1].LengthStr;

			//if (Temp_One_Read.UP_Close.size()) {

			if (Temp_One_Read.MatchedD == Plus)
				{
					Temp_One_Read.LeftMostPos =
						Temp_One_Read.UP_Close[0].AbsLoc + 1 -
						Temp_One_Read.UP_Close[0].LengthStr;
					GetPlus++;
				}
			else
				{
					Temp_One_Read.LeftMostPos =
						Temp_One_Read.UP_Close[0].AbsLoc +
						Temp_One_Read.UP_Close[0].LengthStr - Temp_One_Read.ReadLength;
					GetMinus++;
				}

			//if (isInBin(Temp_One_Read)) {
			LeftReads->push_back (Temp_One_Read);
			//}
			if (LeftReads->size () % 10000 == 0)
				{
					std::cout << LeftReads->size () << std::endl;
				}
		}
	if (NotInVector (Temp_One_Read.Tag, VectorTag))
		{
			VectorTag.push_back (Temp_One_Read.Tag);
		}
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
