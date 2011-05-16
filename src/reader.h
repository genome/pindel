/* 
 * File:   reader.h
 * Author: david
 *
 * Created on 13 mei 2011, 11:24
 */

#ifndef READER_H
#define	READER_H

void GetOneChrSeq (std::ifstream & inf_Seq, std::string & CurrentChr,
									 bool WhetherBuildUp);
bool ReadInBamReads (const char *bam_path, const std::string & FragName,
										 std::string * CurrentChr,
										 std::vector < SPLIT_READ > &LeftReads, int InsertSize,
										 std::string Tag, int binStart, int binEnd);
short ReadInRead (std::ifstream & inf_Seq,
									const std::string & CurrentFragName,
									const std::string & CurrentFrag,
									std::vector < SPLIT_READ > &Reads, const unsigned int lowerBinBorder,
									const unsigned int upperBinBorder);

static int fetch_func (const bam1_t * b1, void *data);

#endif /* READER_H */
