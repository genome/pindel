/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L. -lbam -lz
 */

#ifndef BAM2DEPTH_H_
#define BAM2DEPTH_H_

#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include "bam.h"

struct aux_t {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region has not been specified
	int min_mapQ;    // mapQ filter
    
    void init();
    void destroy();
};

void aux_t::init() 
{
    fp = NULL;
    iter = NULL;
    min_mapQ = 0;
}

void aux_t::destroy() 
{
    bam_close( fp );
    if (iter != NULL) bam_iter_destroy( iter );
}

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps

// This function reads a BAM alignment from one BAM file.
static int read_bam(aux_t * data, bam1_t *b); // read level filters better go here to avoid pileup


int getChromosomeID( bam_header_t *bamHeaderPtr, const std::string & chromosome, const int startPos, const int endPos);

int bam2depth(const std::string & chromosome, 
              const int startPos, 
              const int endPos, 
              const int minBaseQuality, 
              const int minMappingQuality, 
              const std::vector <std::string> & listOfFiles,
              std::vector <double> & averageCoveragePerBam);

int getRelativeCoverage(const std::string& chromosomeName, 
                        const int startPos, 
                        const int endPos, 
                        const int minBaseQuality, 
                        const int minMappingQuality, 
                        const std::vector <std::string> & listOfFiles, 
                        std::vector <double> & standardizedDepthPerBam );

#endif /* BAM2DEPTH_H_ */
