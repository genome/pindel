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
#include "control_state.h"
#include "genotyping.h"

// This function reads a BAM alignment from one BAM file.
//static int read_bam(aux_t * data, bam1_t *b); // read level filters better go here to avoid pileup

void getRelativeCoverage(const std::string & CurrentChrSeq, const int chromosomeID, const ControlState& allGlobalData, Genotyping & OneSV, const Chromosome * chromosome);

void getRelativeCoverageForFiltering(const int chromosomeID, const ControlState& allGlobalData, Genotyping & OneDEL, const Chromosome * chromosome, const std::vector <unsigned> & SampleIDs);

#endif /* BAM2DEPTH_H_ */
