/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L. -lbam -lz
 */
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "bam.h"
#include "bam2depth.h"


// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	/* EW:
	if (region has been specified) {
	   bam_iter_read(bamfile, region, &bam_alignment); <- bam_iter_read reads a specified part of a BAM_file
	}
	else { // no region specified
        bam_read1( bamfile, &bam_alignment ); <- bam_read1 reads entire bam-file
	} */
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	// if the mapping quality in the bam_reads is smaller than the desired mapping quality, set the BAM_FUNMAP-bit in b->core.flag to TRUE
	return ret; // the code returned by bam_read1 or bam_iter_read, probably success or failure.
}

int getChromosomeID( bam_header_t *bamHeaderPtr, const std::string & chromosome, const int startPos, const int endPos )
{
    int dummyBegin, dummyEnd;
    int chromosomeID;
    std::stringstream region;
    region << chromosome << ":" << startPos << "-" << endPos;
    bam_parse_region(bamHeaderPtr, region.str().c_str(), &chromosomeID, &dummyBegin, &dummyEnd);
    return chromosomeID;
}

int bam2depth(const std::string& chromosome, const int startPos, const int endPos, const int minBaseQuality, const int minMappingQuality, const std::vector <std::string> & listOfFiles,
    std::vector< double > & averageCoveragePerBam
)
{
	// initialize the auxiliary data structures
	int numberOfBams = listOfFiles.size(); // the number of BAMs on the command line
    aux_t ** auxOfBamfile = calloc( numberOfBams, sizeof (aux_t *) ); // auxOfBamfiles[i] for the i-th input, one pointer per BAM-file

    // set the default region
	int chromosomeID = -1;

    bam_header_t *h = 0; // BAM header of the 1st input
    // EW: the below would usually be done in a constructor. There are reasons that C++ was invented.
	for (int fileIndex = 0; fileIndex < numberOfBams; ++fileIndex) {
		bam_header_t * h_tmp;
		// per bam-file, create an 'aux_t' structure.
        auxOfBamfile[ fileIndex ] = calloc(1, sizeof( aux_t *));
        auxOfBamfile[ fileIndex ]->init();
		auxOfBamfile[ fileIndex ]->fp = bam_open( listOfFiles[ fileIndex ].c_str(), "r"); // open BAM
		auxOfBamfile[ fileIndex ]->min_mapQ = minMappingQuality;                    // set the mapQ filter
		h_tmp = bam_header_read(auxOfBamfile[ fileIndex ]->fp);         // read the BAM header
		if (fileIndex == 0) {
			h = h_tmp; // keep the header of the 1st BAM
			// if the command line defines a region, parse it and replace the 'global' begin, end and chromosome_id by those.
            chromosomeID = getChromosomeID( h, chromosome, startPos, endPos );
		} else { 
           bam_header_destroy(h_tmp); // if not the 1st BAM, trash the header
        }
        
        // if a region is specified and parsed successfully, set the "iter" of the aux structure
		if (chromosomeID >= 0) {  
			bam_index_t *idx = bam_index_load(listOfFiles[ fileIndex ].c_str());  // load the index of the bamfile
            auxOfBamfile[ fileIndex ]->iter = bam_iter_query(idx, chromosomeID, startPos, endPos); // set the iterator
			bam_index_destroy(idx); // the index is not needed any more; phase out of the memory
		}
	}

	// the core multi-pileup loop
	bam_mplp_t multiPileup = bam_mplp_init(numberOfBams, read_bam, (void**)auxOfBamfile); // initialization [read_bam is function!)
    int * coveragePerBam;
	coveragePerBam = calloc(numberOfBams, sizeof(int)); // coverage_per_bam[i] is the number of covering reads from the i-th BAM // EW: comments are nice, but simpler code and clearer variable names would help a lot already.
    void * supportingReadsPerPosition;
	supportingReadsPerPosition = calloc(numberOfBams, sizeof(void*)); // covering_reads_ptr_per_bam[i] points to the array of covering reads (internal in multiPileup)
    std::vector<int> sumOfReadDepths( numberOfBams, 0);
    int position = 0;
	while (bam_mplp_auto(multiPileup, &chromosomeID, &position, coveragePerBam, supportingReadsPerPosition) > 0) { // go to the next covered position
		if (position < startPos) continue; // out of range; skip
        if (position >= endPos) break; // out of range; skip

		// EW: if you want to check: fputs(h->target_name[chromosome_id], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
		for (int fileIndex = 0; fileIndex < numberOfBams; fileIndex++) { // base level filters have to go here
            sumOfReadDepths[ fileIndex ] += coveragePerBam[ fileIndex ];
			for (int supportingReadIndex = 0; supportingReadIndex < coveragePerBam[ fileIndex ]; ++supportingReadIndex) {
				const bam_pileup1_t *supportingBasePtr = supportingReadsPerPosition[ fileIndex ] + supportingReadIndex; // DON'T modfity covering_reads_ptr_per_bam[][] unless you really know
				if (supportingBasePtr->is_del || supportingBasePtr->is_refskip) sumOfReadDepths[ fileIndex ]--; // having dels or refskips at chromosome_id:pos
				else if (bam1_qual(supportingBasePtr->b)[supportingBasePtr->qpos] < minBaseQuality) sumOfReadDepths[ fileIndex ]--; // low base quality
			}
		}
	}
    averageCoveragePerBam.resize( numberOfBams );
    for (int fileIndex=0; fileIndex< numberOfBams; fileIndex++ ) {
        averageCoveragePerBam[ fileIndex ] = (double)sumOfReadDepths[ fileIndex ] / (endPos - startPos );  
        
        auxOfBamfile[ fileIndex ].destroy();
    }
	// now destroy everything: a good illustration why destructors were invented
	free(coveragePerBam); 
    free(supportingReadsPerPosition);
	bam_mplp_destroy(multiPileup);

	bam_header_destroy(h);

	return 0;
}

/* getRelativeCoverage gets the relative coverage of a SV region relative to the surrounding regions; output is a vector that
    yields, per BAM-file, the SV-region's ploidy relative to the surrounding regions; 2.0 meaning a regular diploid region 
    (or an inversion), 1.0 a heterozygous, deletion, 3.0 a heterozygous duplication, etc. */
int getRelativeCoverage(const std::string & chromosomeName, const int startPos, const int endPos, const int minBaseQuality, const int minMappingQuality, const std::vector <std::string> & listOfFiles, std::vector<double> & standardizedDepthPerBam ) 
{
    const int PLOIDY = 2;
    int numberOfBams = listOfFiles.size();
    standardizedDepthPerBam.resize( numberOfBams );
    std::vector<double> avgCoverageOfRegionBeforeSV( numberOfBams, 0.0 );
    std::vector<double> avgCoverageOfSVRegion( numberOfBams, 0.0 );
    std::vector<double> avgCoverageOfRegionAfterSV( numberOfBams, 0.0 );
    
    int regionLength = endPos - startPos;
    int startOfRegionBeforeSV = (startPos - regionLength >=0) ? startPos - regionLength : 0;
    int chromosomeSize = 100000000; //@@@get this from the chromosome somehow
    int endOfRegionAfterSV = (endPos + regionLength > chromosomeSize ) ? chromosomeSize : endPos + regionLength;
    
    bam2depth( chromosomeName, startOfRegionBeforeSV, startPos, minBaseQuality, minMappingQuality, listOfFiles, avgCoverageOfRegionBeforeSV );
    bam2depth( chromosomeName, startPos, endPos, minBaseQuality, minMappingQuality, listOfFiles, avgCoverageOfSVRegion );
    bam2depth( chromosomeName, endPos, endOfRegionAfterSV, minBaseQuality, minMappingQuality, listOfFiles, avgCoverageOfRegionAfterSV );
    
    for (int fileIndex=0; fileIndex < numberOfBams; fileIndex++ ) {
       if ( avgCoverageOfRegionBeforeSV[ fileIndex ] + avgCoverageOfRegionAfterSV[ fileIndex ] == 0 ) { 
            standardizedDepthPerBam[ fileIndex ] = 2.0; 
       }// if SV fills entire chromosome/contig 
       standardizedDepthPerBam[ fileIndex ] = PLOIDY * (2 * avgCoverageOfSVRegion[ fileIndex] ) 
        / ( avgCoverageOfRegionBeforeSV[ fileIndex ] + avgCoverageOfRegionAfterSV[ fileIndex ] );    
    }
    return 0;
}