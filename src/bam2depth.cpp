/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * Modified from Heng Li's bam2depth.c inside samtools
 *
 */
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "bam.h"
#include "bam2depth.h"
#include "genotyping.h"


/*struct aux_t {     // auxiliary data structure
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
}*/

typedef struct {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ;    // mapQ filter
} aux_t;

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps

std::string Spaces(double input);

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	return ret;
}



// This function reads a BAM alignment from one BAM file.
/*static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	// EW:
	if (region has been specified) {
	   bam_iter_read(bamfile, region, &bam_alignment); <- bam_iter_read reads a specified part of a BAM_file
	}
	else { // no region specified
        bam_read1( bamfile, &bam_alignment ); <- bam_read1 reads entire bam-file
	} */
	/*if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	// if the mapping quality in the bam_reads is smaller than the desired mapping quality, set the BAM_FUNMAP-bit in b->core.flag to TRUE
	return ret; // the code returned by bam_read1 or bam_iter_read, probably success or failure.
}*/

int getChromosomeID( bam_header_t *bamHeaderPtr, const std::string & chromosomeName, const int startPos, const int endPos )
{
    int dummyBegin, dummyEnd;
    int chromosomeID;
    std::stringstream region;
    region << chromosomeName << ":" << startPos << "-" << endPos;
    bam_parse_region(bamHeaderPtr, region.str().c_str(), &chromosomeID, &dummyBegin, &dummyEnd);
    //std::cout << "chromosomeID " << chromosomeID << std::endl; 
    return chromosomeID;
}

int bam2depth(const std::string& chromosomeName, const int startPos, const int endPos, const int minBaseQuality, const int minMappingQuality, const std::vector <std::string> & listOfFiles,
    std::vector< double > & averageCoveragePerBam
)
{
	int i, n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0;
	const bam_pileup1_t **plp;
	char *reg = 0; // specified region
	void *bed = 0; // BED data structure
	bam_header_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	bam_mplp_t mplp;



	// initialize the auxiliary data structures
	n = listOfFiles.size(); // the number of BAMs on the command line
	data = (aux_t **)calloc(n, sizeof(void*)); // data[i] for the i-th input
	beg = 0; end = 1<<30; tid = -1;  // set the default region
	beg = startPos;
	end = endPos;
	for (i = 0; i < n; ++i) {
		bam_header_t *htmp;
		data[i] = (aux_t*)calloc(1, sizeof(aux_t));
		data[i]->fp = bam_open(listOfFiles[i].c_str(), "r"); // open BAM
		data[i]->min_mapQ = mapQ;                    // set the mapQ filter
		htmp = bam_header_read(data[i]->fp);         // read the BAM header
		if (i == 0) {
			h = htmp; // keep the header of the 1st BAM
			//if (reg) bam_parse_region(h, reg, &tid, &beg, &end); // also parse the region
			tid = getChromosomeID( h, chromosomeName, startPos, endPos );
		} else bam_header_destroy(htmp); // if not the 1st BAM, trash the header
		if (tid >= 0) { // if a region is specified and parsed successfully
			bam_index_t *idx = bam_index_load(listOfFiles[i].c_str());  // load the index
			data[i]->iter = bam_iter_query(idx, tid, beg, end); // set the iterator
			bam_index_destroy(idx); // the index is not needed any more; phase out of the memory
		}
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = (int*)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM 
	plp = (const bam_pileup1_t **)calloc(n, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)
   std::vector<int> sumOfReadDepths( n, 0);
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		//fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
		for (i = 0; i < n; ++i) { // base level filters have to go here
			int j, m = 0;
			for (j = 0; j < n_plp[i]; ++j) {
				const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
				if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
				else if (bam1_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
			}
			//printf("\t%d", n_plp[i] - m); // this the depth to output
			sumOfReadDepths[ i ] += n_plp[i] - m;
		}
		//putchar('\n');
	}
  	averageCoveragePerBam.resize( n );
    for (int fileIndex=0; fileIndex< n; fileIndex++ ) {
        averageCoveragePerBam[ fileIndex ] = (double)sumOfReadDepths[ fileIndex ] / (end - beg );  
    }
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);

	bam_header_destroy(h);
	for (i = 0; i < n; ++i) {
		bam_close(data[i]->fp);
		if (data[i]->iter) bam_iter_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data); free(reg);
	if (bed) bed_destroy(bed);
	return 0;
}
/*
int bam2depth(const std::string& chromosome, int chromosomeID, const int startPos, const int endPos, const int minBaseQuality, const int minMappingQuality, const std::vector <std::string> & listOfFiles,
    std::vector< double > & averageCoveragePerBam
)
{   
    std::cout << "bam2depth 1" << std::endl;
	// initialize the auxiliary data structures
	int numberOfBams = listOfFiles.size(); // the number of BAMs on the command line
   aux_t ** auxOfBamfile =(aux_t**) calloc( numberOfBams, sizeof (aux_t *) ); // auxOfBamfiles[i] for the i-th input, one pointer per BAM-file
    std::cout << "bam2depth 2" << std::endl;
    // set the default region
	//int chromosomeID = -1;
    bam_header_t *h = 0; // BAM header of the 1st input
    // EW: the below would usually be done in a constructor. There are reasons that C++ was invented.
	for (int fileIndex = 0; fileIndex < numberOfBams; ++fileIndex) {
		//bam_header_t * h_tmp;
		// per bam-file, create an 'aux_t' structure.
        auxOfBamfile[ fileIndex ] = (aux_t*) calloc(1, sizeof( aux_t *));
        auxOfBamfile[ fileIndex ]->init();
		auxOfBamfile[ fileIndex ]->fp = bam_open( listOfFiles[ fileIndex ].c_str(), "r"); // open BAM
		auxOfBamfile[ fileIndex ]->min_mapQ = minMappingQuality;                    // set the mapQ filter
        
		h_tmp = bam_header_read(auxOfBamfile[ fileIndex ]->fp);  // read the BAM header
        
        std::cout << "bam2depth 2g" << std::endl;
		if (fileIndex == 0) {
			h = h_tmp; // keep the header of the 1st BAM
			// if the command line defines a region, parse it and replace the 'global' begin, end and chromosome_id by those.
            chromosomeID = getChromosomeID( h, chromosome, startPos, endPos );
		} else { 
           bam_header_destroy(h_tmp); // if not the 1st BAM, trash the header
        }
         */
        // if a region is specified and parsed successfully, set the "iter" of the aux structure
	/*	if (chromosomeID >= 0) {  
			bam_index_t *idx = bam_index_load(listOfFiles[ fileIndex ].c_str());  // load the index of the bamfile
            auxOfBamfile[ fileIndex ]->iter = bam_iter_query(idx, chromosomeID, startPos, endPos); // set the iterator
			bam_index_destroy(idx); // the index is not needed any more; phase out of the memory
		}
	}
	// the core multi-pileup loop
    std::cout << "bam2depth 3" << std::endl;
	bam_mplp_t multiPileup = bam_mplp_init(numberOfBams, read_bam, (void**)auxOfBamfile); // initialization [read_bam is function!)
    std::cout << "bam2depth 4" << std::endl;
    int * coveragePerBam =(int *) calloc(numberOfBams, sizeof(int)); // coverage_per_bam[i] is the number of covering reads from the i-th BAM // EW: comments are nice, but simpler code and clearer variable names would help a lot already.
    std::cout << "bam2depth 5" << std::endl;
    const bam_pileup1_t** supportingReadsPerPosition =(const bam_pileup1_t**) calloc(numberOfBams, sizeof(void*)); // covering_reads_ptr_per_bam[i] points to the array of covering reads (internal in multiPileup)
    std::cout << "bam2depth 6" << std::endl;
    std::vector<int> sumOfReadDepths( numberOfBams, 0);
    std::cout << "bam2depth 7" << std::endl;
    int position = 0;
	while (bam_mplp_auto(multiPileup, &chromosomeID, &position, coveragePerBam, supportingReadsPerPosition) > 0) { // go to the next covered position
		if (position < startPos) continue; // out of range; skip
        if (position >= endPos) break; // out of range; skip
        std::cout << "bam_mplp_auto 1" << std::endl;
		// EW: if you want to check: fputs(h->target_name[chromosome_id], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
		for (int fileIndex = 0; fileIndex < numberOfBams; fileIndex++) { // base level filters have to go here
            sumOfReadDepths[ fileIndex ] += coveragePerBam[ fileIndex ];
			for (int supportingReadIndex = 0; supportingReadIndex < coveragePerBam[ fileIndex ]; ++supportingReadIndex) {
				const bam_pileup1_t *supportingBasePtr = supportingReadsPerPosition[ fileIndex ] + supportingReadIndex; // DON'T modfity covering_reads_ptr_per_bam[][] unless you really know
				if (supportingBasePtr->is_del || supportingBasePtr->is_refskip) sumOfReadDepths[ fileIndex ]--; // having dels or refskips at chromosome_id:pos
				else if (bam1_qual(supportingBasePtr->b)[supportingBasePtr->qpos] < minBaseQuality) sumOfReadDepths[ fileIndex ]--; // low base quality
			}
		}
        std::cout << "bam_mplp_auto 2" << std::endl;
	}
std::cout << "bam2depth 8" << std::endl;
    averageCoveragePerBam.resize( numberOfBams );
    for (int fileIndex=0; fileIndex< numberOfBams; fileIndex++ ) {
        averageCoveragePerBam[ fileIndex ] = (double)sumOfReadDepths[ fileIndex ] / (endPos - startPos );  
        std::cout << "Read count in fileIndex " << fileIndex << " " << sumOfReadDepths[ fileIndex ] << std::endl;
    }
    
	// now destroy everything: a good illustration why destructors were invented
	free(coveragePerBam); 
    free(supportingReadsPerPosition);
	bam_mplp_destroy(multiPileup); // [bam_plp_destroy] memory leak: 3. Continue anyway.
     inside samtools pileup
     mp_free(multiPileup->mp, multiPileup->dummy);
     mp_free(multiPileup->mp, multiPileup->head);
     //if (iter->mp->cnt != 0)
     //fprintf(stderr, "[bam_plp_destroy] memory leak: %d. Continue anyway.\n", iter->mp->cnt);
     mp_destroy(multiPileup->mp);
     if (multiPileup->b) bam_destroy1(multiPileup->b);
     free(multiPileup->plp);
     free(multiPileup); 
    */ 
  /*  std::cout << "bam2depth 9" << std::endl;
	bam_header_destroy(h);
    for (int fileIndex=0; fileIndex< numberOfBams; fileIndex++ ) {
        bam_close(auxOfBamfile[fileIndex]->fp);
        if (auxOfBamfile[fileIndex]->iter) bam_iter_destroy(auxOfBamfile[fileIndex]->iter);
        free(auxOfBamfile[fileIndex]);
        //auxOfBamfile[ fileIndex ]->destroy();
    }
    free(auxOfBamfile);
    std::cout << "bam2depth 10" << std::endl;
	return 0;
}*/

/* getRelativeCoverage gets the relative coverage of a SV region relative to the surrounding regions; output is a vector that
    yields, per BAM-file, the SV-region's ploidy relative to the surrounding regions; 2.0 meaning a regular diploid region 
    (or an inversion), 1.0 a heterozygous, deletion, 3.0 a heterozygous duplication, etc. 
	'internal': takes cleaned data set as argument.
*/
void getRelativeCoverageInternal(const std::string & chromosomeName, const int chromosomeSize, const int chromosomeID, const int startPos, const int endPos, const int minBaseQuality, 		const int minMappingQuality, const std::vector <std::string> & listOfFiles, std::vector <double> & standardizedDepthPerBam ) 
{
    const int PLOIDY = 2;
    int numberOfBams = listOfFiles.size();
    standardizedDepthPerBam.resize( numberOfBams );
    std::vector<double> avgCoverageOfRegionBeforeSV( numberOfBams, 0.0 );
    std::vector<double> avgCoverageOfSVRegion( numberOfBams, 0.0 );
    std::vector<double> avgCoverageOfRegionAfterSV( numberOfBams, 0.0 );
    
    int regionLength = endPos - startPos;
    int startOfRegionBeforeSV = (startPos - regionLength >=0) ? startPos - regionLength : 0;
    int endOfRegionAfterSV = (endPos + regionLength > chromosomeSize ) ? chromosomeSize : endPos + regionLength;
    //std::cout << "startOfRegionBeforeSV, Startpos, endpos, endofregionm " << startOfRegionBeforeSV << " " << startPos << " " << endPos << " " << endOfRegionAfterSV << "\n"; 
    //std::cout << "ChrSize " << chromosomeSize << "\n"; 
    //std::cout << "1" << std::endl;
    std::cout.precision(2);
    bam2depth( chromosomeName, startOfRegionBeforeSV, startPos, minBaseQuality, minMappingQuality, listOfFiles, avgCoverageOfRegionBeforeSV );
    std::cout << "before  ";
    for (int fileIndex=0; fileIndex < numberOfBams; fileIndex++ ) {
        std::cout << Spaces(avgCoverageOfRegionBeforeSV[fileIndex]) << std::fixed << avgCoverageOfRegionBeforeSV[fileIndex];
    }
    std::cout << std::endl;
    //std::cout << "2" << std::endl;
    bam2depth( chromosomeName, startPos, endPos, minBaseQuality, minMappingQuality, listOfFiles, avgCoverageOfSVRegion );
    std::cout << "in      "; 
    for (int fileIndex=0; fileIndex < numberOfBams; fileIndex++ ) {
        std::cout << Spaces(avgCoverageOfSVRegion[fileIndex]) << std::fixed << avgCoverageOfSVRegion[fileIndex];
    }
    std::cout << std::endl;
    //std::cout << "3" << std::endl;
    bam2depth( chromosomeName, endPos, endOfRegionAfterSV, minBaseQuality, minMappingQuality, listOfFiles, avgCoverageOfRegionAfterSV );
    std::cout << "after   ";
    for (int fileIndex=0; fileIndex < numberOfBams; fileIndex++ ) {
        std::cout << Spaces(avgCoverageOfRegionAfterSV[fileIndex])  << std::fixed << avgCoverageOfRegionAfterSV[fileIndex];
    }
    std::cout << std::endl;
    //std::cout << "4" << std::endl;
    for (int fileIndex=0; fileIndex < numberOfBams; fileIndex++ ) {
       if ( avgCoverageOfRegionBeforeSV[ fileIndex ] + avgCoverageOfRegionAfterSV[ fileIndex ] == 0 ) { 
            standardizedDepthPerBam[ fileIndex ] = 2.0; 
       } // if SV fills entire chromosome/contig 
       else {
           standardizedDepthPerBam[ fileIndex ] = PLOIDY * (2 * avgCoverageOfSVRegion[ fileIndex] ) 
           / ( avgCoverageOfRegionBeforeSV[ fileIndex ] + avgCoverageOfRegionAfterSV[ fileIndex ] );   
       } 
    }
}

void getRelativeCoverage(const std::string & CurrentChrSeq, const int chromosomeID, const ControlState& allGlobalData, Genotyping & OneSV)
                         //const int startPos, const int endPos, std::vector<double> & standardizedDepthPerBam ) RD_signals
{
    //std::cout << "start " << OneSV.ChrA << "\t" << OneSV.PosA << "\t" << OneSV.CI_A << "\t"
    //<< OneSV.ChrB << "\t" << OneSV.PosB << "\t" << OneSV.CI_B << std::endl;
    std::cout.precision(2);
    //std::cout << "entering getRelativeCoverage" << std::endl;
    const int startPos = OneSV.PosA;
    const int endPos   = OneSV.PosB;
    
    std::string chromosomeName = allGlobalData.CurrentChrName;
	int chromosomeSize = CurrentChrSeq.size() - 2 * g_SpacerBeforeAfter;
	const int MIN_BASE_QUALITY_READDEPTH = 0;
	const int MIN_MAPPING_QUALITY_READDEPTH = 20;
	std::vector<std::string> listOfFiles;
	const std::vector<bam_info>& bamFileData = allGlobalData.bams_to_parse;
	for (unsigned int fileIndex=0; fileIndex<bamFileData.size(); fileIndex++ ) {
		listOfFiles.push_back( bamFileData[ fileIndex ].BamFile );
	}
    //std::cout << "before  getRelativeCoverageInternal " << OneSV.ChrA << "\t" << OneSV.PosA << "\t" << OneSV.CI_A << "\t"
    //<< OneSV.ChrB << "\t" << OneSV.PosB << "\t" << OneSV.CI_B << std::endl;
	getRelativeCoverageInternal( chromosomeName, chromosomeSize, chromosomeID, startPos, endPos, MIN_BASE_QUALITY_READDEPTH, MIN_MAPPING_QUALITY_READDEPTH, listOfFiles, 	
		OneSV.RD_signals);  
    //std::cout << "after  getRelativeCoverageInternal " << OneSV.ChrA << "\t" << OneSV.PosA << "\t" << OneSV.CI_A << "\t"
    //          << OneSV.ChrB << "\t" << OneSV.PosB << "\t" << OneSV.CI_B;
    std::cout << "Genotyping per sample";
    for (unsigned RD_index = 0; RD_index < OneSV.RD_signals.size(); RD_index++) {
        std::cout << "\t" << std::fixed << OneSV.RD_signals[RD_index];
    }
    std::cout << std::endl;
    //std::cout << "leaving getRelativeCoverage" << std::endl;
} 

std::string Spaces(double input) {
   std::string output; 
   if (input > 10000000) output = " ";
   else if (input > 1000000) output = "  ";
   else if (input > 100000) output = "   ";
   else if (input > 10000) output = "    ";
   else if (input > 1000) output = "     ";
   else if (input > 100) output = "      ";
   else if (input > 10) output = "       ";
   else if (input > 1) output = "        ";
   else if (input < 1) output = "        ";
   else output = " "; 
   return output; 
}
