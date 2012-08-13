//
//  search_MEI.cpp
//  pindel_MEI
//
//  Created by Mark Kroon on 20/06/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

// System libraries.
#include <algorithm>
#include <iostream>
#include <math.h>
#include <sstream>

// Other libraries.
#include "pindel.h"
#include "reader.h"
#include "search_MEI.h"


// Maximal distance between two MEI breakpoints found on opposite strands to
// be assumed to refer to the same event.
const int MAX_MEI_BREAKPOINT_DISTANCE = 10;

// Maximal distance between reads for them to be assigned to the same cluster.
const int MAX_DISTANCE_CLUSTER_READS = 100;

// Minimal/maximal expected insert size.
// Todo: get these values from configuration or estimate them from data.
const int INSERT_SIZE_LOWER_BOUND = 450;
const int INSERT_SIZE_UPPER_BOUND = 550;

// Number of extra bases of reference sequence shown in MEI reporting.
const int REFERENCE_REPORT_LENGTH = 10;

// Size of the read buffer when reading split reads locally around potential MEI sites.
const int SPLIT_READ_BUFFER_SIZE = 50000;


// Comparator for sorting connected reads, first on mapped strand, then on position.
static bool comp_simple_read(simple_read* read1, simple_read* read2) {
    if (read1->strand && ! read2->strand) {
        return true;
    } else if (!read1->strand && read2->strand) {
        return false;
    } else {
        return read1->pos < read2->pos;
    }
}


// Returns true for a given read whether it is concordantly mapped together with its mate.
bool is_concordant(const bam1_t* read) {
    if (read->core.flag & BAM_FUNMAP || read->core.flag & BAM_FMUNMAP) {
        // At least one of the pair is unmapped.
        return false;
    }
    if (read->core.tid != read->core.mtid) {
        // Reads mapping to different chromosomes.
        return false;
    }
    if (bam1_strand(read) == bam1_mstrand(read)) {
        // Reads mapping to the same strand.
        return false;
    }
    
    // Return true if insert size is as expected.
    return abs(read->core.isize) > INSERT_SIZE_LOWER_BOUND && abs(read->core.isize) < INSERT_SIZE_UPPER_BOUND;
}


// Return a samtools-formatted sequence as a human-readable string.
std::string get_sequence(uint8_t* sam_seq, int sam_seq_len) {
    std::string sequence;
    for (int i = 0; i < sam_seq_len; ++i) {
        // Append base letter.
        sequence.append(1, bam_nt16_rev_table[bam1_seqi(sam_seq, i)]);
    }
    return sequence;
}


// Cluster reads by their orientation and mapping location.  Each read linked
// in the given vector will be assigned to a single cluster (in the form of a
// vector) in the output.
static std::vector<std::vector<simple_read*>*> cluster_reads(std::vector<simple_read*>& reads) {
    std::vector<std::vector<simple_read*>*> clusters;
    if (reads.size() == 0) {
        return clusters;
    }
    
    // Sort reads first by mapping strand, then on mapping location.
    sort(reads.begin(), reads.end(), comp_simple_read);
    
    // Initialize first cluster.
    std::vector<simple_read*>* current_cluster = new std::vector<simple_read*>;
    simple_read* read_last_analyzed = reads.at(0);
    current_cluster->push_back(read_last_analyzed);
    
    for (size_t i = 1; i < reads.size(); i++) {
        simple_read* current_read = reads.at(i);
        
        if (read_last_analyzed->strand == current_read->strand && 
            (current_read->pos - read_last_analyzed->pos) <= MAX_DISTANCE_CLUSTER_READS) {
            current_cluster->push_back(current_read);
        } else {
            // Close current cluster and initialize a new.
            clusters.push_back(current_cluster);
            current_cluster = new std::vector<simple_read*>;
            current_cluster->push_back(current_read);
        }
        read_last_analyzed = current_read;
    }
    
    // Save remaining cluster.
    if (current_cluster->size() > 0) {
        clusters.push_back(current_cluster);
    }
    return clusters;
}


// Simply print a read to stdout.
static void print_read(simple_read read) {
    std::cout << read.name << " (";
    std::cout << read.pos << ")";
    if (read.strand) {
        std::cout << " <--";
    } else {
        std::cout << " -->";
    }
    std::cout << std::endl;
}


// Utility function to convert an index on a chromosome to an index on the safe
// representation of the chromosome (with N's added at both ends).
//static int get_safe_chr_index(int true_index) {
//    return true_index + g_SpacerBeforeAfter;
//}

// Utility function to convert an index on the safe representation of the chr
// (with N's added at both ends) to the true location on the chromosome.
static int get_true_chr_index(int safety_index) {
    return safety_index - g_SpacerBeforeAfter;
}



// Find split reads around the given cluster edge position (outer_pos), put
// them in split_reads.
static int get_split_reads_for_cluster(std::vector<bam_info>& bam_sources, bool cluster_strand, int outer_pos,
                                       std::string& chr_name, std::string& chr_sequence,
                                       std::vector<SPLIT_READ>& split_reads) {
    
    for (size_t i = 0; i < bam_sources.size(); i++) {
        bam_info source = bam_sources.at(i);
        
        // Determine bounds of the region where mates of split reads may reside 
        // related to the given cluster.
        int upper_bound;
        int lower_bound;
        if (cluster_strand) {
            lower_bound = outer_pos - source.InsertSize;
            upper_bound = outer_pos + 2 * source.InsertSize;
        } else {
            lower_bound = outer_pos - 2 * source.InsertSize;
            upper_bound = outer_pos + source.InsertSize;
        }
        SearchWindow search_window(chr_name, lower_bound, upper_bound);

        // Read split reads in defined window.
        ReadBuffer read_buffer(SPLIT_READ_BUFFER_SIZE, split_reads, chr_sequence);
        int read_result = ReadInBamReads_SR(source.BamFile.c_str(), chr_name, &chr_sequence, split_reads, 
                                            source.InsertSize, source.Tag, search_window, read_buffer);
        if (read_result > 0) {
            return read_result;
        }
    }
    return 0;
}


// Returns a breakpoint for a cluster of connected reads.  If no viable
// breakpoint can be found, it returns a breakpoint with position -1.
// Note: returned pointer must be deleted by caller.
static MEI_breakpoint* get_breakpoint(std::vector<simple_read*>* cluster, std::vector<bam_info>& bam_sources,
                                      bool cluster_strand, std::string& chr_name, std::string& chr_sequence,
                                      std::vector<SPLIT_READ>& associated_split_reads) {
    
    std::vector<SPLIT_READ> split_reads;
    int outer_read_pos = cluster_strand? cluster->at(cluster->size()-1)->pos : cluster->at(0)->pos;
    get_split_reads_for_cluster(bam_sources, cluster_strand, outer_read_pos, chr_name, chr_sequence, split_reads);
        
    // Search for split reads with a mate close to the outer read of the
    // cluster.  Store candidate breakpoints.
    // Todo: speedup by exploiting the fact that both clusters and split reads are sorted
    // by mapping location.
    std::vector<simple_read> simple_split_reads;
    std::vector<int> candidate_breakpoints;
    for (size_t i = 0; i < split_reads.size(); i++) {
        SPLIT_READ read = split_reads.at(i);
        int read_pos = read.MatchedRelPos;
        bool anchor_strand = read.MatchedD == '-';
        if (anchor_strand) {
            read_pos -= read.getReadLength();
        }
        
        if (cluster_strand == anchor_strand) {
            // Anchor mate of current split read is close to cluster, save the breakpoint.
            int candidate_bp = read.getLastAbsLocCloseEnd();
            candidate_breakpoints.push_back(get_true_chr_index(candidate_bp));
            std::cout << "cand. bp: " << candidate_bp << std::endl;

            // Also save the split read contributing to the candidate breakpoint.
            simple_read simple_split_read(read.Name, cluster->at(0)->tid, candidate_bp, !anchor_strand, 
                                          read.getUnmatchedSeq(), "", ""); // todo, fix this!
            simple_split_reads.push_back(simple_split_read);
        }
    }
    
    if (candidate_breakpoints.size() == 0) {
        // No breakpoint was found.
        return new MEI_breakpoint(-1);
    }
    
    // Return the most frequently occurring candidate.
    sort(candidate_breakpoints.begin(), candidate_breakpoints.end());
    int max_count = 0;
    int current_count = 0;
    int current = -1;
    int best = -1;
    for (size_t i = 0; i < candidate_breakpoints.size(); i++) {
        if (candidate_breakpoints.at(i) == current) {
            current_count += 1;
        } else {
            current = candidate_breakpoints.at(i);
            current_count = 1;
        }
        if (current_count > max_count) {
            max_count = current_count;
            best = current;
        }
    }
    
    // Link associated discordant reads (all reads from cluster) and split reads.
    MEI_breakpoint* bp = new MEI_breakpoint(best);
    bp->cluster_strand = cluster_strand;
    bp->associated_split_reads = simple_split_reads;
    for (size_t i = 0; i < cluster->size(); i++) {
        bp->associated_reads.push_back(*cluster->at(i));
    }
    
    return bp;
}


// See documentation in header file.
void searchMEIBreakpoints(MEI_data& currentState, std::vector<bam_info>& bam_sources, std::string& chr_name,
                          std::string& chr_sequence) {
    std::cout << "####################### Searching MEIs #######################" << std::endl;
    std::cout << currentState.discordant_reads.size() << " discordant reads found" << std::endl;
    
    std::vector<std::vector<simple_read*>*> clusters = cluster_reads(currentState.discordant_reads);
    std::cout << "clustering: found " << clusters.size() << " clusters" << std::endl;
    
    // Find breakpoints per cluster.
    for (size_t i = 0; i < clusters.size(); i++) {
        // print cluster debug info
        std::vector<simple_read*>* cluster = clusters.at(i);
        std::cout << "cluster " << i << ": " << cluster->size() << " entries" << std::endl;
        for (size_t j = 0; j < cluster->size(); j++) {
            std::cout << "    ";
            print_read(*cluster->at(j));
        }
        
       if (cluster->size() < 3) {
            // Fluke cluster, skip it. (If there are very few reads in the cluster,
            // we likely won't find enough split reads supporting an insertion)
            continue;
        }
        
        // Find breakpoint for this cluster
        bool cluster_strand = cluster->at(0)->strand;
        std::vector<SPLIT_READ> associated_split_reads;
        MEI_breakpoint* MEI_bp = get_breakpoint(cluster, bam_sources, cluster_strand, chr_name, chr_sequence,
                                                associated_split_reads);
        // Check breakpoint validity.
        if (MEI_bp->breakpoint_pos >= 0) {
            currentState.MEI_breakpoints.push_back(*MEI_bp);
        }
        delete MEI_bp;
    }
    
    // Debug output, to be deleted.
    std::cout << "bps up until now:" << std::endl;
    for (size_t i = 0; i < currentState.MEI_breakpoints.size(); i++) {
        MEI_breakpoint MEI_bp = currentState.MEI_breakpoints.at(i);
        std::cout  << MEI_bp.breakpoint_tid << " " << MEI_bp.breakpoint_pos << " " << MEI_bp.cluster_strand << " names:" << std::endl;
    }
    std::cout << std::endl;
}


// Return the read pairs that are shared between the two given breakpoints.
static std::vector<std::pair<simple_read, simple_read> > get_shared_reads(MEI_breakpoint& bp1, MEI_breakpoint& bp2) {
    std::vector<std::pair<simple_read, simple_read> > output;
    for (size_t i = 0; i < bp1.associated_reads.size(); i++) {
        for (size_t j = 0; j < bp2.associated_reads.size(); j++) {
            // Check if reads from both breakpoints have the same name.
            if (bp1.associated_reads.at(i).name == bp2.associated_reads.at(j).name) {
                std::pair<simple_read, simple_read> share(bp1.associated_reads.at(i), bp2.associated_reads.at(j));
                output.push_back(share);
            }
        }
    }
    return output;
}


// Returns true when the two given breakpoints are assumed to be associated to
// the same element insertion event.
static bool refer_to_same_event(MEI_breakpoint& bp1, MEI_breakpoint& bp2) {
    // Return true if clusters are on different strands and found breakpoints are
    // close together.
    return bp1.cluster_strand != bp2.cluster_strand && 
           bp1.breakpoint_tid == bp2.breakpoint_tid &&
           abs(bp1.breakpoint_pos - bp2.breakpoint_pos) <= MAX_MEI_BREAKPOINT_DISTANCE;
}


// Mobile element insertion event, constructed from two supporting read clusters.
class MEI_event {
public:
    MEI_event(MEI_breakpoint bp1, MEI_breakpoint bp2) : fwd_cluster_bp(bp1), rev_cluster_bp(bp2) {};
    MEI_breakpoint fwd_cluster_bp; // bp estimated from cluster on forward strand.
    MEI_breakpoint rev_cluster_bp; // bp estimated from cluster on reversed strand.

    // Breakpoints from clusters that have shared reads with the two supporting
    // clusters (provided by fwd_cluster_bp and rev_cluster_bp).
    std::vector<MEI_breakpoint> fwd_cluster_connections;
    std::vector<MEI_breakpoint> rev_cluster_connections;
};




static void report_supporting_reads(MEI_breakpoint& MEI_bp, std::vector<MEI_breakpoint>& supporting_bps, 
                                    std::ostream& out) {
    std::vector<MEI_breakpoint>::iterator support_iter;
    for (support_iter = supporting_bps.begin(); support_iter != supporting_bps.end(); ++support_iter) {
//        std::vector<std::string> supporting_reads = get_shared_reads(MEI_bp, *support_iter);
//        std::vector<std::string>::iterator read_iter;
//        for (read_iter = supporting_reads.begin(); read_iter != supporting_reads.end(); ++read_iter) {
//            std::vector<simple_read>::iterator bp_read_iter;
//            for (bp_read_iter = MEI_bp.associated_reads.begin(); bp_read_iter != MEI_bp.associated_reads.end();
//                 ++bp_read_iter) {
//                simple_read read = *bp_read_iter;
//                if (read.name == *read_iter) {
//                    out << read.name << ", (" << read.tid << ", " << read.pos << "): "
//                        << read.sequence << std::endl;
//                }
//            }
//        }
    }
}


// Find contigs from a cluster of reads.
static void get_contigs(MEI_event& event, std::vector<std::pair<MEI_breakpoint, int> >& contigs) {
    for (size_t i = 0; i < event.fwd_cluster_connections.size(); i++) {
        MEI_breakpoint connected_cluster_bp = event.fwd_cluster_connections.at(i);
        // Todo: implement.
    }
}


// Report MEI events.
static void reportMEIevents(MEI_data& mei_data, std::vector<MEI_event>& insertion_events) {
    for (size_t i = 0; i < insertion_events.size(); i++) {
        MEI_event event = insertion_events.at(i);
        std::cout << "####################################################################################################" << std::endl;
        std::cout << "Mobile element insertion found at: " << event.fwd_cluster_bp.breakpoint_tid << ", ";
        std::cout << event.rev_cluster_bp.breakpoint_pos << "--" << event.fwd_cluster_bp.breakpoint_pos;
        std::cout << std::endl;
        
        std::vector<std::pair<MEI_breakpoint, int> > fwd_contigs;
        get_contigs(event, fwd_contigs);
        
        std::cout << "Reads mapping to this insertion (5' end on forward strand):" << std::endl;
        report_supporting_reads(event.fwd_cluster_bp, event.fwd_cluster_connections, std::cout);
        std::cout << "Reads mapping to this insertion (5' end on reverse strand):" << std::endl;
        report_supporting_reads(event.rev_cluster_bp, event.rev_cluster_connections, std::cout);
    }
}


// See documentation in header file.
void searchMEI(MEI_data& finalState) {
    std::vector<MEI_event> insertion_events;
    // Loop over all breakpoint combinations to see if they refer to the same
    // insertion event.
    for (size_t i = 0; i < finalState.MEI_breakpoints.size(); i++) {
        for (size_t j = i+1; j < finalState.MEI_breakpoints.size(); j++) {
            
            MEI_breakpoint& bp1 = finalState.MEI_breakpoints.at(i);
            MEI_breakpoint& bp2 = finalState.MEI_breakpoints.at(j);
            
            if (refer_to_same_event(bp1, bp2)) {
                // Current 2 selected breakpoints refer to same event.
                
                MEI_event event = bp1.cluster_strand? MEI_event(bp2, bp1) : MEI_event(bp1, bp2);
                
                // Find other clusters (with bp) that contain mates of reads in current inserted element.
                for (size_t k = 0; k < finalState.MEI_breakpoints.size(); k++) {
                    MEI_breakpoint& bp_connect = finalState.MEI_breakpoints.at(k);
                    if (event.fwd_cluster_bp.cluster_strand != bp_connect.cluster_strand &&
                        get_shared_reads(event.fwd_cluster_bp, bp_connect).size() > 0) {
                        // Breakpoint is connected when associated cluster shares read names and is on opposite strand.
                        event.fwd_cluster_connections.push_back(bp_connect);
                    } else if (event.rev_cluster_bp.cluster_strand != bp_connect.cluster_strand &&
                               get_shared_reads(event.rev_cluster_bp, bp_connect).size() > 0) {
                        // Breakpoint is connected when associated cluster shares read names and is on opposite strand.
                        event.rev_cluster_connections.push_back(bp_connect);
                    }
                }
                insertion_events.push_back(event);
            }
        }
    }
    
    // Report MEI events constructed from pairing breakpoints.
    reportMEIevents(finalState, insertion_events);
}


//#####################################################################################
//###                                                                               ###
//###      Code below should be integrated with pindel's existing structure.        ###
//###                                                                               ###
//#####################################################################################

#include "control_state.h"

// Return code for not being able to open a BAM file.
const int ERROR_OPENING_ALIGNMENT_FILE = 100;

// Minimal distance between mapped reads to consider them for detection of
// mobile element insertions.
const int MIN_MEI_MAP_DISTANCE = 8000;


static int fetch_disc_read_callback(const bam1_t* alignment, void* data) {
    MEI_data* mei_data = (MEI_data*) data;
    if (!(alignment->core.flag & BAM_FUNMAP || alignment->core.flag & BAM_FMUNMAP) && // Both ends are mapped.
        !is_concordant(alignment) &&                                                  // Ends map discordantly.
        abs(alignment->core.pos - alignment->core.mpos) > MIN_MEI_MAP_DISTANCE) {     // (Very) large mapping distance.

        // Save alignment as simple_read object.
        simple_read* read = new simple_read(bam1_qname(alignment), alignment->core.tid, alignment->core.pos,
                                            bam1_strand(alignment), 
                                            get_sequence(bam1_seq(alignment), alignment->core.l_qseq));
        mei_data->discordant_reads.push_back(read);
    }
    return 0;
}


static int load_discordant_reads(MEI_data& mei_data, std::vector<bam_info>& bam_sources, std::string& chr_name,
                                 const SearchWindow& window) {
    // Loop over associated bam files.
    for (size_t i = 0; i < bam_sources.size(); i++) {
        // Locate file.
        bam_info source = bam_sources.at(i);
        
        // Setup link to bamfile, its index and header.
        bamFile fp = bam_open(source.BamFile.c_str(), "r");
        bam_index_t *idx = bam_index_load (source.BamFile.c_str());
        bam_header_t *header = bam_header_read (fp);
        bam_init_header_hash (header);
        int tid = bam_get_tid (header, chr_name.c_str ());
        
        // Load discordant reads into mei_data.
        bam_fetch(fp, idx, tid, window.getStart(), window.getEnd(), &mei_data, fetch_disc_read_callback);
        bam_index_destroy(idx);
    }
    return 0;
}



// This function is based on Pindel's main function.  Todo: integrate with pindel's main structure.
int searchMEImain(ControlState& current_state, std::ifstream& FastaFile, UserDefinedSettings* userSettings) {

    MEI_data mei_data;
    int result;
    
    bool SpecifiedChrVisited = false;
    while (SpecifiedChrVisited == false && FastaFile >> current_state.CurrentChrName && !FastaFile.eof()) {
        
        std::string emptystr;
        std::getline(FastaFile, emptystr);
        if (userSettings->loopOverAllChromosomes()) {
            GetOneChrSeq(FastaFile, current_state.CurrentChrSeq, true);
        } else if (current_state.CurrentChrName == userSettings->getRegion()->getTargetChromosomeName()) {   
            // just one chr and this is the correct one
            GetOneChrSeq(FastaFile, current_state.CurrentChrSeq, true);
            SpecifiedChrVisited = true;
        } else {   
            // not build up sequence
            GetOneChrSeq(FastaFile, current_state.CurrentChrSeq, false);
            continue;
        }
        
        CONS_Chr_Size = current_state.CurrentChrSeq.size() - 2 * g_SpacerBeforeAfter;
        g_maxPos = 0;
        CurrentChrMask.resize(current_state.CurrentChrSeq.size());
        for (unsigned int i = 0; i < current_state.CurrentChrSeq.size(); i++) {
            CurrentChrMask[i] = 'N';
        }
        
        LoopingSearchWindow currentWindow(userSettings->getRegion(), CONS_Chr_Size, WINDOW_SIZE, 
                                          current_state.CurrentChrName); 
        // loop over one chromosome
        do {
            result = load_discordant_reads(mei_data, current_state.bams_to_parse, current_state.CurrentChrName,
                                           currentWindow);
            if (result) {
                // something went wrong loading the reads, return error code.
                return result;
            }
            
            searchMEIBreakpoints(mei_data, current_state.bams_to_parse, current_state.CurrentChrName,
                                 current_state.CurrentChrSeq);
            mei_data.discordant_reads.clear();
            currentWindow.next();
            
        } while (!currentWindow.finished());
    }
    
    searchMEI(mei_data);
    return 0;
}









