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

// Other libraries.
#include "pindel.h"
#include "reader.h"
#include "search_MEI.h"
#include "smith_waterman_alignment.h"
#include "search_MEI_util.h"


// Maximal distance between two MEI breakpoints found on opposite strands to
// be assumed to refer to the same event.
const int MAX_MEI_BREAKPOINT_DISTANCE = 10;

// Maximal distance between reads for them to be assigned to the same cluster.
const int MAX_DISTANCE_CLUSTER_READS = 100;

// Minimal/maximal expected insert size.
// Todo: get these values from configuration or estimate them from data.
const int INSERT_SIZE_LOWER_BOUND = 450;
const int INSERT_SIZE_UPPER_BOUND = 550;

// Size of the read buffer when reading split reads locally around potential MEI sites.
const int SPLIT_READ_BUFFER_SIZE = 50000;

// String used to characterize comment in text output.
const std::string COMMENT_PREFIX = "# ";

// A contig, being a collection of overlapping reads paired with relative 
// indexes indicating their distance in mapping.
typedef std::vector<std::pair<simple_read, int> > contig;


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


// Find split reads around the given cluster edge position (outer_pos), put
// them in split_reads.
static int get_split_reads_for_cluster(std::vector<bam_info>& bam_sources, bool cluster_strand, int outer_pos,
                                       const Chromosome* chromosome, std::vector<SPLIT_READ>& split_reads) {
    
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
        SearchWindow search_window(chromosome, lower_bound, upper_bound);
        
        std::vector<SPLIT_READ> HangingReads; // Kai: those are reads without closed end mapped, or hanging reads.

        // Read split reads in defined window.
        ReadBuffer read_buffer(SPLIT_READ_BUFFER_SIZE, split_reads, HangingReads, chromosome->getSeq());

        int read_result = ReadInBamReads_SR(source.BamFile.c_str(), chromosome->getName(), &chromosome->getSeq(), 
                                            split_reads, HangingReads, source.InsertSize, source.Tag, search_window, read_buffer);
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
                                      char cluster_strand, const Chromosome* chromosome) {
    
    std::vector<SPLIT_READ> split_reads;
    int outer_read_pos = cluster_strand? cluster->at(cluster->size()-1)->pos : cluster->at(0)->pos;
    get_split_reads_for_cluster(bam_sources, cluster_strand, outer_read_pos, chromosome, split_reads);
        
    // Search for split reads with a mate close to the outer read of the
    // cluster.  Store candidate breakpoints.
    // Todo: speedup by exploiting the fact that both clusters and split reads are sorted
    // by mapping location.
    std::vector<simple_read> simple_split_reads;
    std::vector<int> candidate_breakpoints;
    for (size_t i = 0; i < split_reads.size(); i++) {
        SPLIT_READ read = split_reads.at(i);
        int read_pos = read.MatchedRelPos;
        char anchor_strand = read.MatchedD;
        if (anchor_strand == Minus) {
            read_pos -= read.getReadLength();
        }
        
        if (cluster_strand == anchor_strand) {
            // Anchor mate of current split read is close to cluster, save the breakpoint.
            int candidate_bp = read.getLastAbsLocCloseEnd();
            candidate_breakpoints.push_back(get_true_chr_index(candidate_bp));
//            std::cout << "cand. bp: " << candidate_bp << std::endl;

            // Store the unmatched sequence as it should be matched on the opposite strand of
            // the mapped mate.
            std::string whole_sequence;
            std::string mapped_part;
            std::string unmapped_part;
            if (anchor_strand == Plus) {
                whole_sequence = ReverseComplement(read.getUnmatchedSeq());
                mapped_part = whole_sequence.substr(0, read.CloseEndLength);
                unmapped_part = whole_sequence.substr(read.CloseEndLength, whole_sequence.length());
            } else {
                whole_sequence = read.getUnmatchedSeq();
                mapped_part = whole_sequence.substr(whole_sequence.length() - read.CloseEndLength, 
                                                    whole_sequence.length());
                unmapped_part = whole_sequence.substr(0, whole_sequence.length() - read.CloseEndLength);
            }

            char strand = (anchor_strand == Plus)? Minus : Plus;
            // Save the split read contributing to the candidate breakpoint.
            simple_read simple_split_read(read.Name, cluster->at(0)->tid, candidate_bp, strand, whole_sequence, 
                                          mapped_part, unmapped_part);
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
    bp->chromosome_name = chromosome->getName();
    bp->associated_split_reads = simple_split_reads;
    for (size_t i = 0; i < cluster->size(); i++) {
        bp->associated_reads.push_back(*cluster->at(i));
    }
    
    return bp;
}


// See documentation in header file.
void searchMEIBreakpoints(MEI_data& currentState, std::vector<bam_info>& bam_sources, const Chromosome* chromosome) {
//    std::cout << "####################### Searching MEIs #######################" << std::endl;
//    std::cout << currentState.discordant_reads.size() << " discordant reads found" << std::endl;
    
    std::vector<std::vector<simple_read*>*> clusters = cluster_reads(currentState.discordant_reads);
//    std::cout << "clustering: found " << clusters.size() << " clusters" << std::endl;
    
    // Find breakpoints per cluster.
    for (size_t i = 0; i < clusters.size(); i++) {
        // print cluster debug info
        std::vector<simple_read*>* cluster = clusters.at(i);
//        std::cout << "cluster " << i << ": " << cluster->size() << " entries" << std::endl;
//        for (size_t j = 0; j < cluster->size(); j++) {
//            std::cout << "    ";
//            print_read(*cluster->at(j));
//        }
        
       if (cluster->size() < 3) {
            // Fluke cluster, skip it. (If there are very few reads in the cluster,
            // we likely won't find enough split reads supporting an insertion)
            continue;
        }
        
        // Find breakpoint for this cluster
        char cluster_strand = cluster->at(0)->strand;
        MEI_breakpoint* MEI_bp = get_breakpoint(cluster, bam_sources, cluster_strand, chromosome);
        // Check breakpoint validity.
        if (MEI_bp->breakpoint_pos >= 0) {
            currentState.MEI_breakpoints.push_back(*MEI_bp);
        }
        delete MEI_bp;
    }
    
    // Debug output, to be deleted.
//    std::cout << "bps up until now:" << std::endl;
//    for (size_t i = 0; i < currentState.MEI_breakpoints.size(); i++) {
//        MEI_breakpoint MEI_bp = currentState.MEI_breakpoints.at(i);
//        std::cout  << MEI_bp.breakpoint_tid << " " << MEI_bp.breakpoint_pos << " " << MEI_bp.cluster_strand << " names:" << std::endl;
//    }
//    std::cout << std::endl;
}


// Return the read pairs that are shared between the two given breakpoints.
static std::vector<std::pair<simple_read, simple_read> > get_shared_reads(MEI_breakpoint& bp1, MEI_breakpoint& bp2) {
    std::vector<std::pair<simple_read, simple_read> > output;
    for (size_t i = 0; i < bp1.associated_reads.size(); i++) {
        for (size_t j = 0; j < bp2.associated_reads.size(); j++) {
            // Check if reads from both breakpoints have the same name.
            if (base_read_name(bp1.associated_reads.at(i).name) == base_read_name(bp2.associated_reads.at(j).name)) {
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



static void report_supporting_reads(std::vector<simple_read>& reads, std::ostream& out) {
    out << "# All supporting reads for this insertion:" << std::endl;
    sort(reads.begin(), reads.end(), comp_simple_read_pos);
    std::vector<simple_read>::iterator support_iter;
    for (support_iter = reads.begin(); support_iter != reads.end(); ++support_iter) {
        // Only write the part of the read that falls inside the inserted element.
        if ((*support_iter).is_split) {
            out << "?\t?\t?\t" << (*support_iter).name << "\t" << (*support_iter).unmapped_sequence << std::endl;
        } else {
            out << (*support_iter).tid << "\t" << (*support_iter).pos << "\t" << (*support_iter).strand << "\t" <<
                   (*support_iter).name << "\t" << (*support_iter).sequence << std::endl;
        }
    }
}




// Find contigs from a cluster of reads.
static void get_contigs(std::vector<simple_read> reads, std::vector<contig>& contigs) {
    size_t contig_size;
    bool has_overlap;
    int alignment_index;
    
    // Go on until all reads are assigned to a contig.
    while (!reads.empty()) {
        contig current_contig;
        do {
            contig_size = current_contig.size();
            bool break_read_loop = false;
            
            // Loop over reads to add to current contig.
            for (size_t i = 0; i < reads.size() && !break_read_loop; i++) {
                simple_read current_read = reads.at(i);
                if (contig_size == 0) {
                    // Contig is empty, add current read.
                    std::pair<simple_read, int> init(current_read, 0);
                    current_contig.push_back(init);
                    reads.erase(reads.begin() + i);
                    break;
                }
                for (size_t j = 0; j < contig_size; j++) {
                    std::pair<simple_read, int> contig_indexed_read = current_contig.at(j);
                    alignment_index = get_alignment_pos(contig_indexed_read.first.sequence, current_read.sequence, 
                                                        has_overlap);
//                    std::cout << contig_indexed_read.first.name << "-and-" << current_read.name << " diff: " << alignment_index << " valid:" << has_overlap << std::endl;
                    if (has_overlap) {
                        // Current read overlaps with one in contig, add it.
                        std::pair<simple_read, int> appendum(current_read, 
                                                             alignment_index + contig_indexed_read.second);
                        current_contig.push_back(appendum);
                        reads.erase(reads.begin() + i);
                        break_read_loop = true;
                        break; // Break out of loop over reads as well.
                    }
                }
            }
        } while (current_contig.size() > contig_size); // Continue as long contig size increases.
        contigs.push_back(current_contig);
    }
//    std::cout << "# contigs: " << contigs.size() << std::endl;
}


static void set_reference_highlight(std::string& reference, int breakpoint, bool highlight_until_bp) {
    std::string::iterator ref_iter;
    int counter = 0;
    for (ref_iter = reference.begin(); ref_iter != reference.end(); ++ref_iter) {
        if ((highlight_until_bp && counter < breakpoint) || (!highlight_until_bp && counter >= breakpoint)) {
            *ref_iter = std::toupper((unsigned char) *ref_iter);
        } else {
            *ref_iter = std::tolower((unsigned char) *ref_iter);
        }
        counter++;
    }
}


// Get
static void get_bp_supporting_reads(MEI_breakpoint& breakpoint, std::vector<MEI_breakpoint>& breakpoint_connections,
                                    std::vector<simple_read>& supporting_reads) {
    // Add shared reads of connected clusters.
    for (size_t i = 0; i < breakpoint_connections.size(); i++) {
        MEI_breakpoint connected_cluster_bp = breakpoint_connections.at(i);
        std::vector<std::pair<simple_read, simple_read> > shared_reads = get_shared_reads(breakpoint, 
                                                                                          connected_cluster_bp);
        for (size_t j = 0; j < shared_reads.size(); j++) {
            supporting_reads.push_back(shared_reads.at(j).second);
        }
    }
    for (size_t j = 0; j < breakpoint.associated_split_reads.size(); j++) {
        supporting_reads.push_back(breakpoint.associated_split_reads.at(j));
    }
}


// Report supporting reads for a given breakpoint.
static void report_breakpoint_support(Genome& genome, MEI_breakpoint& breakpoint, std::vector<contig>& contigs,
                                      bool fiveprime_end, std::ostream& out,
                                      std::vector<simple_read>& supporting_reads) {    
    int min_index;
    int max_index;
    int bp_on_screen;
    int offset = 0;
    for (size_t i = 0; i < contigs.size(); i++) {
        contig current_contig = contigs.at(i);
        if (current_contig.size() == 0) {
            // Skip any empty contigs.
            continue;
        }
        if (fiveprime_end) {
            sort(current_contig.begin(), current_contig.end(), comp_indexed_read);
            std::pair<simple_read, int> idx_split_read = *current_contig.begin();
            if (!idx_split_read.first.is_split) {
                // Do not report contig, it is too far inside the inserted element.
                // Note: first read in sorted contig should be a split read if it were to be 
                // at the edge of the element.
                continue;
            }
        
            min_index = idx_split_read.second;
            max_index = (*(current_contig.end()-1)).second + (*(current_contig.end()-1)).first.sequence.length();
//            std::cout << current_contig.at(current_contig.size() - 1).first.name << std::endl;
            bp_on_screen = idx_split_read.first.mapped_sequence.length() + (idx_split_read.second - min_index);
            offset = 1; // todo: resolve this seemingly arbitrary addition for 5' end.
        } else {
            sort(current_contig.rbegin(), current_contig.rend(), comp_indexed_read);
            std::pair<simple_read, int> idx_split_read = *current_contig.begin();
            if (!idx_split_read.first.is_split) {
                // Do not report contig, it is too far inside the inserted element.
                // Note: first read in sorted contig should be a split read if it were to be 
                // at the edge of the element.
                continue;
            }
            
            min_index = (*(current_contig.end()-1)).second;
            max_index = idx_split_read.second + idx_split_read.first.sequence.length();
            bp_on_screen = idx_split_read.first.unmapped_sequence.length() + (idx_split_read.second - min_index);
        }

        std::string ref_prefix = "# Reference: ...";
        bp_on_screen -= (ref_prefix.length() - COMMENT_PREFIX.length());
        std::string reference = get_fasta_subseq(genome, breakpoint.chromosome_name, 
                                                 breakpoint.breakpoint_pos - bp_on_screen + offset, 
                                                 max_index - ref_prefix.length() - COMMENT_PREFIX.length() - min_index);
        set_reference_highlight(reference, bp_on_screen, fiveprime_end);
        out << ref_prefix << reference << std::endl;

        
        for (size_t j = 0; j < current_contig.size(); j++) {
            std::pair<simple_read, int> indexed_read = current_contig.at(j);
            //                out << (indexed_read.second - base) << " ";
            out << "# "; // Comment prefix.
            // This is the c++ way to print a repetition of characters. Really?
            for (int k = 0; k < (indexed_read.second - min_index); k++) {
                out << " ";
            }
            out << indexed_read.first.sequence << " (" << indexed_read.first.name << ") " << std::endl;
        }
    }
}


// Sum and return the spanned width of given contigs.
static int sum_contig_spans(std::vector<contig> contigs, bool fiveprime_end) {
    int summed_spans = 0;
    std::vector<contig>::iterator contig_iter;
    // Compute min and max positions of reads inside each contig.
    for (contig_iter = contigs.begin(); contig_iter != contigs.end(); ++contig_iter) {
        int min = INT_MAX;
        int max_pos = INT_MIN;
        int max = INT_MIN;
        contig::iterator idx_read_iter;
        for (idx_read_iter = (*contig_iter).begin(); idx_read_iter != (*contig_iter).end(); ++idx_read_iter) {
            int pos = (*idx_read_iter).second;
            // Only count the unmapped part of split reads.
            if ((*idx_read_iter).first.is_split) {
                if (fiveprime_end) {
                    pos += (*idx_read_iter).first.mapped_sequence.length();
                } else {
                    pos -= (*idx_read_iter).first.mapped_sequence.length();
                }
            }
            if (pos < min) {
                min = pos;
            }
            if (pos > max_pos) {
                max_pos = pos;
                max = pos + (*idx_read_iter).first.sequence.length();
            }
        }
        if (min > INT_MIN && max_pos < INT_MAX) {
            summed_spans += (max - min);
        }
    }
    return summed_spans;
}

// Report MEI events.
static void reportMEIevents(MEI_data& mei_data, std::vector<MEI_event>& insertion_events, Genome& genome,
                            std::ostream& out) {
    int MEI_counter = 0;
    for (size_t i = 0; i < insertion_events.size(); i++) {
        MEI_event event = insertion_events.at(i);
        
        // Gather all supporting reads for current event.
        std::vector<simple_read> fwd_reads;
        std::vector<simple_read> rev_reads;
        std::vector<simple_read> all_reads;
        get_bp_supporting_reads(event.fwd_cluster_bp, event.fwd_cluster_connections, fwd_reads);
        get_bp_supporting_reads(event.rev_cluster_bp, event.rev_cluster_connections, rev_reads);
        int fwd_split_count = 0;
        for (std::vector<simple_read>::iterator fwd_iter = fwd_reads.begin(); fwd_iter != fwd_reads.end(); ++fwd_iter) {
            all_reads.push_back(*fwd_iter);
            if ((*fwd_iter).is_split) {
                fwd_split_count++;
            }
        }
        int rev_split_count = 0;
        for (std::vector<simple_read>::iterator rev_iter = rev_reads.begin(); rev_iter != rev_reads.end(); ++rev_iter) {
            all_reads.push_back(*rev_iter);
            if ((*rev_iter).is_split) {
                rev_split_count++;
            }
        }
        
        // Assemble supporting reads into contigs.
        std::vector<contig> fwd_contigs;
        get_contigs(fwd_reads, fwd_contigs);
        std::vector<contig> rev_contigs;
        get_contigs(rev_reads, rev_contigs);
        
        // Get lower bound on size of element.
        int summed_contig_span = sum_contig_spans(fwd_contigs, true) + sum_contig_spans(rev_contigs, false);
        
        out << "####################################################################################################" << std::endl;
        // Print machine summary line.
        out << MEI_counter << "\t" << "MEI" << "\t" << event.fwd_cluster_bp.chromosome_name << "\t" << 
               event.rev_cluster_bp.breakpoint_pos << "\t" << event.fwd_cluster_bp.breakpoint_pos << "\t" <<
               summed_contig_span << "\t" << all_reads.size() << "\t" << (fwd_reads.size() - fwd_split_count) << "\t" <<
               fwd_split_count << "\t" << (rev_reads.size() - rev_split_count) << "\t" << rev_split_count << std::endl;
        // Print human-readable summary lines.
        out << COMMENT_PREFIX << "Mobile element insertion (MEI) found on chromosome '" << 
               event.fwd_cluster_bp.chromosome_name << "', between " << event.rev_cluster_bp.breakpoint_pos << "--" << 
               event.fwd_cluster_bp.breakpoint_pos << ", the inserted element is at least " << summed_contig_span << 
               " bp long." << std::endl;
        out << COMMENT_PREFIX << "Found " << all_reads.size() << " supporting reads, of which " <<
               (fwd_reads.size() - fwd_split_count) << " discordant reads and " << fwd_split_count << 
               " split reads at 5' end, " << (rev_reads.size() - rev_split_count) << " discordant reads and " <<
               rev_split_count << " split reads at 3' end." << std::endl;
        
        // Print support for breakpoint at 5' end.
        out << COMMENT_PREFIX << "Supporting reads for insertion location (5' end):" << std::endl;
        report_breakpoint_support(genome, event.fwd_cluster_bp, fwd_contigs, true, out, fwd_reads);
        
        // Print support for breakpoint at 3' end.
        out << COMMENT_PREFIX << "Supporting reads for insertion location (3' end):" << std::endl;
        report_breakpoint_support(genome, event.rev_cluster_bp, rev_contigs, false, out, rev_reads);
        
        // Print all supporting reads and read fragments for the inserted element.
        report_supporting_reads(all_reads, out);
        MEI_counter++;
    }
}


// See documentation in header file.
void searchMEI(MEI_data& finalState, Genome& genome, std::ostream& out) {
    std::vector<MEI_event> insertion_events;
    // Loop over all breakpoint combinations to see if they refer to the same
    // insertion event.
    for (size_t i = 0; i < finalState.MEI_breakpoints.size(); i++) {
        for (size_t j = i+1; j < finalState.MEI_breakpoints.size(); j++) {
            
            MEI_breakpoint& bp1 = finalState.MEI_breakpoints.at(i);
            MEI_breakpoint& bp2 = finalState.MEI_breakpoints.at(j);
            
            if (refer_to_same_event(bp1, bp2)) {
                // Current 2 selected breakpoints refer to same event.
                
                MEI_event event = (bp1.cluster_strand == Plus)? MEI_event(bp1, bp2) : MEI_event(bp2, bp1);
                
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
    reportMEIevents(finalState, insertion_events, genome, out);
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
        std::string read_name = enrich_read_name(bam1_qname(alignment), alignment->core.flag & BAM_FREAD1);
        char strand = bam1_strand(alignment)? Minus : Plus;
        simple_read* read = new simple_read(read_name, alignment->core.tid, alignment->core.pos, strand, 
                                            get_sequence(bam1_seq(alignment), alignment->core.l_qseq));
        mei_data->discordant_reads.push_back(read);
    }
    return 0;
}


static int load_discordant_reads(MEI_data& mei_data, std::vector<bam_info>& bam_sources, const std::string& chr_name,
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
        int tid = bam_get_tid (header, chr_name.c_str());
        
        // Load discordant reads into mei_data.
        bam_fetch(fp, idx, tid, window.getStart(), window.getEnd(), &mei_data, fetch_disc_read_callback);
        bam_index_destroy(idx);
    }
    return 0;
}



// For debugging purposes:
// Class to divert some output stream to multiple destinations.  Author is
// a stackoverflow user named "Loki Astari".
//class ComposeStream: public std::ostream
//{
//    struct ComposeBuffer: public std::streambuf
//    {
//        void addBuffer(std::streambuf* buf)
//        {
//            bufs.push_back(buf);
//        }
//        virtual int overflow(int c)
//        {
//            std::for_each(bufs.begin(),bufs.end(),std::bind2nd(std::mem_fun(&std::streambuf::sputc),c));
//            return c;
//        }
//        
//    private:
//        std::vector<std::streambuf*>    bufs;
//        
//    };  
//    ComposeBuffer myBuffer;
//public: 
//    ComposeStream()
//    :std::ostream(NULL)
//    {
//        std::ostream::rdbuf(&myBuffer);
//    }   
//    void linkStream(std::ostream& out)
//    {
//        out.flush();
//        myBuffer.addBuffer(out.rdbuf());
//    }
//};



// This function is based on Pindel's main function.  Todo: integrate with pindel's main structure.
int searchMEImain(ControlState& current_state, Genome& genome, UserDefinedSettings* userSettings) {

    // Setup output destinations.
//    ComposeStream out;
//    out.linkStream(std::cout);
//    std::ofstream file_output("/Users/mkroon/temp/output_MEI_detection.txt");
//    out.linkStream(file_output);
    
    MEI_data mei_data;
    int result;
    
    do {
        
        const Chromosome* currentChromosome = g_genome.getNextChromosome();
        if (currentChromosome == NULL) {
            break;
        }
        
        g_maxPos = 0;
        CurrentChrMask.resize(currentChromosome->getCompSize());
        for (unsigned int i = 0; i < currentChromosome->getCompSize(); i++) {
            CurrentChrMask[i] = 'N';
        }
        
        LoopingSearchWindow currentWindow(userSettings->getRegion(), currentChromosome, WINDOW_SIZE); 
        // loop over one chromosome
        do {
            result = load_discordant_reads(mei_data, current_state.bams_to_parse, currentChromosome->getName(),
                                           currentWindow);
            if (result) {
                // something went wrong loading the reads, return error code.
                return result;
            }
            
            searchMEIBreakpoints(mei_data, current_state.bams_to_parse, currentChromosome);
            mei_data.discordant_reads.clear();
            currentWindow.next();
            
        } while (!currentWindow.finished());
    } while (true);
    
    searchMEI(mei_data, genome, std::cout);
//    out.flush();
    return 0;
}









