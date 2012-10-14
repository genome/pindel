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
#include "farend_searcher.h"
#include "search_MEI.h"
#include "search_MEI_util.h"

// Logging libraries.
#include "logdef.h"
#include "logstream.h"

// Maximal distance between two MEI breakpoints found on opposite strands to
// be assumed to refer to the same event.
const int MAX_MEI_BREAKPOINT_DISTANCE = 100;

// Maximal distance between reads for them to be assigned to the same cluster.
const int MAX_DISTANCE_CLUSTER_READS = 100;

// Size of the read buffer when reading split reads locally around potential MEI sites.
const int SPLIT_READ_BUFFER_SIZE = 50000;

// String used to characterize comment in text output.
const std::string COMMENT_PREFIX = "# ";

// Minimal size of cluster.
const int MIN_MEI_CLUSTER_SIZE = 3;

// Minimal number of reads needed to support a breakpoint.
const int MIN_MEI_BREAKPOINT_SUPPORT = 3;

// Minimal distance between mapped reads to consider them for detection of
// mobile element insertions.
const int MIN_MEI_MAP_DISTANCE = 8000;

// A contig, being a collection of overlapping reads paired with relative 
// indexes indicating their distance in mapping.
typedef std::vector<std::pair<simple_read, int> > contig;


// Returns true for a given read whether it is concordantly mapped together with its mate.
bool is_concordant(const bam1_t* read, unsigned int insert_size) {
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
    
    // Return true if insert size is as expected.  This definition of concordance is
    // equivalent to read-pair construction elsewhere in Pindel.
    return (unsigned int) abs(read->core.isize) < read->core.l_qseq + 2 * insert_size;
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
static int get_split_reads_for_cluster(std::vector<bam_info>& bam_sources, char cluster_strand, int outer_pos,
                                       const Chromosome* chromosome, std::vector<SPLIT_READ>& split_reads) {
    
    for (size_t i = 0; i < bam_sources.size(); i++) {
        bam_info source = bam_sources.at(i);
        
        // Determine bounds of the region where mates of split reads may reside 
        // related to the given cluster.
        int upper_bound;
        int lower_bound;
        if (cluster_strand == Plus) {
            lower_bound = outer_pos - source.InsertSize;
            upper_bound = outer_pos + 2 * source.InsertSize;
        } else {
            lower_bound = outer_pos - 2 * source.InsertSize;
            upper_bound = outer_pos + source.InsertSize;
        }
        SearchWindow search_window(chromosome, lower_bound, upper_bound);
        
        std::vector<SPLIT_READ> HangingReads; // Kai: those are reads without closed end mapped, or hanging reads.
        std::vector<REF_READ> Ref_Supporting_Reads; // Kai: those reads support the reference allele;
        // Read split reads in defined window.
        ReadBuffer read_buffer(SPLIT_READ_BUFFER_SIZE, split_reads, HangingReads, chromosome->getSeq());

        int read_result = ReadInBamReads_SR(source.BamFile.c_str(), chromosome->getName(), &chromosome->getSeq(), 
                                            split_reads, HangingReads, Ref_Supporting_Reads, source.InsertSize, 
                                            source.Tag, search_window, read_buffer);
        if (read_result > 0) {
            return read_result;
        }
        
        // Remove any split reads for which a far end can be found locally, these are
        // assumed to contribute to some local variants.
        // Todo: determine region that is searched for far end.
        std::vector<SearchWindow> window_holder;
        for (int i = (split_reads.size() - 1); i >= 0; i--) {
            SPLIT_READ current_SR = split_reads.at(i);
            
            int far_end_lower_bound = current_SR.getLastAbsLocCloseEnd();
            int far_end_upper_bound;
            if (cluster_strand == Plus) {
                far_end_upper_bound = far_end_lower_bound + MIN_MEI_MAP_DISTANCE;
            } else {
                far_end_upper_bound = far_end_lower_bound;
                far_end_lower_bound -= MIN_MEI_MAP_DISTANCE;
            }
            
            SearchWindow window(chromosome, far_end_lower_bound, far_end_upper_bound);
            window_holder.clear();
            window_holder.push_back(window);
            SearchFarEndAtPos(chromosome->getSeq(), current_SR, window_holder);
            
            // If far end is found in local window, delete it.
            if (current_SR.goodFarEndFound()) {
                split_reads.erase(split_reads.begin() + i);
            }
        }
    }
    return 0;
}


// Returns a breakpoint for a cluster of connected reads.  If no viable
// breakpoint can be found, it returns a breakpoint with position -1.
// Note: returned pointer must be deleted by caller.
static MEI_breakpoint* get_breakpoint(std::vector<simple_read*>* cluster, std::vector<bam_info>& bam_sources,
                                      int cluster_tid, char cluster_strand, const Chromosome* chromosome) {
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
            simple_read simple_split_read(read.Name, cluster_tid, candidate_bp, strand, read.sample_name, 
                                          whole_sequence, mapped_part, unmapped_part);
            simple_split_reads.push_back(simple_split_read);
        }
    }
    
    LOG_DEBUG(*logStream << time_log() << "Cluster: pos=" << outer_read_pos << " size=" << cluster->size() <<
              " split_reads=" << simple_split_reads.size() << std::endl);
    
    // Return the most frequently occurring candidate.
    std::vector<int> cb_copy = candidate_breakpoints;
    sort(cb_copy.begin(), cb_copy.end());
    int max_count = 0;
    int current_count = 0;
    int current = -1;
    int best = -1;
    for (size_t i = 0; i < cb_copy.size(); i++) {
        if (cb_copy.at(i) == current) {
            current_count += 1;
        } else {
            current = cb_copy.at(i);
            current_count = 1;
        }
        if (current_count > max_count) {
            max_count = current_count;
            best = current;
        }
    }
    
    if (max_count < MIN_MEI_BREAKPOINT_SUPPORT) {
        // No valid breakpoint was found.
        return new MEI_breakpoint(-1);
    }
    
    // Link associated discordant reads (all reads from cluster) and split reads.
    MEI_breakpoint* bp = new MEI_breakpoint(best);
    bp->cluster_strand = cluster_strand;
    bp->chromosome_name = chromosome->getName();
    // Add split reads supporting 'best' candidate breakpoint.
    for (size_t i = 0; i < candidate_breakpoints.size(); i++) {
        if (candidate_breakpoints.at(i) == best) {
            bp->associated_split_reads.push_back(simple_split_reads.at(i));
        }
    }
    for (size_t i = 0; i < cluster->size(); i++) {
        bp->associated_reads.push_back(*cluster->at(i));
    }
    
    return bp;
}


// See documentation in header file.
void searchMEIBreakpoints(MEI_data& currentState, std::vector<bam_info>& bam_sources, const Chromosome* chromosome) {
    LOG_DEBUG(*logStream << time_log() << "Start searching for breakpoints..." << std::endl);
    std::vector<std::vector<simple_read*>*> clusters = cluster_reads(currentState.discordant_reads);
    
    // Find breakpoints per cluster.
    int bp_count = 0;
    for (size_t i = 0; i < clusters.size(); i++) {
        // print cluster debug info
        std::vector<simple_read*>* cluster = clusters.at(i);
        
        if (cluster->size() < ((size_t) MIN_MEI_CLUSTER_SIZE)) {
            // Fluke cluster, skip it. (If there are very few reads in the cluster,
            // we likely won't find enough split reads supporting an insertion)
            continue;
        }
        
        // Find breakpoint for this cluster
        char cluster_strand = cluster->at(0)->strand;
        int cluster_tid = cluster->at(0)->tid;
        MEI_breakpoint* MEI_bp = get_breakpoint(cluster, bam_sources, cluster_tid, cluster_strand, chromosome);
        // Check breakpoint validity.
        if (MEI_bp->breakpoint_pos >= 0) {
            currentState.MEI_breakpoints.push_back(*MEI_bp);
            bp_count += 1;
        }
        delete MEI_bp;
    }
    LOG_DEBUG(*logStream << time_log() << "Found " << bp_count << " breakpoints for " << clusters.size() <<
              " clusters." << std::endl);
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



static void report_supporting_reads(std::vector<simple_read>& reads, std::map<int, std::string>& seq_name_dict,
                                    std::ostream& out) {
    out << "# All supporting sequences for this insertion:" << std::endl;
    sort(reads.begin(), reads.end(), comp_simple_read_pos);
    std::vector<simple_read>::iterator support_iter;
    for (support_iter = reads.begin(); support_iter != reads.end(); ++support_iter) {
        // Only write the part of the read that falls inside the inserted element.
        if ((*support_iter).is_split) {
            out << "?\t?\t?\t" << (*support_iter).name << "\t" << (*support_iter).sample_name << "\t" <<
                   (*support_iter).unmapped_sequence << std::endl;
        } else {
            out << seq_name_dict.at((*support_iter).tid) << "\t" << (*support_iter).pos << "\t" << 
                   (*support_iter).strand << "\t" << (*support_iter).name << "\t" << (*support_iter).sample_name << 
                   "\t" << (*support_iter).sequence << std::endl;
        }
    }
}



// Highlight (capitalize) characters in given string (reference) from/until position defined by
// breakpoint.  If highlight_until_bp is true, all characters up until breakpoint are capitalized
// otherwise all characters from breakpoint to the end of the string are capitalized.
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


// Get all reads (discordant/split) that support a breakpoint.  I.e. all reads that
// either overlap the breakpoint or are associated, but mapped somewhere else.
static void get_bp_supporting_reads(MEI_breakpoint& breakpoint, std::vector<MEI_breakpoint>& breakpoint_connections,
                                    std::vector<simple_read>& supporting_reads) {
    
    for (size_t i = 0; i < breakpoint_connections.size(); i++) {
        MEI_breakpoint connected_cluster_bp = breakpoint_connections.at(i);
        std::vector<std::pair<simple_read, simple_read> > shared_reads = get_shared_reads(breakpoint, 
                                                                                          connected_cluster_bp);
        for (size_t j = 0; j < shared_reads.size(); j++) {
            supporting_reads.push_back(shared_reads.at(j).second);
        }
    }
    
    for (size_t i = 0; i < breakpoint.associated_reads.size(); i++) {
        simple_read read = breakpoint.associated_reads.at(i);
        std::string read_name = base_read_name(read.name);
        bool already_added = false;
        for (size_t j = 0; j < supporting_reads.size(); j++) {
            if (read_name == base_read_name(supporting_reads.at(j).name)) {
                already_added = true;
                break;
            }
        }
        if (!already_added) {
            // Mate not available in breakpoint_connections, switch mate information
            // available in current read.  That way we can at least show the name and
            // position where the mate mapped.
            int temp_pos = read.mate_pos;
            read.mate_pos = read.pos;
            read.pos = temp_pos;
            int temp_tid = read.mate_tid;
            read.mate_tid = read.tid;
            read.tid = temp_tid;
            char temp_strand = read.mate_strand;
            read.mate_strand = read.strand;
            read.strand = temp_strand;
            // Sequence of mate is unknown.
            read.sequence = "?";
            supporting_reads.push_back(read);
        }
    }
    
    for (size_t j = 0; j < breakpoint.associated_split_reads.size(); j++) {
        supporting_reads.push_back(breakpoint.associated_split_reads.at(j));
    }
}


// Output split reads at breakpoint, position reads based on close end mapping.
static void report_split_read_support(Genome& genome, MEI_breakpoint& breakpoint, bool fiveprime_end,
                                      std::ostream& out) {
    
    // Sort on length of mapped part.
    sort(breakpoint.associated_split_reads.begin(), breakpoint.associated_split_reads.end(), comp_simple_read_mapsize);
    
    // Compute on-screen distances.
    int base;
    int end;
    simple_read first = *breakpoint.associated_split_reads.begin();
    simple_read last = *(breakpoint.associated_split_reads.end() - 1);
    if (fiveprime_end) {
        base = first.mapped_sequence.length();
        end = last.unmapped_sequence.length();
    } else {
        base = last.unmapped_sequence.length();
        end = first.mapped_sequence.length();
    }
    
    // Get reference sequence at breakpoint location.
    int offset = (fiveprime_end)? 1 : 0;
    std::string reference = get_fasta_subseq(genome, breakpoint.chromosome_name, 
                                             breakpoint.breakpoint_pos - base + offset, base + end);
    
    // Output local reference sequence.
    set_reference_highlight(reference, base, fiveprime_end);
    std::string REFERENCE_PREFIX = "Reference: ";
    out << COMMENT_PREFIX << REFERENCE_PREFIX << reference << std::endl;

    // Output split reads aligned to reference.
    std::vector<simple_read>::iterator read_iter;
    for (read_iter = breakpoint.associated_split_reads.begin(); read_iter != breakpoint.associated_split_reads.end();
         ++read_iter) {
        simple_read read = (*read_iter);
        int indent = REFERENCE_PREFIX.length();
        indent += (fiveprime_end)? base - read.mapped_sequence.length() : base - read.unmapped_sequence.length();
        out << COMMENT_PREFIX << get_whitespace(indent);
        if (fiveprime_end) {
            out << read.mapped_sequence << read.unmapped_sequence;
        } else {
            out << read.unmapped_sequence << read.mapped_sequence;
        }
        out << " (name: " << read.name << " sample: " << read.sample_name << ") " << std::endl;
    }
}


// Report MEI events.
static void reportMEIevents(MEI_data& mei_data, std::vector<MEI_event>& insertion_events, Genome& genome,
                            std::map<int, std::string>& seq_name_dict, std::ostream& out) {
    int MEI_counter = 0;
    for (size_t i = 0; i < insertion_events.size(); i++) {
        MEI_event event = insertion_events.at(i);
        
        LOG_DEBUG(*logStream << time_log() << 
                  "reporting MEI: #fwd.disc.: " << event.fwd_cluster_bp.associated_reads.size() <<
                  ", #fwd.split: " << event.fwd_cluster_bp.associated_split_reads.size() <<
                  ", #rev.disc.: " << event.rev_cluster_bp.associated_reads.size() <<
                  ", #rev.split: " << event.rev_cluster_bp.associated_split_reads.size() << std::endl);
        
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
        
        // Todo: Get lower bound on size of element.
        int summed_contig_span = 0;
        
        out << "####################################################################################################" << std::endl;
        // Print machine summary line.
        out << MEI_counter << "\t" << "MEI" << "\t" << event.fwd_cluster_bp.chromosome_name << "\t" << 
               event.fwd_cluster_bp.breakpoint_pos << "\t" << event.rev_cluster_bp.breakpoint_pos << "\t" <<
               summed_contig_span << "\t" << all_reads.size() << "\t" << (fwd_reads.size() - fwd_split_count) << "\t" <<
               fwd_split_count << "\t" << (rev_reads.size() - rev_split_count) << "\t" << rev_split_count << std::endl;
        // Print human-readable summary lines.
        out << COMMENT_PREFIX << "Mobile element insertion (MEI) found on chromosome '" << 
               event.fwd_cluster_bp.chromosome_name << "', breakpoint at " << event.fwd_cluster_bp.breakpoint_pos << 
               " (estimated from + strand), " << event.rev_cluster_bp.breakpoint_pos << " (estimated from - strand)," <<
               " the inserted element is at least " << summed_contig_span << " bp long." << std::endl;
        out << COMMENT_PREFIX << "Found " << all_reads.size() << " supporting reads, of which " <<
               (fwd_reads.size() - fwd_split_count) << " discordant reads and " << fwd_split_count << 
               " split reads at 5' end, " << (rev_reads.size() - rev_split_count) << " discordant reads and " <<
               rev_split_count << " split reads at 3' end." << std::endl;
        
        // Print support for breakpoint at 5' end.
        out << COMMENT_PREFIX << "Supporting reads for insertion location (5' end):" << std::endl;
        report_split_read_support(genome, event.fwd_cluster_bp, true, out);
        
        // Print support for breakpoint at 3' end.
        out << COMMENT_PREFIX << "Supporting reads for insertion location (3' end):" << std::endl;
        report_split_read_support(genome, event.rev_cluster_bp, false, out);
        
        // Print all supporting reads and read fragments for the inserted element.
        report_supporting_reads(all_reads, seq_name_dict, out);
        MEI_counter++;
    }
}


// See documentation in header file.
void searchMEI(MEI_data& finalState, Genome& genome, std::map<int, std::string>& seq_name_dict, std::ostream& out) {
    LOG_DEBUG(*logStream << time_log() << "Start calling MEI events from found breakpoints..." << std::endl);
    std::vector<MEI_event> insertion_events;
    
    LOG_DEBUG(*logStream << time_log() << "Examining " << finalState.MEI_breakpoints.size() << 
              " breakpoints in total." << std::endl);
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
                    if (k == i || k == j) {
                        // Avoid connecting a cluster to one of the same event.
                        continue;
                    }
                    MEI_breakpoint& bp_connect = finalState.MEI_breakpoints.at(k);
                    if (get_shared_reads(event.fwd_cluster_bp, bp_connect).size() > 0) {
                        event.fwd_cluster_connections.push_back(bp_connect);
                    } else if (get_shared_reads(event.rev_cluster_bp, bp_connect).size() > 0) {
                        event.rev_cluster_connections.push_back(bp_connect);
                    }
                }
                insertion_events.push_back(event);
            }
        }
    }
    LOG_DEBUG(*logStream << time_log() << "Found " << insertion_events.size() << " MEI events." << std::endl);
    
    // Report MEI events constructed from pairing breakpoints.
    reportMEIevents(finalState, insertion_events, genome, seq_name_dict, out);
}


//#####################################################################################
//###                                                                               ###
//###      Code below should be integrated with pindel's existing structure.        ###
//###                                                                               ###
//#####################################################################################

#include "control_state.h"

// Return code for not being able to open a BAM file.
const int ERROR_OPENING_ALIGNMENT_FILE = 100;

static int fetch_disc_read_callback(const bam1_t* alignment, void* data) {
    MEI_data* mei_data = static_cast<MEI_data*>(data);
    if (!(alignment->core.flag & BAM_FUNMAP || alignment->core.flag & BAM_FMUNMAP) && // Both ends are mapped.
        !is_concordant(alignment, mei_data->current_insert_size) &&                   // Ends map discordantly.
        // Extra check for (very) large mapping distance.  This is done beside the check for read
        // discordance to speed up computation by ignoring signals from small structural variants.
        (alignment->core.tid != alignment->core.mtid || 
         abs(alignment->core.pos - alignment->core.mpos) > MIN_MEI_MAP_DISTANCE)) {

        // Save alignment as simple_read object.
        std::string read_name = enrich_read_name(bam1_qname(alignment), alignment->core.flag & BAM_FREAD1);
        char strand = bam1_strand(alignment)? Minus : Plus;
        char mate_strand = bam1_mstrand(alignment)? Minus : Plus;
        std::string sample_name = get_sample_name(alignment, mei_data->sample_names);    
        
        simple_read* read = new simple_read(read_name, alignment->core.tid, alignment->core.pos, strand, sample_name,
                                            get_sequence(bam1_seq(alignment), alignment->core.l_qseq), 
                                            alignment->core.mtid, alignment->core.mpos, mate_strand);
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
        
        LOG_DEBUG(*logStream << time_log() << "Loading discordant reads from " << source.BamFile << std::endl);
        
        // Setup link to bamfile, its index and header.
        bamFile fp = bam_open(source.BamFile.c_str(), "r");
        bam_index_t *idx = bam_index_load(source.BamFile.c_str());
        bam_header_t *header = bam_header_read(fp);
        bam_init_header_hash(header);
        int tid = bam_get_tid(header, chr_name.c_str());
        
        if (tid < 0) {
            LOG_WARN(*logStream << time_log() << "Could not find sequence in alignment file: '" << chr_name <<
                     "'" << std::endl);
            continue;
        }
        
        mei_data.sample_names = get_sample_dictionary(header);
        
        // Save insert size of current bamfile in data object provided for callback function.
        // Note: the insert size should ideally be separate from the MEI_data object, tried to do
        // this using a std::pair object, which did not work.  Suggestions are welcome here.
        mei_data.current_insert_size = source.InsertSize;
        mei_data.current_chr_name = chr_name;
        
        // Load discordant reads into mei_data.
        bam_fetch(fp, idx, tid, window.getStart(), window.getEnd(), &mei_data, fetch_disc_read_callback);
        bam_index_destroy(idx);
    }
    return 0;
}


// Construct a map linking tid's to sequence names (chromosome names).  This code
// assumes all input bam files have identical set of sequence (identical both in
// name and order).
std::map<int, std::string> get_sequence_name_dictionary(ControlState& state) {
    std::map<int, std::string> dict;
    std::vector<bam_info>::iterator bam_info_iter;
    for (bam_info_iter = state.bams_to_parse.begin(); bam_info_iter != state.bams_to_parse.end(); ++bam_info_iter) {
        bamFile fp = bam_open((*bam_info_iter).BamFile.c_str(), "r");
        bam_header_t *header = bam_header_read(fp);
        bam_init_header_hash(header);
        for (int tid = 0; tid < header->n_targets; tid++) {
            dict.insert(std::make_pair(tid, header->target_name[tid]));
        }
        break; // Skip other bam files, the sequences should be identical.
    }
    return dict;
}


// This function is based on Pindel's main function.  Todo: integrate with pindel's main structure.
int searchMEImain(ControlState& current_state, Genome& genome, UserDefinedSettings* userSettings) {

    std::ofstream file_output(userSettings->getMEIOutputFilename().c_str());
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
            LOG_DEBUG(*logStream << time_log() << "MEI detection current window: " <<
                      currentWindow.getChromosomeName() << " " << currentWindow.getStart() << "--" <<
                      currentWindow.getEnd() << std::endl);
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
    
    std::map<int, std::string> seq_name_dictionary = get_sequence_name_dictionary(current_state);
    
    searchMEI(mei_data, genome, seq_name_dictionary, file_output);
    file_output.close();
    return 0;
}









