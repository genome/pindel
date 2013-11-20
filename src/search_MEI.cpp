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

// Size of the read buffer when reading split reads locally around potential MEI sites.
const int SPLIT_READ_BUFFER_SIZE = 50000;

// String used to characterize comment in text output.
const std::string COMMENT_PREFIX = "# ";

// Minimal fraction of collection of split reads that must have identical bases at each position
// in the sequence to construct a valid consensus sequence.  (e.g. 0.8 means at least 80% of the
// sequences must have the same base at an arbitrary position).
const float MIN_FRACTION_CONSENSUS = 0.8;

// Minimal length of a consensus string (if shorther than this, the strings may be similar
// by coincidence).
const int MIN_CONSENSUS_LENGTH = 15;

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
static void cluster_reads(std::vector<simple_read*>& reads, int insert_size, 
                          std::vector<std::vector<simple_read*> >& clusters, UserDefinedSettings* userSettings) {
    if (reads.size() == 0) {
        return;
    }
    
    // Sort reads first by mapping strand, then on mapping location.
    sort(reads.begin(), reads.end(), comp_simple_read);
    
    // Initialize first cluster.
    std::vector<simple_read*> current_cluster;
    simple_read* read_last_analyzed = reads.at(0);
    simple_read* first_read_in_cluster = read_last_analyzed;
    current_cluster.push_back(read_last_analyzed);
    
    for (size_t i = 1; i < reads.size(); i++) {
        simple_read* current_read = reads.at(i);
        
        if (// Max distance between reads in cluster.
            (current_read->pos - read_last_analyzed->pos) <= userSettings->MAX_DISTANCE_CLUSTER_READS &&
            // Max spanning length of cluster on reference.
            (unsigned) (current_read->pos - first_read_in_cluster->pos) <= 
                (insert_size - first_read_in_cluster->sequence.length()) &&
            // Reads of same cluster must be on same strand.
            read_last_analyzed->strand == current_read->strand) {
            
            current_cluster.push_back(current_read);
        } else {
            // Close current cluster and initialize a new.
            clusters.push_back(current_cluster);
            current_cluster.clear();
            current_cluster.push_back(current_read);
            first_read_in_cluster = current_read;
        }
        read_last_analyzed = current_read;
    }
    
    // Save remaining cluster.
    if (current_cluster.size() > 0) {
        clusters.push_back(current_cluster);
    }
}

int all_sr = 0;
int far_end_found = 0;
int far_end_found_elsewhere = 0;

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
        ReadBuffer read_buffer(SPLIT_READ_BUFFER_SIZE, split_reads, HangingReads, chromosome->getSeq(), *chromosome);

        ReadInBamReads_SR(source.BamFile.c_str(), chromosome->getName(), &chromosome->getSeq(), split_reads, 
                          HangingReads, Ref_Supporting_Reads, source.InsertSize, source.Tag, search_window, 
                          read_buffer, false);
        }
    return 0;
}


// Return a consensus sequence for the mapping parts of reads, if no consensus found return 
// empty string.  strand is the strand to which the reads are mapped -- defining whether the
// sequences should be aligned by left or right side.
static std::string get_consensus_unmapped(std::vector<simple_read>& reads, char strand) {
    
    if (reads.size() == 0) {
        return "";
    }
    
    sort(reads.begin(), reads.end(), comp_simple_read_unmapped_seqsize);
    int max_len = reads.at(0).unmapped_sequence.length();
    
    std::stringstream consensus_ss;
    int index; // Position in the consensus read.
    std::vector<simple_read>::iterator read_iter;
    for (int i = 0; i < max_len; i++) {
        std::map<char, int> char_counts;
        std::map<char, int>::iterator char_iter;
        int read_count = 0; // Nr of reads contributing to current position in consensus read.
        
        // Count character occurrences for current index.
        for (read_iter = reads.begin(); read_iter != reads.end(); ++read_iter) {
            // Get index (either i or i'th position from end of mapped seq)
            index = (strand == Minus)? i : (*read_iter).unmapped_sequence.length() - 1 - i;
            if (index < 0 || index >= (int) (*read_iter).unmapped_sequence.length()) {
                continue;
            }
            read_count++;
            
            char_iter = char_counts.find((*read_iter).unmapped_sequence.at(index));
            if (char_iter != char_counts.end()) {
                char_counts[(*char_iter).first] = (*char_iter).second + 1;
            } else {
                char_counts[(*read_iter).unmapped_sequence.at(index)] = 1;
            }
        }
        
        // Get consensus char from counts.
        char consensus_char = '?';
        int max_count = 0;
        for (char_iter = char_counts.begin(); char_iter != char_counts.end(); ++char_iter) {
            if ((*char_iter).second > max_count) {
                max_count = (*char_iter).second;
                consensus_char = (*char_iter).first;
            }
        }
        
        if (max_count >= MIN_FRACTION_CONSENSUS * read_count) {
            consensus_ss << consensus_char;
        } else {
            // Too much deviating characters at this position.
            break;
        }
    }

    // If needed, reverse the consensus string before returning.
    std::string output = consensus_ss.str();
    if (output.length() < (unsigned) MIN_CONSENSUS_LENGTH) {
        return "";
    }
    if (strand == Plus) {
        std::reverse(output.begin(), output.end());
    }
    return output;
}



// Returns a breakpoint for a cluster of connected reads.  If no viable
// breakpoint can be found, it returns a breakpoint with position -1.
// Note: returned pointer must be deleted by caller.
static void get_breakpoints(std::vector<simple_read*>& cluster, std::vector<bam_info>& bam_sources, int insert_size,
                            int cluster_tid, char cluster_strand, const Chromosome* chromosome, 
                            std::map<std::string, std::string>& sample_dict, std::vector<MEI_breakpoint>& breakpoints,
                            UserDefinedSettings* userSettings) {
    std::vector<SPLIT_READ> split_reads;
    int outer_read_pos = (cluster_strand == Minus)? cluster.at(cluster.size()-1)->pos : cluster.at(0)->pos;
//    int inner_read_pos = (cluster_strand == Minus)? cluster.at(0)->pos : cluster.at(cluster.size()-1)->pos;
    get_split_reads_for_cluster(bam_sources, cluster_strand, outer_read_pos, chromosome, split_reads);
    
    // Search for split reads with a mate close to the outer read of the
    // cluster.  Store candidate breakpoints.
    // Todo: speedup by exploiting the fact that both clusters and split reads are sorted
    // by mapping location.
    std::map<int, std::vector<simple_read> > bio_candidate_breakpoints;
    for (size_t i = 0; i < split_reads.size(); i++) {
        SPLIT_READ read = split_reads.at(i);
        
        char anchor_strand = read.MatchedD;
        if (cluster_strand != anchor_strand) {
            continue;
        }
        
        unsigned int comp_candidate_bp = read.getLastAbsLocCloseEnd();
        unsigned int bio_candidate_bp = get_bio_chr_index(comp_candidate_bp);
        
        if (bio_candidate_breakpoints.find(bio_candidate_bp) == bio_candidate_breakpoints.end()) {
            // New candidate, look ahead to check whether there are enough supporting split reads.
            int SR_support = 1;
            for (size_t j = i + 1; j < split_reads.size(); j++) {
                if (split_reads.at(j).getLastAbsLocCloseEnd() == comp_candidate_bp && 
                    split_reads.at(j).MatchedD == cluster_strand) {
                    SR_support++;
                }
            }
            if (SR_support < userSettings->MIN_DD_BREAKPOINT_SUPPORT) {
                // Not enough support, skip it.
                continue;
            } else {
                std::vector<simple_read> new_bp_split_reads;
                bio_candidate_breakpoints.insert(std::make_pair(bio_candidate_bp, new_bp_split_reads));
            }
        }
        
        // Store the unmatched sequence as it should be matched on the opposite strand of
        // the mapped mate.
        std::string whole_sequence;
        std::string mapped_part;
        std::string unmapped_part;
        if (anchor_strand == Plus) {
            whole_sequence = read.getUnmatchedSeqRev();
            mapped_part = whole_sequence.substr(0, read.CloseEndLength);
            unmapped_part = whole_sequence.substr(read.CloseEndLength, whole_sequence.length());
        } else {
            whole_sequence = read.getUnmatchedSeq();
            mapped_part = whole_sequence.substr(whole_sequence.length() - read.CloseEndLength, 
                                                whole_sequence.length());
            unmapped_part = whole_sequence.substr(0, whole_sequence.length() - read.CloseEndLength);
        }

        std::string sample_name;
        get_sample_name(read.read_group, sample_dict, sample_name);
        simple_read simple_split_read(read.Name, -1, -1, '?', sample_name, whole_sequence, mapped_part, 
                                      unmapped_part);
        (*bio_candidate_breakpoints.find(bio_candidate_bp)).second.push_back(simple_split_read);
    }
    
  
    char SR_mapping_strand = (cluster_strand == Plus)? Minus : Plus;
    
    // Remove any split reads for which a far end can be found locally, these are
    // assumed to contribute to some local variants.
    // Todo: determine region that is searched for far end.
    std::map<int, std::vector<simple_read> >::iterator map_iter;
    for (map_iter = bio_candidate_breakpoints.begin(); map_iter != bio_candidate_breakpoints.end(); ++map_iter) {
        
        std::string mapped_consensus = get_consensus_unmapped((*map_iter).second, SR_mapping_strand);
        std::vector<simple_read> sreads = (*map_iter).second;
        if (mapped_consensus.length() == 0) {
            LOG_DEBUG(*logStream << time_log() << "Consensus building failed for split read mapping ends (" << 
                      map_iter->second.size() << " reads @ " << map_iter->first << ")" << std::endl);
            continue;
        }
        int bio_bp = (*map_iter).first;
                
        // If far end consensus is not found in local window, store breakpoint.
        size_t FE_window_start = std::max(0, get_comp_chr_index(bio_bp) - userSettings->MIN_DD_MAP_DISTANCE);
        size_t FE_window_size = std::min(chromosome->getCompSize() - (unsigned) FE_window_start, 
                                         2 * (unsigned) userSettings->MIN_DD_MAP_DISTANCE);
        if (!contains_subseq_any_strand(mapped_consensus, chromosome->getSeq().substr(FE_window_start, 
                FE_window_size), MIN_CONSENSUS_LENGTH)) {
            MEI_breakpoint bp(cluster_tid, bio_bp, cluster_strand);
            bp.associated_split_reads = (*map_iter).second;

            // Link associated discordant reads (all reads from cluster) and split reads.
            std::vector<simple_read*>::iterator read_iter;
            for (read_iter = cluster.begin(); read_iter != cluster.end(); ++read_iter) {
                bp.associated_reads.push_back(*(*read_iter));
            }
            breakpoints.push_back(bp);
        }
    }
}



// Generate an estimated breakpoint from mapping positions of discordant reads.  New breakpoint is
// added to MEI_bps.  This function assumes that objects in cluster are already ordered by position.
// Another assumption is that there are at least 2 reads in the cluster (i.e. 
// user setting MIN_DD_CLUSTER_SIZE > 1).
void get_breakpoint_estimation(std::vector<simple_read*>& cluster, int insert_size, int cluster_tid, 
                               char cluster_strand, std::vector<MEI_breakpoint>& MEI_bps) {
    // Compute average distance between read mapping positions.
    float dist_mean = 0;
    for (unsigned int i = 0; i < (cluster.size() - 1); i++) {
        dist_mean += (1.0 / (i+1)) * ((cluster.at(i+1)->pos - cluster.at(i)->pos) - dist_mean);
    }

    // Compute mean read lenght in cluster.
    float read_len_mean = 0;
    for (unsigned int i = 0; i < cluster.size(); i++) {
        read_len_mean += (1.0 / (i+1)) * (cluster.at(i)->sequence.length() - read_len_mean);
    }
    
    int outer_pos_high = cluster.at(cluster.size()-1)->pos + cluster.at(cluster.size()-1)->sequence.length();
    int outer_pos_low = cluster.at(0)->pos;
    
    int estimation = (cluster_strand == Plus)? outer_pos_high + dist_mean : outer_pos_low - dist_mean;

    MEI_breakpoint new_bp(cluster_tid, estimation, cluster_strand);
    // Link associated discordant reads (all reads from cluster) and split reads.
    std::vector<simple_read*>::iterator read_iter;
    for (read_iter = cluster.begin(); read_iter != cluster.end(); ++read_iter) {
        new_bp.associated_reads.push_back(*(*read_iter));
    }
    MEI_bps.push_back(new_bp);
}



// See documentation in header file.
void searchMEIBreakpoints(MEI_data& currentState, std::vector<bam_info>& bam_sources, const Chromosome* chromosome,
                          UserDefinedSettings* userSettings) {
    LOG_DEBUG(*logStream << time_log() << "Start searching for breakpoints..." << std::endl);
    std::vector<std::vector<simple_read*> > clusters;
    cluster_reads(currentState.discordant_reads, currentState.current_insert_size, clusters, userSettings);
    
    // Find breakpoints per cluster.
    int bp_count = 0;
    for (size_t i = 0; i < clusters.size(); i++) {
        // print cluster debug info
        std::vector<simple_read*> cluster = clusters.at(i);
   
        if (cluster.size() < ((size_t) userSettings->MIN_DD_CLUSTER_SIZE)) {
            // Fluke cluster, skip it. (If there are very few reads in the cluster,
            // we likely won't find enough split reads supporting an insertion)
            continue;
        }


        
        // Find breakpoint for this cluster
        char cluster_strand = cluster.at(0)->strand;
        int cluster_tid = cluster.at(0)->tid;
        std::vector<MEI_breakpoint> MEI_bps;
        std::vector<MEI_breakpoint>::iterator MEI_iter;
        get_breakpoints(cluster, bam_sources, currentState.current_insert_size, cluster_tid, cluster_strand, chromosome,
                        currentState.sample_names, MEI_bps, userSettings);
        if (MEI_bps.size() > 1) {
            // More than one breakpoints found for current cluster.  Select only the one with the
            // most split reads supporting it.
            size_t best_support = 0;
            MEI_breakpoint best_bp;
            for (MEI_iter = MEI_bps.begin(); MEI_iter != MEI_bps.end(); ++MEI_iter) {
                if (MEI_iter->associated_split_reads.size() > best_support) {
                    best_bp = *MEI_iter;
                    best_support = MEI_iter->associated_split_reads.size();
                }
            }
            MEI_bps.clear();
            MEI_bps.push_back(best_bp);
        } else if (MEI_bps.size() == 0) {
            // No breakpoint found with split read support.  Estimate breakpoint from cluster reads.
            get_breakpoint_estimation(cluster, currentState.current_insert_size, cluster_tid, cluster_strand, MEI_bps);
        }
        // Check breakpoint validity.
        for (MEI_iter = MEI_bps.begin(); MEI_iter != MEI_bps.end(); ++MEI_iter) {
            currentState.MEI_breakpoints.push_back(*MEI_iter);
            bp_count += 1;
            LOG_INFO(*logStream << "Found potential DD breakpoint: " << (*MEI_iter).breakpoint_tid << ", " <<
                     (*MEI_iter).breakpoint_pos << ", " << (*MEI_iter).cluster_strand << ", " <<
                     (*MEI_iter).associated_reads.size() << ", " <<
                     (*MEI_iter).associated_split_reads.size() << std::endl);
        }
    }
    LOG_DEBUG(*logStream << time_log() << "Found " << bp_count << " breakpoints for " << clusters.size() <<
              " clusters." << std::endl);
}


// Mobile element insertion event, constructed from two supporting read clusters.
class MEI_event {
public:
    MEI_event(){};
    MEI_event(MEI_breakpoint bp1, MEI_breakpoint bp2) : fwd_cluster_bp(bp1), rev_cluster_bp(bp2) {};
    
    MEI_breakpoint fwd_cluster_bp; // bp estimated from cluster on forward strand.
    MEI_breakpoint rev_cluster_bp; // bp estimated from cluster on reversed strand.

    // Breakpoints from clusters that have shared reads with the two supporting
    // clusters (provided by fwd_cluster_bp and rev_cluster_bp).
    std::vector<MEI_breakpoint> fwd_cluster_connections;
    std::vector<MEI_breakpoint> rev_cluster_connections;
    
    // Reads mapping inside event.
    std::vector<simple_read> fwd_mapping_reads;
    std::vector<simple_read> rev_mapping_reads;
};



static void report_supporting_reads(std::vector<simple_read>& reads, std::map<int, std::string>& seq_name_dict,
                                    std::ostream& out) {
    out << "# All supporting sequences for this insertion (i.e. sequences that map inside the inserted element):" << 
            std::endl;
    sort(reads.begin(), reads.end(), comp_simple_read_pos);
    std::vector<simple_read>::iterator support_iter;
    for (support_iter = reads.begin(); support_iter != reads.end(); ++support_iter) {
        // Only write the part of the read that falls inside the inserted element.
        if ((*support_iter).is_split) {
            out << "?\t?\t?\t" << (*support_iter).name << "\t" << (*support_iter).sample_name << "\t" <<
                   (*support_iter).evidence_strand << "\t" << (*support_iter).unmapped_sequence << std::endl;
        } else {
            out << seq_name_dict.at((*support_iter).tid) << "\t" << (*support_iter).pos << "\t" << 
                   (*support_iter).strand << "\t" << (*support_iter).name << "\t" << (*support_iter).sample_name << 
                   "\t" << (*support_iter).evidence_strand << "\t" << (*support_iter).sequence << std::endl;
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


// Add evidence strand information to reads of given event.
static void set_evidence_strands(MEI_event& event) {
    for (size_t i=0; i < event.fwd_cluster_bp.associated_reads.size(); i++) {
        event.fwd_cluster_bp.associated_reads.at(i).evidence_strand = Plus;
    }
    for (size_t i=0; i < event.fwd_cluster_bp.associated_split_reads.size(); i++) {
        event.fwd_cluster_bp.associated_split_reads.at(i).evidence_strand = Plus;
    }
    for (size_t i=0; i < event.fwd_mapping_reads.size(); i++) {
        event.fwd_mapping_reads.at(i).evidence_strand = Plus;
    }
    for (size_t i=0; i < event.rev_cluster_bp.associated_reads.size(); i++) {
        event.rev_cluster_bp.associated_reads.at(i).evidence_strand = Minus;
    }
    for (size_t i=0; i < event.rev_cluster_bp.associated_split_reads.size(); i++) {
        event.rev_cluster_bp.associated_split_reads.at(i).evidence_strand = Minus;
    }
    for (size_t i=0; i < event.rev_mapping_reads.size(); i++) {
        event.rev_mapping_reads.at(i).evidence_strand = Minus;
    }
}


// Get all reads (discordant/split) that support a breakpoint.  I.e. all reads that
// either overlap the breakpoint or are associated, but mapped somewhere else.
static void get_event_supporting_reads(MEI_event& event, std::vector<simple_read>& supporting_reads) {
    
    // Add all split reads and discordant reads if available.
    supporting_reads.insert(supporting_reads.end(), event.fwd_mapping_reads.begin(), event.fwd_mapping_reads.end());
    supporting_reads.insert(supporting_reads.end(), event.fwd_cluster_bp.associated_split_reads.begin(),
                            event.fwd_cluster_bp.associated_split_reads.end());
    supporting_reads.insert(supporting_reads.end(), event.rev_mapping_reads.begin(), event.rev_mapping_reads.end());
    supporting_reads.insert(supporting_reads.end(), event.rev_cluster_bp.associated_split_reads.begin(),
                            event.rev_cluster_bp.associated_split_reads.end());
    
    // For non-available discordant reads, add known info as pseudo read.
    std::vector<simple_read> tmp_all_associated_reads = event.fwd_cluster_bp.associated_reads;
    tmp_all_associated_reads.insert(tmp_all_associated_reads.end(), event.rev_cluster_bp.associated_reads.begin(),
                                    event.rev_cluster_bp.associated_reads.end());
    for (size_t i = 0; i < tmp_all_associated_reads.size(); i++) {
        simple_read read = tmp_all_associated_reads.at(i);
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
}


// Output split reads at breakpoint, position reads based on close end mapping.
static void report_split_read_support(Genome& genome, MEI_breakpoint& breakpoint, bool fiveprime_end,
                                      std::map<int, std::string>& seq_name_dict, std::ostream& out) {
    if (breakpoint.associated_split_reads.size() == 0) {
        // No split reads to report.
        return;
    }
    
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
    std::string reference = get_fasta_subseq(genome, seq_name_dict.at(breakpoint.breakpoint_tid), 
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
static void reportMEIevent(MEI_data& mei_data, MEI_event& event, int MEI_count, Genome& genome,
                           std::map<int, std::string>& seq_name_dict, std::ostream& out) {
    
    // Set evidence strand for event's reads (they'll be reported).
    set_evidence_strands(event);
    
    // List all read info that needs to be reported.
    std::vector<simple_read> all_reads;
    get_event_supporting_reads(event, all_reads);
    
    LOG_DEBUG(*logStream << time_log()
                         << "reporting DD: #fwd.disc.: " << event.fwd_cluster_bp.associated_reads.size()
                         << ", #fwd.split: " << event.fwd_cluster_bp.associated_split_reads.size()
                         << ", #rev.disc.: " << event.rev_cluster_bp.associated_reads.size()
                         << ", #rev.split: " << event.rev_cluster_bp.associated_split_reads.size() << std::endl);
    
    size_t all_read_count = event.fwd_cluster_bp.associated_reads.size() +
                            event.fwd_cluster_bp.associated_split_reads.size() +
                            event.rev_cluster_bp.associated_reads.size() +
                            event.rev_cluster_bp.associated_split_reads.size();
    
    out << "####################################################################################################" << std::endl;
    // Print machine summary line.
    out << MEI_count << "\t" << "DD" << "\t";
    out << seq_name_dict.at(event.fwd_cluster_bp.breakpoint_tid) << "\t" <<
           event.fwd_cluster_bp.breakpoint_pos << "\t" << event.rev_cluster_bp.breakpoint_pos;
    
    out << "\t" << all_read_count << "\t" << event.fwd_cluster_bp.associated_reads.size() << "\t"
        << event.fwd_cluster_bp.associated_split_reads.size();
    out << "\t" << event.rev_cluster_bp.associated_reads.size() << "\t"
        << event.rev_cluster_bp.associated_split_reads.size() << std::endl;
    
    // Print human-readable summary lines.
    out << COMMENT_PREFIX << "Dispersed Duplication insertion (DD) found on chromosome '" <<
           seq_name_dict.at(event.fwd_cluster_bp.breakpoint_tid) << "', breakpoint at " <<
           event.fwd_cluster_bp.breakpoint_pos << " (estimated from + strand), " <<
           event.rev_cluster_bp.breakpoint_pos << " (estimated from - strand)" << std::endl;
    out << COMMENT_PREFIX << "Found " << all_read_count << " supporting reads, of which " <<
           event.fwd_cluster_bp.associated_reads.size() << " discordant reads and " <<
           event.fwd_cluster_bp.associated_split_reads.size() << " split reads at 5' end, " <<
           event.rev_cluster_bp.associated_reads.size() << " discordant reads and " <<
           event.rev_cluster_bp.associated_split_reads.size() << " split reads at 3' end." << std::endl;

    // Print support for breakpoint at 5' end.
    out << COMMENT_PREFIX << "Supporting reads for insertion location (5' end):" << std::endl;
    report_split_read_support(genome, event.fwd_cluster_bp, true, seq_name_dict, out);
    // Print support for breakpoint at 3' end.
    out << COMMENT_PREFIX << "Supporting reads for insertion location (3' end):" << std::endl;
    report_split_read_support(genome, event.rev_cluster_bp, false, seq_name_dict, out);
        
    // Print all supporting reads and read fragments for the inserted element.
    report_supporting_reads(all_reads, seq_name_dict, out);
}


bool comp_breakpoint_pos(const MEI_breakpoint& bp1, const MEI_breakpoint& bp2) {
    return bp1.breakpoint_tid < bp2.breakpoint_tid || 
           (bp1.breakpoint_tid == bp2.breakpoint_tid && bp1.breakpoint_pos < bp2.breakpoint_pos);
}



#include "control_state.h"

// Return code for not being able to open a BAM file.
const int ERROR_OPENING_ALIGNMENT_FILE = 100;

static int fetch_disc_read_callback(const bam1_t* alignment, void* data) {
    //    MEI_data* mei_data = static_cast<MEI_data*>(data);
    std::pair<MEI_data*, UserDefinedSettings*>* env = static_cast<std::pair<MEI_data*, UserDefinedSettings*>*>(data);
    MEI_data* mei_data = env->first;
    UserDefinedSettings* userSettings = env->second;
    if (!(alignment->core.flag & BAM_FUNMAP || alignment->core.flag & BAM_FMUNMAP) && // Both ends are mapped.
        !is_concordant(alignment, mei_data->current_insert_size) &&                   // Ends map discordantly.
        // Extra check for (very) large mapping distance.  This is done beside the check for read
        // discordance to speed up computation by ignoring signals from small structural variants.
        (alignment->core.tid != alignment->core.mtid ||
         abs(alignment->core.pos - alignment->core.mpos) > userSettings->MIN_DD_MAP_DISTANCE)) {
            
            // Save alignment as simple_read object.
            std::string read_name = enrich_read_name(bam1_qname(alignment), alignment->core.flag & BAM_FREAD1);
            char strand = bam1_strand(alignment)? Minus : Plus;
            char mate_strand = bam1_mstrand(alignment)? Minus : Plus;
            std::string read_group;
            get_read_group(alignment, read_group);
            std::string sample_name;
            get_sample_name(read_group, mei_data->sample_names, sample_name);
            
            simple_read* read = new simple_read(read_name, alignment->core.tid, alignment->core.pos, strand, sample_name,
                                                get_sequence(bam1_seq(alignment), alignment->core.l_qseq),
                                                alignment->core.mtid, alignment->core.mpos, mate_strand);
            mei_data->discordant_reads.push_back(read);
        }
    return 0;
}


static int load_discordant_reads(MEI_data& mei_data, std::vector<bam_info>& bam_sources, const std::string& chr_name,
                                 const SearchWindow& window, UserDefinedSettings* userSettings) {
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
        
        // Set up environment variable for callback function.
        std::pair<MEI_data*, UserDefinedSettings*> env = std::make_pair(&mei_data, userSettings);
        
        // Load discordant reads into mei_data.
        bam_fetch(fp, idx, tid, window.getStart(), window.getEnd(), &env, fetch_disc_read_callback);
        bam_index_destroy(idx);
    }
    return 0;
}





static int append_cluster_connections(std::vector<MEI_event>& insertion_events, ControlState& current_state,
                                      UserDefinedSettings* userSettings) {
    
    // Setup maps for base read names of mates we need to collect.  Also setup 'exclude_names' holding
    // the original read names (we don't want those, only their mates, which fall inside the event).
    std::map<std::string, size_t> fwd_name_links, rev_name_links, exclude_names;
    std::string tmp_basename;
    for (size_t i = 0; i < insertion_events.size(); i++) {
        MEI_event event = insertion_events.at(i);
        for (size_t j = 0; j < event.fwd_cluster_bp.associated_reads.size(); j++) {
            tmp_basename = base_read_name(event.fwd_cluster_bp.associated_reads.at(j).name);
            fwd_name_links.insert(std::make_pair(tmp_basename, i));
            exclude_names.insert(std::make_pair(event.fwd_cluster_bp.associated_reads.at(j).name, i));
        }
        for (size_t j = 0; j < event.rev_cluster_bp.associated_reads.size(); j++) {
            tmp_basename = base_read_name(event.rev_cluster_bp.associated_reads.at(j).name);
            rev_name_links.insert(std::make_pair(tmp_basename, i));
            exclude_names.insert(std::make_pair(event.rev_cluster_bp.associated_reads.at(j).name, i));
        }
    }
    
    
    // Loop over whole genome to find mates of discordant reads near DD breakpoints.
    g_genome.reset();
    MEI_data mei_data;
    SearchRegion region = SearchRegion("ALL");
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
        
        LoopingSearchWindow currentWindow(&region, currentChromosome, WINDOW_SIZE);
        // loop over one chromosome
        do {
            LOG_DEBUG(*logStream << time_log() << "Dispersed Duplication detection current window: " <<
                     currentWindow.getChromosomeName() << " " << currentWindow.getStart() << "--" <<
                     currentWindow.getEnd() << std::endl);
            // Load discordant reads.
            int result = load_discordant_reads(mei_data, current_state.bams_to_parse, currentChromosome->getName(),
                                               currentWindow, userSettings);
            if (result) {
                // something went wrong loading the reads, return error code.
                return result;
            }
            
            // todo: store reads here!!!
            std::map<std::string, size_t>::iterator name_match;
            size_t disc_read_count = mei_data.discordant_reads.size();
            for (size_t i = 0; i < disc_read_count; i++) {
                // Determine event and strand for which mate is evidence.
                tmp_basename = base_read_name(mei_data.discordant_reads.at(i)->name);
                int event_idx = -1;
                char strand = Plus;
                name_match = fwd_name_links.find(tmp_basename);
                if (name_match != fwd_name_links.end()) {
                    // Current read referenced by a DD event near bp on fwd strand.
                    event_idx = (*name_match).second;
                } else {
                    name_match = rev_name_links.find(tmp_basename);
                    if (name_match != rev_name_links.end()) {
                        // Current read referenced by a DD event near bp on rev strand.
                        event_idx = (*name_match).second;
                        strand = Minus;
                    }
                }
                
                if (event_idx == -1) {
                    // No match found, this read is not related to an event.
                    continue;
                }
                
                if (exclude_names.find(mei_data.discordant_reads.at(i)->name) != exclude_names.end()) {
                    // read name in exlude list, this is one of the reads we used for calling
                    // the breakpoint, skip it!
                    continue;
                }

                if (strand == Plus) {
                    insertion_events.at(event_idx).fwd_mapping_reads.push_back(*(mei_data.discordant_reads.at(i)));
                } else {
                    insertion_events.at(event_idx).rev_mapping_reads.push_back(*(mei_data.discordant_reads.at(i)));
                }
            }
            
            cleanup_reads(mei_data.discordant_reads);
            currentWindow.next();
            
        } while (!currentWindow.finished());
    } while (true);
    return 0;
}


void searchMEI(MEI_data& finalState, Genome& genome, std::map<int, std::string>& seq_name_dict,
               UserDefinedSettings* userSettings, ControlState& current_state, std::ostream& out) {
    LOG_INFO(*logStream << time_log() << "Start calling dispersed duplication events from found breakpoints..." << std::endl);
    std::vector<MEI_event> insertion_events;
    
    size_t bp_amount = finalState.MEI_breakpoints.size();
    LOG_INFO(*logStream << time_log() << "Examining " << bp_amount << 
             " breakpoints in total." << std::endl);
    
    std::sort(finalState.MEI_breakpoints.begin(), finalState.MEI_breakpoints.end(), comp_breakpoint_pos);
    LOG_DEBUG(*logStream << time_log() << "Sorted breakpoints." << std::endl);
        
    for (size_t i = 0; i < (bp_amount-1); i++) {
        if (finalState.MEI_breakpoints.at(i).cluster_strand == finalState.MEI_breakpoints.at(i+1).cluster_strand ||
            (finalState.MEI_breakpoints.at(i+1).breakpoint_pos - finalState.MEI_breakpoints.at(i).breakpoint_pos) > 
                userSettings->MAX_DD_BREAKPOINT_DISTANCE ||
            finalState.MEI_breakpoints.at(i).breakpoint_tid != finalState.MEI_breakpoints.at(i+1).breakpoint_tid) {
            // Current two consecutive breakpoints cannot be combined into an event.
            continue;
        }
        
        MEI_event event;
        if (finalState.MEI_breakpoints.at(i).cluster_strand == Plus) {
            event = MEI_event(finalState.MEI_breakpoints.at(i), finalState.MEI_breakpoints.at(i+1));
        } else {
            event = MEI_event(finalState.MEI_breakpoints.at(i+1), finalState.MEI_breakpoints.at(i));
        }
        
        insertion_events.push_back(event);
    }
    
    if (userSettings->DD_REPORT_DUPLICATION_READS) {
        // Append information about reads mapping inside DDs.
        LOG_INFO(*logStream << time_log() << "Collecting discordant read information for dispersed duplication "
                 << "events." << std::endl);
        append_cluster_connections(insertion_events, current_state, userSettings);
    }
    
    LOG_INFO(*logStream << time_log() << "Reporting " << insertion_events.size() << " dispersed duplication events to "
             << userSettings->getMEIOutputFilename().c_str() << std::endl);
    
    // Report events.
    for (size_t i = 0; i < insertion_events.size(); i++) {
        reportMEIevent(finalState, insertion_events.at(i), i + 1, genome, seq_name_dict, out);
    }
    
    LOG_INFO(*logStream << time_log() << "Found " << insertion_events.size() << " dispersed duplication events."
             << std::endl);
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
    
    // Reset genome before traversal.
    g_genome.reset();
    
    std::ofstream file_output(userSettings->getMEIOutputFilename().c_str());
    MEI_data mei_data;
    int result;
    
    do {
        const Chromosome* currentChromosome = g_genome.getNextChromosome();
        if (currentChromosome == NULL) {
            break;
        }
        // Skip current chromosome if another one is specifically targeted.
        if (!userSettings->loopOverAllChromosomes() && 
            currentChromosome->getName() != userSettings->getRegion()->getTargetChromosomeName()) {
            continue;
        }
        
        g_maxPos = 0;
        CurrentChrMask.resize(currentChromosome->getCompSize());
        for (unsigned int i = 0; i < currentChromosome->getCompSize(); i++) {
            CurrentChrMask[i] = 'N';
        }
        
        LoopingSearchWindow currentWindow(userSettings->getRegion(), currentChromosome, WINDOW_SIZE);
        // loop over one chromosome
        do {
            LOG_INFO(*logStream << time_log() << "Dispersed Duplication detection current window: " <<
                     currentWindow.getChromosomeName() << " " << currentWindow.getStart() << "--" <<
                     currentWindow.getEnd() << std::endl);
            result = load_discordant_reads(mei_data, current_state.bams_to_parse, currentChromosome->getName(),
                                           currentWindow, userSettings);
            if (result) {
                // something went wrong loading the reads, return error code.
                return result;
            }
            
            searchMEIBreakpoints(mei_data, current_state.bams_to_parse, currentChromosome, userSettings);
            cleanup_reads(mei_data.discordant_reads);
            currentWindow.next();
            
        } while (!currentWindow.finished());
    } while (true);
    
    // Reset genome for subsequent traversals.
    g_genome.reset();
    
    std::map<int, std::string> seq_name_dictionary = get_sequence_name_dictionary(current_state);
   
    searchMEI(mei_data, genome, seq_name_dictionary, userSettings, current_state, file_output);
    file_output.close();
    return 0;
}









