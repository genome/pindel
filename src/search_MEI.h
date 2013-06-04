//
//  search_MEI.h
//  pindel_MEI
//
//  Created by Mark Kroon on 20/06/12.
//  Copyright (c) 2012 LUMC. All rights reserved.
//

#ifndef pindel_MEI_search_MEI_h
#define pindel_MEI_search_MEI_h

#include <vector>

#include "sam.h"

// Forward declarations.
class ControlState;
class flags_hit;

// Class storing a read with a link to its mate.
class simple_read {
public:
    simple_read(std::string name, int32_t tid, int32_t pos, char strand, std::string sample_name, std::string sequence,
                int32_t mate_tid, int32_t mate_pos, char mate_strand) : name(name), tid(tid), pos(pos), strand(strand), 
            sample_name(sample_name), sequence(sequence), mate_tid(mate_tid), mate_pos(mate_pos), 
            mate_strand(mate_strand), is_split(false) {};
    simple_read(std::string name, int32_t tid, int32_t pos, char strand, std::string sample_name,
                std::string sequence, std::string mapped_sequence, std::string unmapped_sequence) : name(name), 
            tid(tid), pos(pos), strand(strand), sample_name(sample_name), sequence(sequence), mate_tid(-1),
            mate_pos(-1), mate_strand('?'), is_split(true), mapped_sequence(mapped_sequence), 
            unmapped_sequence(unmapped_sequence) {};
    // Read name.
    std::string name;

    // Mapped chromosome ID.
    int32_t tid;
    
    // Mapped position.
    int32_t pos;
    
    // Mapped strand.
    char strand;
    
    // Evidence strand (what side of DD event is this read evidence for, used for reporting)
    char evidence_strand;
    
    // Name of sample.
    std::string sample_name;
    
    // DNA sequence.
    std::string sequence;
    
    // Mate's mapped chromosome ID.
    int32_t mate_tid;
    
    // Mate's mapped position.
    int32_t mate_pos;
    
    // Mate's mapped strand.
    char mate_strand;
    
    // Extra fields for split reads (should be subclassed in an ideal world).
    bool is_split;
    std::string mapped_sequence;
    std::string unmapped_sequence;
};



class MEI_breakpoint {
public:
    MEI_breakpoint() {};
    MEI_breakpoint(int breakpoint_pos) : breakpoint_pos(breakpoint_pos) {};
    MEI_breakpoint(int breakpoint_tid, int breakpoint_pos, char cluster_strand) : 
    breakpoint_tid(breakpoint_tid), breakpoint_pos(breakpoint_pos), cluster_strand(cluster_strand) {};
    int breakpoint_tid;                                 // breakpoint chromosome identifier
    int breakpoint_pos;                                 // breakpoint location on chromosome
    char cluster_strand;                                // strand to which reads of associated cluster are mapped to
    std::string chromosome_name;                        // display name for chromosome identified by breakpoint_tid
    std::vector<simple_read> associated_reads;          // reads from cluster
    std::vector<simple_read> associated_split_reads;    // split reads around breakpoint
};


// Class containing relevant data while gathering discordant reads and finding breakpoints.
class MEI_data {
public:
    // Split reads needed for precise breakpoint estimation.
    std::vector<SPLIT_READ> split_reads;
    // Reads that are discordantly mapped.
    std::vector<simple_read*> discordant_reads;
    // Potential breakpoints found for MEI variants.
    std::vector<MEI_breakpoint> MEI_breakpoints;
    
    // Bamfile-specific environment values
    unsigned int current_insert_size;
    std::string current_chr_name;
    std::map<std::string, std::string> sample_names;
};


bool is_concordant(const bam1_t* read);

std::string get_sequence(uint8_t* sam_seq, int sam_seq_len);

// Finds all breakpoints indicating mobile element insertions inside the
// region spanned by the given ControlState (currentState).  This method
// assumes that currentState already contains the split reads found in
// this region.
void searchMEIBreakpoints(MEI_data& currentState, std::vector<bam_info>& bam_sources, std::string& chr_name,
                          std::string& chr_sequence);

// Generates mobile element insertion calls MEI-breakpoints contained by
// the given MEI data storage object (mei_data).  This method assumes that
// searchMEIBreakpoints has already been called previously for finalState.
void searchMEI(MEI_data& mei_data, Genome& genome, std::map<int, std::string>& seq_name_dict, std::ostream& out);

// Main entry point for MEI detection.
int searchMEImain(ControlState& current_state, Genome& genome, UserDefinedSettings* user_settings);

#endif
