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
    simple_read(std::string name, int32_t tid, int32_t pos, bool strand, std::string sequence) : 
            name(name), tid(tid), pos(pos), strand(strand), sequence(sequence) {};
    simple_read(std::string name, int32_t tid, int32_t pos, bool strand, std::string sequence, 
                std::string mapped_sequence, std::string unmapped_sequence) : name(name), tid(tid), pos(pos), 
            strand(strand), sequence(sequence), mapped_sequence(mapped_sequence), 
            unmapped_sequence(unmapped_sequence) {};
    std::string name;
    int32_t tid;
    int32_t pos;
    bool strand;
    std::string sequence;
    
    // Extra fields for split reads (should be subclassed in an ideal world).
    std::string mapped_sequence;
    std::string unmapped_sequence;
};



class MEI_breakpoint {
public:
    MEI_breakpoint(int breakpoint_pos) : breakpoint_pos(breakpoint_pos) {};
    MEI_breakpoint(int breakpoint_tid, int breakpoint_pos, bool cluster_strand) : 
    breakpoint_tid(breakpoint_tid), breakpoint_pos(breakpoint_pos), cluster_strand(cluster_strand) {};
    int breakpoint_tid;
    int breakpoint_pos;
    bool cluster_strand;
    std::vector<simple_read> associated_reads;          // reads from cluster
    std::vector<simple_read> associated_split_reads;    // split reads around breakpoint
};


class MEI_data {
public:
    // Split reads needed for precise breakpoint estimation.
    std::vector<SPLIT_READ> split_reads;
    // Reads that are discordantly mapped.
    std::vector<simple_read*> discordant_reads;
    // Potential breakpoints found for MEI variants.
    std::vector<MEI_breakpoint> MEI_breakpoints;
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
// the given ControlState (finalState).  This method assumes that
// searchMEIBreakpoints has already been called previously for finalState.
void searchMEI(MEI_data& finalState);

// Main entry point for MEI detection.
int searchMEImain(ControlState& current_state, std::ifstream&, UserDefinedSettings* user_settings);

#endif
