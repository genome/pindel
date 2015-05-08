/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "search_MEI_util.h"
#include "sam_header.h"

// Logging libraries.
#include "logdef.h"
#include "logstream.h"


// Comparator for sorting connected reads, first on mapped strand, then on position.
bool comp_simple_read(simple_read* read1, simple_read* read2) {
    if (read1->strand == '+' && read2->strand != '+') {
        return true;
    } else if (read1->strand != '+' && read2->strand == '+') {
        return false;
    } else {
        return read1->pos < read2->pos;
    }
}


// Compare two reads on position.
bool comp_simple_read_pos(const simple_read& read1, const simple_read& read2) {
    if (read1.tid != read2.tid) {
        // By chromosome.
        return read1.tid < read2.tid;
    }
    if (read1.pos != read2.pos) {
        // By position on chromosome.
        return read1.pos < read2.pos;
    }
    // By strand.
    return read1.strand < read2.strand;
}


// Compare two split reads on size of mapped part.
bool comp_simple_read_mapsize(const simple_read& read1, const simple_read& read2) {
    return read1.mapped_sequence.length() > read2.mapped_sequence.length();
}

// Compare two split reads on size of unmapped part.
bool comp_simple_read_unmapped_seqsize(const simple_read& read1, const simple_read& read2) {
    return read1.unmapped_sequence.length() > read2.unmapped_sequence.length();
}


// Simply print a read to stdout.
void print_read(simple_read read) {
    std::cout << read.name << " (";
    std::cout << read.pos << ")";
    if (read.strand == '+') {
        std::cout << " <--";
    } else {
        std::cout << " -->";
    }
    std::cout << std::endl;
}


// Return an 'enriched' read name given it's base name and whether 
std::string enrich_read_name(char* base_name, bool is_first_read) {
    std::stringstream output;
    output << "@" << base_name;
    if (is_first_read) {
        output << "/1";
    } else {
        output << "/2";
    }
    return output.str();
}


// Inversed of enrich_read_name(), if it cannot be done, simply return the given string unchanged.
std::string base_read_name(std::string& enriched_read_name) {
    size_t found = enriched_read_name.find("/", 0);
    if (found != std::string::npos && found > 0) {
        return enriched_read_name.substr(1, found - 1);
    }
    return enriched_read_name;
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


// Utility function to convert an index on a chromosome to an index on the safe
// representation of the chromosome (with N's added at both ends).
int get_comp_chr_index(int true_index) {
    return true_index + g_SpacerBeforeAfter;
}


// Utility function to convert an index on the safe representation of the chr
// (with N's added at both ends) to the true location on the chromosome.
int get_bio_chr_index(int safety_index) {
    return safety_index - g_SpacerBeforeAfter;
}


// Look into the given fasta input (fasta_file) for a subsequence of the given
// chromosome (chr_name), starting at given position (chr_position) with given 
// length (subsequence_length).
std::string get_fasta_subseq(Genome& genome, std::string& chr_name, int chr_position, int subsequence_length) {
    
    // Loop over chromosomes.
    genome.reset();
    const Chromosome* chr = genome.getNextChromosome();
    while (chr != NULL) {
        
        if (chr->getName() == chr_name) {
            return chr->getSeq().substr(get_comp_chr_index(chr_position), subsequence_length);
        }
        chr = genome.getNextChromosome();
    }
    // Requested chromosome not encountered.  Todo: throw exception.
    return "";
}


// Comparator for sorting indexed reads from a contig, by their relative positions.
bool comp_indexed_read(std::pair<simple_read, int> read1, std::pair<simple_read, int> read2) {
    return read1.second < read2.second;
}


// Get a given amount of spaces in a string.
std::string get_whitespace(unsigned int amount) {
    std::stringstream ss;
    for (unsigned int i = 0; i < amount; i++) {
        ss << " ";
    }
    return ss.str();
}


std::string time_log() {
    std::stringstream ss;
    ss.precision(2);
    ss << std::fixed << ((float) clock()) / CLOCKS_PER_SEC << ": ";
    return ss.str();
}


// Returns true if a subsequence of query of at least min_match_length in length matches
// some subsequence in db with a score of at least min_score.  Score is computed by adding
// 2 for each match, -1 for each gap or substitution.
bool contains_subseq(const std::string& query, const std::string& db, int min_length) {
    
    size_t db_length = db.length();
    size_t query_length = query.length();
    
    // Minimum nr of matches in valid alignment.
    int min_match_length = min_length - (int) g_maxMismatch[min_length];
    
    int MATCH_SCORE = 1;
    int MISMATCH_SCORE = -2;
    
    // Setup 2 rows of the SW alignment matrix.
    int mismatch_countA[db_length];
    int mismatch_countB[db_length];
    int* prev_mc = mismatch_countA;
    int* current_mc = mismatch_countB;
    int* temp;
    
    // Setup array for storing current match length. (i.e. the number of matches in the
    // optimal alignment)
    int alignment_lengthsA[db_length];
    int alignment_lengthsB[db_length];
    int* prev_al = alignment_lengthsA;
    int* current_al = alignment_lengthsB;
    
    for (size_t i = 0; i < db_length; i++) {
        *(prev_mc + i) = 0;
        *(prev_al + i) = 0;
    }
    
    // Counters for current entry in matrix and totals for current row.
    int score;
    int max_score;
    char action;
    int max_alignment_length_current_row;
    int min_mismatch_count_current_row;

    // Loop over query.
    for (size_t i = 0; i < query_length; i++) {
        min_mismatch_count_current_row = 0;
        max_alignment_length_current_row = 0;
        
        // Set value of first column.
        *current_mc = 0;
        *current_al = 0;
        if (db[0] == query[i]) {
            *current_mc = 0;
            *current_al = 1;
        } else {
        }
        
        // Loop over db.
        for (size_t j = 1; j < db_length; j++) {
            
            max_score = 0;
            action = 'n';
            
            // Get score for alignment with match of current base on query and db.
            score = (*(prev_al + (j - 1)) + 1) * MATCH_SCORE + *(prev_mc + (j - 1)) * MISMATCH_SCORE;
            if (query[i] == db[j] && max_score < score) {
                max_score = score;
                action = 'm';
            } else {
                // Get score for alignment with substitution of current bases.
                score = *(prev_al + (j - 1)) * MATCH_SCORE + (*(prev_mc + (j - 1)) + 1) * MISMATCH_SCORE;
                if (max_score < score) {
                    max_score = score;
                    action = 'M';
                }
            }
            
            // Get score for alignment with gap on query.
            score = (*(current_al + (j - 1))) * MATCH_SCORE + (*(current_mc + (j - 1)) + 1) * MISMATCH_SCORE;
            if (max_score < score) {
                max_score = score;
                action = 'g';
            }
            
            // Get score for alignment with gap on db.
            score = (*(prev_al + j) + 1) * MATCH_SCORE + (*(prev_mc + j) + 1) * MISMATCH_SCORE;
            if (max_score < score) {
                max_score = score;
                action = 'G';
            }
                        
            // Set values for current mismatch count and alignment length.
            switch (action) {
                case 'g':
                    // Gap on query.
                    *(current_mc + j) = *(current_mc + (j - 1)) + 1;
                    *(current_al + j) = *(current_al + (j - 1));
                    break;
                case 'G':
                    // Gap on db.
                    *(current_mc + j) = *(prev_mc + j) + 1;
                    *(current_al + j) = *(prev_al + j) + 1;
                    break;
                case 'm':
                    // Match.
                    *(current_mc + j) = *(prev_mc + (j - 1));
                    *(current_al + j) = *(prev_al + (j - 1)) + 1;
                    break;
                case 'M':
                    // Mismatch (substitution).
                    *(current_mc + j) = *(prev_mc + (j - 1)) + 1;
                    *(current_al + j) = *(prev_al + (j - 1)) + 1;
                    break;
                default:
                    // New alignment (action == 'n').
                    if (query[i] == db[j]) {
                        *(current_mc + j) = 0;
                        *(current_al + j) = 1;
                    } else {
                        *(current_mc + j) = 1;
                        *(current_al + j) = 1;
                    }
                    break;
            }
            
            if (*(current_al + j) >= min_length && *(current_mc + j) <= (int) g_maxMismatch[*(current_al + j)]) {
                // Valid alignment found! No need to look further.
                return true;
            }
            
            // Keep track of best score for alignments at this point in query.
            if (*(current_mc + j) < min_mismatch_count_current_row) {
                min_mismatch_count_current_row = *(current_mc + j);
            }
            if (*(current_al + j) > max_alignment_length_current_row) {
                max_alignment_length_current_row = *(current_al + j);
            }
        }

        if ((int) (query_length - i - 1) + (max_alignment_length_current_row - min_mismatch_count_current_row) < 
            min_match_length) {
            // With ideal conditions (lowest number of mismatches, perfect matching in rest of query),
            // the current longest alignment is still too short to reach the required quality settings 
            // (min_match_length).
            return false;
        }
        
        // Swap arrays, such that current becomes prev for next iteration.
        temp = prev_al;
        prev_al = current_al;
        current_al = temp;
        temp = prev_mc;
        prev_mc = current_mc;
        current_mc = temp;
    }
    
    // No local alignment found with score above threshold.
    // Note: this statement will probably not be reached, but previous 'return false' statement would.
    return false;
}


// Returns true if query matches a subequence of db or a subsequence of db reverse
// complimented, where the match has a minimal score of min_score.  Score is computed
// by adding 2 for each matching base, -1 for each gap or mismatch.
bool contains_subseq_any_strand(const std::string& query, const std::string& db, int min_score) {
    return contains_subseq(query, db, min_score) || contains_subseq(ReverseComplement(query), db, min_score);
}



// Set up map linking read group ids with sample names.
std::map<std::string, std::string> get_sample_dictionary(bam_header_t* header) {
    std::map<std::string, std::string> sample_dict;
    int num_ids;
    int num_names;
    char** tmp_ids;
    if (header->dict == 0) header->dict = sam_header_parse2(header->text);
    // Convert string literals to char* to feed into samtools api.
    char* rg_tag = new char[3];
    strncpy(rg_tag, "RG", 2); 
    rg_tag[2] = '\0';
    char* id_tag = new char[3];
    strncpy(id_tag, "ID", 2);
    id_tag[2] = '\0';
    char* sm_tag = new char[3];
    strncpy(sm_tag, "SM", 2);
    id_tag[2] = '\0';
    tmp_ids = sam_header2list(header->dict, rg_tag, id_tag, &num_ids);
    char** tmp_names;
    tmp_names = sam_header2list(header->dict, rg_tag, sm_tag, &num_names);
    if (num_ids > 0 && num_ids == num_names) {
        for (int i = 0; i < num_ids; i++) {
            sample_dict.insert(std::make_pair(tmp_ids[i], tmp_names[i]));
        }
    }
    free(tmp_ids);
    free(tmp_names);
    delete rg_tag;
    delete id_tag;
    delete sm_tag;
    return sample_dict;
}



void get_sample_name(std::string& read_group, std::map<std::string, std::string>& sample_dictionary, std::string& sample_name) {
    try {
        sample_name = sample_dictionary.at(read_group);
    } catch (std::exception& e) {
        // Read group id not found in mapping. todo: increase severity of warning?
        LOG_DEBUG(*logStream << "Could not find sample name for read group: " << read_group << std::endl);
        if (g_sampleNames.size() == 1) {
            // Only one sample name known in pindel's input.  Assume current read group is
            // from this sample.
            sample_name = *(g_sampleNames.begin());
        }
    }
}

void cleanup_reads(std::vector<simple_read*>& read_vector) {
    size_t vec_size = read_vector.size();
    for (size_t i = 0; i < vec_size; i++) {
        delete read_vector.at(i);
    }
    read_vector.clear();
}



