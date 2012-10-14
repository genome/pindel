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
    return read1.pos < read2.pos;
}


// Compare two split reads on size of mapped part.
bool comp_simple_read_mapsize(const simple_read& read1, const simple_read& read2) {
    return read1.mapped_sequence.length() > read2.mapped_sequence.length();
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
int get_safe_chr_index(int true_index) {
    return true_index + g_SpacerBeforeAfter;
}


// Utility function to convert an index on the safe representation of the chr
// (with N's added at both ends) to the true location on the chromosome.
int get_true_chr_index(int safety_index) {
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
            return chr->getSeq().substr(get_safe_chr_index(chr_position), subsequence_length);
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

