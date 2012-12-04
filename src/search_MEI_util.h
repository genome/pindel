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

#ifndef MEI_detection_search_MEI_util_h
#define MEI_detection_search_MEI_util_h

#include <sstream>

#include "pindel.h"
#include "search_MEI.h"

// Return current time in string format as used for logging.
std::string time_log();

bool comp_simple_read(simple_read* read1, simple_read* read2);

bool comp_simple_read_pos(const simple_read& read1, const simple_read& read2);

bool comp_simple_read_mapsize(const simple_read& read1, const simple_read& read2);

bool comp_simple_read_unmapped_seqsize(const simple_read& read1, const simple_read& read2);

void print_read(simple_read read);

std::string enrich_read_name(char* base_name, bool is_first_read);

std::string base_read_name(std::string& enriched_read_name);

std::string get_sequence(uint8_t* sam_seq, int sam_seq_len);

int get_comp_chr_index(int true_index);

int get_bio_chr_index(int safety_index);

std::string get_fasta_subseq(Genome& genome, std::string& chr_name, int chr_position, int subsequence_length);

bool comp_indexed_read(std::pair<simple_read, int> read1, std::pair<simple_read, int> read2);

std::string get_whitespace(unsigned int amount);

bool contains_subseq(const std::string& query, const std::string& db, int min_length);

bool contains_subseq_any_strand(const std::string& query, const std::string& db, int min_score=12);

// Construct a mapping from read group ID to sample name given a BAM header.
std::map<std::string, std::string> get_sample_dictionary(bam_header_t* header);

void get_sample_name(std::string& read_group, std::map<std::string, std::string>& sample_dictionary, 
                     std::string& sample_name);

void cleanup_reads(std::vector<simple_read*>& read_vector);

#endif
