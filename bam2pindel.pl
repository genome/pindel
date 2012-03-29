#!/usr/local/bin/perl

########################################################################################################
#                                                                                                      #
# CGP Software License                                                                                 #
#                                                                                                      #
# Copyright (c) 2010 Genome Research Ltd.                                                              #
# Author: Cancer Genome Project, cgpit@sanger.ac.uk                                                    #
#                                                                                                      #
# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT   #
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND               #
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES #
# OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN  #
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                           #
#                                                                                                      #
# This code is free software; you can redistribute it and/or modify it under the terms of the BSD      #
# License.                                                                                             #
#                                                                                                      #
# Any redistribution or derivation in whole or in part including any substantial portion of this code  #
# must include this copyright and permission notice.                                                   #
#                                                                                                      #
########################################################################################################

use strict;
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use Getopt::Long 'GetOptions';
use Data::Dumper;
use Carp;
use File::Path qw(make_path);

BEGIN {
  use Cwd qw(abs_path);
	my $prog_path = abs_path($0);
	$prog_path =~ s/bam2pindel\.pl$/Adaptor.pm/;
  my $adaptor_path;
  if(-e $prog_path) {
    $adaptor_path = $prog_path;
  }
  if($ENV{'BAM_2_PINDEL_ADAPT'}) {
    $adaptor_path = $ENV{'BAM_2_PINDEL_ADAPT'};
  }

	if(!$adaptor_path) {
		die "\nAdaptor.pm was not found in same location as bam2pindel.pl, nor was the Environment variable of BAM_2_PINDEL_ADAPT set with the path to your Adaptor.pm\n\n";
	}
	require $adaptor_path;
};

my $proper_pm_mask = 2; #0x0002
my $query_unmap_mask = 4; #0x0004
my $mate_unmap_mask = 8; #0x0008
my $strand_mask = 16; #0x0010 (set when reverse)
my $first_in_pair_mask = 64;
my $second_in_pair_mask = 128;

my %file_handles;
my %cent_coords;
my @record_store;
my $record_count = 0;
my %total_records; # used to check output file is complete

my $loops = 0; # for profiling
my $no_rec_in_write = 0;

my $file_count_loc;
my $bam_position = 0;
my $resume_pos = 0;

my $debug = 0;

{
  my $options = option_builder();
  rg_parse($options);
  if($options->{'c'}) {
  	centromere_parser($options->{'c'});
  }
  if($options->{'q'}) {
    qual_profile($options);
  }
  resume_run();
  run_parser($options);
  unlink $file_count_loc;
}

sub resume_run {
  if(-e $file_count_loc && -s $file_count_loc) {
    # first stat time file was created
    my @stat_arr = stat($file_count_loc);
    my $mtime = $stat_arr[9];
    open my $RESUME, '<', $file_count_loc or croak "Unable to open $file_count_loc for reading";
    while (my $line = <$RESUME>) {
      chomp $line;
      if($line =~ m/^(.*)\t(.*)\t([0-9]+)$/) {
        # is a file, chr_arm and count
        my $file = $1;
        my $chr_arm = $2;
        my $count = $3;
        @stat_arr = stat($file);
        if($stat_arr[9] > $mtime) {
          # this indicates that a partial write occurred
          # unable to resume from this
          undef %file_handles;
          undef %total_records;
          $resume_pos = 0;
          last;
        }
        $file_handles{$chr_arm}{'fn'} = $file;
        $total_records{$file} = $count;
      }
      elsif($line =~ m/^([0-9]+)$/) {
        # this is the bam position
        $resume_pos = $1;
      }
    }
    close $RESUME;
    if($resume_pos != 0) {
      print STDERR 'Able to resume at record ',$resume_pos,"\n";
    }
    else {
      print STDERR "Unable to resume, starting from scratch\n";
    }
  }
  return 1;
}

sub write_records {
  my ($options, $input) = @_;

  if($input) {
    #print STDERR "start cache\n";
    foreach my $rec(@{$input}) {
      if($rec && @{$rec} > 0) {
        my $chr_arm = get_chr_arm($options, $rec->[1], $rec->[2]);
        if($file_handles{$chr_arm}{'data'}) {
          push @{$file_handles{$chr_arm}{'data'}}, $rec->[0];
        }
        else {
          $file_handles{$chr_arm}{'data'} = [$rec->[0]];
        }
        $record_count++;
      }
    }
    #print STDERR "end cache\n";
  }

  if($record_count >= $options->{'mr'} || !$input) {
    my $file_count_str = q{};
    print STDERR 'Starting write of ',$record_count,' records',"\n";
    foreach my $chr_arm(keys %file_handles) {
      if($file_handles{$chr_arm}{'data'} && @{$file_handles{$chr_arm}{'data'}} > 0) {
        my ($CHR_FH, $file_path) = get_chr_fh($options, $chr_arm);
        foreach my $rec(@{$file_handles{$chr_arm}{'data'}}) {
          print $CHR_FH $rec;
          # keeps track of number of records that should be in each output file
          if($total_records{$file_path}) {
            $total_records{$file_path}++;
          }
          else {
            $total_records{$file_path} = 1;
          }
        }
				close $CHR_FH;
				$file_count_str .= $file_path."\t".$chr_arm."\t".$total_records{$file_path}."\n";
				$file_handles{$chr_arm}{'data'} = [];
      }
    }
    $record_count = 0;

    open my $FILE_COUNT, '>', $file_count_loc or croak "Unable to open for writing $file_count_loc";
    print $FILE_COUNT $file_count_str or croak "Unable to write output to $file_count_loc";
    print $FILE_COUNT $bam_position,"\n" or croak "Unable to write output to $file_count_loc";
    close $FILE_COUNT or croak "Failed to close $file_count_loc";

    print STDERR 'Finished writing records',"\n";
  }
  return;
}

sub centromere_parser {
	my ($cent_file) = @_;
	open my $CENTS, '<', $cent_file or croak "Could not open file $cent_file";
	while (my $line=<$CENTS>) {
		if(index($line,'#') == 0) {
			next;
		}
		my @gff_bits = split /\t/, $line;
		my $cent_chr = $gff_bits[0];
		my $cent_start = $gff_bits[3];
		my $cent_end = $gff_bits[4];
		$cent_coords{$cent_chr} = [$cent_start, $cent_end];
	}
	close $CENTS;
	return;
}

sub rg_parse {
  my ($options) = @_;
  # get the header only
  my $command = find_prog('samtools').' view -H '.$options->{'i'};
  my $pid = open my $PROC, '-|', $command or croak "Could not fork: $OS_ERROR";

  while( my $tmp = <$PROC> ) {
    if($tmp !~ /^\@RG/xs) {
      next;
    }
    chomp $tmp;
    my @rg_data = split /\t/xms, $tmp;
    # have to find tag each time as merged data may have different tag positions
    my $id_pos = tag_finder('ID', \@rg_data, 1); # no point parsing 0, as just @RG tag
    my $pi_pos = tag_finder('PI', \@rg_data, 1); # no point parsing 0, as just @RG tag
    if($pi_pos) {
    	if(!$options->{'rgs'}) {
				$options->{'rgs'} = {};
			}
			my (undef, $id) = split /:/xms, $rg_data[$id_pos];
			my (undef, $pi) = split /:/xms, $rg_data[$pi_pos];
			next if($options->{'rg'} && $id != $options->{'rg'});
			$options->{'rgs'}->{$id} = $pi;
    }
  }
  close $PROC or croak 'Failed to close process: '.$command."\nERROR CODE: ".$CHILD_ERROR;

  if(!$options->{'rgs'} && !$options->{'pi'}) {
  	croak q{No @RG header lines with PI tags detected, you must define '-pi'};
  }
}

# will return false if not a name sorted bam
sub name_sorted {
  my ($bam) = @_;
  my $command = find_prog('samtools').' view -H '.$bam.q{ | grep -F '@HD' | grep -F 'SO:queryname'};
  my $res = `$command`;
  return length($res);
}

sub find_prog {
  my ($prog) = @_;

  my $prog_path = abs_path($0);
	$prog_path =~ s/bam2pindel\.pl$/$prog/;
  if(!-e $prog_path) {
    undef $prog_path;
    my @possible_locs = split ':', $ENV{'PATH'};

    foreach my $path(@possible_locs) {
      $prog_path = $path.'/'.$prog;
      if(-e $prog_path) {
        last;
      }
      undef $prog_path;
    }
  }
  return $prog_path;
}

sub run_parser {
  my ($options) = @_;

  my $command = find_prog('samtools').' view ';
  if($options->{'om'}) { # use oldmethod for parsing of bam file - SLOW
  	$command = find_prog('samtools').' view ';
  }
  elsif($options->{'nm'}) {
  	$command = find_prog('samgroupbyname').' -p ';
  	if($options->{'rg'}) {
  	  $command .= '-r '.$options->{'rg'}.' ';
  	}
  }
  else {
    # no method supplied - check bam is name sorted
    if(!name_sorted($options->{'i'})) {
      bam_sort_order_error($options->{'i'});
    }
  }
  $command .= $options->{'i'};

  print STDERR "Forked pipe to run: ",$command,"\n";
  my $pid = open my $PROC, '-|', $command or croak "Could not fork: $OS_ERROR";

  my @pair;
  eval {
		MAIN: while( my $tmp = <$PROC> ) {
		  last MAIN unless(defined $tmp);
		  if(index($tmp, '@') == 0) {
				next;
			}
		  $bam_position++;
		  next if($bam_position < $resume_pos);
			my @records;
			chomp $tmp;
			my @comps = split /\t/, $tmp;
			my $read_1 = Pindel::Bam::Adaptor->new(\@comps);
			@pair = ($read_1);

			# this forces a paired set of reads into @pair
			my $paired = 0;

			while($paired == 0) {
				$tmp = <$PROC>;
				last MAIN unless(defined $tmp);
				$bam_position++;
				chomp $tmp;
				my @comps2 = split /\t/, $tmp;
				my $read_2 = Pindel::Bam::Adaptor->new(\@comps2);
				if($pair[0]->name() ne $read_2->name()) {
					# at this point take the existing [0]
					# read and run through single end recovery
					if($pair[0]->mapped()) {
						#debug('would run recovery on mapped unpaired read');
						push @records, @{read_recovery($pair[0], $options)};
					}

					$pair[0] = $read_2;
				}
				else {
					push @pair, $read_2;
					$paired = 1;
				}
			}

			my ($used_1, $used_2, $proc_res) = process_data(\@pair, $options);
			if($used_1 == 1 || $used_2 == 1) {
				push @records, @{$proc_res};
			}
			# run recovery over any read not used
			if($used_1 == 0) {
				push @records, @{read_recovery($pair[0], $options)};
			}
			if($used_2 == 0) {
				push @records, @{read_recovery($pair[1], $options)};
			}
			write_records($options, \@records) if(@records > 0);
		}
  };
	if($EVAL_ERROR) {
		my $error =  'Error while processing: '.$command."\nEVAL ERROR: ".$EVAL_ERROR."\nCHILD ERROR: ".$CHILD_ERROR."\n\nReads in buffer:\n";
		$error .= Dumper(\@pair);
		$error .= "\nBroke at record $bam_position\n";
		croak $error;
	}
  # close samtools view
  close $PROC or croak 'Failed to close process: '.$command."\nERROR CODE: ".$CHILD_ERROR;

  write_records($options);
  check_output_files();
}

sub check_output_files {
  my @files = keys %total_records;
  foreach my $file(@files) {
    my $res = `wc -l $file`;
    chomp $res;
    my ($count,undef) = split / /, $res;
    if(! defined $count) {
      croak 'failed to parse wc -l';
    }
    if($total_records{$file} != ($count/3)) {
      croak 'Output file is corrupt, please re-run';
      unlink $file_count_loc;
    }
  }
  return 1;
}

sub get_chr_arm {
  my ($options, $chr, $pos) = @_;
  my $file_chr = $chr;
	if($options->{'c'} && $chr ne 'MT') {
		if($pos > $cent_coords{$chr}->[0]) {
			$file_chr .= '-q';
		}
		else {
			$file_chr .= '-p';
		}
	}
	return $file_chr;
}

sub get_chr_fh {
	my ($options, $file_chr) = @_;

	my $CHR_FH;
	my $file_path;
	if(!$file_handles{$file_chr}{'fn'}) {
	  $file_path = $options->{'o'}.'_'.$file_chr.'.txt';
		open $CHR_FH, '>', $file_path  or croak 'Could not create file: '.$file_path;
		$file_handles{$file_chr}{'fn'} = $file_path;
    #$file_handles{$file_chr}{'count'} = 1;
    print STDERR 'First open: ',$file_path," (1)\n";
	}
	else {
  	$file_path = $file_handles{$file_chr}{'fn'};
    open $CHR_FH, '>>', $file_path  or croak 'Could not append to file: '.$file_path;
    #$file_handles{$file_chr}{'count'}++;
    #print STDERR 'Reopened: ',$file_path,' (',$file_handles{$file_chr}{'count'},")\n";
	}
	return ($CHR_FH, $file_path);
}

sub process_data {
  my ($read_pair, $options) = @_;

	my ($used_1, $used_2) = (0,0);

	# setup for basic short circuits to skip readpairs
	{
		if($options->{'r'} && $read_pair->[0]->rname() ne $options->{'r'} && $read_pair->[1]->rname() ne $options->{'r'}) {
			return ($used_1, $used_2);
		}
	}

	my $read_1 = $read_pair->[0];
	my $read_2 = $read_pair->[1];

  #if($debug) {
  	#debug("\n");
		##debug(join "\t", @{$read_1->sam()});
		##debug(join "\t", @{$read_2->sam()});
		#debug(Dumper($read_1));
		#debug(Dumper($read_2));
  #}

  if(!$read_1->mapped() && !$read_2->mapped()) { # both reads unmapped
  	#debug($read_pair->[0]->name().' Both reads unmapped');
    return ($used_1, $used_2);
  }

  # skip if both are SW
  if($read_1->sw() && $read_2->sw()) {
  	#debug($read_pair->[0]->name().' Both reads are SW');
    return ($used_1, $used_2);
  }

  # skip if both are unique with no edit distanct
  # NOTE: added NM check as hides all small indels otherwise
  if($read_1->unique() && $read_2->unique() && $read_1->edits() == 0 && $read_2->edits() == 0) {
  	#debug($read_pair->[0]->name().' Both reads are unique with no edits');
    return ($used_1, $used_2);
  }

  # skip if both are suboptimal
  if($read_1->suboptimal() && $read_2->suboptimal()) {
	  #debug($read_pair->[0]->name().' Both reads suboptimal');
    return ($used_1, $used_2);
  }

  my $max_edit = int ($read_1->read_length() / $options->{'e'}) + 1;

  # to be useful for Pindel at least one read must have full match with max of 2 edits
  if($read_1->edits() > $max_edit && $read_2->edits() > $max_edit) {
  	#debug($read_pair->[0]->name().' Both reads have more than '.$max_edit.' edits (allowing 1 edit per '.$options->{'e'}.' bps + 1)');
  	return ($used_1, $used_2);
  }

  my $uniq;

	if($read_1->mapped() && $read_2->mapped()) {
		if($read_1->unique() && $read_1->edits() <= $max_edit && $read_2->edits() > 0) {
			$uniq = 1;
		}
		elsif($read_2->unique() && $read_2->edits() <= $max_edit && $read_1->edits() > 0) {
			$uniq = 2;
		}
		else {
			#debug($read_pair->[0]->name().' Both reads mapped at sufficient accuracy, edit distance not useful');
			return ($used_1, $used_2);
		}
	}
	elsif(!$read_1->mapped() && $read_2->edits() <= $max_edit) {
		$uniq = 2;
	}
	elsif(!$read_2->mapped() && $read_1->edits() <= $max_edit) {
		$uniq = 1;
	}
	else {
		#debug($read_pair->[0]->name().' One read mapped at sufficient accuracy, but nm for anchor is >'.$max_edit.' edits (allowing 1 edit per '.$options->{'e'}.' bps + 1)');
		return ($used_1, $used_2);
	}

	my @records;
	if($uniq == 1)  {
		#debug('Read 1 is anchor');
		my $ref_a = build_record($read_1, $read_2, $options);
		if($ref_a) {
			push @records, $ref_a;
			$used_2 = 1;
		}
		if($read_2->unique() && $read_2->edits() <= $max_edit && $read_1->edits() > 0) {
			#debug('Read 2 is also anchor');
			my $ref = build_record($read_2, $read_1, $options);
			if($ref) {
				push @records, $ref;
				$used_1 = 1;
			}
		}
	}
	elsif($uniq == 2) {
		#debug('Read 2 is anchor');
		my $ref = build_record($read_2, $read_1, $options);
		if($ref) {
			push @records, $ref;
			$used_1 = 1;
		}
	}

  return ($used_1, $used_2, \@records);
}

sub read_recovery {
	my ($read_1, $options) = @_;
	my @results;

	if($options->{'r'} && $read_1->rname() ne $options->{'r'}) {
		return \@results;
	}

	# read_1 should be the read not previously used as a query
	#debug(Dumper($read_1));

	# try to use read_1 without ref to read 2
	if($read_1->mapped()) {
		if(cigar_is_clean($read_1->cigar())) {
			#debug('clean cigar in read recovery, no use to pindel');
			return \@results;
		}

		### new addition to try and reduce false positives
		if($read_1->sw()) {
			#debug('sw in read recovery, no use to pindel');
			return \@results;
		}
		if($read_1->suboptimal()) {
			#debug('suboptimal in recovery, no use to pindel');
			return \@results;
		}
		my @tmp_sam = @{$read_1->sam()}; # no longer a reference as would break stuff
		my $pi;
		if($options->{'pi'}) {
			$pi = $options->{'pi'};
		}
		else {
			$pi = $options->{'rgs'}{$read_1->read_group()};
		}
		if($read_1->strand()) { # +ve
			# make negative
			$tmp_sam[1] += 16;
			# as now negative + pi *1.2;
			$tmp_sam[3] += int($pi*1.2);
		}
		else {
			$tmp_sam[1] -= 16;
			# as now positive - pi *1.2;
			$tmp_sam[3] -= int($pi*1.2);
		}

		my $fake_anchor = Pindel::Bam::Adaptor->new(\@tmp_sam);
		#debug(Dumper($fake_anchor));
		#debug('fake anchor used');
		push @results, build_record($fake_anchor, $read_1, $options);
		return \@results;
	}
	else {
		#debug('Unmapped in recovery, discard');
		return \@results;
	}
	# not actually possible to get here
	return \@results;
}

# meaning all M
sub cigar_is_clean {
	my $cigar = shift;
	$cigar =~ s/[0-9]//g;
	if($cigar eq 'M') {
		return 1;
	}
	return 0;
}

sub trim_on_qual {
  my ($seq, $quals_in, $thresh, $flip) = @_;
  my $trim_at = 0;
  my @quals;
  if($flip) {
    $seq = reverse $seq;
    @quals = reverse @{$quals_in};
  }
  my $i = 0;
  for(; $i<@quals-1; $i++) {
    #print "$i - ",$quals[$i],"\n";
    if($quals[$i] > $thresh && $quals[$i+1] > $thresh) {
      #print "Threshold found";
      last;
    }
  }
  #print $i,"\n";
  #print $seq,"\n";

  my $new_seq;
  #if($i >= 4) {
    $new_seq = substr $seq, $i;
    if(length $new_seq < ((length $seq) * 0.4) ) {
      return undef;
    }
  #}
  #else {
  #  $new_seq = $seq;
  #}

  #print $new_seq,"\n";
  if($flip) {
    $new_seq = reverse $new_seq;
    #print $new_seq,"\n\n";
  }
  return $new_seq;
}

sub unmapped_seq {
  my ($options, $unmap_read) = @_;
  my $seq = $unmap_read->sam()->[9];

  # this component added to apply quality trimming
  if($options->{'q'}) {
    my $quals = qual_to_vals($unmap_read->sam()->[10]);
    my $qual_thresh = $options->{'rg_profile'}->{$unmap_read->read_group}->{'qual_thresh'};
    if($quals->[0] < $qual_thresh) {
      $seq = trim_on_qual($seq, $quals, $qual_thresh, 0);
    }
    elsif($quals->[-1] < $qual_thresh) {
      $seq = trim_on_qual($seq, $quals, $qual_thresh, 1);
    }

    if(!$seq) {
      return $seq;
    }
  }

  if(!$unmap_read->strand()) {
    $seq =~ tr/ACGT/TGCA/;
    $seq = reverse $seq;
  }

  if(!$unmap_read->mapped() || (index($unmap_read->sam()->[5], 'S') > -1)) {
  	$seq =~ s/^N+//;
  	$seq =~ s/N+$//;
  }

  return $seq;
}

sub length_correct {
  my ($map_read) = @_;
  my $length = 0;
  my ($seq, $cigar) = ($map_read->sam()->[9], $map_read->sam()->[5]);
  my @cigar_lengths = split /[A-Z]+/, $cigar;
  my @cigar_operations = split /[0-9]+/, $cigar;
  shift @cigar_operations; # will always have empty first element;
  if(@cigar_lengths != @cigar_operations) {
  	print STDERR Dumper($map_read);
    croak 'Cigar string has unmatched lengths and operators';
  }
  if(@cigar_lengths == 1) {
    if($cigar_lengths[0] == length $seq) {
      $length = $cigar_lengths[0];
    }
    else {
    	print STDERR Dumper($map_read);
      croak 'Cigar string has no variation yet does not match sequence length';
    }
  }
  else {
    for my $i(0..@cigar_lengths-1) {
      if($cigar_operations[$i] eq 'M') {
        $length += $cigar_lengths[$i];
      }
      elsif($cigar_operations[$i] eq 'I') {
        $length += $cigar_lengths[$i];
      }
      elsif($cigar_operations[$i] eq 'S') {
        $length += $cigar_lengths[$i];
      }
      elsif($cigar_operations[$i] eq 'D') {
        $length -= $cigar_lengths[$i];
      }
      else {
        $debug = 1;
        #debug('CIGAR: '.$cigar.' seq length: '.length $seq);
        #debug('Operation Length:', Dumper(\@cigar_lengths));
        #debug('Operation Type:', Dumper(\@cigar_operations));
        croak 'I only understand M,I,S,D for length correction';
      }
    }
  }
  return $length;
}

sub trim_n_from_ends {
	my ($seq) = @_;
	$seq =~ s/^N+//;
	$seq =~ s/N+$//;
	my $seq_len = length $seq; # seq withouth runs of N's at each end
	my $n_in_seq = $seq =~ tr/N//; # returns count of N's within remaning seq
	return ($seq_len, $n_in_seq);
}

sub build_record {
  my ($map_read, $unmap_read, $options) = @_;
  my $pindel_rec;
  my @extended_rec;

  my $sequence = unmapped_seq($options, $unmap_read);
  if(!$sequence) {
    return \@extended_rec;
  }

  my $read_name = $map_read->name();

	my ($seq_len_unmap, $n_in_seq_unmap) = trim_n_from_ends($unmap_read->sam()->[9]);
	my ($seq_len_map, $n_in_seq_map) = trim_n_from_ends($map_read->sam()->[9]);

	my $max_n_unmap = int ($seq_len_unmap * $options->{'mn'});
	my $max_n_map = int ($seq_len_map * $options->{'mn'});

	if($seq_len_unmap >= $options->{'ms'} && $n_in_seq_unmap <= $max_n_unmap
	&& $seq_len_map >= $options->{'ms'} && $n_in_seq_map <= $max_n_map) {
		my $length = -1; # pindel is 0 based not 1 based so coords of +ve mappings will be -1 to SAM record
		if(!$map_read->strand()) {
			# -1 so return the location of the last base not the one after
			$length = length_correct($map_read) - 1;
		}

		my $pi_val;
		if(!$map_read->read_group() && !$options->{'pi'}) {
			print STDERR Dumper($map_read);
			my $error = "Readpair does not have a RG field and no default PI has been specified\n\n";
			$error .= join "\t", @{$map_read->sam()};
			$error .= "\n";
			$error .= join "\t", @{$unmap_read->sam()};
			$error .= "\n";
			croak $error;
		}
		else {
			$pi_val = $options->{'pi'};
		}

		if(!$pi_val) {
			if($options->{'rgs'}->{$map_read->read_group()}) {
				$pi_val = $options->{'rgs'}->{$map_read->read_group()};
			}
			elsif($options->{'rgs'}) {
				my $error = 'Header @RG line does not include a PI value for this readpairs readgroup'."\n\n";
				$error .= join "\t", @{$map_read->sam()};
				$error .= "\n";
				$error .= join "\t", @{$unmap_read->sam()};
				$error .= "\n";
				croak $error;
			}
		}

		if($map_read->sam()->[2] eq $unmap_read->sam()->[2]
		&& $map_read->sam()->[3] != $unmap_read->sam()->[3]
		&& abs($map_read->sam()->[8]) < $pi_val) {
			# fix for reads that map closely together
			if($map_read->strand()) {
				$length -= $pi_val;
			}
			else {
				$length += $pi_val;
			}
		}

		my $read_num;
		if(($unmap_read->sam()->[1] | $first_in_pair_mask) == $unmap_read->sam()->[1]) {
			$read_num = 1;
		}
		elsif(($unmap_read->sam()->[1] | $second_in_pair_mask) == $unmap_read->sam()->[1]) {
			$read_num = 2;
		}
		if($unmap_read->sam()->[1] == $map_read->sam()->[1]) {
			if($read_num == 1) {
				$read_num = 2;
			}
			else {
				$read_num = 1;
			}
		}
		$pindel_rec = '@'.$unmap_read->sam()->[0].'/'.$read_num."\n"; # the readname of unmapped record
		$pindel_rec .= $sequence."\n"; # sequence of unmapped read
		if($map_read->strand()) {
			$pindel_rec .= '+';
		}
		else {
			$pindel_rec .= '-';
		}
		$pindel_rec .= "\t".$map_read->sam()->[2]."\t"; # refseq name

		my $calc_pos = ($map_read->sam()->[3] + $length);
		if($calc_pos < 1) {
			$calc_pos = 1;
		}

		$pindel_rec .= $calc_pos."\t".$map_read->sam()->[4]."\t"; # position, quality

		$pindel_rec .= $pi_val."\t"; # predicted insert size
		$pindel_rec .= $options->{'s'}."\n"; # sample tag
	}
	#else {
		#debug($unmap_read->sam()->[0].' more than allowed Ns in trimmed sequence (allowing max '.($options->{'mn'}*100).'%) or hard limit of remaning sequence reached - '.$options->{'ms'});
	#}

	if($pindel_rec) {
		@extended_rec = ($pindel_rec, $map_read->sam()->[2], $map_read->sam()->[3]);
	}

  return \@extended_rec;
}

sub tag_finder {
  my ($tag, $comp_ref, $start_at) = @_;
  my $tag_pos;
  for my $i($start_at..@{$comp_ref}-1) {
    my @tag = split ":", $comp_ref->[$i];
    if($tag[0] eq $tag) {
      $tag_pos = $i;
      last;
    }
  }
  return $tag_pos;
}

sub qual_profile {
  my ($options) = @_;
  my %group_profiles;
  #$options->{'rgs'}->{$id};
  #map{ord($_)- 33 }

  foreach my $rg_id(keys %{$options->{'rgs'}}) {
    $group_profiles{$rg_id}->{'qual_sum'} = 0;
    $group_profiles{$rg_id}->{'base_count'} = 0;

    # must be mapped and in specified rg
    my $command = find_prog('samtools').' view -F 4 -r '.$rg_id.' '.$options->{'i'}; #.' | cut -f 11 | head -n 500000';
    my $pid = open my $PROC, '-|', $command or croak "Could not fork: $OS_ERROR";
    my $count = 0;
    while( my $sam = <$PROC> ) {
      chomp $sam;
      my $tmp = (split /\t/,$sam)[10];
      #print $tmp,"\n";
      my $quals = qual_to_vals($tmp);

      #print join ',',@quals;
      #print "\n\n";

      foreach my $qual(@{$quals}) {
        $group_profiles{$rg_id}->{'qual_sum'} += $qual;
        $group_profiles{$rg_id}->{'base_count'}++;
      }

      $count++;
      last if($count >= 500_000);
    }
    close $PROC;# or croak 'Failed to close process, '.$command;
    if($count < 250_000) {
      warn 'Insufficient reads in RG '.$rg_id.' to create a profile for quality trimming, defaulting to 8';
      $group_profiles{$rg_id}->{'qual_thresh'} = 8;
    }
    else {
      $group_profiles{$rg_id}->{'qual_thresh'} = ($group_profiles{$rg_id}->{'qual_sum'} / $group_profiles{$rg_id}->{'base_count'}) * 0.3;
    }
  }

  print Dumper(\%group_profiles);
  $options->{'rg_profile'} = \%group_profiles;
  return 1;
}

sub qual_to_vals {
  my ($qual) = @_;
  my @vals = map{ord($_)- 33 } split //, $qual;
  return \@vals;
}

sub option_builder {
	my ($factory) = @_;

	my %opts = ();

  my $result = GetOptions (
    'i|input=s'   => \$opts{'i'},
    'o|output=s'  => \$opts{'o'},
    's|sample=s'  => \$opts{'s'},
    'c|centromeres=s'  => \$opts{'c'},
    'om|oldmethod'  => \$opts{'om'},
    'nm|newmethod'  => \$opts{'nm'},
    'pi|insert=n'  => \$opts{'pi'},
    'r|restrict=s'  => \$opts{'r'},
    'rg|readgroup=s'  => \$opts{'rg'},
    'd|debug'   => \$opts{'d'},
    'h|help'    => \$opts{'h'},
    'mn|max_n=f'    => \$opts{'mn'},
    'ms|min_seq=n'    => \$opts{'ms'},
    'e|edits=n'    => \$opts{'e'},
    'mr|max_rec=n'    => \$opts{'mr'},
    'q|qual_trim' => \$opts{'q'},
	);
	if(!$result || !$opts{'i'} || !$opts{'o'} || !$opts{'s'} || $opts{'h'}) {
		if($opts{'h'}) {
			licence();
		}
  	usage();
	}
  if($opts{'d'}) {
    $debug = 1;
  }

  if(!$opts{'mn'}) {
  	$opts{'mn'} = 0.1;
  }
  if(!$opts{'ms'}) {
  	$opts{'ms'} = 22;
  }
  if(!$opts{'e'}) {
  	$opts{'e'} = 20;
  }

  if(!$opts{'mr'}) {
  	$opts{'mr'} = 200000;
  }

  if($opts{'rg'}) {
    # need to modify the output path
    my @bits = split /\//, $opts{'o'};
    unshift @bits, '.' if(@bits == 1);
    my $stub = pop @bits;
    my $full_path = join '/', @bits;
    $full_path .= '/'.$opts{'rg'}.'/';
    make_path($full_path);
    $full_path .= $stub;
    $opts{'o'} = $full_path;
  }

  $file_count_loc = $opts{'o'}.'_write.count';

	return \%opts;
}

## UTILITIES
sub debug {
  my @data = @_;
  if($debug) {
    foreach my $item(@data) {
      warn $item."\n";
    }
  }
  return;
}

sub bam_sort_order_error {
  my ($file) = @_;
  print STDERR "Error while processing input: $file\n";
  print STDERR <<'BAM_SORT_ORDER_MESSAGE';

Input file is not 'queryname' sorted, explained below.

No preference for parsing method has been defined so the default of 'samtools view'
was applied.  Inorder to use this method your input file should be sorted by queryname.

If you are confident that the file is appropriately sorted (but missing the appropriate
header tag) you can suppress this error by adding '-om' to the command.

  ------------------
  ALTERNATIVE OPTION
  ------------------
  '-nm' option can be applied if you have 'samgroupbyname' installed.
  At time of writing this is not publically available, however please check with the
  authors if this has been released in the interim (-h for contact details).

BAM_SORT_ORDER_MESSAGE
  exit 1;
}

sub licence {
	print <<'LICENCE_DOC';

########################################################################################################
#                                                                                                      #
# CGP Software License                                                                                 #
#                                                                                                      #
# Copyright (c) 2010 Genome Research Ltd.                                                              #
# Author: Cancer Genome Project, cgpit@sanger.ac.uk                                                    #
#                                                                                                      #
# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT   #
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND               #
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES #
# OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN  #
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                           #
#                                                                                                      #
# This code is free software; you can redistribute it and/or modify it under the terms of the BSD      #
# License.                                                                                             #
#                                                                                                      #
# Any redistribution or derivation in whole or in part including any substantial portion of this code  #
# must include this copyright and permission notice.                                                   #
#                                                                                                      #
########################################################################################################
LICENCE_DOC
}

sub usage {
	print <<'USAGE_DOC';

bam2pindel.pl

Description - Extract readpairs for use by Pindel

WARNING 1. This tool was designed for BWA based BAM/SAM data
        2. You must prepare a name sorted bam file
             samtools sort -n -m 4000000000 ORIGINAL.bam NEW_namesorted
              - NOTE: the above command will need 8GB of RAM

Required inputs...
    -i|input        : Input BAM file (req)
    -o|output       : Output filtered ready for pindel
                       (prefix only, _<ref>.txt will be appended, based on anchor read) (req)
    -s|sample       : Sample or label (sampA,sampB...) (req)

Optional inputs...
    -om|oldmethod   : Switch - Use samtools as converter
    -nm|newmethod   : Switch - use samgroupbyname to group automatically without need for namesorted bam
    -c|centromeres  : Location of file listing centromere coordinates (gff3 format)
                       - creates a p and q arm output file for none MT chromosomes
    -pi|insert      : Required if BAM file does not have PI tag in header RG record
                       - See SAM format specification
    -q|qual_trim    : Trim back reads where quality drops and continues to end below 30% of the average base quality for the readgroup.
                       - Especially useful for dealing with sequence interaction effects.
    -r|restrict     : Restrict to chromosome xx
                       - You may get output for other references when each end of a readpair maps to a different reference

Tuning inputs (default when not specified)...
    -mn|max_n       : Max fraction of read seq that can be N after runs of N trimmed from each end (0.1).
    -ms|min_seq     : Min bases of read seq allowed after runs of N trimmed from each end (22).
    -e|edits        : 1 edit allowed for every x bps in total read length (20).

Other...
    -rg|readgroup   : Valid only with '-nm', instruct samgroupbyname to only parse a particular readgroup
		-mr|max_rec     : Max records to cache before writing to file (default is 200,000 - requires ~150MB RAM)
    -d|debug        : The reason why each readpair is discarded will be printed to STDERR
                       - Output is very large, recommend pairing with '-r'
    -h|help         : Print this message

Outputs - A file for each reference (split by centromere if requested) with fields relevant for pindel

Example...
  ./bam2pindel_bwa.pl -i namesorted.bam -o ./pindel_input/output_prefix -s tumour -om

Author : Keiran M Raine (cgpit@sanger.ac.uk)

USAGE_DOC

  exit 0;
}
