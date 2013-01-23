package Pindel::Bam::Adaptor;

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
use warnings FATAL => 'all';
use Carp;
use English qw( -no_match_vars );

my $q_unmap_mask = 4; #0x0004 (set when not mapped)
my $strand_mask = 16; #0x0010 (set when reverse)

our $regexps;

sub new {
    my ($class, $read_ref) = @_;
    my $self = {
    							'edits' => -1, # must start at -1 as 0 has meaning
    							'mapped' => 0,
									'unique' => 0,
									'repeat' => 0,
									#'insert' => undef,
									'sw' => 0,
									'suboptimal' => 1,
									'rg' => undef,
									'strand' => 1, # as used in tests set 1 for + and 0 for -
    						};
    bless $self, $class;
    if(!$read_ref) {
      croak 'Adaptor should be instantiated with array reference of sam read split on tabs';
    }
    if(!$regexps) {
    	# pre-compiles regular expressions for _tag_finder
    	$regexps = {'RG' => qr/^RG\:/,
									'NM' => qr/^NM\:/,
									'XT' => qr/^XT\:/,
									'X0' => qr/^X0\:/,
									'X1' => qr/^X1\:/,
									};
    }
    $self->_populate($read_ref);
    return $self;
}

sub _populate {
	my ($self, $read_ref) = @_;
	$self->{'sam'} = $read_ref;
	
	$self->{'name'} = $read_ref->[0];
	$self->{'rname'} = $read_ref->[2];
	$self->{'cigar'} = $read_ref->[5];
	
	$self->{'read_length'} = length $read_ref->[9];
	if(($read_ref->[1] | $strand_mask) == $read_ref->[1]) {
		$self->{'strand'} = 0;
	}
	
	my $rg_pos = _tag_finder($regexps->{'RG'}, $read_ref, 11); # fixed tags are first 11 items so skip
	if($rg_pos) {
		my (undef, undef, $rg) = split /:/, $read_ref->[$rg_pos];
		$self->{'rg'} = $rg;
	}
	
	
	if(($read_ref->[1] | $q_unmap_mask) != $read_ref->[1]) {
		$self->{'mapped'} = 1;
		my ($xt, $nm);
		my ($best, $subopt) = (0,0);
		
		(undef, undef, $nm) = split /:/, $read_ref->[_tag_finder($regexps->{'NM'}, $read_ref, 11)];
    $self->{'edits'} = 0 + $nm;
		
		(undef, undef, $xt) = split /:/, $read_ref->[_tag_finder($regexps->{'XT'}, $read_ref, 11)];
		if($xt eq 'U') {
			$self->{'unique'} = 1;
		}
		elsif($xt eq 'M') {
			$self->{'sw'} = 1;
			if($self->{'edits'} <= 2) {
    		# don't think it is realistically possible for a SW alignment to be considered
    		# optimal however this will allow those that only have soft-clipping and minimal
    		# edits to be used (Kai recommended 2)
    		$self->{'suboptimal'} = 0;
    	}
		}
		
		if($xt ne 'M' && $xt ne 'N') {
    	(undef, undef, $best) = split /:/, $read_ref->[_tag_finder($regexps->{'X0'}, $read_ref, 11)];
    	#my $x1_pos = _tag_finder($regexps->{'X1'}, $read_ref, 11);
    	#if($x1_pos) {
	    #	(undef, undef, $subopt) = split /:/, $read_ref->[$x1_pos];
    	#}
			if($self->{'unique'} && $best == 1) {
				# this deals with reads where a 100% match and several <100% matches were found
				$self->{'suboptimal'} = 0;
			}
    }
	}
	return;
}

#sub repeat {
#	return shift->{'repeat'};
#}

#sub insert {
#	return shift->{'insert'};
#}

sub cigar {
	return shift->{'cigar'};
}

sub rname {
	return shift->{'rname'};
}

sub name {
	return shift->{'name'};
}

# 1 if positive, 0 if neg
sub strand {
	return shift->{'strand'};
}

sub read_group {
	return shift->{'rg'};
}

sub read_length {
	return shift->{'read_length'};
}

sub mapped {
	return shift->{'mapped'};
}

sub unique {
	return shift->{'unique'};
}

sub sw {
	return shift->{'sw'};
}

sub edits {
	return shift->{'edits'};
}

sub suboptimal {
	return shift->{'suboptimal'};
}

sub sam {
	return shift->{'sam'};
}

sub _tag_finder {
  my ($reg_ex, $comp_ref, $start_at) = @_;
  for my $i($start_at..@{$comp_ref}-1) {
  	if($comp_ref->[$i] =~ $reg_ex) {
      return $i;
    }
  }
  return undef;
}

1;