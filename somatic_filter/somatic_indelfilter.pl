#!/usr/bin/perl
# INDEL filter for GenomeVIP 
#----------------------------------
# $Authors: Beifang Niu $ Jay Mash 
# $Date: Wed Jan 22 12:59:27 CST 2014 ( The Jan 22 12:59:27 CST 2014 ) $ 
# $Revision:  $
# $URL: $
#----------------------------------
use strict;
use warnings;

use Cwd;
use Carp;
use FileHandle;
use IO::File;
use Getopt::Long;
use POSIX qw( WIFEXITED );
use File::Temp qw/ tempfile /;

# get paras from config file
my %paras; 
map { chomp; @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } (<>);
#map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;
#
# file exist ? 
#
unless ( -e $paras{'input'} ) { die "input indels not exist !!! \n"; }
my $input_fh = IO::File->new( $paras{'input'} ) or die " could not open input indel file for reading $! ";

# create temp file
#
# no complex temp file
#
my ( undef, $nocomplex_output ) = tempfile();
my $nocomplex_output_fh = IO::File->new( $nocomplex_output, ">" ) or die "Temporary file could not be created. $!";
map { chomp; my @t = split;
    if ($t[32] + $t[34] + $t[36] >= $paras{'cov'} && $t[33] + $t[34] + $t[36] >= $paras{'cov'} && $t[39] + $t[41] + $t[43] >= $paras{'cov'} && $t[40] + $t[41] + $t[43] >= $paras{'cov'}) {
        if ( ($t[34] + $t[36])  == 0 && ($t[41] + $t[43])/($t[39] + $t[41] + $t[43]) >= $paras{'vaf'} && (($t[41] + $t[43])/($t[40] + $t[41] + $t[43]) >= $paras{'vaf'} ) ) {
            #print; print "\n";
            # no complex
            if ( $t[1] eq "I" || ($t[1] eq "D" && $t[4] == 0) ) {
                $nocomplex_output_fh->print($_."\n");
            }
        }
    }
} <$input_fh>;
$input_fh->close;
$nocomplex_output_fh->close; 
#
# pindel2vcf output temp file
#
my ( undef, $pindel2vcf_output ) = tempfile();
#my $pindel2vcf_output_fh = IO::File->new( $pindel2vcf_output, ">" ) or die "Temporary file could not be created. $!";
#
# run pindel2vcf
#
my $pindel2vcf_command = "";
$pindel2vcf_command = "$paras{'pindel2vcf'} -R $paras{'referencename'} -r $paras{'reference'} -p $nocomplex_output -d $paras{'referencedate'} -v $pindel2vcf_output";
print $pindel2vcf_command."\n";
system( $pindel2vcf_command );

my $filter_output_fh = IO::File->new( $paras{'output'}, ">" ) or die "Filtered output file could not be created. $!";
my $pindel2vcf_output_fh = IO::File->new( $pindel2vcf_output ) or die "Temporary file could not be opened. $!";
while ( <$pindel2vcf_output_fh> ) {
    print;
    if ( /^#/ ) { $filter_output_fh->print($_); next };
    my @a= split /\t/; 
    my @b = split/\;/, $a[7]; 
    for ( my $i=0; $i<scalar(@b); $i++) { 
        if ( $b[$i]=~/^HOMLEN/ ) { 
            my @c = split/=/, $b[$i]; 
            if ( $c[1] <= $paras{'hom'} ) { $filter_output_fh->print($_); last; } 
        } 
    }
} <$pindel2vcf_output_fh>;
$pindel2vcf_output_fh->close;
$filter_output_fh->close;

