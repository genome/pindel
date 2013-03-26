# Quick Start

<style type="text/css">
.scrollable {
  width: 585px;
  overflow: auto;
  white-space: pre;
}
</style>

***

This is a very brief introduction on how to run Pindel.

We will use the fasta file hs_ref_chr20.fa as the reference genome. The reads that we want to compare against the reference are in a Pindel format file called COLO-829_20-p_ok. We want Pindel to write the output files to a separate directory called output (so before calling Pindel, we can create the directory: <code>mkdir output</code>).

To run Pindel, we can use the following command on the data in the demo folder:

<pre class="terminal scrollable">
mkdir output
./pindel -f demo/hs_ref_chr20.fa -p demo/COLO-829_20-p_ok.txt -c 20 -o output/ref
./pindel -f demo/simulated_reference.fa -i demo/simulated_config.txt -c ALL -o output/simulated
</pre>

The output directory will now look like:
<pre class="terminal scrollable">
-rw-r--r-- 1 root root  17961311 2011-05-23 09:54 ref_D
-rw-r--r-- 1 root root     40203 2011-05-23 09:54 ref_TD
-rw-r--r-- 1 root root     38546 2011-05-23 09:54 ref_INV
-rw-r--r-- 1 root root  16364718 2011-05-23 09:54 ref_SI
-rw-r--r-- 1 root root     22946 2011-05-23 09:54 ref_LI
-rw-r--r-- 1 root root    285266 2011-05-23 09:54 ref_BP
</pre>

In general, pindel is used as follows:

<pre class="terminal scrollable">
  ./pindel -f &lt;reference.fa&gt; -p &lt;pindel_input&gt; [and/or -i bam_configuration_file] -c &lt;chromosome_name&gt; -o &lt;prefix_for_output_files&gt;
</pre>
If you wish to use Pindel directly on a BAM-file, instead of first converting the BAM-file into a Pindel input file with bam2pindel, you need to take the following steps:

1) create a bam-configuration file, consisting of one line per BAM-file that you want Pindel to process. Every line should contain the name of a bam-file, the expected average insert size, and a label to indicate the identity of the sample. For example:

<pre class="terminal scrollable">
tumor_sample_1222.bam 250 TUMOR_1222
somatic_sample_1222.bam 250 HEALTHY_1222
</pre>

2) run pindel with the -i option, so

<pre class="terminal scrollable">
./pindel -f hs_ref_GRCh37.fa -i 1222config.txt -c ALL -o sample_1222
</pre>

The output files should then contain all detected indels and SVs relative to the reference, and the labels will indicate in which samples each indel and SV occurred.

Pindel has many parameters that you can set to increase speed (multithreading) or change the balance between sensitivity and specificity. The parameters are:

#####Required parameters

<pre class="terminal scrollable">
-f/--fasta               the reference genome sequences in fasta format
-p/--pindel-file         the Pindel input file; (either this or a bam configuration file is required).
-i/--config-file         the bam config file; either this or a pindel input file is required. Per line: path and file name of bam, insert 
                            size and sample tag. For example: /data/tumour.bam  400  tumour
-o/--output-prefix       Output prefix
-c/--chromosome          Which chr/fragment. Pindel will process reads for one chromosome each time. ChrName must be the same as in reference 
                            sequence and in read file. '-c ALL' will make Pindel loop 
                            over all chromosomes. The search for indels and SVs can also be limited to a specific region; 
                                -c 20:10,000,000 will only look for indels and SVs after position 10,000,000 == [10M, end], 
                                -c 20:5,000,000-15,000,000 will report indels in the range between and including the bases at 
                                  position 5,000,000 and 15,000,000 = [5M, 15M]	
</pre>

#####Parameters affecting runtime and memory usage

<pre class="terminal scrollable">
-T/--number_of_threads   the number of threads Pindel will use (default 1). More threads assures lower runtime, but requires 
                            multiple processors
-w/--window_size         for saving RAM, divides the reference in bins of X million bases and only analyzes the reads per bin 
                            (default 10 (=10 million)). A smaller bin size will reduce memory but will increase runtime slightly.
</pre>

#####Parameters affecting which structural variants are reported

<pre class="terminal scrollable">
-x/--max_range_index             the maximum size of structural variations to be detected; the higher this number, the greater the 
                                    number of SVs reported, but the computational cost and memory requirements increase, as does the 
                                    rate of false positives. 1=128, 2=512, 3=2,048, 4=8,092, 5=32,368, 6=129,472, 7=517,888, 8=2,071,552, 
                                    9=8,286,208 (maximum 9, default 5)
-r/--report_inversions           report inversions (default true)
-t/--report_duplications         report tandem duplications (default true)
-l/--report_long_insertions      report insertions of which the full sequence cannot be deduced because of their length (default true)
-k/--report_breakpoints          report breakpoints (default true)
-s/--report_close_mapped_reads   report reads of which only one end (the one closest to the mapped read of the paired-end read) could 
                                    be mapped (default false)
-n/--min_NT_size                 only report inserted (NT) sequences in deletions greater than this size (default 50)
-v/--min_inversion_size          only report inversions greater than this number of bases (default 50)
</pre>

#####Parameters affecting sensitivity and selectivity

<pre class="terminal scrollable">
-d/--min_num_matched_bases           only consider reads as evidence if they map with more than this number of bases to the reference (default 30)
-a/--additional_mismatch             Pindel will only map part of a read to the reference genome if there are no other candidate positions 
                                        with no more than the specified number of mismatches position. The bigger this value, the more accurate 
                                        but less sensitive. (default value 1)
-m/--min_perfect_match_around_BP     at the point where the read is split into two, there should at least be this number of perfectly matching bases 
                                        between read and reference (default value 3)
-e/--sequencing_error_rate           the expected fraction of sequencing errors (default 0.05)
-u/--maximum_allowed_mismatch_rate   only reads with fewer mismatches with the reference genome than this fraction will be considered (default 0.1)
</pre>

#####Miscellaneous parameters
<pre class="terminal scrollable">
-b/--breakdancer         [file name]. Pindel is able to use calls from other SV methods such as BreakDancer to further increase sensitivity and specificity.  
                            BreakDancer result or calls from any methods must in the format:   ChrA LocA stringA ChrB LocB stringB other
-Q                       [file name] The list of BreakDancer calls with Pindel support information. 
                            Format: chr   Loc_left   Loc_right   size   type   index
                            For example, "1	72766323 	72811840 	45516	D	11970" means the deletion event chr1:72766323-72811840 of size 45516 is 
                            reported as an event with index 11970 in Pindel report of deletion. 
-h/--help                show the command line options of Pindel
</pre>
