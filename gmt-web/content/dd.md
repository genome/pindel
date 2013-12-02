# Detection of dispersed duplications (beta)

As of version 0.2.5, pindel supports calling dispersed duplication insertions (DD).  DD includes all duplications like mobile element insertions and segmental duplications as long as there exists a copy in the provided reference genome.  This module of pindel focuses on _dispersed_ duplications, tandem duplications are typically not included, but the locality of which events are called can be controlled with the --MIN_DD_MAP_DISTANCE parameter.

## Running pindel with dispersed duplication detection

One can run pindel with the usual mandatory parameters and use "-q" to turn on detection of DD.  Pindel then collects evidence for DD events in the form of discordant read pairs (where one end maps inside the duplicated segment and one just outside) and split reads (where one read end covers the breakpoint of the insertion on the reference).

All parameters related to DD detection:

-q/--detect_DD
           Flag indicating whether to detect dispersed duplications. (default: 
           false) 

/--MAX_DD_BREAKPOINT_DISTANCE
           Maximum distance between dispersed duplication breakpoints to assume 
           they refer to the same event. (default: 350) 

/--MAX_DISTANCE_CLUSTER_READS
           Maximum distance between reads for them to provide evidence for a 
           single breakpoint for dispersed duplications. (default: 100) 

/--MIN_DD_CLUSTER_SIZE
           Minimum number of reads needed for calling a breakpoint for dispersed 
           duplications. (default: 3) 

/--MIN_DD_BREAKPOINT_SUPPORT
           Minimum number of split reads for calling an exact breakpoint for 
           dispersed duplications. (default: 3) 

/--MIN_DD_MAP_DISTANCE
           Minimum mapping distance of read pairs for them to be considered 
           discordant. (default: 8000) 

/--DD_REPORT_DUPLICATION_READS
           Report discordant sequences and positions for mates of reads mapping 
           inside dispersed duplications. (default: false)


## Output file format

DD events are reported in a file post-fixed with "_DD".  An example output for DD detection looks like this:

    ####################################################################################################
    1	DD	reference	7068	7069	48	19	5	19	5
    # Dispersed Duplication insertion (DD) found on chromosome 'reference', breakpoint at 7068 (estimated from + strand), 7069 (estimated from - strand)
    # Found 48 supporting reads, of which 19 discordant reads and 5 split reads at 5' end, 19 discordant reads and 5 split reads at 3' end.
    # Supporting reads for insertion location (5' end):
    # Reference: TCGCCTATCTCACGATCGCCTCAATGCACCCGACGATAGGGCTCCCGTTGACCTTCAACAGCTTCGGTGGCTACTAGATACTCtattaaagggtcattggcgaaaaggcatagttgccgagggctcatggaagccagattcttcgtagattacacgacacagttcgc
    #            TCGCCTATCTCACGATCGCCTCAATGCACCCGACGATAGGGCTCCCGTTGACCTTCAACAGCTTCGGTGGCTACTAGATACTCCAATCCTGGCTAATCTC (name: @read_393/2 sample: sample1) 
    #                             GCCTCAATGCACCCGACGATAGGGCTCCCGTTGACCTTCAACAGCTTCGGTGGCTACTAGATACTCCAATCCTGGCTAATCTCTCATACCGGCACCGCTC (name: @read_394/2 sample: sample1) 
    #                                              GATAGGGCTCCCGTTGACCTTCAACAGCTTCGGTGGCTACTAGATACTCCAATCCTGGCTAATCTCTCATACCGGCACCGCTCTGTCGGTCGCGAAATGC (name: @read_395/2 sample: sample1) 
    #                                                              ACCTTCAACAGCTTCGGTGGCTACTAGATACTCCAATCCTGGCTAATCTCTCATACCGGCACCGCTCTGTCGGTCGCGAAATGCAACGCCCACGTTATGG (name: @read_396/2 sample: sample1) 
    #                                                                               TGGCTACTAGATACTCCAATCCTGGCTAATCTCTCATACCGGCACCGCTCTGTCGGTCGCGAAATGCAACGCCCACGTTATGGTGGGAGGCTTCCGCAGC (name: @read_397/2 sample: sample1) 
    # Supporting reads for insertion location (3' end):
    # Reference: atctcacgatcgcctcaatgcacccgacgatagggctcccgttgaccttcaacagcttcggtggctactagatactcTATTAAAGGGTCATTGGCGAAAAGGCATAGTTGCCGAGGGCTCATGGAAGCCAGATTCTTCGTAGATTACACGACACAGTTCGCCACAGC
    #                                                                               TCGCGGCATTTATTAAAGGGTCATTGGCGAAAAGGCATAGTTGCCGAGGGCTCATGGAAGCCAGATTCTTCGTAGATTACACGACACAGTTCGCCACAGC (name: @read_457/1 sample: sample1) 
    #                                                              TGTTCCCCACACAGCGCTCGCGGCATTTATTAAAGGGTCATTGGCGAAAAGGCATAGTTGCCGAGGGCTCATGGAAGCCAGATTCTTCGTAGATTACACG (name: @read_456/1 sample: sample1) 
    #                                             ATAGGATTGGCTCAAACTGTTCCCCACACAGCGCTCGCGGCATTTATTAAAGGGTCATTGGCGAAAAGGCATAGTTGCCGAGGGCTCATGGAAGCCAGAT (name: @read_455/1 sample: sample1) 
    #                            ATCCAGCTGGTGTTAATATAGGATTGGCTCAAACTGTTCCCCACACAGCGCTCGCGGCATTTATTAAAGGGTCATTGGCGAAAAGGCATAGTTGCCGAGG (name: @read_454/1 sample: sample1) 
    #            TGACCCTCTATCTCAAATCCAGCTGGTGTTAATATAGGATTGGCTCAAACTGTTCCCCACACAGCGCTCGCGGCATTTATTAAAGGGTCATTGGCGAAAA (name: @read_453/1 sample: sample1) 
    # All supporting sequences for this insertion (i.e. sequences that map inside the inserted element):
    ?	?	?	@read_457/1	sample1	-	TCGCGGCATT
    ?	?	?	@read_456/1	sample1	-	TGTTCCCCACACAGCGCTCGCGGCATT
    ?	?	?	@read_455/1	sample1	-	ATAGGATTGGCTCAAACTGTTCCCCACACAGCGCTCGCGGCATT
    ?	?	?	@read_454/1	sample1	-	ATCCAGCTGGTGTTAATATAGGATTGGCTCAAACTGTTCCCCACACAGCGCTCGCGGCATT
    ?	?	?	@read_453/1	sample1	-	TGACCCTCTATCTCAAATCCAGCTGGTGTTAATATAGGATTGGCTCAAACTGTTCCCCACACAGCGCTCGCGGCATT
    ?	?	?	@read_393/2	sample1	+	CAATCCTGGCTAATCTC
    ?	?	?	@read_394/2	sample1	+	CAATCCTGGCTAATCTCTCATACCGGCACCGCTC
    ?	?	?	@read_395/2	sample1	+	CAATCCTGGCTAATCTCTCATACCGGCACCGCTCTGTCGGTCGCGAAATGC
    ?	?	?	@read_396/2	sample1	+	CAATCCTGGCTAATCTCTCATACCGGCACCGCTCTGTCGGTCGCGAAATGCAACGCCCACGTTATGG
    ?	?	?	@read_397/2	sample1	+	CAATCCTGGCTAATCTCTCATACCGGCACCGCTCTGTCGGTCGCGAAATGCAACGCCCACGTTATGGTGGGAGGCTTCCGCAGC
    reference	136603	-	@read_452/1	sample1	-	TTAATAAATGCCGCGAGCGCTGTGTGGGGAACAGTTTGAGCCAATCCTATATTAACACCAGCTGGATTTGAGATAGAGGGTCAATCGGGTGCCCTGTGAC
    reference	136620	-	@read_451/1	sample1	-	CGCTGTGTGGGGAACAGTTTGAGCCAATCCTATATTAACACCAGCTGGATTTGAGATAGAGGGTCAATCGGGTGCCCTGTGACCCCGTAGCATGGGCATA
    reference	136637	-	@read_450/1	sample1	-	TTTGAGCCAATCCTATATTAACACCAGCTGGATTTGAGATAGAGGGTCAATCGGGTGCCCTGTGACCCCGTAGCATGGGCATAGGTAAGCTGAGCCTCAT
    reference	136654	-	@read_449/1	sample1	-	TTAACACCAGCTGGATTTGAGATAGAGGGTCAATCGGGTGCCCTGTGACCCCGTAGCATGGGCATAGGTAAGCTGAGCCTCATCGTCCGAACTTCCGTCA
    reference	136670	-	@read_448/1	sample1	-	TTGAGATAGAGGGTCAATCGGGTGCCCTGTGACCCCGTAGCATGGGCATAGGTAAGCTGAGCCTCATCGTCCGAACTTCCGTCAGGATAAAGGCTGGAAG
    reference	136687	-	@read_447/1	sample1	-	TCGGGTGCCCTGTGACCCCGTAGCATGGGCATAGGTAAGCTGAGCCTCATCGTCCGAACTTCCGTCAGGATAAAGGCTGGAAGAAGTTCAGGTTCGCTAG
    reference	136704	-	@read_446/1	sample1	-	CCGTAGCATGGGCATAGGTAAGCTGAGCCTCATCGTCCGAACTTCCGTCAGGATAAAGGCTGGAAGAAGTTCAGGTTCGCTAGTGCGGGGAGAAGCGTTC
    reference	136721	-	@read_445/1	sample1	-	GTAAGCTGAGCCTCATCGTCCGAACTTCCGTCAGGATAAAGGCTGGAAGAAGTTCAGGTTCGCTAGTGCGGGGAGAAGCGTTCTTCGGCCCAACTAGGAC
    reference	136737	-	@read_444/1	sample1	-	CGTCCGAACTTCCGTCAGGATAAAGGCTGGAAGAAGTTCAGGTTCGCTAGTGCGGGGAGAAGCGTTCTTCGGCCCAACTAGGACTCCTCGTTAACTGCCG
    reference	136754	-	@read_443/1	sample1	-	GGATAAAGGCTGGAAGAAGTTCAGGTTCGCTAGTGCGGGGAGAAGCGTTCTTCGGCCCAACTAGGACTCCTCGTTAACTGCCGTGCCTCTTTGATTTTTA
    reference	136771	-	@read_442/1	sample1	-	AGTTCAGGTTCGCTAGTGCGGGGAGAAGCGTTCTTCGGCCCAACTAGGACTCCTCGTTAACTGCCGTGCCTCTTTGATTTTTATGACGCTGAGAGGCTCG
    reference	136788	-	@read_441/1	sample1	-	GCGGGGAGAAGCGTTCTTCGGCCCAACTAGGACTCCTCGTTAACTGCCGTGCCTCTTTGATTTTTATGACGCTGAGAGGCTCGATGATCACTCATATGTC
    reference	136804	-	@read_440/1	sample1	-	TTCGGCCCAACTAGGACTCCTCGTTAACTGCCGTGCCTCTTTGATTTTTATGACGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGG
    reference	136807	+	@read_416/2	sample1	+	GGCCCAACTAGGACTCCTCGTTAACTGCCGTGCCTCTTTGATTTTTATGACGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGGTGG
    reference	136821	-	@read_439/1	sample1	-	TCCTCGTTAACTGCCGTGCCTCTTTGATTTTTATGACGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCC
    reference	136823	+	@read_415/2	sample1	+	CTCGTTAACTGCCGTGCCTCTTTGATTTTTATGACGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCCGC
    reference	136838	-	@read_438/1	sample1	-	GCCTCTTTGATTTTTATGACGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTG
    reference	136840	+	@read_414/2	sample1	+	CTCTTTGATTTTTATGACGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTGCA
    reference	136855	-	@read_437/1	sample1	-	GACGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGT
    reference	136857	+	@read_413/2	sample1	+	CGCTGAGAGGCTCGATGATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGTGT
    reference	136871	-	@read_436/1	sample1	-	ATGATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGTGTGGCGTATGGCTCGC
    reference	136874	+	@read_412/2	sample1	+	ATCACTCATATGTCCGACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGTGTGGCGTATGGCTCGCTTC
    reference	136888	-	@read_435/1	sample1	-	CGACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGTGTGGCGTATGGCTCGCTTCAGGCCTGAGCAAGC
    reference	136890	+	@read_411/2	sample1	+	ACGTTGCCACAAGGTGGCTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGTGTGGCGTATGGCTCGCTTCAGGCCTGAGCAAGCCG
    reference	136905	-	@read_434/1	sample1	-	GGCTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGTGTGGCGTATGGCTCGCTTCAGGCCTGAGCAAGCCGAGCACCGTCACAATC
    reference	136907	+	@read_410/2	sample1	+	CTAGATCATTTCCCGCACGCAGGTCATATTGCATCGTGTGCCAGTAGTGTGGCGTATGGCTCGCTTCAGGCCTGAGCAAGCCGAGCACCGTCACAATCAA
    reference	136924	+	@read_409/2	sample1	+	CGCAGGTCATATTGCATCGTGTGCCAGTAGTGTGGCGTATGGCTCGCTTCAGGCCTGAGCAAGCCGAGCACCGTCACAATCAATTGCAGTACAAAATTCG
    reference	136941	+	@read_408/2	sample1	+	CGTGTGCCAGTAGTGTGGCGTATGGCTCGCTTCAGGCCTGAGCAAGCCGAGCACCGTCACAATCAATTGCAGTACAAAATTCGTGACCGGTCGTCGTATC
    reference	136957	+	@read_407/2	sample1	+	GGCGTATGGCTCGCTTCAGGCCTGAGCAAGCCGAGCACCGTCACAATCAATTGCAGTACAAAATTCGTGACCGGTCGTCGTATCACATGGAGCTGTAATG
    reference	136974	+	@read_406/2	sample1	+	AGGCCTGAGCAAGCCGAGCACCGTCACAATCAATTGCAGTACAAAATTCGTGACCGGTCGTCGTATCACATGGAGCTGTAATGAGCCGAATCGGTAGCAG
    reference	136991	+	@read_405/2	sample1	+	GCACCGTCACAATCAATTGCAGTACAAAATTCGTGACCGGTCGTCGTATCACATGGAGCTGTAATGAGCCGAATCGGTAGCAGTAGCGCTATCCAGGGTC
    reference	137008	+	@read_404/2	sample1	+	TGCAGTACAAAATTCGTGACCGGTCGTCGTATCACATGGAGCTGTAATGAGCCGAATCGGTAGCAGTAGCGCTATCCAGGGTCTCAGACGACCCCACAAC
    reference	137024	+	@read_403/2	sample1	+	TGACCGGTCGTCGTATCACATGGAGCTGTAATGAGCCGAATCGGTAGCAGTAGCGCTATCCAGGGTCTCAGACGACCCCACAACACTCAACGACGACTGA
    reference	137041	+	@read_402/2	sample1	+	ACATGGAGCTGTAATGAGCCGAATCGGTAGCAGTAGCGCTATCCAGGGTCTCAGACGACCCCACAACACTCAACGACGACTGATGCTGCGGAAGCCTCCC
    reference	137058	+	@read_401/2	sample1	+	GCCGAATCGGTAGCAGTAGCGCTATCCAGGGTCTCAGACGACCCCACAACACTCAACGACGACTGATGCTGCGGAAGCCTCCCACCATAACGTGGGCGTT
    reference	137075	+	@read_400/2	sample1	+	AGCGCTATCCAGGGTCTCAGACGACCCCACAACACTCAACGACGACTGATGCTGCGGAAGCCTCCCACCATAACGTGGGCGTTGCATTTCGCGACCGACA
    reference	137091	+	@read_399/2	sample1	+	TCAGACGACCCCACAACACTCAACGACGACTGATGCTGCGGAAGCCTCCCACCATAACGTGGGCGTTGCATTTCGCGACCGACAGAGCGGTGCCGGTATG
    reference	137108	+	@read_398/2	sample1	+	ACTCAACGACGACTGATGCTGCGGAAGCCTCCCACCATAACGTGGGCGTTGCATTTCGCGACCGACAGAGCGGTGCCGGTATGAGAGATTAGCCAGGATT

DD calls are separated with a string of hash characters (#).  Each DD call starts with a line of tab-separated values summarizing the event.  The values are as follows (in this order):

1. Event identification number (integer).
2. Type of event (currently simply "DD" for all events).
3. Sequence name on which DD event is located.
4. Location of DD event as estimated from evidence based on the forward strand.
5. Location of DD event as estimated from evidence based on the reverse strand.
6. Total number of reads (both split and discordant reads) supporting the event.
7. Number of supporting discordant reads on forward strand.
8. Number of supporting split reads on forward strand.
9. Number of supporting discordant reads on reverse strand.
10. Number of supporting split reads on reverse strand.

The next few lines for each event are prefixed with a hash character (#) and show the event summary in a human readable way.  If possible, the supporting split reads are also shown here aligned to the local reference sequence depicting the possible exact breakpoint position of the event.

Finally a number of lines per event are printed that give information on the read ends that (partly) map the duplicated segment's sequence.  These are the following tab-separated values:

1. Sequence name to which the read end was alternatively mapped ("?" for split reads).
2. Location on sequence to which the read end was alternatively mapped ("?" for split reads).
3. Strand to which the read end was alternatively mapped (forward "+", reverse "-", again "?" in case of a split read).
4. The read name.
5. Name of sample where read originated from.
6. Strand to which the mate of this read end mapped (forward "+" or reverse "-").
7. (Part of) the sequence that maps inside the duplicated segment.



