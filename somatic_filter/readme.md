
This is the initial implementation of somatic indel filtering for Pindel indel calls based on ICGC DREAM challenge round 3 test result. Please advise. contact <kaiye@xjtu.edu.cn>

1. Extract indel summary lines using `grep ChrID [pindel output file] > [output file]`, for example:

       grep ChrID pindel_output_D > D.head
       grep ChrID pindel_output_SI > SI.head
          
   You may optionally use `cat` to combine the files.

1. Generate the configuration file for the run. The file requires the following entries, with values and entries separated by an equal sign. An example is provided here as somatic.indel.filter.config.
 
     | Field name | Definition |
     | ---------- | ---------- |
     | `input` | The indel summary file generated in the previous step. |
     | `vaf` | The minimum variant allele frequency at the _tumor_ sample to pass the filter. |
     | `cov` | The minimum coverage at _both_ samples to pass the filter. |
     | `hom` | The maximum homopolymer length at breakpoint to pass the filter. |
     | `pindel2vcf` | The location of pindel's `pindel2vcf` binary. |
     | `reference` | The reference genome. Passed to the `-r` parameter in `pindel2vcf`. |
     | `referencename` | The name of the reference genome. Passed to the `-R` parameter in `pindel2vcf`. |
     | `referencedate` | The date of the reference genome. Passed to the `-d` parameter in `pindel2vcf`. |
     | `output` | The name of the output VCF file. Passed to the `-v` parameter in `pindel2vcf`. |
1. Execute `somatic_indelfilter.pl` with the configuration file
    perl somatic_indelfilter.pl somatic.indel.filter.config

Somatic indel calls will be stored as a VCF file with the name specified in the `output` field in the configuration file.

We assume normal-tumor paired genome runs by `pindel`, and that the normal sample comes first in the pindel output file.