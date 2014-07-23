# Frequently Asked Questions

<style type="text/css">
.scrollable {
  width: 585px;
  overflow: auto;
  white-space: pre;
}
</style>

***

##Q: When trying to use bam2pindel.pl to convert bam files into Pindel input I get the error message:
<pre class="scrollable terminal">"make_path" is not exported by the File::Path module
Can't continue after import errors at bam2pindel.pl line 31
BEGIN failed--compilation aborted at bam2pindel.pl line 31.
</pre>

*A:* This error message indicates that you have an old version of Perl, specifically of Perl's File::Path module (File::Path version before 2.07). Installing a newer version of File::Path would help; an alternative is modifying two lines in the bam2pindel.pl file:

Modify "use File::Path qw(make_path);" into "use File::Path qw(mkpath);" and "make_path($full_path);" into "mkpath($full_path);"

## Q:Receiving error messages when compiling Pindel

*A:* First of all, we'd recommend you to download and use the build script. If that fails, please try the following steps:
1. First, we'd recommend that you download the newest [SAMtools source](http://sourceforge.net/projects/samtools/files/samtools/).
2. If, when compiling Pindel, you get the error message "./samtools-0.1.13/bgzf.h:30: fatal error: zlib.h: No such file or directory" then download [zlib](http://zlib.net/). Unpack it, go to the zlib directory, and give the commands "./configure" ; "make test", and "make install". 
3. If when compiling Pindel, you get the error message "/usr/bin/ld: cannot find -lbam" then please go to the SAMtools directory, and give the command "make".
4. If, when making SAMtools, you receive the message "bam_tview.c:5: fatal error: curses.h: No such file or directory" then: either install curses, or, following the instructions in the SAMtools README file, edit the Makefile of SAMtools to comment the LIBCURSES line out (put # before that line), and set CURSES_LIB to 0.

## Q: I get a segmentation fault or bad::alloc error while running Pindel: what's happening?

*A:* Especially in older versions of Pindel (before 0.2.4d) certain parts of Pindel would try allocate big consecutive blocks of memory (say 1G); this can be a problem, even on a computer with 8 GB or more RAM, as other software can 'clutter' parts of the memory, which has the consequence that while the total amount of free RAM is many GB, it is splintered into small pieces all smaller than 1G. Newer versions of Pindel circumvent this problem. The problem can basically be solved in two ways:

* for older versions of pindel: deactivate the search for breakpoints and long insertions by adding the parameters "-l false -k false" to the command line. (in some cases, this will also help for the newer versions)
* for newer versions of pindel (who should generally not have this problem): check your window size; if window sizes are too big (say >100) Pindel may have trouble allocating the memory blocks; try smaller settings, like window size of 10 or so. In general, whenever you get memory problems, setting the window size smaller will help; it will degrade speed a little (10-20%), but in most cases that is not too big of a problem, especially not compared to your program aborting prematurely. 

## Q: What is insert size? 

*A:* Insert size is the length of sequence between the paired-end adapters in paired-end sequence. It is normally longer than the sum of lengths of paired-end reads, if reads are not overlapping. 

See [this SeqAnswers thread for more details](http://seqanswers.com/forums/showthread.php?t=13098).

##Q: Cannot compile Pindel source code 

*A:* If it is caused by warnings, you need to remove -Werror in the make file. You need to specify the samtools FOLDER, not samtools binary, as Pindel needs libraries of samtools.

The correct way to compile is

    ./INSTALL /path/to/samtools_FOLDER/
