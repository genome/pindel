# Manual Install From Source

<style type="text/css">
.scrollable {
  width: 585px;
  overflow: auto;
  white-space: pre;
}
</style>

***

##step 1: get the code
<pre class="scrollable terminal"> 
git clone git://github.com/genome/pindel.git
</pre>

##step 2: compile source code
<pre class="scrollable terminal"> 
cd pindel 
./INSTALL /path/to/samtools_FOLDER/
</pre>

##step 3: test the code
<pre class="scrollable terminal"> 
cd demo
../pindel -f simulated_reference.fa -i simulated_config.txt -o output
</pre>






