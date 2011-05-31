# Makefile for pindel and pindel-debug
# Include the Makefile.local from the main project for configuration
include ../Makefile.local

SOURCES=pindel.cpp reader.cpp reporter.cpp searcher.cpp parameter.cpp \
    refreader.cpp ControlState.cpp search_deletions.cpp search_indels.cpp \
    search_inversions.cpp search_inversions_NT.cpp search_short_insertions.cpp \
    search_tandem_duplications.cpp search_tandem_duplications_NT.cpp \
    outputsorter.cpp FarEndSearcher.cpp searchvariant.cpp searchshortinsertions.cpp searchdeletions.cpp
 
all:	pindel pindel-debug
include $(SOURCES:.cpp=.d)
include $(SOURCES:.cpp=-debug.d)
test:	pindel $(addprefix cppcheck-result-,$(SOURCES:.cpp=.xml))

# Debug version of pindel
pindel-debug: $(SOURCES:.cpp=-debug.o)
	$(CXX) $(LDFLAGS) $^ -fopenmp -g -pg -fprofile-arcs -ftest-coverage \
	  -lm -lz -L$(SAMTOOLS) -lbam -o $@

# Normal version of pindel	
pindel:	$(SOURCES:.cpp=.o)
	$(CXX) $(LDFLAGS) $^ -O3 -fopenmp -lm -lz -L$(SAMTOOLS) -lbam -o $@

clean:
	$(RM) pindel pindel-debug $(OBJECTS) $(OBJECTS_DEBUG) $(CPPCHECKS) \
	    *.gcov *.gcda *.gcno gmon.out *.d *.o

# Pattern rule to create normal .o files	
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -Wall -Werror -g -c -O3 -fopenmp -I$(SAMTOOLS) $< -o $@

# Pattern rule to create debug version of .o files
%-debug.o: %.cpp
	$(CXX) $(CXXFLAGS) -Wall -Werror -g -c -fopenmp -g -pg \
	  -fprofile-arcs -ftest-coverage -I$(SAMTOOLS) $< -o $@

# Pattern rule to create cppcheck files	
cppcheck-result-%.xml: %.cpp
	cppcheck --quiet --enable=all --xml -I$(SAMTOOLS) $(realpath $<) 2> $@

# Pattern rule to create the .d dependency files
%.d: %.cpp
	@set -e; rm -f $@; \
	    $(CXX) -I$(SAMTOOLS) -MM $(CPPFLAGS) $< > $@.$$$$; \
	    sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	    $(RM) $@.$$$$

# Pattern rule to create the .d dependency files for the debug version
%-debug.d: %.cpp
	@set -e; rm -f $@; \
	    $(CXX) -I$(SAMTOOLS) -MM $(CPPFLAGS) $< > $@.$$$$; \
	    sed 's,\($*\)\.o[ :]*,\1-debug.o $@ : ,g' < $@.$$$$ > $@; \
	    $(RM) $@.$$$$
	
# Pseudo targets for configuration
.PHONY: all clean test
.IGNORE: clean
