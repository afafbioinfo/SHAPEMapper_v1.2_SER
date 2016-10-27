all:  parseAlignment.cpp countMutations.cpp parsed_reads.h ref_seq.h string_funcs.h
	g++ parseAlignment.cpp -o parseAlignment -O3
	g++ countMutations.cpp -o countMutations -O3

