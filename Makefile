
CFLAGS=-Wall -Wextra -I libs/seq_file/ -I libs/string_buffer/ -I libs/htslib/htslib -L libs/htslib/htslib
LINKING=-lhts -lpthread
SRCS=global.c libs/string_buffer/string_buffer.c

ifdef DEBUG
	OPT=-O0 -g -ggdb
else
	OPT=-O2
endif

all: bin/vcfref bin/vcfcombine

bin/vcfref: vcf_ref.c bin Makefile
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfref vcf_ref.c $(SRCS) $(LINKING) -lz

bin/vcfcombine: vcf_combine.c bin Makefile
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfcombine vcf_combine.c $(SRCS) $(LINKING) -lz

bin:
	mkdir -p bin

clean:
	rm -rf bin/* *.greg *.dSYM

.PHONY: all clean
