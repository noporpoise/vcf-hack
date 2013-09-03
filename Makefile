
CFLAGS=-Wall -Wextra -I libs/seq_file/ -I libs/string_buffer/ -I libs/htslib/htslib -L libs/htslib/htslib
LINKING=-lhts -lpthread

ifdef DEBUG
	OPT=-O0 -g -ggdb
else
	OPT=-O2
endif

all: vcfref

vcfref: vcf_ref.c bin Makefile
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfref vcf_ref.c libs/string_buffer/string_buffer.c $(LINKING) -lz

bin:
	mkdir -p bin

clean:
	rm -rf bin/* *.greg *.dSYM

.PHONY: all clean
