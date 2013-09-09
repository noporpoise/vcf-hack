
CFLAGS=-Wall -Wextra -I libs/seq_file/ -I libs/string_buffer/ -I libs/htslib/htslib/ -L libs/htslib/htslib
LINKING=-lhts -lpthread
SRCS=global.c libs/string_buffer/string_buffer.c

LIBS=libs/string_buffer/string_buffer.c libs/htslib/htslib/libhts.a libs/seq_file/seq_file.h

REQ=$(LIBS) bin Makefile

ifdef DEBUG
	OPT=-O0 -g -ggdb
else
	OPT=-O2
endif

all: bin/vcfref bin/vcfcombine

bin/vcfref: vcf_ref.c $(SRCS) $(REQ)
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfref vcf_ref.c $(SRCS) $(LINKING) -lz

bin/vcfcombine: vcf_combine.c $(SRCS) $(REQ)
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfcombine vcf_combine.c $(SRCS) $(LINKING) -lz

bin:
	mkdir -p bin

$(LIBS):
	cd libs; make

clean:
	rm -rf bin/* *.greg *.dSYM

.PHONY: all clean
