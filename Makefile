
CFLAGS=-Wall -Wextra -I libs/seq_file/ -I libs/string_buffer/ \
       -I libs/htslib/htslib/ -L libs/htslib/htslib \
       -I libs/bit_array/ -L libs/bit_array/
LINKING=-lhts -lpthread
SRCS=global.c libs/string_buffer/string_buffer.c libs/bit_array/libbitarr.a

LIBS=libs/bit_array/libbitarr.a \
     libs/string_buffer/string_buffer.c libs/string_buffer/libstrbuf.a \
     libs/htslib/htslib/libhts.a libs/seq_file/seq_file.h

REQ=$(LIBS) bin Makefile

ifdef DEBUG
	OPT=-O0 -g -ggdb -DDEBUG=1
else
	OPT=-O2
endif

all: bin/vcfref bin/vcfcombine bin/vcfcombo

bin/vcfref: vcf_ref.c $(SRCS) | $(REQ)
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfref vcf_ref.c $(SRCS) $(LINKING) -lz

bin/vcfcombine: vcf_combine.c $(SRCS) | $(REQ)
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfcombine vcf_combine.c $(SRCS) $(LINKING) -lz

bin/vcfcombo: vcf_combo.c $(SRCS) | $(REQ)
	$(CC) $(CFLAGS) $(OPT) -o bin/vcfcombo vcf_combo.c $(SRCS) $(LINKING) -lz

bin/mask2vcf: mask2vcf.c $(SRCS) | $(REQ)
	$(CC) $(CFLAGS) $(OPT) -o bin/mask2vcf mask2vcf.c $(SRCS) $(LINKING) -lz

bin:
	mkdir -p bin

$(LIBS):
	cd libs; make

clean:
	rm -rf bin/* *.greg *.dSYM

.PHONY: all clean
