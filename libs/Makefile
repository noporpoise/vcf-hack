

all: htslib string_buffer seq_file bit_array

htslib/Makefile:
	git clone git://github.com/samtools/htslib.git

string_buffer/Makefile:
	git clone git://github.com/noporpoise/string_buffer.git

seq_file/Makefile:
	git clone git://github.com/noporpoise/seq_file.git

bit_array/Makefile:
	git clone https://github.com/noporpoise/BitArray.git
	mv BitArray bit_array

pull:
	cd htslib; git pull
	cd string_buffer; git pull
	cd seq_file; git pull
	cd bit_array; git pull

common: pull htslib string_buffer seq_file bit_array

htslib: htslib/Makefile
	cd htslib; make

string_buffer: string_buffer/Makefile
	cd string_buffer; make

seq_file: seq_file/Makefile
	cd seq_file; make HTSLIB=../htslib

bit_array: bit_array/Makefile
	cd bit_array; make

clean:
	cd htslib; make clean
	cd string_buffer; make clean
	cd seq_file; make clean

.PHONY: all clean pull htslib string_buffer seq_file
