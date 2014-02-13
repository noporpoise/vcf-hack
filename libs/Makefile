
all: htslib string_buffer seq_file bit_array
	cd htslib && git pull && make
	cd bit_array &&  git pull && make
	cd string_buffer && git pull && make
	cd seq_file && git pull && make

clean:
	cd htslib && make clean
	cd bit_array && make clean
	cd string_buffer && make clean
	cd seq_file && make clean
	rm -rf BitArray

.PHONY: all clean

htslib/Makefile:
	git clone git://github.com/samtools/htslib.git

string_buffer/Makefile:
	git clone git://github.com/noporpoise/string_buffer.git

seq_file/Makefile:
	git clone git://github.com/noporpoise/seq_file.git

bit_array/Makefile:
	git clone https://github.com/noporpoise/BitArray.git
	mv BitArray bit_array
