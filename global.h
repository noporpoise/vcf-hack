#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "seq_file.h"

#define SWAP(x,y,tmp) ((tmp) = (x), (x) = (y), (y) = (tmp))
#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) (0x1UL << (64 - __builtin_clzl(x)))
#endif

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)
void call_die(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)))
__attribute__((noreturn));

#define warn(fmt, ...) call_warn(__FILE__, __LINE__, fmt, ##__VA_ARGS__)
void call_warn(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)));

void print_usage(const char *usage, const char *errfmt,  ...);

char parse_entire_int(char *str, int *result);

size_t count_char(const char *str, char c);

// VCF: CHROM-POS-ID-REF-ALT-QUAL-FILTER-INFO-FORMAT[-SAMPLE0...] '-' is '\t'
#define VCHR  0
#define VPOS  1
#define VID   2
#define VREF  3
#define VALT  4
#define VQUAL 5
#define VFLTR 6
#define VINFO 7
#define VFRMT 8

void vcf_columns(char *vcfline, char *fields[9]);

void load_reads(const char *path, read_t **reads, size_t *capcty, size_t *nchroms);

#endif /* GLOBAL_H_ */
