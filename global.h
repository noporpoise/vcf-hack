#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "seq_file.h"

#define SWAP(x,y,tmp) ((tmp) = (x), (x) = (y), (y) = (tmp))
#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)
void call_die(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)))
__attribute__((noreturn));

#define warn(fmt, ...) call_warn(__FILE__, __LINE__, fmt, ##__VA_ARGS__)
void call_warn(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)));

void print_usage(const char *usage, const char *errfmt,  ...);

char parse_entire_int(char *str, int *result);

void load_reads(const char *path, read_t **reads, size_t *capcty, size_t *nchroms);

#endif /* GLOBAL_H_ */
