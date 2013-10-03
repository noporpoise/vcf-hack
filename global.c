#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "global.h"

void print_usage(const char* usage, const char *errfmt,  ...)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);
    if(errfmt[strlen(errfmt)-1] != '\n') fputc('\n', stderr);
  }
  fputs(usage, stderr);
  exit(EXIT_FAILURE);
}

void call_die(const char *file, int line, const char *fmt, ...)
{
  va_list argptr;
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Error: ", file, line);
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  if(*(fmt + strlen(fmt) - 1) != '\n') fputc('\n', stderr);
  exit(EXIT_FAILURE);
}

void call_warn(const char *file, int line, const char *fmt, ...)
{
  va_list argptr;
  fflush(stdout);
  fprintf(stderr, "[%s:%i] Warning: ", file, line);
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  if(*(fmt + strlen(fmt) - 1) != '\n') fputc('\n', stderr);
  fflush(stderr);
}

char parse_entire_int(char *str, int *result)
{
  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);
  if(tmp > INT_MAX || tmp < INT_MIN || *strtol_last_char_ptr != '\0') return 0;
  *result = (int)tmp;
  return 1;
}

size_t count_char(const char *str, char c)
{
  const char *tmp = str;
  size_t count = 0;
  while((tmp = strchr(tmp, c)) != NULL) { count++; tmp++; }
  return count;
}

void vcf_columns(char *vcfline, char *fields[9])
{
  size_t i;
  fields[0] = vcfline;
  for(i = 1; i < 9; i++) {
    fields[i] = strchr(fields[i-1]+1, '\t');
    if(fields[i] == NULL) die("Invalid VCF line: %s", vcfline);
    fields[i]++;
  }
}

void load_reads(const char *path, read_t **reads, size_t *capcty, size_t *nchroms)
{
  seq_file_t *sf;
  if((sf = seq_open(path)) == NULL) die("Cannot open file: %s\n", path);
  while(1)
  {
    if(*nchroms == *capcty && (*reads = realloc(*reads, *capcty *= 2)) == NULL)
      die("Out of memory");
    read_t *r = *reads+*nchroms;
    if(seq_read_alloc(r) == NULL) die("Out of memory");
    if(seq_read(sf, r) <= 0) { seq_read_dealloc(r); break; }
    seq_read_truncate_name(r);
    (*nchroms)++;
  }
  seq_close(sf);
}
