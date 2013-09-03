#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <zlib.h>

#include "seq_file.h"
#include "string_buffer.h"
#include "khash.h"

static const char usage[] =
"usage: vcfref [-s] <in.vcf[.gz]> [in.fa ...]\n"
"  Remove VCF entries that do not match the reference. Biallelic only.\n"
"  -s swaps alleles if it fixes ref mismatch\n";

void print_usage(const char *errfmt,  ...)
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

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)
void call_die(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)))
__attribute__((noreturn));

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

KHASH_MAP_INIT_STR(ghash, read_t*)

void load_reads(const char *path, read_t **reads, size_t *capcty, size_t *nchroms)
{
  fprintf(stderr, "Loading %s\n", path);
  seq_file_t *sf;
  if((sf = seq_open(path)) == NULL) die("Cannot open file: %s\n", path);
  while(1)
  {
    if(*nchroms == *capcty && (*reads = realloc(*reads, *capcty *= 2)) == NULL)
      die("Out of memory");
    if(seq_read_alloc(*reads+*nchroms) == NULL) die("Out of memory");
    if(seq_read(sf,*reads+*nchroms) <= 0) {
      seq_read_dealloc(*reads+*nchroms);
      break;
    }
    (*nchroms)++;
  }
  seq_close(sf);
}

int main(int argc, char **argv)
{
  if(argc < 2) print_usage(NULL);

  char swap_alleles = 0;

  int c;
  while((c = getopt(argc, argv, "s")) >= 0) {
    switch (c) {
      case 's': swap_alleles = 1; break;
      default: die("Unknown option: %c", c);
    }
  }

  if(optind == argc) print_usage("Not enough arguments");

  char *inputpath = argv[optind];
  char **refpaths = argv + optind + 1;
  size_t num_refs = argc - optind - 1;

  gzFile gzin = gzopen(inputpath, "r");
  if(gzin == NULL) die("Cannot read file: %s", inputpath);

  size_t i, nchroms = 0, capacity = 1024;
  khash_t(ghash) *genome = kh_init(ghash);
  read_t *reads = malloc(capacity * sizeof(read_t)), *r;
  int hret;
  khiter_t k;

  for(i = 0; i < num_refs; i++)
    load_reads(refpaths[i], &reads, &capacity, &nchroms);

  if(num_refs == 0)
    load_reads("-", &reads, &capacity, &nchroms);

  if(nchroms == 0) die("No chromosomes loaded");

  for(i = 0; i < nchroms; i++) {
    r = reads + i;
    k = kh_put(ghash, genome, r->name.b, &hret);
    if(hret == 0) fprintf(stderr, "Warn: dup read name (taking first): %s\n", r->name.b);
    else kh_value(genome, k) = r;
  }

  // Now read VCF
  StrBuf line;
  strbuf_alloc(&line, 1024);
  char *sep0, *sep1, *sep2, *sep3, *sep4, *chr;
  int pos, reflen, altlen;

  while(strbuf_reset_gzreadline(&line, gzin) > 0)
  {
    if(line.buff[0] == '#') fputs(line.buff, stdout);
    else if((sep0 = strchr(line.buff, '\t')) != NULL &&
            (sep1 = strchr(sep0+1, '\t')) != NULL &&
            (sep2 = strchr(sep1+1, '\t')) != NULL &&
            (sep3 = strchr(sep2+1, '\t')) != NULL &&
            (sep4 = strchr(sep3+1, '\t')) != NULL)
    {
      *sep0 = *sep1 = '\0';
      chr = line.buff;
      pos = atoi(sep0+1)-1;
      k = kh_get(ghash, genome, chr);
      r = kh_value(genome, k);
      *sep0 = *sep1 = '\t';
      reflen = sep3 - sep2 - 1;
      altlen = sep4 - sep3 - 1;
      if(k == kh_end(genome)) fprintf(stderr, "Cannot find chrom: %s", chr);
      else if(pos < 0) fprintf(stderr, "Bad line: %s\n", line.buff);
      else if((reflen == 1 && altlen == 1) || sep2[1] == sep3[1])
      {
        if((unsigned)pos + reflen <= r->seq.end &&
           strncasecmp(r->seq.b+pos,sep2+1,reflen) == 0)
        {
          fputs(line.buff, stdout);
        }
        else if(swap_alleles && (unsigned)pos + altlen <= r->seq.end &&
                strncasecmp(r->seq.b+pos,sep3+1,altlen) == 0)
        {
          // swap alleles
          char tmp[altlen], *ref = sep2+1, *alt = sep3+1;
          memcpy(tmp, alt, altlen);
          memmove(ref+altlen+1, ref, reflen);
          memcpy(ref, tmp, altlen);
          ref[altlen] = '\t';
          fputs(line.buff, stdout);
        }
        // else printf("FAIL0\n");
      }
      // else printf("FAIL1\n");
    }
    else { fprintf(stderr, "Error: bad line\n%s\n", line.buff); exit(-1); }
  }

  strbuf_dealloc(&line);
  gzclose(gzin);

  for(i = 0; i < nchroms; i++) seq_read_dealloc(reads+i);
  free(reads);

  return 0;
}
