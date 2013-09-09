#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <zlib.h>

#include "global.h"
#include "seq_file.h"
#include "string_buffer.h"
#include "khash.h"

static const char usage[] =
"usage: vcfref [-s] <in.vcf[.gz]> [in.fa ...]\n"
"  Remove VCF entries that do not match the reference. Biallelic only.\n"
"  -s swaps alleles if it fixes ref mismatch\n";

KHASH_MAP_INIT_STR(ghash, read_t*)

int main(int argc, char **argv)
{
  if(argc < 2) print_usage(usage, NULL);

  char swap_alleles = 0;

  int c;
  while((c = getopt(argc, argv, "s")) >= 0) {
    switch (c) {
      case 's': swap_alleles = 1; break;
      default: die("Unknown option: %c", c);
    }
  }

  if(optind == argc) print_usage(usage, "Not enough arguments");

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

  for(i = 0; i < num_refs; i++) {
    fprintf(stderr, "Loading %s\n", refpaths[i]);
    load_reads(refpaths[i], &reads, &capacity, &nchroms);
  }

  if(num_refs == 0) {
    fprintf(stderr, "Loading from stdin\n");
    load_reads("-", &reads, &capacity, &nchroms);
  }

  if(nchroms == 0) die("No chromosomes loaded");

  for(i = 0; i < nchroms; i++) {
    r = reads + i;
    fprintf(stderr, "Loaded: '%s'\n", r->name.b);
    k = kh_put(ghash, genome, r->name.b, &hret);
    if(hret == 0) warn("Duplicate read name (taking first): %s", r->name.b);
    else kh_value(genome, k) = r;
  }

  // Now read VCF
  StrBuf line;
  strbuf_alloc(&line, 1024);
  char *fields[9];
  char *chr;
  int pos, reflen, altlen;

  while(strbuf_reset_gzreadline(&line, gzin) > 0)
  {
    if(line.buff[0] == '#') fputs(line.buff, stdout);
    else
    {
      strbuf_chomp(&line);
      vcf_columns(line.buff, fields);
      fields[1][-1] = fields[2][-1] = '\0';
      chr = line.buff;
      pos = atoi(fields[1])-1;
      k = kh_get(ghash, genome, chr);
      r = kh_value(genome, k);
      fields[1][-1] = fields[2][-1] = '\t';
      reflen = fields[4] - fields[3] - 1;
      altlen = fields[5] - fields[4] - 1;
      if(k == kh_end(genome)) warn("Cannot find chrom: %s", chr);
      else if(pos < 0) warn("Bad line: %s\n", line.buff);
      else if((reflen == 1 && altlen == 1) || fields[3][0] == fields[4][0])
      {
        if((unsigned)pos + reflen <= r->seq.end &&
           strncasecmp(r->seq.b+pos,fields[3],reflen) == 0)
        {
          fputs(line.buff, stdout);
          fputc('\n', stdout);
        }
        else if(swap_alleles && (unsigned)pos + altlen <= r->seq.end &&
                strncasecmp(r->seq.b+pos,fields[4],altlen) == 0)
        {
          // swap alleles
          char tmp[altlen], *ref = fields[3], *alt = fields[4];
          memcpy(tmp, alt, altlen);
          memmove(ref+altlen+1, ref, reflen);
          memcpy(ref, tmp, altlen);
          ref[altlen] = '\t';
          fputs(line.buff, stdout);
          fputc('\n', stdout);
        }
        // else printf("FAIL0\n");
      }
      // else printf("FAIL1\n");
    }
  }

  strbuf_dealloc(&line);
  gzclose(gzin);

  for(i = 0; i < nchroms; i++) seq_read_dealloc(reads+i);
  free(reads);

  fprintf(stderr, " Done.\n");

  return 0;
}
