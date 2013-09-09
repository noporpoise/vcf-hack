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
"usage: vcfcombine <k> <in.vcf[.gz]> [in.fa ...]\n"
"  Combine variants within k bases of each other\n";

KHASH_MAP_INIT_STR(ghash, read_t*)

// ref:"G" alts:"A,T" offset: 1; rlen:1; ref: "TGA"; out: "TAA,TTA"
// rlen is the number of bases in the ref
static void merge_alts(char *alts, size_t offset, size_t rlen,
                       const char *ref, size_t mergelen, StrBuf *out)
{
  char *str = strtok(alts, ",");
  size_t len;

  do {
    len = strlen(str);
    strbuf_append_char(out, ',');
    strbuf_append_strn(out, ref, offset);
    strbuf_append_strn(out, str, len);
    strbuf_append_strn(out, ref+offset+rlen, mergelen-offset-rlen);
  } while((str = strtok(NULL, ",")) != NULL);
}

static int strptrcmp(const void *a, const void *b) {
  return strcasecmp(a, b);
}

// Remove duplicate alternative alleles
static void reduce_alts(char *alts, StrBuf *out)
{
  char *str;
  size_t i, num = 1;

  str = alts;
  while((str = strchr(str+1, ',')) != NULL) num++;

  char *alleles[num];
  alleles[0] = str = alts;
  for(i = 1; i < num; i++) {
    str = strchr(str+1, ',');
    *str = '\0';
    alleles[i] = str+1;
  }
  qsort(alleles, num, sizeof(char*), strptrcmp);

  strbuf_append_str(out, alleles[0]);

  for(i = 1; i < num; i++) {
    if(strcmp(alleles[i],alleles[i-1]) != 0) {
      strbuf_append_char(out, ',');
      strbuf_append_str(out, alleles[i]);
    }
  }
}

// Merge line0, line1 into out. Returns reflen [ e.g. strlen(REF) ]
// VCF: CHROM-POS-ID-REF-ALT-QUAL-FILTER-INFO-FORMAT[-SAMPLE0...] '-' is '\t'
static int merge_vcf_lines(StrBuf *line0, StrBuf *line1,
                           StrBuf *tmpbuf, StrBuf *out, read_t *r)
{
  strbuf_reset(out);

  char *fields0[8], *fields1[8];
  int i, pos0, pos1, reflen, reflen0, reflen1, altlen0, altlen1;
  fields0[0] = line0->buff;
  fields1[0] = line1->buff;
  for(i = 1; i < 8 && (fields0[i] = strchr(fields0[i-1]+1, '\t')) != NULL; i++) fields0[i]++;
  if(i != 8) die("Invalid line: %s", line0->buff);
  for(i = 1; i < 8 && (fields1[i] = strchr(fields1[i-1]+1, '\t')) != NULL; i++) fields1[i]++;
  if(i != 8) die("Invalid line: %s", line1->buff);

  // Copy "CHROM-POS-ID-"
  out->len = 0; out->buff[0] = '\0';
  strbuf_append_strn(out, line0->buff, fields0[3] - fields0[0]);

  // Convert tabs to NUL
  for(i = 1; i < 8; i++) fields0[i][-1] = fields1[i][-1] = '\0';

  if(!parse_entire_int(fields0[1],&pos0))
    die("Invalid entry: %s:%s", fields0[0], fields0[1]);
  if(!parse_entire_int(fields1[1],&pos1))
    die("Invalid entry: %s:%s", fields1[0], fields1[1]);

  reflen0 = fields0[4]-fields0[3]-1;
  reflen1 = fields1[4]-fields1[3]-1;
  altlen0 = fields0[5]-fields0[4]-1;
  altlen1 = fields1[5]-fields1[4]-1;
  reflen = MAX2(pos1+reflen1-pos0, reflen0);

  if(pos1 + reflen1 > (signed)r->seq.end)
    die("Out of bounds: %s:%s", fields1[0], fields1[1]);

  const char *ref = r->seq.b+(pos0-1);

  // Print "REF"
  strbuf_append_str(out, fields0[3]);
  strbuf_append_strn(out, ref+reflen0, reflen-reflen0);
  strbuf_append_char(out, '\t');

  // Merge alt alleles
  strbuf_reset(tmpbuf);
  merge_alts(fields0[4], 0, reflen0, ref, reflen, tmpbuf);
  merge_alts(fields1[4], pos1-pos0, reflen1, ref, reflen, tmpbuf);
  reduce_alts(tmpbuf->buff+1, out);
  strbuf_append_char(out, '\t');

  // Revert
  for(i = 1; i < 8; i++) fields0[i][-1] = fields1[i][-1] = '\t';  

  // Append remaining
  strbuf_append_str(out, fields0[5]);

  return reflen;
}

int main(int argc, char **argv)
{
  char *inputpath, **refpaths;
  gzFile gzin;
  size_t i, nchroms = 0, capacity = 1024, num_refs;
  int hret, overlap = 0;
  khiter_t hpos;

  if(argc < 3) print_usage(usage, NULL);

  int c;
  while((c = getopt(argc, argv, "")) >= 0) {
    switch (c) {
      default: die("Unknown option: %c", c);
    }
  }

  if(optind == argc) print_usage(usage, "Not enough arguments");

  if(!parse_entire_int(argv[optind], &overlap) || overlap < 0)
    die("Invalid <overlap> value: %s %i", argv[optind], overlap);

  inputpath = argv[optind+1];
  refpaths = argv + optind + 2;
  num_refs = argc - optind - 2;

  gzin = gzopen(inputpath, "r");
  if(gzin == NULL) die("Cannot read file: %s", inputpath);

  khash_t(ghash) *genome = kh_init(ghash);
  read_t *reads = malloc(capacity * sizeof(read_t)), *r;

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
    hpos = kh_put(ghash, genome, r->name.b, &hret);
    if(hret == 0) fprintf(stderr, "Warn: dup read name (taking first): %s\n", r->name.b);
    else kh_value(genome, hpos) = r;
  }

  // Now read VCF
  StrBuf sbuf0, sbuf1, sbuftmp0, sbuftmp1; // line and next line
  StrBuf *line, *nline, *tmpbuf, *tmpout, *swap_buf;
  char *sep0, *sep1, *sep2, *sep3, *nchr;
  int pos, npos, reflen, nreflen, same_chr;
  size_t chrlen, nchrlen, print;

  strbuf_alloc(&sbuf0, 1024);
  strbuf_alloc(&sbuf1, 1024);
  strbuf_alloc(&sbuftmp0, 1024);
  strbuf_alloc(&sbuftmp1, 1024);
  line = &sbuf0;
  nline = &sbuf1;
  tmpbuf = &sbuftmp0;
  tmpout= &sbuftmp1;

  #define prntbf(sbuf) ({ fputs((sbuf)->buff, stdout); fputc('\n', stdout); })

  while(strbuf_reset_gzreadline(line, gzin) > 0) {
    strbuf_chomp(line);
    if(line->buff[0] == '#') prntbf(line);
    else if(line->len > 0) break;
  }

  if(line->len == 0) die("Empty VCF");

  // Parse first VCF entry
  if((sep0 = strchr(line->buff, '\t')) == NULL ||
     (sep1 = strchr(sep0+1, '\t')) == NULL ||
     (sep2 = strchr(sep1+1, '\t')) == NULL ||
     (sep3 = strchr(sep2+1, '\t')) == NULL)
  {
    die("Bad VCF line\n%s\n", nline->buff);
  }

  *sep0 = *sep1 = '\0';
  pos = atoi(sep0+1)-1;
  chrlen = strlen(line->buff);
  *sep0 = *sep1 = '\t';
  reflen = sep3 - sep2 - 1;

  // VCF fields: CHROM POS ID REF ALT ...
  while(strbuf_reset_gzreadline(nline, gzin) > 0)
  {
    print = 0;
    strbuf_chomp(nline);

    if((sep0 = strchr(nline->buff, '\t')) == NULL ||
       (sep1 = strchr(sep0+1, '\t')) == NULL ||
       (sep2 = strchr(sep1+1, '\t')) == NULL ||
       (sep3 = strchr(sep2+1, '\t')) == NULL)
    {
      die("Bad VCF line\n%s\n", nline->buff);
    }

    *sep0 = *sep1 = '\0';
    nchr = nline->buff;
    npos = atoi(sep0+1)-1;
    nchrlen = strlen(nchr);
    hpos = kh_get(ghash, genome, nchr);
    *sep0 = *sep1 = '\t';
    nreflen = sep3 - sep2 - 1;

    if(hpos == kh_end(genome)) { warn("Cannot find chr: %s", nchr); print = 1; }
    else if(npos < 0) { warn("Bad line: %s", nline->buff); print = 1; }
    else
    {
      same_chr = (strncmp(nchr, line->buff, nchrlen) == 0);
      if(same_chr && pos > npos) die("VCF not sorted: %s", nline->buff);
      if(same_chr && npos - (pos+reflen-1) <= overlap) {
        // Overlap - merge
        r = kh_value(genome, hpos);
        reflen = merge_vcf_lines(line, nline, tmpbuf, tmpout, r);
        SWAP(line, tmpout, swap_buf);
      }
      else {
        // No overlap
        print = 1;
      }
    }

    if(print) {
      prntbf(line);
      // next line become current line
      SWAP(line,nline,swap_buf);
      pos = npos;
      reflen = nreflen;
    }
  }

  // Print last line
  prntbf(line);

  strbuf_dealloc(&sbuf0);
  strbuf_dealloc(&sbuf1);
  strbuf_dealloc(&sbuftmp0);
  strbuf_dealloc(&sbuftmp1);
  gzclose(gzin);

  for(i = 0; i < nchroms; i++) seq_read_dealloc(reads+i);
  free(reads);

  fprintf(stderr, " Done.\n");

  return 0;
}
