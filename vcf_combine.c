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
static int merge_vcf_lines(StrBuf *line0, char *fields1[9],
                           StrBuf *tmpbuf, StrBuf *out, read_t *r)
{
  char *fields0[9];
  int i, pos0, pos1, reflen, reflen0, reflen1;

  strbuf_reset(out);
  vcf_columns(line0->b, fields0);

  // Copy "CHROM-POS-ID-"
  strbuf_append_strn(out, line0->b, fields0[3] - fields0[0]);

  // Convert tabs to NUL
  for(i = 1; i < 9; i++) fields0[i][-1] = fields1[i][-1] = '\0';

  if(!parse_entire_int(fields0[1],&pos0))
    die("Invalid entry: %s:%s", fields0[0], fields0[1]);
  if(!parse_entire_int(fields1[1],&pos1))
    die("Invalid entry: %s:%s", fields1[0], fields1[1]);

  reflen0 = fields0[4]-fields0[3]-1;
  reflen1 = fields1[4]-fields1[3]-1;
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
  reduce_alts(tmpbuf->b+1, out);
  strbuf_append_char(out, '\t');

  // Revert
  for(i = 1; i < 9; i++) fields0[i][-1] = fields1[i][-1] = '\t';

  // Append remaining
  strbuf_append_str(out, fields0[5]);

  return reflen;
}

int main(int argc, char **argv)
{
  // compiler complains about unused function without these linese
  (void)kh_clear_ghash;
  (void)kh_del_ghash;

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
    if(hret == 0) warn("Duplicate read name (taking first): %s", r->name.b);
    else kh_value(genome, hpos) = r;
  }

  // Now read VCF
  StrBuf sbuf0, sbuf1, sbuftmp0, sbuftmp1; // line and next line
  StrBuf *line, *nline, *tmpbuf, *tmpout, *swap_buf;
  char *fields[9];
  char *nchr, *trm;
  int pos, npos, reflen, nreflen, same_chr;
  size_t chrlen, nchrlen, print;

  strbuf_alloc(&sbuf0, 1024);
  strbuf_alloc(&sbuf1, 1024);
  strbuf_alloc(&sbuftmp0, 1024);
  strbuf_alloc(&sbuftmp1, 1024);
  line = &sbuf0;
  nline = &sbuf1;
  tmpbuf = &sbuftmp0;
  tmpout = &sbuftmp1;

  #define prntbf(sbuf) ({ fputs((sbuf)->b, stdout); fputc('\n', stdout); })

  while(strbuf_reset_gzreadline(line, gzin) > 0) {
    strbuf_chomp(line);
    if(strncmp(line->b, "##", 2) == 0) prntbf(line);
    else if(line->end > 0) break;
  }

  if(strncmp(line->b,"#CHROM",6) != 0)
    die("Expected header: '%s'", line->b);

  // Drop sample information from #CHROM POS ... header line
  vcf_columns(line->b, fields);
  if((trm = strchr(fields[8], '\t')) != NULL) strbuf_shrink(line, trm-line->b);
  prntbf(line);

  strbuf_reset_gzreadline(line, gzin);
  if(line->end == 0) die("Empty VCF");

  // Parse first VCF entry
  vcf_columns(line->b, fields);

  fields[1][-1] = fields[2][-1] = '\0';
  pos = atoi(fields[1])-1;
  chrlen = strlen(line->b);
  fields[1][-1] = fields[2][-1] = '\t';
  reflen = fields[4] - fields[3] - 1;

  // Drop sample information
  if((trm = strchr(fields[8], '\t')) != NULL) strbuf_shrink(line, trm-line->b);

  // VCF fields: CHROM POS ID REF ALT ...
  while(strbuf_reset_gzreadline(nline, gzin) > 0)
  {
    print = 0;
    strbuf_chomp(nline);
    vcf_columns(nline->b, fields);

    fields[1][-1] = fields[2][-1] = '\0';
    nchr = nline->b;
    npos = atoi(fields[1])-1;
    nchrlen = strlen(nchr);
    hpos = kh_get(ghash, genome, nchr);
    fields[1][-1] = fields[2][-1] = '\t';
    nreflen = fields[4] - fields[3] - 1;
    
    // Drop sample information
    if((trm = strchr(fields[8], '\t')) != NULL)
      strbuf_shrink(nline, trm-nline->b);

    if(hpos == kh_end(genome)) { warn("Cannot find chr: %s", nchr); print = 1; }
    else if(npos < 0) { warn("Bad line: %s", nline->b); print = 1; }
    else
    {
      same_chr = (chrlen == nchrlen && strncmp(nchr, line->b, nchrlen) == 0);
      if(same_chr && pos > npos) die("VCF not sorted: %s", nline->b);
      if(same_chr && npos - (pos+reflen-1) <= overlap) {
        // Overlap - merge
        r = kh_value(genome, hpos);
        reflen = merge_vcf_lines(line, fields, tmpbuf, tmpout, r);
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
      chrlen = nchrlen;
      pos = npos;
      reflen = nreflen;
    }
  }

  // Print last line
  prntbf(line);

  kh_destroy(ghash, genome);
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
