#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <zlib.h>

#include "global.h"
#include "seq_file.h"
#include "string_buffer.h"
#include "bit_array.h"
#include "khash.h"

static const char usage[] =
"usage: vcfcombo <k> <in.vcf[.gz]> [in.fa ...]\n"
"  Combine variants within k bases of each other\n";

KHASH_MAP_INIT_STR(ghash, read_t*)

#define prntbf(sbuf) ({ fputs((sbuf)->buff, stdout); fputc('\n', stdout); })

typedef struct {
  char *ref;
  char **alts;
  size_t pos, reflen, num_alts, cap_alts;
} Var;

#define var_is_ins(var) ((var)->ref[0] == '\0')

static inline char var_is_del(Var *var) {
  size_t i;
  for(i = 0; i < var->num_alts; i++)
    if(var->alts[i][0] == '\0')
      return 1;
  return 0;
}

#define var_is_indel(var) (var_is_ins(var) || var_is_del(var))

static inline void var_alloc(Var *var) {
  var->cap_alts = 16;
  var->alts = malloc(var->cap_alts * sizeof(char*));
  var->num_alts = 0;
}

static inline void var_dealloc(Var *var) {
  free(var->alts);
}

static inline void var_alt_capacity(Var *var, size_t len) {
  if(len > var->cap_alts) {
    var->cap_alts = ROUNDUP2POW(len);
    var->alts = realloc(var->alts, var->cap_alts * sizeof(char*));
  }
}

// static void var_print(const Var *var) {
//   size_t i;
//   printf("ref: %s; alts: %s", var->ref, var->alts[0]);
//   for(i = 1; i < var->num_alts; i++) printf(", %s", var->alts[i]);
//   printf("; pos: %zu; reflen: %zu; num_alts: %zu\n",
//          var->pos, var->reflen, var->num_alts);
// }

static int varcmp(const void *a, const void *b) {
  const Var *v1 = (const Var*)a, *v2 = (const Var*)b;
  int cmp = v1->pos - v2->pos;
  return cmp == 0 ? v1->reflen - v2->reflen : cmp;
}

static void vars_sort(Var *vars, size_t nvars)
{
  qsort(vars, nvars, sizeof(Var), varcmp);
}

// Check if two variants are compatible (v1 must be <= v2)
int vars_compatible(const Var *v1, const Var *v2)
{
  return (v1->pos + v1->reflen <= v2->pos) &&
         (!var_is_ins(v1) || !var_is_ins(v2) || v1->pos != v2->pos);
}

void construct_genotype(const Var **vars, size_t nvars,
                        const size_t *alleles, const char *ref, size_t reflen,
                        StrBuf *out)
{
  // printf("genotype: ref: '%s' reflen: %zu\n", ref, reflen);
  size_t i, end = 0;
  strbuf_append_char(out, ',');
  for(i = 0; i < nvars; i++) {
    // if(vars[i]->pos > end) printf("{%.*s}", (int)(vars[i]->pos - end), ref + end);
    // printf("-%s-", vars[i]->alts[alleles[i]]);
    if(vars[i]->pos > end) strbuf_append_strn(out, ref+end, vars[i]->pos-end);
    strbuf_append_str(out, vars[i]->alts[alleles[i]]);
    end = vars[i]->pos + vars[i]->reflen;
  }
  // printf("%.*s\n", (int)(reflen - end), ref + end);
  strbuf_append_strn(out, ref + end, reflen - end);
}

// Returns number of genotypes added
static size_t print_genotypes(const Var **vars, size_t nvars,
                              const char *ref, size_t reflen, StrBuf *out)
{
  size_t gt, i, num_genotypes = 1, alleles[nvars];
  memset(alleles, 0, nvars * sizeof(size_t));

  for(i = 0; i < nvars; i++)
    num_genotypes *= vars[i]->num_alts;

  // printf("print_genotypes: [reflen: %zu]\n", reflen);
  // for(i = 0; i < nvars; i++)
  //   var_print(vars[i]);

  for(gt = 0; gt < num_genotypes; gt++) {
    for(i = nvars-1; i < SIZE_MAX; i--) {
      alleles[i]++;
      if(alleles[i] == vars[i]->num_alts) alleles[i] = 0;
      else break;
    }
    construct_genotype(vars, nvars, alleles, ref, reflen, out);
  }

  return num_genotypes;
}

// Try combining a given set of variants
int try_var_combination(const Var *vars, size_t nvars, BIT_ARRAY *bitset,
                        const char *ref, size_t reflen, StrBuf *out, size_t *gtcount)
{
  size_t i, j, num = 0;
  const Var *set[nvars], *prevvar = NULL;

  for(i = 0, j = nvars-1; i < nvars; i++, j--) {
    if(bit_array_get(bitset, j)) {
      if(prevvar == NULL || vars_compatible(prevvar, &vars[i])) {
        set[num++] = &vars[i];
        prevvar = &vars[i];
      }
      else return j;
    }
  }

  *gtcount += print_genotypes(set, num, ref, reflen, out);

  return -1;
}

static size_t generate_var_combinations(const Var *vars, size_t nvars,
                                        const char *ref, size_t reflen,
                                        BIT_ARRAY *bitset, StrBuf *out)
{
  size_t i, num_var_gt = 0, max = 1UL<<nvars;
  int ret;

  bit_array_resize(bitset, nvars);
  bit_array_clear_all(bitset);
  bit_array_set_bit(bitset, 0);

  for(i = 1; i < max; )
  {
    // bit_array_print_substr(bitset, 0, bit_array_length(bitset), stdout, '1', '0', 0);
    // printf("\n");
    ret = try_var_combination(vars, nvars, bitset, ref, reflen, out, &num_var_gt);
    if(ret != -1) {
      // printf("incompatible %i\n", ret);
      bit_array_add_word(bitset, ret, 1); // add 1<<ret
      i += 1UL<<ret;
    }
    else { bit_array_add(bitset, 1); i++; }
  }

  return num_var_gt;
}

static int strptrcmp(const void *a, const void *b) {
  const char *x = *(const char**)a, *y = *(const char**)b;
  return strcasecmp(x, y);
}

static void reduce_alts(char **alts, size_t num, StrBuf *out)
{
  size_t i;
  qsort(alts, num, sizeof(char**), strptrcmp);
  strbuf_append_str(out, alts[0]);
  for(i = 1; i < num; i++) {
    if(strcmp(alts[i],alts[i-1]) != 0) {
      strbuf_append_char(out, ',');
      strbuf_append_str(out, alts[i]);
    }
  }
}

// Trim matching start bases
static void var_trim_alts_starts(Var *var)
{
  size_t i, offset;
  char c;
  for(offset = 0; (c = var->ref[offset]) != '\0'; offset++) {
    for(i = 0; i < var->num_alts && var->alts[i][offset] == c; i++);
    if(i < var->num_alts) break;
  }
  if(offset > 0) {
    var->pos += offset;
    var->ref += offset;
    var->reflen -= offset;
    for(i = 0; i < var->num_alts; i++) var->alts[i] += offset;
  }
}

// Trim matching end bases
static void var_trim_alts_ends(Var *var)
{
  size_t i, trim, minlen = var->reflen, lens[var->num_alts];
  char c;
  for(i = 0; i < var->num_alts; i++) {
    lens[i] = strlen(var->alts[i]);
    minlen = MIN2(minlen, lens[i]);
  }
  for(trim = 0; trim < minlen; trim++) {
    c = var->ref[var->reflen-trim-1];
    for(i = 0; i < var->num_alts && var->alts[i][lens[i]-trim-1] == c; i++);
    if(i < var->num_alts) break;
  }
  var->reflen -= trim;
  for(i = 0; i < var->num_alts; i++) var->alts[i][lens[i]-trim] = '\0';
}

static void vcf_to_var(char *fields[9], Var *var)
{
  size_t i;
  var->num_alts = 1 + count_char(fields[VALT], ',');
  var_alt_capacity(var, var->num_alts);
  var->alts[0] = strtok(fields[VALT], ",");
  for(i = 1; i < var->num_alts; i++) var->alts[i] = strtok(NULL, ",");
  fields[VALT+1][-1] = '\0';
  var->pos = atoi(fields[VPOS]);
  var->ref = fields[VREF];
  var->reflen = fields[VALT] - fields[VREF] - 1;
}

typedef struct {
  StrBuf *lines;
  Var *vars;
  size_t nvars, cap_vars;
} VarSet;

static inline void varset_alloc(VarSet *vset)
{
  size_t i;
  vset->cap_vars = 16;
  vset->vars = malloc(vset->cap_vars * sizeof(Var));
  vset->lines = malloc(vset->cap_vars * sizeof(StrBuf));
  vset->nvars = 0;
  for(i = 0; i < vset->cap_vars; i++) {
    var_alloc(&vset->vars[i]);
    strbuf_alloc(&vset->lines[i], 1024);
  }
}

static inline void varset_dealloc(VarSet *vset) {
  size_t i;
  for(i = 0; i < vset->cap_vars; i++) {
    var_dealloc(&vset->vars[i]);
    strbuf_dealloc(&vset->lines[i]);
  }
  free(vset->vars);
  free(vset->lines);
}

static inline void varset_capacity(VarSet *vset, size_t len) {
  if(len > vset->cap_vars) {
    size_t i, oldcap = vset->cap_vars;
    vset->cap_vars = ROUNDUP2POW(len);
    vset->vars = realloc(vset->vars, vset->cap_vars * sizeof(Var));
    for(i = oldcap; i < vset->cap_vars; i++) {
      var_alloc(&vset->vars[i]);
      strbuf_alloc(&vset->lines[i], 1024);
    }
  }
}

static inline void vset_merge(VarSet *vset, BIT_ARRAY *bitset, const char *ref,
                              StrBuf *tmp, StrBuf *out)
{
  Var *var = &vset->vars[0];
  size_t i, pos, reflen, num_alts, minstart = SIZE_MAX, maxend = 0;
  char *fields[9];

  // printf(" MERGE!\n");

  for(i = 0; i < vset->nvars; i++) {
    var = &vset->vars[i];
    vcf_columns(vset->lines[i].buff, fields);
    vcf_to_var(fields, var);
    var_trim_alts_starts(var);
    var_trim_alts_ends(var);
    var->pos--; // convert to 0-based
    pos = var->pos; reflen = var->reflen;
    if(var_is_indel(var)) { pos--; reflen++; }
    minstart = MIN2(minstart, pos);
    maxend = MAX2(maxend, pos + reflen);
    // printf(" pos: %zu reflen: %zu\n", var->pos, var->reflen);
  }
  vars_sort(vset->vars, vset->nvars);

  for(i = 0; i < vset->nvars; i++) vset->vars[i].pos -= minstart;

  strbuf_reset(tmp);

  num_alts = generate_var_combinations(vset->vars, vset->nvars,
                                       ref+minstart, maxend-minstart,
                                       bitset, tmp);

  // printf("BUF: '%s'\n", tmp->buff);

  char *alts[num_alts];
  alts[0] = tmp->buff+1;
  for(i = 1; i < num_alts; i++) {
    alts[i] = strchr(alts[i-1], ',')+1;
    alts[i][-1] = '\0';
  }

  StrBuf *line = &vset->lines[vset->nvars-1];
  var = &vset->vars[vset->nvars-1];
  for(i = 1; i < var->num_alts; i++) var->alts[i][-1] = ',';
  fields[VALT+1][-1] = '\t';

  // vcf_columns(line->buff, fields);
  strbuf_reset(out);

  // Copy "CHROM-POS-ID-"
  strbuf_append_strn(out, line->buff, fields[3] - fields[0]);
  // Copy "REF-"
  strbuf_append_strn(out, ref+minstart, maxend-minstart);
  strbuf_append_char(out, '\t');
  // ALT
  reduce_alts(alts, num_alts, out);
  strbuf_append_char(out, '\t');
  // Append remaining
  strbuf_append_str(out, fields[5]);

  // printf("OUT: %s\n", out->buff);
  prntbf(out);
}

static inline void vcf_print(VarSet *vset, BIT_ARRAY *bitset,
                             khash_t(ghash) *genome, StrBuf *tmpbuf, StrBuf *out)
{
  khiter_t hpos;
  StrBuf *line = &vset->lines[0];
  if(vset->nvars > 1) {
    char *chr = line->buff;
    char *tmp = strchr(chr, '\t');
    *tmp = '\0';
    hpos = kh_get(ghash, genome, chr);
    *tmp = '\t';
    if(hpos == kh_end(genome)) {
      *tmp = '\0'; warn("Cannot find chr: %s", chr); *tmp = '\t';
      prntbf(line);
    } else {
      read_t *r = kh_value(genome, hpos);
      vset_merge(vset, bitset, r->seq.b, tmpbuf, out);
    }
  }
  else prntbf(line);
}

// ACCAT
// 1 A T
// 1 AC A
// 2 CCA C
// 4 A C

// 0 'A' 'T'
// 1 'C' ''
// 2 'CA' ''
// 3 'A' 'C','T'

// [A|T][C|]C[A|C]
// ACCA 000 ref
// ACCC 001 var2
// A-CA 010 var1
// A-CC 011 var1+var2
// TCCA 100 var0
// TCCC 101 var0+var2
// T-CA 110 var0+var1
// T-CC 111 var0+var1+var2

// ACCA 0000 ref
// ACCC 0001 var3
// AC-- 0010 var2
// xxxx 0011 var2+var3
// A-CA 0100 var1
// A-CC 0101 var1+var3
// A--- 0110 var1+var2
// xxxx 0111 var1+var2+var3
// TCCA 1000 var0
// TCCC 1001 var0+var3
// TC-- 1010 var0+var2
// xxxx 1011 var0+var2+var3
// T-CA 1100 var0+var1
// T-CC 1101 var0+var1+var3
// T--- 1110 var0+var1+var2
// xxxx 1111 var0+var1+var2+var3

void test()
{
  BIT_ARRAY bitset;
  StrBuf tmpbuf;

  bit_array_alloc(&bitset, 64);
  strbuf_alloc(&tmpbuf, 1024);

  char empty[] = "", baseA[] = "A", baseC[] = "C", baseG[] = "G", baseT[] = "T";
  char baseCA[] = "CA";
  char *alts0[1] = {baseT}, *alts1[1] = {empty};
  char *alts2[1] = {empty}, *alts3[2] = {baseC, baseT};
  (void)baseG;

  Var vars[4] = {{.ref = baseA, .alts = alts0, .reflen = 1, .pos = 0, .num_alts = 1},
                 {.ref = baseC, .alts = alts1, .reflen = 1, .pos = 1, .num_alts = 1},
                 {.ref = baseA, .alts = alts3, .reflen = 1, .pos = 3, .num_alts = 2},
                 {.ref = baseCA, .alts = alts2, .reflen = 2, .pos = 2, .num_alts = 1}};

  vars_sort(vars, 4);
  generate_var_combinations(vars, 4, "ACCA", 4, &bitset, &tmpbuf);

  printf(" ALTS: '%s'\n", tmpbuf.buff);

  bit_array_dealloc(&bitset);
  strbuf_dealloc(&tmpbuf);
}

int main(int argc, char **argv)
{
  // test();

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
  char *fields[9];
  char *nchr, *trm;
  int pos, npos, reflen, nreflen, same_chr;
  size_t chrlen, nchrlen, print;

  StrBuf tmpbuf, outbuf;
  strbuf_alloc(&tmpbuf, 1024);
  strbuf_alloc(&outbuf, 1024);

  VarSet vset;
  varset_alloc(&vset);

  BIT_ARRAY bitset;
  bit_array_alloc(&bitset, 64);

  StrBuf *line = &vset.lines[0], *nline;

  while(strbuf_reset_gzreadline(line, gzin) > 0) {
    strbuf_chomp(line);
    if(strncmp(line->buff, "##", 2) == 0) prntbf(line);
    else if(line->len > 0) break;
  }

  if(strncmp(line->buff,"#CHROM",6) != 0)
    die("Expected header: '%s'", line->buff);

  // Drop sample information from #CHROM POS ... header line
  vcf_columns(line->buff, fields);
  if((trm = strchr(fields[8], '\t')) != NULL) strbuf_shrink(line, trm-line->buff);
  prntbf(line);

  strbuf_reset_gzreadline(line, gzin);
  if(line->len == 0) die("Empty VCF");

  // Parse first VCF entry
  vcf_columns(line->buff, fields);

  fields[1][-1] = fields[2][-1] = '\0';
  pos = atoi(fields[1])-1;
  chrlen = strlen(line->buff);
  reflen = fields[4] - fields[3] - 1;
  fields[1][-1] = fields[2][-1] = '\t';

  // Drop sample information
  if((trm = strchr(fields[8], '\t')) != NULL)
    strbuf_shrink(line, trm-line->buff);

  vset.nvars = 1;

  // VCF fields: CHROM POS ID REF ALT ...
  while(1)
  {
    varset_capacity(&vset, vset.nvars+1);
    line = &vset.lines[0];
    nline = &vset.lines[vset.nvars];

    if(strbuf_reset_gzreadline(nline, gzin) <= 0) break;

    print = 0;
    strbuf_chomp(nline);
    vcf_columns(nline->buff, fields);

    fields[1][-1] = fields[2][-1] = '\0';
    nchr = nline->buff;
    npos = atoi(fields[1])-1;
    nchrlen = strlen(nchr);
    nreflen = fields[4] - fields[3] - 1;
    fields[1][-1] = fields[2][-1] = '\t';
    
    // Drop sample information
    if((trm = strchr(fields[8], '\t')) != NULL)
      strbuf_shrink(nline, trm-nline->buff);

    if(npos < 0) die("Bad line: %s\n", nline->buff);

    same_chr = (chrlen == nchrlen && strncmp(nchr, line->buff, nchrlen) == 0);
    if(same_chr && pos > npos) die("VCF not sorted: %s", nline->buff);
    if(same_chr && npos - (pos+reflen-1) <= overlap) {
      // Overlap - merge
      // printf("OVERLAP\n");
      vset.nvars++;
    }
    else
    {
      // No overlap -> print buffered lines
      // printf("PRINT\n");

      vcf_print(&vset, &bitset, genome, &tmpbuf, &outbuf);

      // next line become current line
      StrBuf swapbuf;
      SWAP(vset.lines[0], vset.lines[vset.nvars], swapbuf);
      vset.nvars = 1;
      chrlen = nchrlen;
      pos = npos;
      reflen = nreflen;
    }
  }

  // Print last line
  vcf_print(&vset, &bitset, genome, &tmpbuf, &outbuf);

  gzclose(gzin);

  kh_destroy(ghash, genome);
  bit_array_dealloc(&bitset);
  varset_dealloc(&vset);

  strbuf_dealloc(&tmpbuf);
  strbuf_dealloc(&outbuf);

  for(i = 0; i < nchroms; i++) seq_read_dealloc(reads+i);
  free(reads);

  fprintf(stderr, " Done.\n");

  return 0;
}
