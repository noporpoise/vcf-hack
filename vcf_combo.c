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
  StrBuf line;
  char *fields[9], *ref, **alts;
  size_t pos, reflen, num_alts, cap_alts;
} Var;

typedef struct {
  Var *vars;
  size_t nvars, cap_vars;
} VarSet;

#define var_is_ins(var) ((var)->ref[0] == '\0')

static inline char var_is_del(Var *var) {
  size_t i;
  for(i = 0; i < var->num_alts; i++)
    if(var->alts[i][0] == '\0')
      return 1;
  return 0;
}

#define var_is_indel(var) (var_is_ins(var) || var_is_del(var))

// To be a SNP all alleles must be exactly one base long
static inline int alts_are_snps(char **alts, size_t num_alts)
{
  size_t i;
  for(i = 0; i < num_alts; i++)
    if(alts[i][0] == '\0' || alts[i][1] != '\0') return 0;
  return 1;
}

static inline void var_alloc(Var *var) {
  strbuf_alloc(&var->line, 512);
  var->cap_alts = 16;
  var->alts = malloc(var->cap_alts * sizeof(char*));
  var->num_alts = 0;
}

static inline void var_dealloc(Var *var) {
  strbuf_dealloc(&var->line);
  free(var->alts);
}

static inline void var_alt_capacity(Var *var, size_t len) {
  if(len > var->cap_alts) {
    var->cap_alts = ROUNDUP2POW(len);
    var->alts = realloc(var->alts, var->cap_alts * sizeof(char*));
  }
}

static void var_construct(Var *var)
{
  size_t i; char *trm, **fields = var->fields;
  strbuf_chomp(&var->line);
  // printf("READ: %s\n", var->line.buff);
  vcf_columns(var->line.buff, fields);
  for(i = 1; i < VFRMT; i++) fields[i][-1] = '\0';

  var->num_alts = 1 + count_char(fields[VALT], ',');
  var_alt_capacity(var, var->num_alts);
  var->alts[0] = strtok(fields[VALT], ",");
  for(i = 1; i < var->num_alts; i++) var->alts[i] = strtok(NULL, ",");

  int pos = atoi(fields[VPOS]);
  var->pos = pos - 1;
  var->ref = fields[VREF];
  var->reflen = fields[VALT] - fields[VREF] - 1;

  if(pos == 0) {
    for(i = 1; i < VFRMT; i++) fields[i][-1] = '\t';
    die("Bad line: %s\n", var->line.buff);
  }

  // Drop sample information
  if((trm = strchr(fields[8], '\t')) != NULL)
    strbuf_shrink(&var->line, trm - var->line.buff);
}

// returns 1 or 0
static int vars_overlap(Var *v0, Var *v1, size_t overlap)
{
  size_t len0 = v0->fields[VPOS] - v0->fields[VCHR] - 1;
  size_t len1 = v1->fields[VPOS] - v1->fields[VCHR] - 1;
  int same_chr = (len0 == len1 && !strncmp(v0->fields[VCHR], v1->fields[VCHR], len0));
  if(same_chr && v0->pos > v1->pos) die("VCF not sorted: %s", v1->line.buff);
  return (same_chr && v0->pos + v0->reflen + overlap - 1 >= v1->pos);
}

#ifdef DEBUG
static void var_print(const Var *var) {
  size_t i;
  printf("ref: %s; alts: %s", var->ref, var->alts[0]);
  for(i = 1; i < var->num_alts; i++) printf(", %s", var->alts[i]);
  printf("; pos: %zu; reflen: %zu; num_alts: %zu\n",
         var->pos, var->reflen, var->num_alts);
}
#endif

static int varcmp(const Var *v1, const Var *v2) {
  int cmp = (long)v1->pos - v2->pos;
  return cmp == 0 ? (long)v1->reflen - (long)v2->reflen : cmp;
}

// Order by ref position then by ref length
static int varcmp2(const void *a, const void *b)
{
  return varcmp((const Var*)a, (const Var*)b);
}

static void vars_sort(Var *vars, size_t nvars)
{
  qsort(vars, nvars, sizeof(Var), varcmp2);
}

// Check if two variants are compatible (v1 must be <= v2)
int vars_compatible(const Var *v1, const Var *v2)
{
  return (v1->pos + v1->reflen <= v2->pos) &&
         (!var_is_ins(v1) || !var_is_ins(v2) || v1->pos != v2->pos);
}

// Returns 1 if var contains allele, 0 otherwise
int var_contains_allele(const Var *var, char *allele)
{
  size_t i;
  if(strcmp(var->ref, allele) == 0) return 1;
  for(i = 0; i < var->num_alts; i++)
    if(strcmp(var->alts[i], allele) == 0)
      return 1;
  return 0;
}

void vars_merge(Var *dst, const Var *src)
{
  size_t i;
  for(i = 0; i < src->num_alts; i++) {
    if(!var_contains_allele(dst, src->alts[i])) {
      var_alt_capacity(dst, dst->num_alts+1);
      dst->alts[dst->num_alts++] = src->alts[i];
    }
  }
}

static inline void copy_from_ref(StrBuf *sbuf, const char *ref, size_t len)
{
  strbuf_append_strn(sbuf, ref, len);
  char *end = sbuf->buff + sbuf->len, *ptr = end - len;
  for(; ptr < end; ptr++) *ptr = toupper(*ptr);
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
    if(vars[i]->pos > end) copy_from_ref(out, ref+end, vars[i]->pos-end);
    strbuf_append_str(out, vars[i]->alts[alleles[i]]);
    end = vars[i]->pos + vars[i]->reflen;
  }
  // printf("%.*s\n", (int)(reflen - end), ref + end);
  copy_from_ref(out, ref + end, reflen - end);
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
                        const char *ref, size_t reflen,
                        StrBuf *out, size_t *gtcount)
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

  strbuf_reset(out);

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
  return strcmp(x, y);
}

// Padding base is -1 or base char
static void reduce_alt_strings(char **alts, size_t num,
                               int padding_base, StrBuf *out)
{
  size_t i;
  qsort(alts, num, sizeof(char**), strptrcmp);
  if(padding_base != -1) strbuf_append_char(out, padding_base);
  strbuf_append_str(out, alts[0]);
  for(i = 1; i < num; i++) {
    if(strcmp(alts[i],alts[i-1]) != 0) {
      strbuf_append_char(out, ',');
      if(padding_base != -1) strbuf_append_char(out, padding_base);
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

static void var_sort_alts(Var *var) {
  qsort(var->alts, var->num_alts, sizeof(char*), strptrcmp);
}

static void var_remove_dup_alts(Var *var)
{
  size_t i; char *tmp;
  for(i = 0; i < var->num_alts; ) {
    if(strcmp(var->alts[i], var->ref) == 0 ||
       (i > 1 && strcmp(var->alts[i], var->alts[i-1]) == 0))
    {
      var->num_alts--;
      SWAP(var->alts[i], var->alts[var->num_alts], tmp);
    }
    else i++;
  }
}

//
// VarSet set of Vars
//

static inline void varset_alloc(VarSet *vset)
{
  size_t i;
  vset->cap_vars = 16;
  vset->vars = malloc(vset->cap_vars * sizeof(Var));
  vset->nvars = 0;
  for(i = 0; i < vset->cap_vars; i++) var_alloc(&vset->vars[i]);
}

static inline void varset_dealloc(VarSet *vset) {
  size_t i;
  for(i = 0; i < vset->cap_vars; i++) var_dealloc(&vset->vars[i]);
  free(vset->vars);
}

static inline void varset_capacity(VarSet *vset, size_t len) {
  if(len > vset->cap_vars) {
    size_t i, oldcap = vset->cap_vars;
    vset->cap_vars = ROUNDUP2POW(len);
    vset->vars = realloc(vset->vars, vset->cap_vars * sizeof(Var));
    for(i = oldcap; i < vset->cap_vars; i++) var_alloc(&vset->vars[i]);
  }
}

// Remove duplicates: same pos, same alts
// vset->vars should be sorted first with:
//   vars_sort(vset->vars, vset->nvars);
// result is vset is merged duplicate variants
static inline void varset_remove_duplicates(VarSet *vset)
{
  size_t i; Var tmp;
  for(i = 1; i < vset->nvars; ) {
    if(varcmp(&vset->vars[i], &vset->vars[i-1]) == 0)
    {
      vset->nvars--;
      SWAP(vset->vars[i], vset->vars[vset->nvars], tmp);
      vars_merge(&vset->vars[i], &vset->vars[vset->nvars]);
    }
    else i++;
  }
}

static inline void varset_dump(const VarSet *vset)
{
  size_t i, v; Var *var;
  for(v = 0; v < vset->nvars; v++) {
    var = &vset->vars[v];
    for(i = 1; i < VFRMT; i++) var->fields[i][-1] = '\t';
    prntbf(&var->line);
  }
}

static inline void varset_print(VarSet *vset, khash_t(ghash) *genome,
                                BIT_ARRAY *bitset, StrBuf *tmp, StrBuf *out)
{
  size_t i, num_alts, minstart = SIZE_MAX, maxend = 0;
  char *ref;
  Var *var = &vset->vars[0];
  khiter_t hpos;

  if(vset->nvars == 1) {
    varset_dump(vset);
    return;
  }

  // Find reference chromosome
  if((hpos = kh_get(ghash, genome, var->fields[VCHR])) == kh_end(genome))
  {
    warn("Cannot find chr: %s", var->fields[VCHR]);
    varset_dump(vset);
    return;
  }
  else {
    read_t *r = kh_value(genome, hpos);
    ref = r->seq.b;
  }

  #ifdef DEBUG
  printf(" MERGE! [nvars=%zu]\n", vset->nvars);
  #endif

  for(i = 0; i < vset->nvars; i++)
  {
    var = &vset->vars[i];
    var_trim_alts_starts(var);
    var_trim_alts_ends(var);
    var_sort_alts(var);
    var_remove_dup_alts(var);
    minstart = MIN2(minstart, var->pos);
    maxend = MAX2(maxend, var->pos + var->reflen);
    #ifdef DEBUG
      var_print(var);
    #endif
  }

  vars_sort(vset->vars, vset->nvars);
  varset_remove_duplicates(vset);

  for(i = 0; i < vset->nvars; i++) vset->vars[i].pos -= minstart;

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

  int padding_base = -1;
  if(minstart+1 != maxend || !alts_are_snps(alts, num_alts)) {
    padding_base = minstart == 0 ? 'N' : toupper(ref[minstart-1]);
    #ifdef DEBUG
      printf("pad: %c\n", padding_base);
    #endif
  }

  strbuf_reset(out);
  var = &vset->vars[0];

  // Copy "CHROM-POS-ID-"
  int pos = minstart + 1 - (padding_base != -1);
  strbuf_sprintf(out, "%s\t%i\t%s\t", var->fields[VCHR], pos, var->fields[VID]);
  // Copy "REF-"
  if(padding_base != -1) strbuf_append_char(out, padding_base);
  copy_from_ref(out, ref+minstart, maxend-minstart);
  strbuf_append_char(out, '\t');
  // ALT
  reduce_alt_strings(alts, num_alts, padding_base, out);
  strbuf_append_char(out, '\t');
  // Append remaining
  for(i = 1; i < VFRMT; i++) var->fields[i][-1] = '\t';
  strbuf_append_str(out, var->fields[VQUAL]);

  prntbf(out);
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
  StrBuf tmpbuf, outbuf;
  strbuf_alloc(&tmpbuf, 1024);
  strbuf_alloc(&outbuf, 1024);

  VarSet vset;
  varset_alloc(&vset);

  BIT_ARRAY bitset;
  bit_array_alloc(&bitset, 64);

  while(strbuf_reset_gzreadline(&tmpbuf, gzin) > 0) {
    strbuf_chomp(&tmpbuf);
    if(strncmp(tmpbuf.buff, "##", 2) == 0) prntbf(&tmpbuf);
    else if(tmpbuf.len > 0) break;
  }

  if(strncmp(tmpbuf.buff,"#CHROM",6) != 0)
    die("Expected header: '%s'", tmpbuf.buff);

  // Drop sample information from #CHROM POS ... header tmpbuf
  char *fields[9], *trm;
  vcf_columns(tmpbuf.buff, fields);
  if((trm = strchr(fields[8], '\t')) != NULL) strbuf_shrink(&tmpbuf, trm-tmpbuf.buff);
  prntbf(&tmpbuf);

  // Parse first VCF entry
  if(strbuf_reset_gzreadline(&vset.vars[0].line, gzin) == 0) die("Empty VCF");
  var_construct(&vset.vars[0]);
  vset.nvars = 1;

  // VCF fields: CHROM POS ID REF ALT ...
  while(1)
  {
    varset_capacity(&vset, vset.nvars+1);
    Var *var = &vset.vars[0], *nvar = &vset.vars[vset.nvars];
    if(strbuf_reset_gzreadline(&nvar->line, gzin) <= 0) break;
    var_construct(nvar);

    if(vars_overlap(var, nvar, overlap)) {
      // Overlap - merge
      vset.nvars++;
    }
    else
    {
      // No overlap -> print buffered lines
      varset_print(&vset, genome, &bitset, &tmpbuf, &outbuf);

      // next line become current line
      Var swap_var;
      SWAP(*var, *nvar, swap_var);
      vset.nvars = 1;
    }
  }

  // Print last line
  varset_print(&vset, genome, &bitset, &tmpbuf, &outbuf);

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
