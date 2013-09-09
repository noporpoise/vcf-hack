
# Correct mismatches to the reference
./bin/vcfref -s tests/calls.vcf tests/ref.fa > tests/refcorrect.vcf

# Combine calls within 10bp of each other
./bin/vcfcombine 10 tests/refcorrect.vcf tests/ref.fa > tests/combined.vcf
