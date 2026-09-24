#!/bin/bash
# Convert raw dog inputs into the forms pheweb reads, for a dataset folder that
# already has a config.py with species = 'dog' (see example_dog_datadir/).
# Reads only the originals; writes only new files -- nothing is deleted or
# overwritten in place, so re-running is harmless.
#
# Does four things:
#   1. gene BED   4-col -> 5-col (dups gene name into the ensg slot), drops X
#   2. FASTA      one multi-FASTA -> per-chrom headerless, newline-stripped
#                 files named reference-canFam4-chrom-{1..38}.fa
#   3. sumstats   *.loco.mlma -> *.loco.csv (ref=A2, alt=A1,
#                 maf=min(Freq,1-Freq)), sorted by chrom,pos
#   4. pheno-list.json from the csv filenames (PLACEHOLDER labels -- edit them)
#
# Usage: ./prep_dog_dataset.sh <data_dir> <raw_genome.fa> <raw_genes.bed>
#
# Example (what was run for dog_compulsive_disorder_pheweb):
#   ./prep_dog_dataset.sh dog_compulsive_disorder_pheweb \
#       dog_compulsive_disorder_pheweb/UU_Cfam_GSD_1.0_ROSY.fa \
#       dog_compulsive_disorder_pheweb/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.genes.validchr.bed
set -euo pipefail

D="${1:?Usage: $0 <data_dir> <raw_genome.fa> <raw_genes.bed>}"
RAW_FA="${2:?Usage: $0 <data_dir> <raw_genome.fa> <raw_genes.bed>}"
RAW_BED="${3:?Usage: $0 <data_dir> <raw_genome.fa> <raw_genes.bed>}"
D="$(cd "$D" && pwd)"
RES="$D/generated-by-pheweb/resources"
mkdir -p "$RES" "$D/data"
export LC_ALL=C

echo "### 1/4  gene BED: 4-col -> 5-col, drop X"
# pheweb wants: chrom start end gene_name ensg (the 5th duplicates the name,
# same as make_bed.sh does for dog). Autosomes only -- these mlma files have no X.
awk -F'\t' -v OFS='\t' '$1 ~ /^[0-9]+$/ { print $1,$2,$3,$4,$4 }' \
    "$RAW_BED" \
    > "$D/UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed"
wc -l < "$D/UU_CFAM_GSD_1.0_rosy.refseq.ensformat.bed" | xargs echo "  genes kept:"

echo "### 2/4  reference FASTA: split chr1-chr38, strip headers + newlines"
# pheweb seek()s to a byte offset, so each file must be pure sequence.
awk -v res="$RES" '
    /^>/ {
        if (out != "") { close(out); out = "" }
        name = substr($1, 2)
        if (name ~ /^chr([1-9]|[1-2][0-9]|3[0-8])$/) {
            n = name; sub(/^chr/, "", n)
            out = res "/reference-canFam4-chrom-" n ".fa"
            print "  writing chrom " n > "/dev/stderr"
        }
        next
    }
    out != "" { printf "%s", $0 > out }
' "$RAW_FA"
ls "$RES" | grep -c '\.fa$' | xargs echo "  fasta files written:"

echo "### 3/4  mlma -> csv (ref=A2, alt=A1, maf=min(Freq,1-Freq)), sorted"
# Same mapping as mlma_to_csv.py, in awk+sort so it streams instead of loading
# a 9.9M-row DataFrame (and needs no pandas, which isn't a pheweb dep).
for f in "$D"/data/*.loco.mlma; do
    out="${f%.mlma}.csv"
    echo "  $(basename "$f") -> $(basename "$out")"
    { echo 'chrom,pos,ref,alt,pval,beta,maf'
      awk -F'\t' 'NR > 1 && $1 ~ /^[0-9]+$/ {
              maf = ($6 < 1 - $6) ? $6 : 1 - $6
              print $1 "," $3 "," $5 "," $4 "," $9 "," $7 "," maf
          }' "$f" \
      | sort -t, -k1,1n -k2,2n -S 2G
    } > "$out"
    wc -l < "$out" | xargs echo "    rows (incl. header):"
done

echo "### 4/4  pheno-list.json"
python3 - "$D" <<'PY'
import json, os, re, sys
d = sys.argv[1]
labels = {'F1': 'Factor 1', 'F2': 'Factor 2', 'F3': 'Factor 3'}
plist = []
for name in sorted(os.listdir(os.path.join(d, 'data'))):
    if not name.endswith('.loco.csv'):
        continue
    code = name.split('_')[0]                       # F1 / F2 / F3
    n = int(re.search(r'_N-(\d+)', name).group(1))  # sample size from filename
    plist.append({'assoc_files': ['data/' + name], 'phenocode': code,
                  'phenostring': labels.get(code, code),
                  'category': 'Compulsive Disorder', 'num_samples': n})
with open(os.path.join(d, 'pheno-list.json'), 'w') as f:
    json.dump(plist, f, indent=4)
print('  ' + json.dumps(plist, indent=4).replace('\n', '\n  '))
PY

echo "### done"
