#/bin/bash
set -beEu -o pipefail

mkdir -p data results auspice

# Follow the steps in
# https://docs.nextstrain.org/en/latest/tutorials/creating-a-phylogenetic-workflow.html
# but start with viral_usher-generated genbank.fa.xz and metadata.tsv.gz,
# skipping the augur index and augur filter steps.

# Decompress ../rsvA/genbank.fasta.xz
xzcat ../rsvA/genbank.fasta.xz > data/sequences.fasta
# Prepare metadata: use accession for strain column to match fasta,
# doctor incomplete dates to TreeTime's expected -XX
grep '^>' data/sequences.fasta | sed -e 's/^>//' > tmp
set +o pipefail
gzip -dc ../rsvA/metadata.tsv.gz \
| head -1 \
    > data/metadata.tsv
set -o pipefail
gzip -dc ../rsvA/metadata.tsv.gz \
| grep -Fwf tmp \
| awk -F$'\t' -v OFS=$'\t' '{$1 = $2; if (length($6) == 4) {$6 = $6 "-XX-XX";} else if (length($6) == 7) {$6 = $6 "-XX";} print;}' \
    >> data/metadata.tsv
rm tmp

# Copy files from ../rsvA because nextstrain docker container can see only current dir
cp -p ../rsvA/refseq.gbff data/

sequences_fasta=data/sequences.fasta
reference_name=NC_038235.1
reference_gbff=data/refseq.gbff
metadata_tsv=data/metadata.tsv
config_dir=zika-tutorial/config
auspice_json=auspice/rsv.json

time nextstrain shell . <<EOF

# Changed --reference-sequence $reference_gbff to --reference-name $reference_name because
# augur align said it could not recognize $reference_gbff as GenBank or fasta.
# Added --nthreads auto because the default is 1 thread!
time augur align \
  --sequences $sequences_fasta \
  --reference-name $reference_name \
  --output results/aligned.fasta \
  --nthreads auto \
  --fill-gaps
#MacARM: ran for 255m46.847s then got auto-killed or something, no output
#linux: real    160m58.216s
#linux: real    164m21.261s

# Added --nthreads auto
time augur tree \
  --alignment results/aligned.fasta \
  --nthreads auto \
  --output results/tree_raw.nwk
#linux: real    147m17.485s
#linux: real    167m27.089s

time augur refine \
  --tree results/tree_raw.nwk \
  --alignment results/aligned.fasta \
  --metadata $metadata_tsv \
  --output-tree results/tree.nwk \
  --output-node-data results/branch_lengths.json \
  --timetree \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4
#linux: real    211m12.688s
#linux: real    246m7.767s

# Only country, not region because viral_usher doesn't make that metadata column
time augur traits \
  --tree results/tree.nwk \
  --metadata $metadata_tsv \
  --output-node-data results/traits.json \
  --columns country \
  --confidence
# linux: real    5m2.708s
# linux: real    6m10.003s

time augur ancestral \
  --tree results/tree.nwk \
  --alignment results/aligned.fasta \
  --output-node-data results/nt_muts.json \
  --inference joint
#linux: real    14m21.355s
#linux: real    16m34.040s

time augur translate \
  --tree results/tree.nwk \
  --ancestral-sequences results/nt_muts.json \
  --reference-sequence $reference_gbff \
  --output-node-data results/aa_muts.json
#linux: real    1m1.978s
#linux: real    1m12.171s

time augur export v2 \
  --tree results/tree.nwk \
  --metadata $metadata_tsv \
  --node-data results/branch_lengths.json \
              results/traits.json \
              results/nt_muts.json \
              results/aa_muts.json \
  --colors $config_dir/colors.tsv \
  --lat-longs $config_dir/lat_longs.tsv \
  --auspice-config $config_dir/auspice_config.json \
  --output $auspice_json
#linux: real    0m58.549s
#linux: real    1m2.688s

exit
EOF
#linux: real    540m55.724s
#linux: real    602m59.173s
