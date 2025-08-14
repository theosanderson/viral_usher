#/bin/bash
set -beEu -o pipefail

mkdir -p data results auspice

# Follow the steps in
# https://docs.nextstrain.org/en/latest/tutorials/creating-a-phylogenetic-workflow.html
# but start with viral_usher-generated genbank.fa.xz and metadata.tsv.gz,
# skipping the augur index and augur filter steps.

# Decompress ../measles/genbank.fasta.xz
xzcat ../measles/genbank.fasta.xz > data/sequences.fasta
# Prepare metadata: use accession for strain column to match fasta,
# doctor incomplete dates to TreeTime's expected -XX
grep '^>' data/sequences.fasta | sed -e 's/^>//' > tmp
set +o pipefail
gzip -dc ../measles/metadata.tsv.gz \
| head -1 \
    > data/metadata.tsv
set -o pipefail
gzip -dc ../measles/metadata.tsv.gz \
| grep -Fwf tmp \
| awk -F$'\t' -v OFS=$'\t' '{$1 = $2; if (length($6) == 4) {$6 = $6 "-XX-XX";} else if (length($6) == 7) {$6 = $6 "-XX";} print;}' \
    >> data/metadata.tsv
rm tmp

# Copy files from ../measles because nextstrain docker container can see only current dir
cp -p ../measles/refseq.gbff data/

sequences_fasta=data/sequences.fasta
reference_name=NC_001498.1
reference_gbff=data/refseq.gbff
metadata_tsv=data/metadata.tsv
config_dir=zika-tutorial/config
auspice_json=auspice/measles.json

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
#MacARM: real	5m41.506s

# Added --nthreads auto
time augur tree \
  --alignment results/aligned.fasta \
  --nthreads auto \
  --output results/tree_raw.nwk
#MacARM: real	4m54.421s

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
#MacARM: real	4m21.997s

# Only country, not region because viral_usher doesn't make that metadata column
time augur traits \
  --tree results/tree.nwk \
  --metadata $metadata_tsv \
  --output-node-data results/traits.json \
  --columns country \
  --confidence
#MacARM: real	0m4.766s

time augur ancestral \
  --tree results/tree.nwk \
  --alignment results/aligned.fasta \
  --output-node-data results/nt_muts.json \
  --inference joint
#MacARM: real	0m17.491s

time augur translate \
  --tree results/tree.nwk \
  --ancestral-sequences results/nt_muts.json \
  --reference-sequence $reference_gbff \
  --output-node-data results/aa_muts.json
#MacARM: real	0m2.212s

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
#	WARNING:  "Viet Nam", a value of the geographic resolution "country", appears in the tree but not in the metadata.
#		This will cause transmissions & demes involving this location not to be displayed in Auspice
#	WARNING:  The filter "region" does not appear as a property on any tree nodes.
#	WARNING:  The filter "author" does not appear as a property on any tree nodes.
#Validation of 'auspice/measles.json' succeeded, but there were warnings you may want to resolve.
#MacARM: real	0m1.601s

exit
EOF
#15m5.037s