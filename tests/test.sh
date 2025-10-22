#!/bin/bash

set -beEu -o pipefail

cd ~/github/viral_usher
# Build a local docker image
time docker build --platform linux/amd64 -t viral_usher_test .

# Quick test of simple build from scratch, with segment 4 of epizootic hemorrhagic disease virus
# which is small and doesn't have very many sequences:
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -w $testdir \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# Now run build on the same config, but with --update
time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml \
    -u

# Now start over but use viral_usher_trees as a base
rm -rf $testdir
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -w $testdir \
    --use_viral_usher_trees \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# Again, with --update
time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml \
    -u

# Now build from scratch with extra_fasta (no metadata)
rm -rf $testdir
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -f tests/test_data/ehdv_extra.fa \
    -w $testdir \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# Again, with --update
time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml \
    -u

# Build from scratch with extra_fasta and extra_metadata
rm -rf $testdir
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -f tests/test_data/ehdv_extra.fa \
    -m tests/test_data/ehdv_extra.metadata.tsv \
    -w $testdir \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# TODO: make sure the generated metadata has the expected additional columns and
# metadata rows for user-added sequences.

# Again, with --update
time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml \
    -u

# TODO: make sure the generated metadata has the expected additional columns and
# metadata rows for user-added sequences.

# Build from scratch: PRRSV2 with nextclade annotations, no extra fasta
rm -rf $testdir
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_038291.1 \
    -s "PRRSV2" \
    -t 2499685 \
    -x community/isuvdl/mazeller/prrsv2/orf5/yimim2023 \
    -w $testdir \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# Again, with --update
time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml \
    -u

# Build from viral_usher_trees: PRRSV2 with nextclade annotations and extra fasta
rm -rf $testdir
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_038291.1 \
    -s "PRRSV2" \
    -t 2499685 \
    -x community/isuvdl/mazeller/prrsv2/orf5/yimim2023 \
    --use_viral_usher_trees \
    -f tests/test_data/prrsv2_extra.fa \
    -w $testdir \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# Again, with --update
time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml \
    -u

rm -rf $testdir
echo Success
