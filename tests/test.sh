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

# Now build from scratch with extra_fasta (no metadata)
rm -rf $testdir
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -f tests/test_data/ehdv_extra.fasta \
    -w $testdir \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# Build from scratch with extra_fasta and extra_metadata
rm -rf $testdir
testdir=$(mktemp -d -t viral_usher_test.XXXXXX)
time viral_usher init \
    -r NC_013399.1 \
    -s "epizootic hemorrhagic disease virus" \
    -t 3431252 \
    -f tests/test_data/ehdv_extra.fasta \
    -m tests/test_data/ehdv_extra.metadata.tsv \
    -w $testdir \
    -c $testdir/config.toml

time viral_usher build \
    -d viral_usher_test \
    -c $testdir/config.toml

# TODO: make sure the generated metadata has the expected additional columns and
# metadata rows for user-added sequences.

rm -rf $testdir
echo Success
