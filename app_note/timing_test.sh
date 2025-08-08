#!/bin/bash
source ~/.bashrc
conda activate viral_usher
set -beEu -o pipefail

# Make config files non-interactively
viral_usher init -s mumps     -r NC_002200.1 -t 2560602 -x '' -w mumps -c mumps/config.toml
viral_usher init -s measles   -r NC_001498.1 -t 11234   -x 'nextstrain/measles/N450/WHO-2012' -w measles -c measles/config.toml
viral_usher init -s oropouche -r NC_005775.1 -t 118655  -x '' -w oropouche_m -c oropouche_m/config.toml
viral_usher init -s 'RSV A'   -r NC_038235.1 -t 11250   -x 'nextstrain/rsv/a/EPI_ISL_412866' -w rsvA -c rsvA/config.toml
viral_usher init -s mpox      -r NC_003310.1 -t 10244   -x 'nextstrain/mpox/all-clades' -w mpox -c mpox/config.toml

for species in mumps measles oropouche_m rsvA mpox; do
    for i in 1 2 3 4 5; do
        echo starting $species run $i
        time viral_usher build --config $species/config.toml >& $species.$i.log
        echo done $species $i
        echo ""
    done
done
