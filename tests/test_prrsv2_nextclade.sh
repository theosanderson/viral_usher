# Build from scratch: PRRSV2 with nextclade annotations, no extra fasta
set -beEu -o pipefail

workdir=$1

testdir=$(realpath $(dirname "${BASH_SOURCE[0]}"))
docker_image=${VIRAL_USHER_DOCKER_IMAGE:-ghcr.io/theosanderson/viral_usher:development}

mkdir -p $workdir
workdir=$(realpath $workdir)

time viral_usher init \
    -r NC_038291.1 \
    -s "PRRSV2" \
    -t 2499685 \
    -x community/isuvdl/mazeller/prrsv2/orf5/yimim2023 \
    -w $workdir \
    -c $workdir/config.toml

time viral_usher build \
    -d $docker_image \
    -c $workdir/config.toml

# Again, with --update
time viral_usher build \
    -d $docker_image \
    -c $workdir/config.toml \
    -u

# Check outputs exist
cd $workdir
source $testdir/check_output_files.sh

