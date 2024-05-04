# Script to build (and push) docker images for ThermoCodegen.

# Note that this assumes that an appropriate multiplatform docker builder is available.
# i.e. it may be necessary to run:
# > docker buildx create --use
# first.

usage() {
    echo "Usage:" 1>&2
    echo "bash $0 [-h] [-t tag<string>] [-p platform<string>,platform<string>] [-d]" 1>&2
    echo "  -h: print this help message and exit" 1>&2
    echo "  -b: specify a base image (defaults to registry.gitlab.com/cianwilson/tcg_integrations:tcg[-debug])" 1>&2
    echo "  -t: specify a tag (defaults to registry.gitlab.com/cianwilson/tcg_slb_database:<branch>-integrations[-debug]-<platform>)" 1>&2
    echo "  -p: comma separated list of platforms (defaults to native platform and appends default tag name)" 1>&2
    echo "  -d: enable debugging (if default, suffixes tag and base image with -debug)" 1>&2
}

error() {
    usage
    exit 1
}

# realpath not available by default on macs so define it here
realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

full_path=$(realpath $0)
script_path=$(dirname $full_path)
repo_path=$(dirname $script_path)
repo=$(basename $repo_path)
build_path=$(dirname $repo_path)

# parse the arguments
TAG=''
DEBUG=false
FLAGS=''
PLATFORMS=''
BASE='registry.gitlab.com/cianwilson/tcg_integrations:tcg'
BRANCH=`git rev-parse --abbrev-ref HEAD`
TARGET='db'

while getopts "t:b:p:dh" opt; do
    case $opt in
        h )
           usage
           exit 0
           ;;
        b )
           BASE=${OPTARG}
           ;;
        t )
           TAG=${OPTARG}
           ;;
        p )
           PLATFORMS="--platform ${OPTARG}"
           ;;
        d )
           DEBUG=true
           FLAGS=' --debug'
           BASE='registry.gitlab.com/cianwilson/tcg_integrations:tcg-debug'
           TARGET='intermediate'
           ;;
        * )
           echo "ERROR: unknown option -${OPTARG}." 1>&2
           error
           ;;
    esac
done

PTAG='-userplatform'
if [ -z "$PLATFORMS" ]; then
    PROC=`uname -m`
    if [ "$PROC" == "x86_64" ]; then
        PTAG="-amd64"
    elif [ "$PROC" == "arm64" ]; then
        PTAG="-arm64"
    fi
else
    BUILDER='buildx'
fi

# if no tag is specified use default
if [ -z "$TAG" ]; then
    TAG="registry.gitlab.com/cianwilson/tcg_slb_database:${BRANCH}"
    if [ "$DEBUG" = true ]; then
        if [ "$BASE" = 'registry.gitlab.com/cianwilson/tcg_integrations:tcg-debug' ]; then
          TAG="${TAG}-integrations-debug${PTAG}"
        else
          TAG="${TAG}-userbase-debug${PTAG}"
        fi
    else
        if [ "$BASE" = 'registry.gitlab.com/cianwilson/tcg_integrations:tcg' ]; then
          TAG="${TAG}-integrations${PTAG}"
        else
          TAG="${TAG}-userbase${PTAG}"
        fi
    fi
fi

cd $build_path
docker $BUILDER build --target "$TARGET" --file $script_path/Dockerfile \
                      --build-arg BASEIMAGE="$BASE" \
                      --build-arg DIR="$repo" \
                      --build-arg FLAGS="$FLAGS" \
                      --tag "${TAG}" $PLATFORMS --push .

