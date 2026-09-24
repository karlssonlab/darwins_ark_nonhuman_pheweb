#!/bin/bash
# Run a PheWeb dataset locally (macOS/arm64 friendly) via the container in
# docker/Dockerfile. The dataset folder is bind-mounted at /data, so the same
# image serves any dog or cat dataset -- `species` in the folder's config.py
# decides which (see pheweb/species.py).
#
# Usage:
#   ./run_local.sh build                          # build/rebuild the image
#   ./run_local.sh process <data_dir>             # pheweb process
#   ./run_local.sh serve   <data_dir> [port]      # pheweb serve  (default 8000)
#   ./run_local.sh run     <data_dir> <args...>   # any pheweb subcommand
#   ./run_local.sh recomb  <data_dir> <map_dir>   # build the recombination map
#                                                 # (only needed on a Mac; on a
#                                                 #  cluster run prep_recomb_map.py
#                                                 #  directly in the conda env)
#   ./run_local.sh shell   [data_dir]             # bash in the container
#   ./run_local.sh test                           # pytest, no dataset needed
#
# Examples:
#   ./run_local.sh build
#   ./run_local.sh process example_dog_datadir
#   ./run_local.sh serve   example_dog_datadir 8000     # dog  -> localhost:8000
#   ./run_local.sh serve   example_cat_datadir 8001     # cat  -> localhost:8001
#
# Note: the pheweb source is baked into the image (an editable install with a
# compiled cffi extension), so after editing pheweb/ re-run `./run_local.sh
# build` -- only the last two layers rebuild, it takes seconds.

#bash strict mode
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="dap-pheweb:py38"
PLATFORM="linux/amd64"

die() { echo "error: $*" >&2; exit 1; }

# `docker run -t` fails with "the input device is not a TTY" when stdout/stdin
# aren't a terminal (piped output, nohup, sbatch), so only ask for one when
# there is one. TTY_FLAGS is intentionally word-split at the call sites.
if [[ -t 0 && -t 1 ]]; then TTY_FLAGS='-it'; else TTY_FLAGS=''; fi

abs_datadir() {
    local d="${1:?a data dir is required}"
    [[ -d "$d" ]] || die "no such data dir: $d"
    (cd "$d" && pwd)
}

build() {
    docker build --platform "$PLATFORM" -f "$REPO/docker/Dockerfile" -t "$IMAGE" "$REPO"
}

have_image() { docker image inspect "$IMAGE" >/dev/null 2>&1; }

ensure_image() { have_image || { echo "image $IMAGE not found; building it first"; build; }; }

cmd="${1:-}"; shift || true

case "$cmd" in
build)
    build
    ;;

process)
    ensure_image
    datadir="$(abs_datadir "${1:-}")"
    docker run --rm $TTY_FLAGS --platform "$PLATFORM" \
        -v "$datadir:/data" -e PHEWEB_DATADIR=/data -w /data \
        "$IMAGE" pheweb process
    ;;

serve)
    ensure_image
    datadir="$(abs_datadir "${1:-}")"
    port="${2:-8000}"
    echo "serving $datadir at http://localhost:$port"
    docker run --rm $TTY_FLAGS --platform "$PLATFORM" \
        -v "$datadir:/data" -e PHEWEB_DATADIR=/data -w /data \
        -p "$port:8000" \
        "$IMAGE" pheweb serve --port 8000 --no-reloader
    ;;

run)
    ensure_image
    datadir="$(abs_datadir "${1:-}")"; shift
    [[ $# -gt 0 ]] || die "give a pheweb subcommand, e.g. 'run <data_dir> phenolist verify'"
    docker run --rm $TTY_FLAGS --platform "$PLATFORM" \
        -v "$datadir:/data" -e PHEWEB_DATADIR=/data -w /data \
        "$IMAGE" pheweb "$@"
    ;;

recomb)
    # Build the tabixed recombination map a data dir's region view serves.
    # <map_dir> holds the per-chromosome *_chr<N>_map.txt files; it is mounted
    # read-only since this only ever reads them.
    ensure_image
    datadir="$(abs_datadir "${1:-}")"
    mapdir="$(abs_datadir "${2:-}")"
    docker run --rm $TTY_FLAGS --platform "$PLATFORM" \
        -v "$datadir:/data" -v "$mapdir:/map:ro" -e PHEWEB_DATADIR=/data -w /app \
        "$IMAGE" python /app/prep_recomb_map.py "${@:3}" /map /data
    ;;

shell)
    ensure_image
    if [[ -n "${1:-}" ]]; then
        datadir="$(abs_datadir "$1")"
        docker run --rm $TTY_FLAGS --platform "$PLATFORM" \
            -v "$datadir:/data" -e PHEWEB_DATADIR=/data -w /data "$IMAGE" bash
    else
        docker run --rm $TTY_FLAGS --platform "$PLATFORM" "$IMAGE" bash
    fi
    ;;

test)
    ensure_image
    docker run --rm $TTY_FLAGS --platform "$PLATFORM" "$IMAGE" python -m pytest -q
    ;;

*)
    sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 1
    ;;
esac
