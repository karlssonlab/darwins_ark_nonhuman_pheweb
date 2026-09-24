#!/bin/sh
# Startup script for hosted deployment (SciLifeLab Serve requires one at the
# working directory).
#
# PORT defaults to 8000. Serve only permits 3000-9999, so anything outside that
# range is rejected here rather than failing obscurely at bind time.
set -eu

PORT="${PORT:-8000}"
if [ "$PORT" -lt 3000 ] || [ "$PORT" -gt 9999 ]; then
    echo "error: PORT=$PORT is outside the 3000-9999 range the host allows" >&2
    exit 1
fi

# PHEWEB_DATADIR is baked into the image (see Dockerfile.serve). Fail loudly if
# the dataset is missing rather than serving an empty site.
: "${PHEWEB_DATADIR:?PHEWEB_DATADIR is not set}"
if [ ! -f "$PHEWEB_DATADIR/pheno-list.json" ]; then
    echo "error: no pheno-list.json in $PHEWEB_DATADIR" >&2
    exit 1
fi
if [ ! -d "$PHEWEB_DATADIR/generated-by-pheweb" ]; then
    echo "error: no generated-by-pheweb/ in $PHEWEB_DATADIR -- the image has no results" >&2
    exit 1
fi

echo "serving $PHEWEB_DATADIR on port $PORT as uid $(id -u)"

# `cd` into the data dir: pheno-list.json holds relative assoc_files paths, which
# pheweb resolves against the working directory (same reason run_local.sh uses -w
# /data). --no-reloader because the reloader forks and is pointless in a container.
cd "$PHEWEB_DATADIR"
exec pheweb serve --port "$PORT" --no-reloader
