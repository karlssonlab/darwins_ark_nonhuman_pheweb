# Local dev/run environment for DAP_pheweb (dog + cat).
#
# Why a container: the pinned deps (Flask 1.1, numpy 1.19, pysam 0.16,
# gevent 21.1) have no wheels for modern Python, nor for macOS arm64 at any
# Python. Pinning the platform to linux/amd64 + cp38 gets manylinux wheels for
# all of them, so the install is reproducible. On Apple Silicon this runs
# emulated -- slower than native, fine for loading/serving a dataset locally.
FROM --platform=linux/amd64 python:3.8-slim-bullseye

# Install build tools from a PINNED snapshot of the Debian archive rather than the
# live mirrors.
#
# Why: bullseye is oldstable, and the live archive (both deb.debian.org and
# security.debian.org -- they share a Fastly CDN) has been observed serving a
# Packages index that advertises versions whose .deb files 404 in the pool. The
# failing package differs from attempt to attempt (xz-utils, libc-dev-bin,
# linux-libc-dev, libdpkg-perl), so it is not one bad file and retrying does not
# help. snapshot.debian.org serves a consistent, frozen point in time.
#
# This lay hidden for months because the layer was cached and the packages were
# never actually fetched; it surfaced the first time the build cache was cleared.
# Pinning also makes the image reproducible, which the live archive never was.
#
# Check-Valid-Until is disabled because a snapshot's Release file is, by design,
# older than apt's freshness window.
RUN printf '%s\n' \
      'deb http://snapshot.debian.org/archive/debian/20250601T000000Z bullseye main' \
      'deb http://snapshot.debian.org/archive/debian-security/20250601T000000Z bullseye-security main' \
      'deb http://snapshot.debian.org/archive/debian/20250601T000000Z bullseye-updates main' \
      > /etc/apt/sources.list \
    && echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99snapshot \
    && apt-get update && apt-get install -y --no-install-recommends \
        build-essential zlib1g-dev perl gzip \
    && rm -rf /var/lib/apt/lists/*

# setup.py is a legacy setup_requires/cffi_modules build; newer setuptools
# refuses it, so hold setuptools back.
RUN pip install --no-cache-dir 'setuptools<60' 'wheel<0.38' 'cffi~=1.15' \
        'pytest-runner~=5.2' 'pytest~=6.2'

# deps in their own layer so editing pheweb source doesn't reinstall them
COPY docker/requirements-py38.txt /tmp/requirements-py38.txt
RUN pip install --no-cache-dir -r /tmp/requirements-py38.txt

WORKDIR /app
COPY . /app
# --no-deps: everything is already installed at the pins above
RUN pip install --no-cache-dir --no-deps -e .

# The dataset folder is bind-mounted here; `species` in its config.py picks
# dog vs cat (see pheweb/species.py).
ENV PHEWEB_DATADIR=/data
EXPOSE 8000
CMD ["pheweb", "serve", "--port", "8000", "--no-reloader"]
