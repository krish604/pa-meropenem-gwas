#!/usr/bin/env bash
#
# Build the gubbins core binary from source into tools/gubbins/.
#
# WHY THIS EXISTS INSTEAD OF `conda install gubbins`
# ==================================================
# The bioconda `osx-arm64` build of gubbins 3.4.3 (build string
# `py310hdfa5cb7_1`, the newest published osx-arm64 build) ships a
# `libgubbins.0.dylib` that references the zlib symbols `_gzopen`, `_gzread`
# and `_gzclose` but links no zlib at all:
#
#     $ otool -L libgubbins.0.dylib
#     libgubbins.0.dylib:
#             @rpath/libgubbins.0.dylib
#             /usr/lib/libSystem.B.dylib          # <- no zlib
#
# libSystem does not export `_gzopen` on macOS, so the lazy-binding stub jumps
# to a null pointer and the process dies with SIGSEGV the moment it is handed an
# alignment file that exists. `run_gubbins.py` hides this behind a misleading
# "Gubbins crashed, please ensure you have enough free memory" message
# (gubbins/common.py:246 swallows the real error).
#
# gubbins' own `configure.ac` requires zlib - `PKG_CHECK_MODULES([zlib], [zlib])`
# at configure.ac:40 - so the defect is in the conda recipe, not the C++ sources.
# A source build on this machine links zlib correctly and produces byte-identical
# output to the conda binary once zlib is supplied to it.
#
# Filed upstream: https://github.com/nickjcroucher/gubbins/issues/452
#   "gubbins 3.4.3 segfaults immediately on macOS/arm64 when installed from
#    bioconda (libgubbins is not linked against zlib)"
#   OPEN as of 2026-10-02, no maintainer response yet.
#
# This script is the laptop's workaround ONLY. config/machines/bigmachine.yaml is
# Linux, where `libz` is a normal system library and the same link succeeds - the
# conda package is correct there. Do not "fix" bigmachine to match this file.
#
# The defect is macOS-specific; this build is not needed anywhere else.
#
# ---------------------------------------------------------------------------
# Homebrew build dependencies (macOS only)
# ---------------------------------------------------------------------------
#     brew install autoconf automake libtool check autoconf-archive
#
#   autoconf, automake, libtool, check
#       The four the gubbins README lists under "OSX/Linux - from source"
#       ("autoconf, libtool, gcc, check, etc...").
#   autoconf-archive
#       NOT listed in that section, and its absence is a documentation bug.
#       configure.ac:13 calls `AX_ADD_AM_MACRO_STATIC([])`, an autoconf-archive
#       macro that *generates* `aminclude_static.am`. The file is not vendored in
#       the repo (m4/ tracks only ax_pthread.m4 and ax_python_devel.m4, on
#       v3.4.3 and on master), so without this package `autoreconf -i` fails with:
#           automake: error: cannot open < aminclude_static.am
#       autoconf-archive appears only in the README's sibling
#       "installing from the repository" section. See the filed issue.
#
# Note `/usr/bin/libtool` is Apple's, not GNU libtool, and cannot drive
# autoreconf; the GNU one from Homebrew is required, and its commands live under
# an `gnubin` directory that is not on PATH by default. Handled below.
#
# `check` is keg-only, so configure never finds it and `make check` fails with
# `'check.h' file not found` unless CPPFLAGS/LDFLAGS point at it. Handled below.
#
# ---------------------------------------------------------------------------
# What this installs, and what it deliberately does not
# ---------------------------------------------------------------------------
# It builds and installs the C++ core only (`bin/gubbins` + `lib/libgubbins.*`),
# via `make -C src`. It does NOT run the top-level `make`, because
# `python/Makefile.am` hardcodes a bare `pip install .` while the README
# documents `pip3`; on macOS only /usr/bin/pip3 exists, so the top-level target
# fails with `pip: No such file or directory` *after* the C++ has already built.
# That is an upstream inconsistency and is not this script's to paper over - the
# Python wrapper is the concern of whoever wires stage 8, and is not needed to
# run the binary.
#
# The `gubbins` name on PATH is NOT shadowed: this installs into
# tools/gubbins/bin, which the laptop overlay declares under
# `paths.tool_search_dirs`. PATH stays authoritative and is searched first.

set -euo pipefail

# --- pinned version -------------------------------------------------------
# v3.4.3 to match the bioconda build the laptop config was verified against, so
# the two are the same code and the zlib linkage is the only difference. Bump
# this and re-verify output equality before changing it.
GUBBINS_TAG="v3.4.3"
GUBBINS_REPO="https://github.com/nickjcroucher/gubbins.git"

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${REPO_ROOT}/tools/gubbins"

# Homebrew's GNU libtool puts its commands in an `gnubin` dir that is not on
# PATH, and autoreconf needs `glibtoolize` from it. Derived, never hard-coded:
# this is a machine fact, and hard rule 2 says read it, never type it.
BREW_PREFIX="${BREW_PREFIX:-$(brew --prefix 2>/dev/null || true)}"
LIBTOOL_GNUBIN="${BREW_PREFIX}/opt/libtool/libexec/gnubin"
CHECK_PREFIX="${CHECK_PREFIX:-$(brew --prefix check 2>/dev/null || true)}"

die() { printf '\nERROR: %s\n' "$*" >&2; exit 1; }
step() { printf '\n=== %s\n' "$*"; }

# --- preflight ------------------------------------------------------------
[ -n "${BREW_PREFIX}" ] || die "brew not found; install the build deps listed in this script's header"

for tool in autoreconf automake make; do
    command -v "$tool" >/dev/null 2>&1 || die "$tool missing; see the Homebrew deps in this script's header"
done

if [ -d "${LIBTOOL_GNUBIN}" ]; then
    PATH="${LIBTOOL_GNUBIN}:${PATH}"
    export PATH
fi
# autoreconf's LT_INIT invokes `libtoolize`, which Homebrew installs under
# libexec/gnubin (not on PATH by default). Apple's /usr/bin/libtool is a
# different program entirely and cannot drive autoreconf, so resolve by name
# rather than assuming gnubin is laid out the same way.
command -v libtoolize >/dev/null 2>&1 \
    || die "GNU libtoolize not found; brew install libtool (Apple's /usr/bin/libtool cannot drive autoreconf)"

# autoconf-archive is easy to omit and its absence produces a confusing
# automake error about a file the repo does not ship. Fail with the real reason.
if [ -z "${CHECK_PREFIX}" ] || [ ! -d "${CHECK_PREFIX}" ]; then
    die "Homebrew 'check' not found; needed by 'make check'. brew install check"
fi
if ! ls "${BREW_PREFIX}/opt/autoconf-archive/share/aclocal/ax_add_am_macro_static.m4" >/dev/null 2>&1; then
    die "autoconf-archive not found; 'autoreconf -i' cannot generate aminclude_static.am without it.
       brew install autoconf-archive   (see the filed issue - the README omits this from the from-source deps)"
fi

# --- build in a throwaway tree -------------------------------------------
# A fresh checkout every run, so this genuinely reproduces from nothing and
# never inherits a stale object file. `--depth 1 --branch <tag>` pins to exactly
# the tag with no history and no network cost.
BUILD_DIR="$(mktemp -d "${TMPDIR:-/tmp}/gubbins-build.XXXXXXXX")"
cleanup() { rm -rf "${BUILD_DIR}"; }
trap cleanup EXIT

step "cloning ${GUBBINS_REPO} at ${GUBBINS_TAG}"
git clone --quiet --depth 1 --branch "${GUBBINS_TAG}" "${GUBBINS_REPO}" "${BUILD_DIR}/gubbins"
cd "${BUILD_DIR}/gubbins"

# Confirm we really are on the pinned tag. A moved or retagged tag would
# silently change what "v3.4.3" means, which is the one thing this script
# exists to keep honest.
ACTUAL_TAG="$(git describe --tags --exact-match 2>/dev/null || echo unknown)"
[ "${ACTUAL_TAG}" = "${GUBBINS_TAG}" ] || die "expected tag ${GUBBINS_TAG}, got '${ACTUAL_TAG}'"

step "autoreconf -i"
autoreconf -i || die "autoreconf failed. If the message is 'cannot open < aminclude_static.am', autoconf-archive is missing."

step "configure --prefix=${PREFIX}"
CONFIGURE_LOG="${BUILD_DIR}/configure.log"
./configure --prefix="${PREFIX}" >"${CONFIGURE_LOG}" 2>&1 || { cat "${CONFIGURE_LOG}" >&2; die "configure failed"; }

# Belt and braces: configure.ac:40 requires zlib, but confirm it actually found
# one rather than trusting the exit code, since a silent zlib-less configure is
# precisely the failure being worked around. Checked against configure's own
# output: `pkg-config --cflags zlib` is legitimately empty on macOS (the .pc
# file carries only `-lz`), so automake drops ZLIB_CFLAGS from the Makefile
# entirely and its absence there proves nothing.
grep -q '^checking for zlib\.\.\. yes$' "${CONFIGURE_LOG}" \
    || { cat "${CONFIGURE_LOG}" >&2; die "configure did not report 'checking for zlib... yes'; zlib was not found (see header)"; }

step "make -C src"
make -C src || die "C++ build failed"

# Fail loudly, as required. The C suite is in src/; the top-level `make check`
# additionally invokes the broken `pip install .` in python/ and would report a
# failure that has nothing to do with the code.
step "make -C src check"
CHECK_LOG="${BUILD_DIR}/make-check.log"
if ! make -C src check \
        CPPFLAGS="-I${CHECK_PREFIX}/include" \
        LDFLAGS="-L${CHECK_PREFIX}/lib -pthread" >"${CHECK_LOG}" 2>&1; then
    cat "${CHECK_LOG}" >&2
    die "make check FAILED - refusing to install a binary whose test suite does not pass"
fi
# automake writes the summary to stdout and can exit 0 in odd states; require an
# explicit pass line so a silent no-op cannot masquerade as success.
if ! grep -qE '^PASS: run_all_tests' "${CHECK_LOG}"; then
    cat "${CHECK_LOG}" >&2
    die "make check produced no 'PASS: run_all_tests'; refusing to install"
fi
sed -n '/Testsuite summary/,/^# ERROR/p' "${CHECK_LOG}" | sed 's/^/  /'

step "make -C src install -> ${PREFIX}"
rm -rf "${PREFIX}"
make -C src install || die "install failed"

# --- verify the thing we actually came here for --------------------------
# The whole point is the zlib linkage. Assert it on the installed binary so a
# future toolchain change that reintroduces the conda defect fails here, loudly,
# instead of segfaulting three stages into a run.
step "verifying installed binary"
BIN="${PREFIX}/bin/gubbins"
[ -x "${BIN}" ] || die "expected ${BIN} to exist and be executable"

if ! otool -L "${BIN}" | grep -qE 'libz'; then
    die "installed gubbins does not link zlib - this is the exact defect filed upstream.
       ${BIN} would SIGSEGV on any existing alignment file."
fi

printf '\nzlib linkage: %s\n' "$(otool -L "${BIN}" | grep -E 'libz' | sed 's/^\s*//')"
printf 'installed:    %s\n' "${BIN}"
printf 'library:      %s\n' "$(ls "${PREFIX}"/lib/libgubbins.* 2>/dev/null | tr '\n' ' ')"
printf '\nOK. Add %s/bin to the machine overlay tool_search_dirs if not already there.\n' "${PREFIX}"

cat <<'NOTE'

NOTE  The installed binary records an absolute install_name for libgubbins, so
      tools/gubbins must stay where it is, and the repo must not be moved after
      building. Re-run this script after moving either.
NOTE