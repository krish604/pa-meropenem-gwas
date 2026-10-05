#!/usr/bin/env python
# =============================================================================
# VENDORED FROM pyseer - DO NOT EDIT THE SECTION BELOW.
#
# Upstream : https://github.com/mgalardini/pyseer/blob/master/scripts/count_patterns.py
# Retrieved: 2026-10-02, from the `master` branch.
# Licence  : Apache License 2.0, the licence pyseer 1.1.2 itself declares.
#
# sha256 of the upstream file as retrieved, whole, shebang included:
#     1081bd1c1a2409cc115facd639af0d876ba0b01ab10ccafd1cabdf926e817643
#
# sha256 of the UPSTREAM LOGIC in this file - everything from the copyright
# line below to the end, trailing newline normalised, which is exactly what
# `tests/unit/test_gwas_unique_patterns.py` recomputes:
#     2d1898a029095c82a1ce6aa095b4fa190a1e4d646291fac1aecfae4e9a34b1a7
#
# The second digest is the one CI enforces, and it is recorded here rather than
# only in the test so that the two cannot drift: the test recomputes this
# section's digest and fails if it does not match the line above. Editing the
# logic without updating this digest fails; updating the digest without editing
# the logic also fails. The first digest is for a human to check the vendoring
# once against the file it claims to come from.
#
# WHY IT IS HERE RATHER THAN INSTALLED
# -----------------------------------
# pyseer 1.1.2 is the newest version conda can solve on any channel we use;
# 1.2.0+ require `glmnet_py`, which is PyPI-only. The 1.1.2 wheel declares
# exactly five console scripts - `pyseer`, `annotate_hits_pyseer`,
# `scree_plot_pyseer`, `square_mash` and `phandango_mapper` (see
# `entry_points.txt` in the installed dist-info). This script is not among
# them: it lives in the GitHub source tree and is never packaged, and there is
# no `unique_patterns` module for it to be re-exported from either.
#
# Per pyseer's own documentation this script is where the multiple-testing
# threshold comes from, and its whole logic is
#
#     LC_ALL=C sort -u <patterns> | wc -l      # unique patterns
#     threshold = alpha / count               # Bonferroni
#
# fed by pyseer's own `--output-patterns` flag, which 1.1.2 does have. Vendoring
# it is the only way to compute the threshold pyseer documents, and
# substituting a different correction method would be a silent change of the
# science. See docs/environment-arm64.md section 4 and ticket 16.
#
# WHY THE BODY IS BYTE-FOR-BYTE UPSTREAM
# --------------------------------------
# A vendored helper that has been "tidied" is no longer the helper the upstream
# project documents, and nothing in the test suite could tell the difference.
# The only changes here are the shebang and this comment block. The pipeline
# wraps this file rather than importing it, so it runs exactly as a user who
# cloned pyseer would run it.
#
# HOW THE PIPELINE CALLS IT, AND THREE THINGS THE CALLER OWNS
# --------------------------------------------------------
# 1. `--memory`, `--cores` and `--temp` are defaults upstream
#    (`1024`, `1`, `/tmp`). They are environmental facts, so the pipeline passes
#    them from the machine overlay rather than inheriting guesses. Nothing here
#    reads the environment.
# 2. Upstream interpolates BOTH `options.patterns` (line 52 below) and
#    `options.temp` into a `shell=True` string without quoting. A path
#    containing a space silently changes the command, and one containing a
#    shell metacharacter executes it. The caller cannot fix this without editing
#    upstream, so it validates both resolved paths and refuses anything outside
#    a conservative character set. Do not relax that check here instead - fix
#    the caller's validation.
# 3. Upstream passes `--memory - mem_adjust` to `sort`, and adds
#    `--parallel=<n>` when `--cores > 1`. A small `--memory` therefore reaches
#    `sort` as a negative buffer size, and `--parallel` is not in POSIX, so not
#    every `sort` has it. The caller checks both before invoking.
# =============================================================================

# Copyright 2017 Marco Galardini and John Lees

'''Count unique patterns'''

mem_adjust = 10


def get_options():
    import argparse

    description = 'Calculate p-value threshold using Bonferroni correction'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('patterns',
                        help='File of patterns from pyseer')

    parser.add_argument('--threshold',
                        default=False,
                        action='store_true',
                        help='Only print p-value threshold')
    parser.add_argument('--alpha',
                        default=0.05,
                        type=float,
                        help='Family-wise error rate')
    parser.add_argument('--cores',
                        default=1,
                        help='Number of cores to use')
    parser.add_argument('--memory',
                        default=1024,
                        help='Maximum memory to use (in Mb)')
    parser.add_argument('--temp',
                        default='/tmp',
                        help='Directory to write tmp files to')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import sys
    import subprocess
    from decimal import Decimal

    command = "LC_ALL=C sort -u "
    if int(options.cores) > 1:
        command +=  "--parallel=" + str(options.cores)
    command += (" -S " + str(int(options.memory) - mem_adjust) + "M" +
               " -T " + options.temp +
               " " + options.patterns +
               " | wc -l")

    p = subprocess.check_output(command, shell=True, universal_newlines=True)
    if not options.threshold:
        print("Patterns:\t" + p.rstrip())
        print("Threshold:\t" + '%.2E' % Decimal(options.alpha/float(p.rstrip())))
    else:
        print('%.2E' % Decimal(options.alpha/float(p.rstrip())))
