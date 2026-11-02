#!/bin/sh
# run_bend_on_vm.sh -- check the .bend files on a Linux machine over ssh and write receipt_bend.txt.
# The logic is in run_bend.py (tested by validation/proofs/test_proofs.py):
#   - meta_laws.bend must print "All terms check." with exit status 0;
#   - a mutant passes only when Bend reports a type error at the mutant's own law AND
#     meta_laws.bend holds a refute_mutant_<X> law whose hypothesis is that law;
#   - scp, ssh, wrapper, parse, import and unknown-name failures are errors, never a
#     rejection (exit 2 for infrastructure).
# Setup (once): a Linux machine you can ssh to, with a folder holding an executable
# bend.sh that runs Bend2 v2.0.20 (https://github.com/bendlang/bend). Then:
#   BEND_HOST=<ssh host> BEND_DIR=<that folder> BEND_HOME=<HOME for bend> sh run_bend_on_vm.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
exec "$ROOT/.venv/bin/python" "$HERE/run_bend.py"
