#!/bin/bash
set -euxo pipefail
LD_LIBRARY_PATH=/local/miniconda3/lib ./scov1  $@
