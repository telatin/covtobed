#!/bin/bash

set -euxo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd "$SCRIPT_DIR/.."

c++ -std=c++11 *.cpp  -I/usr/include/bamtools /usr/lib/x86_64-linux-gnu/libbamtools.a  /usr/lib/gcc/x86_64-linux-gnu/5/libstdc++.a \
 -o binaries/covtobed -lz
cd -
