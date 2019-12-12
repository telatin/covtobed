#!/bin/bash

set -euxo pipefail

c++ -std=c++11 *.cpp  -I/usr/include/bamtools /usr/lib/x86_64-linux-gnu/libbamtools.a   \
 -o binaries/covtobed -lz
