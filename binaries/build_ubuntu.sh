#!/bin/bash

if [ ! -e "/usr/lib/x86_64-linux-gnu/libbamtools.a" ]; then
	echo "Bamtools not found at: /usr/lib/x86_64-linux-gnu/libbamtools.a"
	exit 1
fi

if [ ! -e "/usr/lib/gcc/x86_64-linux-gnu/5/libstdc++.a" ]; then
	echo "Stdlib not found at: /usr/lib/gcc/x86_64-linux-gnu/5/libstdc++.a"
	exit 1
fi

set -euxo pipefail
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd "$SCRIPT_DIR/.."


# Make DEBUG
sed -i 's/#define debug if(.*)/#define debug if(true)/' base.cpp
c++ -std=c++11 ./*.cpp  -I/usr/include/bamtools /usr/lib/x86_64-linux-gnu/libbamtools.a  /usr/lib/gcc/x86_64-linux-gnu/5/libstdc++.a \
 -o binaries/covtobed_debug -lz


# Make NORMAL
sed -i 's/#define debug if(.*)/#define debug if(false)/' base.cpp
c++ -std=c++11 ./*.cpp  -I/usr/include/bamtools /usr/lib/x86_64-linux-gnu/libbamtools.a  /usr/lib/gcc/x86_64-linux-gnu/5/libstdc++.a \
 -o binaries/covtobed -lz -O3
cd -
