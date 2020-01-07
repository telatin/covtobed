#!/bin/bash
if [ ! -d "$HOME/miniconda3" ]; then
	echo "MINICONDA not found in $HOME/miniconda3"
	exit 1;
fi
if [ ! -e "$HOME/miniconda3/lib/libbamtools.a" ]; then
	echo "'$HOME/miniconda3/lib/libbamtools.a' not found"
	exit 2;
fi

set -euxo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd "$SCRIPT_DIR/.."

sed -i.bak 's/#define debug if(.*)/#define debug if(true)/' base.cpp
c++ -std=c++11 *.cpp -I${HOME}/miniconda3/include/bamtools/ -L${HOME}/miniconda3/lib/ $HOME/miniconda3/lib/libbamtools.a  \
 -o binaries/covtobed_mac_debug -lz

sed -i.bak 's/#define debug if(.*)/#define debug if(false)/' base.cpp
c++ -std=c++11 *.cpp -I${HOME}/miniconda3/include/bamtools/ -L${HOME}/miniconda3/lib/ $HOME/miniconda3/lib/libbamtools.a  \
 -o binaries/covtobed_mac -lz
cd -
