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

c++ -std=c++11 *.cpp -I${HOME}/miniconda3/include/bamtools/ -L${HOME}/miniconda3/lib/ $HOME/miniconda3/lib/libbamtools.a  \
 -o binaries/covtobed_mac -lz
