# Automatic test

[![Build Status](https://travis-ci.org/telatin/covtobed.svg?branch=master)](https://travis-ci.org/telatin/covtobed)

This repository is continuosly tested with [TravisCI](https://travis-ci.org/telatin/covtobed), using a `test.sh` script 
that is also callable by the user. 

The `test.sh` script will use a binary called "covtobed" in the root of the repository (that is automatically built
by TravisCI, or that you should produce if manually compiling the source), 
if this is not found it will copy a pre-compiled binary with statically linked libraries from the "binaries" directory. 

[Read more in the Wiki](https://github.com/telatin/covtobed/wiki/)
