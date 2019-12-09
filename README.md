# covtobed

Convert

## Compile

This tool requires **libbamtools**.

### With Conda
Dynamic library:
```
c++ -std=c++11 *.cpp -I${HOME}/miniconda3/include/bamtools/ -L${HOME}/miniconda3/lib/ -lbamtools -o covtobed
```

Static library:
```
c++ -std=c++11 *.cpp -I${HOME}/miniconda3/include/bamt  path/to/libbamtools.a -o covtobed -lz
```
