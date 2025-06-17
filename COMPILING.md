# Compiling

this project requires libbamtools in addition to the more commonly available zlib.

```bash
c++ -std=c++11 *.cpp  -I/path/to/libbamtools/  -lbamtools \
	 -o covtobed -lz
```
