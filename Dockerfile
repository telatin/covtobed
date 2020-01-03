FROM ubuntu:14.04
MAINTAINER Andrea Telatin <andrea@telatin.com>
COPY Docker/coverage /bin/coverage
COPY Docker/v2.3.0.tar.gz /
RUN mkdir /covtobed
COPY *.cpp *.h /covtobed/
RUN apt-get update && apt-get install -y software-properties-common

RUN add-apt-repository universe && apt-get update
RUN apt-get install -y build-essential cmake wget libz-dev
RUN tar xvfz "v2.3.0.tar.gz"
RUN mkdir /bamtools-2.3.0/build && cd /bamtools-2.3.0/build && cmake .. && make
RUN cp /bamtools-2.3.0/lib/libbamtools* /usr/lib/x86_64-linux-gnu/

RUN cd /covtobed && c++ -std=c++11 *.cpp  -I/bamtools-2.3.0/include/  /usr/lib/x86_64-linux-gnu/libbamtools.a   \
	 -o /bin/covtobed -lz;
